#!/usr/bin/env python3
"""Summarize (and optionally run) the LAPACK, BLAS and CBLAS test suites.

This script analyzes the ``.out`` files written by the LAPACK testing
drivers (``xlintst*``, ``xeigtst*`` and ``xdmdeigtst*``) and prints a
summary table of the number of tests run and the number of failures per
precision (s/d/c/z).  With ``--run`` it executes the testing drivers
first and then analyzes their output.

The BLAS (``xblat[123]?``) and CBLAS (``x?cblat[123]``) test drivers are
analyzed too, from their own testing directories, and are reported in
their own summary sections.  Those drivers report their test counts in
lines of the form::

     SGEMV      COMPUTATIONAL TESTS:     3456 RUN,        0 FAILED
     SGEMV      ERROR-EXIT TESTS:           6 RUN,        0 FAILED

Computational failures are counted as numerical errors and error-exit
failures as other errors.  Output produced by a build whose drivers do
not report counts is still summarized, by counting one test per verdict.

When index-64 extended API outputs (``*_64.out``, produced by CMake
builds with ``BUILD_INDEX64_EXT_API=ON``) are present, they are analyzed
as well and reported in a separate "extended API" section so that the
default-API totals remain comparable across builds.  With
``--merge-apis`` a library whose two API variants report the same errors
is summarized in one combined section instead.

Examples:
    ./lapack_testing.py -n
        Print the numbers of failed tests by analyzing the LAPACK output.

    ./lapack_testing.py -n -r -p s
        Run the REAL precision tests, then print the numbers of failures.

    ./lapack_testing.py -n -p s -t eig
        Print the numbers of failures in REAL precision by analyzing only
        the eigenproblem test output.

    ./lapack_testing.py -t blas
        Summarize only the BLAS test output.

    ./lapack_testing.py -s --junit-xml results.xml
        Print only the summary table and also write a JUnit XML report
        of the analyzed output files, e.g. for GitLab CI test reports.
"""

from __future__ import annotations

import argparse
import io
import math
import os
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Dict, List, Optional, Sequence, TextIO, Tuple

# Precision letters and the labels used in the summary table.
PRECISIONS: "Tuple[Tuple[str, str], ...]" = (
    ("s", "REAL"),
    ("d", "DOUBLE PRECISION"),
    ("c", "COMPLEX"),
    ("z", "COMPLEX16"),
)

# Summary table label of each precision letter, e.g. "s" -> "REAL".
PRECISION_NAMES: "Dict[str, str]" = dict(PRECISIONS)

# Second precision letter of the mixed-precision linear equation tests.
MIXED_PARTNER: "Dict[str, str]" = {"d": "s", "z": "c"}

# Eigenproblem test sets using the classic ``alasum``/``alasvm`` summary
# format: (name, has shared input file, description).  Sets with a shared
# input read e.g. ``nep.in``; the others read e.g. ``sec.in``/``dec.in``.
EIG_STANDARD_SETS: "Tuple[Tuple[str, bool, str], ...]" = (
    ("nep", True, "Nonsymmetric Eigenvalue Problem"),
    ("sep", True, "Symmetric Eigenvalue Problem"),
    ("se2", True, "Symmetric Eigenvalue Problem 2-stage"),
    ("svd", True, "Singular Value Decomposition"),
    ("ec", False, "Eigen Condition"),
    ("ed", False, "Nonsymmetric Eigenvalue"),
    ("gg", False, "Nonsymmetric Generalized Eigenvalue Problem"),
    ("gd", False, "Nonsymmetric Generalized Eigenvalue Problem driver"),
    ("sb", False, "Symmetric Eigenvalue Problem"),
    ("sg", False, "Symmetric Eigenvalue Generalized Problem"),
    ("bb", False, "Banded Singular Value Decomposition routines"),
    ("glm", True, "Generalized Linear Regression Model routines"),
    ("gqr", True, "Generalized QR and RQ factorization routines"),
    ("gsv", True, "Generalized Singular Value Decomposition routines"),
    ("csd", True, "CS Decomposition routines"),
    ("lse", True, "Constrained Linear Least Squares routines"),
)

# Balancing/backtransformation test sets, which use the ``schkbl``-style
# "total number of examples tested" summary format.
EIG_BALANCE_SETS: "Tuple[Tuple[str, str], ...]" = (
    ("bal", "Matrix Balancing"),
    ("bak", "Balancing Backtransformation"),
    ("gbal", "Generalized Matrix Balancing"),
    ("gbak", "Generalized Balancing Backtransformation"),
)

# Linear equation test sets: (family, name stem, executable prefix,
# description).  Names and executables are built from the precision
# letter plus the stem, e.g. stest.in/stest.out/xlintsts; the mixed
# precision sets use the letter of both precisions (dstest, xlintstds)
# and exist only for the precisions in MIXED_PARTNER.
LIN_SETS: "Tuple[Tuple[str, str, str, str], ...]" = (
    ("lin", "test", "xlintst", "Linear Equation routines"),
    ("mixed", "test", "xlintst", "Mixed Precision linear equation routines"),
    ("rfp", "test_rfp", "xlintstrf", "RFP linear equation routines"),
)

# BLAS and CBLAS test sets, one per BLAS level: (level, description).
# The Level 1 drivers read no input file and write to standard output.
BLAS_LEVELS: "Tuple[Tuple[int, str], ...]" = (
    (1, "Level 1 BLAS routines"),
    (2, "Level 2 BLAS routines"),
    (3, "Level 3 BLAS routines"),
)

# Libraries, in reporting order.  Each has its own testing directory and
# its own section in the summary table.
LIBRARY_LAPACK = "LAPACK"
LIBRARY_BLAS = "BLAS"
LIBRARY_CBLAS = "CBLAS"
LIBRARIES: "Tuple[str, ...]" = (LIBRARY_LAPACK, LIBRARY_BLAS, LIBRARY_CBLAS)

# LAPACK test families, in reporting order per precision.
LAPACK_FAMILIES: "Tuple[str, ...]" = ("eig",) + tuple(s[0] for s in LIN_SETS) + ("dmd",)

# All test families, in reporting order per precision.
ALL_FAMILIES: "Tuple[str, ...]" = LAPACK_FAMILIES + ("blas", "cblas")

# Which library each family belongs to.
FAMILY_LIBRARY: "Dict[str, str]" = dict(
    [(family, LIBRARY_LAPACK) for family in LAPACK_FAMILIES]
    + [("blas", LIBRARY_BLAS), ("cblas", LIBRARY_CBLAS)]
)

# API suffixes that may exist: default API and index-64 extended API.
KNOWN_SUFFIXES: "Tuple[str, ...]" = ("", "_64")

RESULTS_FILENAME = "testing_results.txt"

# Classic summary lines printed by alasum.f/alasvm.f:
#   "  All tests for XYZ routines passed the threshold ( ddd tests run)"
#   "  XYZ:  ddd out of  ddd tests failed to pass the threshold"
RE_TESTS_RUN = re.compile(r"(\d+)\s+tests run\)")
RE_TESTS_FAILED = re.compile(r"(\d+)\s+out of\s+(\d+)")

# Failure records printed by the eigencondition checkers (schkec.f and
# friends), e.g. " Error in STRSYL: RMAX =..." — one per failing routine.
RE_EC_ERROR = re.compile(r"^ ?Error in \w+")

# Summary lines printed by the balancing checkers (schkbl.f and friends).
# The complex generalized checkers use slightly different wording
# ("ratio of largest test error", "ILO or IHI is wrong").
RE_EXAMPLES_TESTED = re.compile(r"total number of examples tested\s*=\s*(\d+)")
RE_INFO_NOT_ZERO = re.compile(r"number of examples where info is not 0\s*=\s*(\d+)")
RE_ILO_IHI_WRONG = re.compile(
    r"example number where ILO or IHI (?:is )?wrong\s*=\s*(\d+)"
)
RE_LARGEST_ERROR = re.compile(r"(?:value|ratio) of largest test error\s*=\s*(\S+)")

# Per-test verdict lines printed by the DMD checkers (schkdmd.f90 and
# friends), e.g. ">>>> Z - U*V test PASSED.".  The word boundary keeps
# aggregate lines such as "SGEDMD :: ALL TESTS PASSED." from matching.
RE_DMD_VERDICT = re.compile(r"\btest\s+(PASSED|FAILED)\b", re.IGNORECASE)

# Test counts reported by the BLAS/CBLAS drivers, e.g.
#   " SGEMV      COMPUTATIONAL TESTS:     3456 RUN,        0 FAILED"
#   " cblas_sgemv      ROW-MAJOR    COMPUTATIONAL TESTS:  3456 RUN, ..."
# The routine name field width differs per driver, so never match on
# column positions.
RE_BLAS_COUNTS = re.compile(
    r"^\s*\S+\s+(?:(?:COLUMN-MAJOR|ROW-MAJOR)\s+)?"
    r"(COMPUTATIONAL|ERROR-EXIT) TESTS:\s*(\d+) RUN,\s*(\d+) FAILED\s*$"
)

# Per-routine verdicts.  These are only counted when the driver did not
# report counts (output from a build without the counting instrumentation).
RE_BLAS_PASSED = re.compile(
    r"^\s*\S+\s+PASSED THE (?:(?:COLUMN-MAJOR|ROW-MAJOR)\s+)?"
    r"(?:COMPUTATIONAL TESTS|TESTS OF ERROR-EXITS)\b"
)
RE_BLAS_SUSPECT = re.compile(
    r"\bCOMPLETED THE (?:(?:COLUMN-MAJOR|ROW-MAJOR)\s+)?COMPUTATIONAL TESTS\b"
)
RE_BLAS_FAILED_COMPUTATIONAL = re.compile(r"\bFAILED ON CALL NUMBER:")
RE_BLAS_FAILED_ERROR_EXIT = re.compile(r"\bFAILED THE TESTS OF ERROR-EXITS\b")

# Driver-level breakage that the per-routine counts cannot express: the
# run was abandoned, misconfigured, or never reached its footer.
RE_BLAS_ABANDONED = re.compile(r"\*{5,7} (?:FATAL ERROR - )?TESTS ABANDONED \*{5,7}")
RE_BLAS_NOT_RECOGNIZED = re.compile(r"^\s*SUBPROGRAM NAME .* NOT RECOGNIZED")
RE_BLAS_DOT_PRODUCTS = re.compile(r"^\s*ERROR IN [SDCZ]M[VM]T?CH\b")
RE_BLAS_INTERNAL = re.compile(r"Shouldn't be here in CHECK")
RE_BLAS_INPUT_ERROR = re.compile(
    r"^\s*(?:NUMBER OF VALUES OF |VALUE OF [NK] IS LESS THAN"
    r"|ABSOLUTE VALUE OF INCX OR INCY )"
)

# Detail lines that sit behind a verdict which is already counted.  They
# are worth showing but must not be counted: the ``cblat2_64.out`` fixture
# has 91 of them behind just 17 failing routines.
RE_BLAS_DETAIL = re.compile(
    r"XERBLA WAS CALLED WITH"
    r"|ILLEGAL VALUE OF PARAMETER NUMBER"
    r"|FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE"
    r"|FATAL ERROR - PARAMETER NUMBER"
    r"|FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL"
    r"|BUT WITH MAXIMUM TEST RATIO"
    r"|WARNING: Skipping xerbla tests"
)
RE_BLAS_NOT_TESTED = re.compile(r"^\s*\S+\s+WAS NOT TESTED\s*$")

# Level 2/3 footer.  Note this is printed even when routines failed, so
# it means "not truncated", not "passed"; its absence means the driver
# died part way through.
RE_BLAS_END_OF_TESTS = re.compile(r"^\s*END OF TESTS\s*$")

# Level 1 drivers have no counts of their own in an uninstrumented build
# and no footer at all; one "Test of subprogram number" block is one
# subprogram, followed by either a PASS line or FAIL detail.
RE_BLAS_L1_CASE = re.compile(r"^\s*Test of subprogram number\s*\d+")
RE_BLAS_L1_PASS = re.compile(r"^\s*-{5} PASS -{5}\s*$")
RE_BLAS_L1_FAIL = re.compile(r"^\s*FAIL\s*$")

# Parser kinds, used by TestCase.parser.
PARSER_STANDARD = "standard"
PARSER_BALANCE = "balance"
PARSER_DMD = "dmd"
PARSER_BLAS1 = "blas1"
PARSER_BLAS23 = "blas23"


@dataclass
class Counts:
    """Accumulated test counts for one or more test output files."""

    runs: int = 0
    numerical: int = 0
    illegal: int = 0
    info: int = 0

    @property
    def other(self) -> int:
        """Return the number of non-numerical errors (illegal + info).

        Returns:
            The combined number of "illegal value" and INFO errors.
        """
        return self.illegal + self.info

    @property
    def errors(self) -> int:
        """Return the total number of errors of any kind.

        Returns:
            The combined number of numerical and other errors.
        """
        return self.numerical + self.other

    def add(self, other: "Counts") -> None:
        """Accumulate another set of counts into this one.

        Args:
            other: The counts to add in place.
        """
        self.runs += other.runs
        self.numerical += other.numerical
        self.illegal += other.illegal
        self.info += other.info


@dataclass
class FileReport:
    """Parsing result for a single test output file."""

    counts: Counts = field(default_factory=Counts)
    notable_lines: "List[str]" = field(default_factory=list)


@dataclass
class SectionResult:
    """Accumulated counts of one library/API section of the summary."""

    # Per-precision rows of the summary table, in reporting order.
    precisions: "List[Tuple[str, Counts]]" = field(default_factory=list)
    total: Counts = field(default_factory=Counts)
    # Counts per output file, keyed by the API-independent output name.
    # Used to compare one API variant against another; a file that was
    # missing has no entry, so a partial run never compares equal.
    case_counts: "Dict[str, Counts]" = field(default_factory=dict)

    def error_map(self) -> "Dict[str, Tuple[int, int, int]]":
        """Return the per-file error counts, ignoring the run counts.

        Returns:
            The (numerical, illegal, info) triple of every analyzed
            output file, keyed by its API-independent name.
        """
        return {
            name: (counts.numerical, counts.illegal, counts.info)
            for name, counts in self.case_counts.items()
        }


@dataclass(frozen=True)
class TestCase:
    """One test driver invocation and its expected output file."""

    precision: str
    family: str
    description: str
    input_name: "Optional[str]"
    output_name: str
    executable: str
    parser: str
    library: str = LIBRARY_LAPACK
    # True when the API suffix also applies to the input file name.  The
    # BLAS Level 2/3 inputs name the output file on their first line, so
    # the _64 run needs the generated _64 input; every other driver takes
    # the same input for both APIs.
    input_suffixed: bool = False
    # False when the driver opens its own output file, so the harness must
    # not also redirect standard output onto it.
    redirect_stdout: bool = True

    def suffixed_output(self, suffix: str) -> str:
        """Return the output file name for an API suffix.

        Args:
            suffix: The API suffix, either ``""`` or ``"_64"``.

        Returns:
            The output file name, e.g. ``snep_64.out`` for suffix
            ``"_64"`` and base output name ``snep.out``.
        """
        stem = self.output_name[: -len(".out")]
        return "{}{}.out".format(stem, suffix)

    def suffixed_input(self, suffix: str) -> "Optional[str]":
        """Return the input file name for an API suffix.

        Args:
            suffix: The API suffix, either ``""`` or ``"_64"``.

        Returns:
            The input file name, or None for the drivers that read no
            input at all.
        """
        if self.input_name is None or not self.input_suffixed or not suffix:
            return self.input_name
        stem, _, extension = self.input_name.rpartition(".")
        return "{}{}.{}".format(stem, suffix, extension)

    def suffixed_executable(self, suffix: str) -> str:
        """Return the test driver name for an API suffix.

        Args:
            suffix: The API suffix, either ``""`` or ``"_64"``.

        Returns:
            The executable name, e.g. ``xeigtsts_64``.
        """
        return self.executable + suffix


@dataclass
class CaseOutcome:
    """Analysis outcome of one test case in one API variant.

    Collected in analysis order for the JUnit XML report: the parsing
    result of the output file (or None when the file was missing), the
    error message of a driver run that failed under ``--run`` or of an
    output file that could not be read, and the wall-clock duration of
    the driver run when ``--run`` was given.
    """

    case: TestCase
    suffix: str
    run_error: "Optional[str]" = None
    report: "Optional[FileReport]" = None
    duration: "Optional[float]" = None


def build_test_cases(letters: str, families: "Sequence[str]") -> "List[TestCase]":
    """Build the list of test cases for the selected precisions/families.

    Args:
        letters: Precision letters to include, in order (subset of
            ``"sdcz"``).
        families: Test families to include (subset of ``lin``, ``eig``,
            ``mixed``, ``rfp``, ``dmd``).

    Returns:
        The test cases in reporting order: for each precision, the
        eigenproblem sets, then the linear equation, mixed precision,
        RFP and DMD sets.
    """
    cases: "List[TestCase]" = []
    for letter in letters:
        if "eig" in families:
            for name, shared_input, description in EIG_STANDARD_SETS:
                cases.append(
                    TestCase(
                        precision=letter,
                        family="eig",
                        description=description,
                        input_name=(name if shared_input else letter + name) + ".in",
                        output_name=letter + name + ".out",
                        executable="xeigtst" + letter,
                        parser=PARSER_STANDARD,
                    )
                )
            for name, description in EIG_BALANCE_SETS:
                cases.append(
                    TestCase(
                        precision=letter,
                        family="eig",
                        description=description,
                        input_name=letter + name + ".in",
                        output_name=letter + name + ".out",
                        executable="xeigtst" + letter,
                        parser=PARSER_BALANCE,
                    )
                )
        for family, stem, executable_prefix, description in LIN_SETS:
            if family not in families:
                continue
            if family == "mixed":
                if letter not in MIXED_PARTNER:
                    continue
                letters_part = letter + MIXED_PARTNER[letter]
            else:
                letters_part = letter
            cases.append(
                TestCase(
                    precision=letter,
                    family=family,
                    description=description,
                    input_name=letters_part + stem + ".in",
                    output_name=letters_part + stem + ".out",
                    executable=executable_prefix + letters_part,
                    parser=PARSER_STANDARD,
                )
            )
        if "dmd" in families:
            cases.append(
                TestCase(
                    precision=letter,
                    family="dmd",
                    description="Dynamic Mode Decomposition",
                    input_name=letter + "dmd.in",
                    output_name=letter + "dmd.out",
                    executable="xdmdeigtst" + letter,
                    parser=PARSER_DMD,
                )
            )
        if "blas" in families:
            for level, description in BLAS_LEVELS:
                # Level 1 reads no input; Level 2/3 read e.g. sblat2.in,
                # whose first line names the output file, so the _64 run
                # needs the generated sblat2_64.in.
                cases.append(
                    TestCase(
                        precision=letter,
                        family="blas",
                        description=description,
                        input_name=(
                            None if level == 1 else "{}blat{}.in".format(letter, level)
                        ),
                        output_name="{}blat{}.out".format(letter, level),
                        executable="xblat{}{}".format(level, letter),
                        parser=PARSER_BLAS1 if level == 1 else PARSER_BLAS23,
                        library=LIBRARY_BLAS,
                        input_suffixed=level != 1,
                        redirect_stdout=level == 1,
                    )
                )
        if "cblas" in families:
            for level, description in BLAS_LEVELS:
                # The CBLAS inputs carry no output file name, so the same
                # input serves both APIs and the harness does the
                # redirection for every level.
                cases.append(
                    TestCase(
                        precision=letter,
                        family="cblas",
                        description="C interface to " + description,
                        input_name=(
                            None if level == 1 else "{}in{}".format(letter, level)
                        ),
                        output_name="{}test{}.out".format(letter, level),
                        executable="x{}cblat{}".format(letter, level),
                        parser=PARSER_BLAS1 if level == 1 else PARSER_BLAS23,
                        library=LIBRARY_CBLAS,
                    )
                )
    return cases


def parse_standard(lines: "Sequence[str]") -> FileReport:
    """Parse a test output file in the classic alasum/alasvm format.

    Counts runs from both the passing summary lines (``... tests run)``)
    and the failing summary lines (``N out of M tests failed ...``), so
    that failing test sets contribute to the run total as well.  The
    eigencondition checkers report failures as ``Error in <routine>``
    records instead; each such record counts as one numerical failure.

    Args:
        lines: The lines of the output file.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    report = FileReport()
    for line in lines:
        match = RE_TESTS_RUN.search(line)
        if match:
            report.counts.runs += int(match.group(1))
            continue
        match = RE_TESTS_FAILED.search(line)
        if match:
            report.counts.numerical += int(match.group(1))
            report.counts.runs += int(match.group(2))
            report.notable_lines.append(line)
            continue
        if RE_EC_ERROR.match(line):
            report.counts.numerical += 1
            report.notable_lines.append(line)
            continue
        if "illegal" in line or "Illegal" in line:
            report.counts.illegal += 1
            report.notable_lines.append(line)
            continue
        if " INFO" in line:
            report.counts.info += 1
            report.notable_lines.append(line)
    return report


def parse_balance(lines: "Sequence[str]") -> FileReport:
    """Parse a balancing/backtransformation test output file.

    These checkers (``schkbl.f`` and friends) do not use the alasum
    summary format.  Runs are taken from the ``total number of examples
    tested`` line, INFO errors from the ``number of examples where info
    is not 0`` line.  A non-finite ``value of largest test error`` or a
    nonzero ``example number where ILO or IHI wrong`` is counted as one
    numerical failure.

    Args:
        lines: The lines of the output file.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    report = FileReport()
    for line in lines:
        match = RE_EXAMPLES_TESTED.search(line)
        if match:
            report.counts.runs += int(match.group(1))
            continue
        match = RE_INFO_NOT_ZERO.search(line)
        if match:
            info_errors = int(match.group(1))
            report.counts.info += info_errors
            if info_errors > 0:
                report.notable_lines.append(line)
            continue
        match = RE_ILO_IHI_WRONG.search(line)
        if match:
            if int(match.group(1)) != 0:
                report.counts.numerical += 1
                report.notable_lines.append(line)
            continue
        match = RE_LARGEST_ERROR.search(line)
        if match:
            # Fortran prints double precision exponents as 0.1D+01.
            token = match.group(1).replace("D", "E").replace("d", "e")
            try:
                value = float(token)
            except ValueError:
                value = math.inf
            if not math.isfinite(value):
                report.counts.numerical += 1
                report.notable_lines.append(line)
    return report


def parse_dmd(lines: "Sequence[str]") -> FileReport:
    """Parse a dynamic mode decomposition test output file.

    Each per-test verdict line (``... test PASSED.`` or ``... test
    FAILED ...``) counts as one test run; each FAILED verdict counts as
    one numerical failure (the line itself reports how many individual
    cases failed).

    Args:
        lines: The lines of the output file.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    report = FileReport()
    for line in lines:
        match = RE_DMD_VERDICT.search(line)
        if match:
            report.counts.runs += 1
            if match.group(1).upper() == "FAILED":
                report.counts.numerical += 1
                report.notable_lines.append(line)
    return report


def parse_blas(lines: "Sequence[str]", level_one: bool) -> FileReport:
    """Parse a BLAS or CBLAS test output file.

    Test counts come from the ``... TESTS: n RUN, m FAILED`` lines the
    drivers report per routine: computational failures are numerical
    errors, error-exit failures are other errors.  Output from a build
    whose drivers do not report counts is still summarized, by falling
    back to one test per verdict.

    The drivers exit with status 0 even when they abandon the run, and
    print ``END OF TESTS`` even when routines failed, so breakage is
    detected from the text: an abandoned or misconfigured run, and a
    Level 2/3 file that never reached its footer, each count as one other
    error.

    Args:
        lines: The lines of the output file.
        level_one: True for the Level 1 drivers, which have no footer.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    report = FileReport()
    reported_counts = False
    saw_footer = False
    abandoned = False
    # Verdict tallies, used only if the driver reported no counts.
    verdicts = Counts()
    cases = 0
    case_failed = False

    for line in lines:
        match = RE_BLAS_COUNTS.match(line)
        if match:
            reported_counts = True
            report.counts.runs += int(match.group(2))
            failures = int(match.group(3))
            if match.group(1) == "COMPUTATIONAL":
                report.counts.numerical += failures
            else:
                report.counts.illegal += failures
            continue

        if RE_BLAS_END_OF_TESTS.match(line):
            saw_footer = True
            continue

        # Driver-level breakage, which no per-routine count can express.
        if (
            RE_BLAS_ABANDONED.search(line)
            or RE_BLAS_NOT_RECOGNIZED.match(line)
            or RE_BLAS_DOT_PRODUCTS.match(line)
            or RE_BLAS_INTERNAL.search(line)
            or RE_BLAS_INPUT_ERROR.match(line)
        ):
            abandoned = True
            report.counts.info += 1
            report.notable_lines.append(line)
            continue

        if RE_BLAS_FAILED_ERROR_EXIT.search(line):
            verdicts.runs += 1
            verdicts.illegal += 1
            report.notable_lines.append(line)
            continue
        if RE_BLAS_FAILED_COMPUTATIONAL.search(line):
            verdicts.runs += 1
            verdicts.numerical += 1
            report.notable_lines.append(line)
            continue
        if RE_BLAS_SUSPECT.search(line):
            verdicts.runs += 1
            verdicts.numerical += 1
            report.notable_lines.append(line)
            continue
        if RE_BLAS_PASSED.match(line):
            verdicts.runs += 1
            continue

        if RE_BLAS_DETAIL.search(line) or RE_BLAS_NOT_TESTED.match(line):
            report.notable_lines.append(line)
            continue

        if level_one:
            if RE_BLAS_L1_CASE.match(line):
                cases += 1
                case_failed = False
                continue
            # A NRM2 stress failure prints FAIL without clearing PASS, so
            # a case can report both; treat any FAIL as a failure.
            if RE_BLAS_L1_FAIL.match(line) and not case_failed:
                case_failed = True
                verdicts.numerical += 1
                report.notable_lines.append(line)

    if not reported_counts:
        if level_one:
            verdicts.runs += cases
        report.counts.add(verdicts)

    # A run that reported why it stopped has already been counted.
    if not level_one and not saw_footer and not abandoned:
        report.counts.info += 1
        report.notable_lines.append(
            "output ends without 'END OF TESTS': the driver did not finish\n"
        )
    return report


def parse_lines(parser: str, lines: "Sequence[str]") -> FileReport:
    """Parse test output lines with the parser kind of a test case.

    Args:
        parser: One of ``PARSER_STANDARD``, ``PARSER_BALANCE``,
            ``PARSER_DMD``, ``PARSER_BLAS1`` and ``PARSER_BLAS23``.
        lines: The lines of the output file.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    if parser == PARSER_BALANCE:
        return parse_balance(lines)
    if parser == PARSER_DMD:
        return parse_dmd(lines)
    if parser in (PARSER_BLAS1, PARSER_BLAS23):
        return parse_blas(lines, level_one=parser == PARSER_BLAS1)
    return parse_standard(lines)


def find_unrecognized_outputs(directories: "Dict[str, Path]") -> "List[str]":
    """Find ``.out`` files in the test directories this script cannot analyze.

    A file is unrecognized if its name matches no known test case of the
    library that owns its directory, in any precision, family or API
    variant — typically a test that was added to the harness without
    extending this script's test tables, or a renamed output such as
    those of ``make variants_testing``.  The current ``-p``/``-t``
    selection is deliberately ignored: a deselected file is not an
    unrecognized one.

    Args:
        directories: The existing testing directory of each library.

    Returns:
        The unrecognized file names, prefixed by their directory when
        more than one directory was scanned, sorted alphabetically.
    """
    all_cases = build_test_cases("sdcz", ALL_FAMILIES)
    unrecognized: "List[str]" = []
    for library, directory in directories.items():
        known = {
            case.suffixed_output(suffix)
            for case in all_cases
            if case.library == library
            for suffix in KNOWN_SUFFIXES
        }
        known.add(RESULTS_FILENAME)
        for path in directory.glob("*.out"):
            if path.name in known:
                continue
            unrecognized.append(
                path.name if len(directories) == 1 else str(directory / path.name)
            )
    return sorted(unrecognized)


def discover_suffixes(cases: "Sequence[TestCase]", directory: Path) -> "List[str]":
    """Detect which API variants have output files in a testing directory.

    Args:
        cases: The selected test cases of one library.
        directory: That library's testing directory.

    Returns:
        The suffixes (out of ``""`` and ``"_64"``) for which at least one
        expected output file exists; ``[""]`` if none exist at all.
    """
    suffixes = [
        suffix
        for suffix in KNOWN_SUFFIXES
        if any((directory / case.suffixed_output(suffix)).is_file() for case in cases)
    ]
    return suffixes or [""]


def find_executable(name: str, bin_dir: "Optional[str]") -> "Optional[Path]":
    """Locate a test driver executable.

    Args:
        name: The executable name without platform suffix, e.g.
            ``xlintsts``.
        bin_dir: The directory passed via ``--bin``, or None to probe the
            usual locations of CMake and Makefile builds relative to the
            current working directory.

    Returns:
        The absolute path of the executable, or None if it was not found.
    """
    if bin_dir is not None:
        directories = [Path(bin_dir)]
    else:
        directories = [
            Path("bin"),
            Path("bin") / "Release",
            Path("bin") / "Debug",
            Path("TESTING") / "LIN",
            Path("TESTING") / "EIG",
            Path("BLAS") / "TESTING",
            Path("CBLAS") / "testing",
        ]
    for directory in directories:
        for filename in (name, name + ".exe"):
            candidate = directory / filename
            if candidate.is_file():
                return candidate.resolve()
    return None


SOURCE_INPUT_DIRS: "Dict[str, str]" = {
    LIBRARY_LAPACK: "TESTING",
    LIBRARY_BLAS: "BLAS/TESTING",
    LIBRARY_CBLAS: "CBLAS/testing",
}


def run_test_case(
    case: TestCase, suffix: str, test_dir: Path, bin_dir: "Optional[str]"
) -> "Optional[str]":
    """Run one test driver, capturing its output in the ``.out`` file.

    Args:
        case: The test case to run.
        suffix: The API suffix, either ``""`` or ``"_64"``.
        test_dir: The directory containing the ``.in`` files; the driver
            runs there and the ``.out`` file is written there.
        bin_dir: The directory containing the test drivers, or None to
            probe the usual locations.

    Returns:
        An error message if the driver could not be run or exited with a
        nonzero status, otherwise None.
    """
    executable_name = case.suffixed_executable(suffix)
    executable = find_executable(executable_name, bin_dir)
    if executable is None:
        return "executable {} not found".format(executable_name)

    input_name = case.suffixed_input(suffix)
    input_path: "Optional[Path]" = None
    if input_name is not None:
        input_path = test_dir / input_name
        if not input_path.is_file():
            # CMake build trees hold only the .out files; the .in files
            # live in the source tree next to this script.  The _64 input
            # of a BLAS Level 2/3 driver is generated into the build tree
            # and has no source-tree counterpart.
            source_input = (
                Path(__file__).resolve().parent
                / SOURCE_INPUT_DIRS[case.library]
                / input_name
            )
            if source_input.is_file():
                input_path = source_input
            elif case.input_suffixed and suffix:
                # Falling back to the default-API input would make the
                # driver overwrite the default-API output file.
                return (
                    "input file {} not found (needed to keep the {} output "
                    "separate)".format(input_path, suffix)
                )
            else:
                return "input file {} not found".format(input_path)

    output_path = test_dir / case.suffixed_output(suffix)
    # Write to a temporary file first so that a driver that cannot even
    # start does not clobber the results of an earlier run.  Drivers that
    # open their own output file (the BLAS Level 2/3 testers, which take
    # its name from the first line of their input) write it directly.
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    try:
        with (
            open(str(input_path), "rb")
            if input_path is not None
            else open(os.devnull, "rb")
        ) as stdin, open(
            str(temporary_path) if case.redirect_stdout else os.devnull, "wb"
        ) as stdout:
            process = subprocess.run(
                [str(executable)],
                stdin=stdin,
                stdout=stdout,
                stderr=subprocess.STDOUT,
                cwd=str(test_dir),
            )
    except OSError as error:
        if case.redirect_stdout:
            try:
                temporary_path.unlink()
            except OSError:
                pass
        return "{} could not be run: {}".format(executable_name, error)
    if case.redirect_stdout:
        temporary_path.replace(output_path)
    if process.returncode != 0:
        return "{} exited with status {}".format(executable_name, process.returncode)
    return None


def read_output_file(path: Path) -> "Optional[List[str]]":
    """Read a test output file.

    Args:
        path: The path of the ``.out`` file.

    Returns:
        The lines of the file, or None if the file does not exist or
        cannot be read (the cause is then reported on standard error).
    """
    if not path.is_file():
        return None
    try:
        with open(str(path), encoding="utf-8", errors="replace") as handle:
            return handle.readlines()
    except OSError as error:
        print(
            "lapack_testing.py: cannot read {}: {}".format(path, error),
            file=sys.stderr,
        )
        return None


# Fixed-width summary table columns: label, run count and one block of
# "count (percent)" per error kind.  Values are right-aligned so that
# they line up under the ==== rule of their column.  The 73-column
# table is indented to sit centered under the 80-column headings.
SUMMARY_INDENT = "  "
SUMMARY_HEADER = SUMMARY_INDENT + "{:<18}   {:>13}   {:>18}   {:>18}".format(
    "SUMMARY", "nb test run", "numerical error", "other error"
)
SUMMARY_RULE = SUMMARY_INDENT + "   ".join(("=" * 18, "=" * 13, "=" * 18, "=" * 18))


def format_summary_row(label: str, counts: Counts) -> str:
    """Format one row of the summary table.

    Args:
        label: The row label, e.g. a precision name.
        counts: The counts to report in the row.

    Returns:
        The formatted table row without trailing newline.
    """
    if counts.runs > 0:
        numerical_percent = 100.0 * counts.numerical / counts.runs
        other_percent = 100.0 * counts.other / counts.runs
    else:
        numerical_percent = 0.0
        other_percent = 0.0
    numerical_percent_str = "({:.1f}%)".format(numerical_percent)
    other_percent_str = "({:.1f}%)".format(other_percent)
    return SUMMARY_INDENT + "{:<18}   {:>13}   {:>9} {:>8}   {:>9} {:>8}".format(
        label,
        counts.runs,
        counts.numerical,
        numerical_percent_str,
        counts.other,
        other_percent_str,
    )


def api_name(suffix: str) -> str:
    """Return the human-readable name of one API variant.

    Args:
        suffix: The API suffix, either ``""`` or ``"_64"``.

    Returns:
        The API name, e.g. ``"Extended API (_64)"``.
    """
    if suffix:
        return "Extended API ({})".format(suffix)
    return "Default API"


def section_title(library: str, suffixes: "Sequence[str]") -> str:
    """Return the human-readable name of a summary section.

    Args:
        library: The library name, e.g. ``"BLAS"``.
        suffixes: The API suffixes the section covers.  More than one
            means the section reports them together in a single table.

    Returns:
        The section name used in headings and messages, e.g.
        ``"BLAS: Extended API (_64)"`` for one API variant, or
        ``"BLAS: Default API and Extended API (_64)"`` for two.
    """
    names = [api_name(suffix) for suffix in suffixes]
    if len(names) == 1:
        return "{}: {}".format(library, names[0])
    return "{}: {} and {}".format(library, ", ".join(names[:-1]), names[-1])


def section_heading(title: str) -> str:
    """Return an 80-column dashed heading with a centered title.

    Args:
        title: The heading text.

    Returns:
        str: A line of exactly 80 characters — dashes running to both
        edges (one space at the beginning and end) with the title
        centered.
    """
    return " {} \n".format(" {} ".format(title).center(78, "-"))


class SummaryLog:
    """Collector for the detailed results file (``testing_results.txt``)."""

    def __init__(self, handle: "Optional[TextIO]") -> None:
        """Initialize the collector.

        Args:
            handle: The open results file, or None if it could not be
                opened (details are then discarded).
        """
        self._handle = handle

    def record(self, header: str, lines: "Sequence[str]") -> None:
        """Append one analyzed output file to the results file.

        Args:
            header: A short description of the file (its name).
            lines: The lines of the file.
        """
        if self._handle is None:
            return
        self._handle.write("==== {} ====\n".format(header))
        self._handle.writelines(lines)
        self._handle.flush()

    def close(self) -> None:
        """Close the results file if it was open."""
        if self._handle is not None:
            self._handle.close()


# Characters that must not appear in XML 1.0 text (the complement of its
# Char production).  Fortran test output can contain control characters,
# which would make consumers reject the whole report.
RE_XML_FORBIDDEN = re.compile("[^\t\n\r\x20-\ud7ff\ue000-\ufffd\U00010000-\U0010ffff]")

# Cap on the text of one JUnit failure/error element.  GitLab CI shows
# the text in the test details; an abandoned BLAS run can log megabytes.
JUNIT_TEXT_LIMIT = 16 * 1024


def sanitize_xml_text(text: str) -> str:
    """Replace characters that must not appear in XML 1.0 text.

    Args:
        text: The text to sanitize.

    Returns:
        The text with each forbidden character replaced by U+FFFD, the
        same replacement character used when decoding the output files.
    """
    return RE_XML_FORBIDDEN.sub("\ufffd", text)


def counts_message(counts: Counts) -> str:
    """Format the one-line counts summary of an analyzed output file.

    Args:
        counts: The counts of the file.

    Returns:
        A summary such as ``"17 numerical error(s), 2 other error(s)
        (illegal: 2, info: 0), 1298 test(s) run"``.
    """
    parts: "List[str]" = []
    if counts.numerical:
        parts.append("{} numerical error(s)".format(counts.numerical))
    if counts.other:
        parts.append(
            "{} other error(s) (illegal: {}, info: {})".format(
                counts.other, counts.illegal, counts.info
            )
        )
    parts.append("{} test(s) run".format(counts.runs))
    return ", ".join(parts)


def junit_testcase(outcome: CaseOutcome) -> "ET.Element":
    """Build the JUnit ``<testcase>`` element of one analyzed test case.

    The element carries at most one status child: an ``<error>`` when
    the driver could not be run under ``--run`` or the output file
    could not be read, a ``<skipped>`` when the output file was
    missing, a ``<failure>`` for numerical errors,
    an ``<error>`` for other (illegal value / INFO) errors, and none
    when everything passed.  GitLab CI displays the element text of the
    status child and ignores its ``message`` attribute, so the text
    always carries the full story: the message first, then the notable
    lines of the output file.

    Args:
        outcome: The analysis outcome of the test case.

    Returns:
        The ``<testcase>`` element.
    """
    case = outcome.case
    element = ET.Element(
        "testcase",
        {
            "classname": "{}{}.{}".format(case.library, outcome.suffix, case.family),
            "name": "{} ({} {})".format(
                case.suffixed_output(outcome.suffix),
                PRECISION_NAMES[case.precision],
                case.description,
            ),
        },
    )
    if case.input_name is not None:
        # The source-tree input file, as a repository-relative path.
        element.set(
            "file", "{}/{}".format(SOURCE_INPUT_DIRS[case.library], case.input_name)
        )
    report = outcome.report
    if report is not None:
        element.set("assertions", str(report.counts.runs))
    if outcome.duration is not None:
        element.set("time", "{:.3f}".format(outcome.duration))

    details: "List[str]" = []
    if outcome.run_error is not None:
        status = "error"
        message = outcome.run_error
        if report is not None:
            details.append(counts_message(report.counts))
    elif report is None:
        status = "skipped"
        message = "expected output file {} was missing".format(
            case.suffixed_output(outcome.suffix)
        )
    elif report.counts.errors > 0:
        status = "failure" if report.counts.numerical > 0 else "error"
        message = counts_message(report.counts)
    else:
        return element
    if report is not None:
        details.extend(line.rstrip("\n") for line in report.notable_lines)

    child = ET.SubElement(element, status)
    child.set("message", sanitize_xml_text(message))
    text = "\n".join([message] + details)
    if len(text) > JUNIT_TEXT_LIMIT:
        text = text[:JUNIT_TEXT_LIMIT] + "\n... [output truncated]"
    child.text = sanitize_xml_text(text)
    return element


def build_junit_tree(
    outcomes: "Sequence[CaseOutcome]", unrecognized: "Sequence[str]"
) -> "ET.ElementTree":
    """Build the JUnit XML document for the analyzed test cases.

    The document holds one ``<testsuite>`` per (library, API) section,
    in analysis order, with one ``<testcase>`` per output file.  Output
    files this script does not recognize are reported as one extra
    failing test case in a synthetic ``lapack_testing.py`` suite, so
    that the report does not look clean while ``--fail-on-unrecognized``
    fails the run.

    Args:
        outcomes: The analysis outcomes, in analysis order.
        unrecognized: The names of the unrecognized ``.out`` files.

    Returns:
        The document; its root is a ``<testsuites>`` element.
    """
    root = ET.Element("testsuites", {"name": "lapack_testing"})
    grouped: "Dict[Tuple[str, str], List[CaseOutcome]]" = {}
    for outcome in outcomes:
        grouped.setdefault((outcome.case.library, outcome.suffix), []).append(outcome)

    total_tests = 0
    total_failures = 0
    total_errors = 0
    total_skipped = 0
    for (library, suffix), suite_outcomes in grouped.items():
        suite = ET.SubElement(
            root, "testsuite", {"name": section_title(library, [suffix])}
        )
        failures = 0
        errors = 0
        skipped = 0
        suite_time = 0.0
        timed = False
        for outcome in suite_outcomes:
            element = junit_testcase(outcome)
            suite.append(element)
            if element.find("failure") is not None:
                failures += 1
            elif element.find("error") is not None:
                errors += 1
            elif element.find("skipped") is not None:
                skipped += 1
            if outcome.duration is not None:
                suite_time += outcome.duration
                timed = True
        suite.set("tests", str(len(suite_outcomes)))
        suite.set("failures", str(failures))
        suite.set("errors", str(errors))
        suite.set("skipped", str(skipped))
        if timed:
            suite.set("time", "{:.3f}".format(suite_time))
        total_tests += len(suite_outcomes)
        total_failures += failures
        total_errors += errors
        total_skipped += skipped

    if unrecognized:
        message = (
            "{} .out file(s) in the testing directories are not known to "
            "this script and were not analyzed".format(len(unrecognized))
        )
        suite = ET.SubElement(
            root,
            "testsuite",
            {
                "name": "lapack_testing.py",
                "tests": "1",
                "failures": "1",
                "errors": "0",
                "skipped": "0",
            },
        )
        testcase = ET.SubElement(
            suite,
            "testcase",
            {"classname": "lapack_testing", "name": "unrecognized .out files"},
        )
        failure = ET.SubElement(testcase, "failure")
        failure.set("message", sanitize_xml_text(message))
        failure.text = sanitize_xml_text("\n".join([message] + list(unrecognized)))
        total_tests += 1
        total_failures += 1

    root.set("tests", str(total_tests))
    root.set("failures", str(total_failures))
    root.set("errors", str(total_errors))
    root.set("skipped", str(total_skipped))
    return ET.ElementTree(root)


def write_junit_xml(
    path: Path, outcomes: "Sequence[CaseOutcome]", unrecognized: "Sequence[str]"
) -> "Optional[str]":
    """Write the JUnit XML report requested via ``--junit-xml``.

    The file is written atomically: first to a temporary file next to
    the target, which is renamed over it only when complete, so an
    aborted run does not leave a truncated report behind.

    Args:
        path: The target path of the report; missing parent directories
            are created.
        outcomes: The analysis outcomes, in analysis order.
        unrecognized: The names of the unrecognized ``.out`` files.

    Returns:
        An error message if the report could not be written, otherwise
        None.
    """
    # Path.with_name below would raise ValueError for such a path.
    if not path.name:
        return "cannot write {}: the path has no file name".format(path)
    tree = build_junit_tree(outcomes, unrecognized)
    # ET.indent is Python 3.9+; without it the report is one long line,
    # which every consumer accepts just the same.
    indent = getattr(ET, "indent", None)
    if indent is not None:
        indent(tree)
    temporary_path = path.with_name(path.name + ".tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tree.write(str(temporary_path), encoding="UTF-8", xml_declaration=True)
        temporary_path.replace(path)
    except OSError as error:
        try:
            temporary_path.unlink()
        except OSError:
            pass
        return "cannot write {}: {}".format(path, error)
    return None


def parse_args(argv: "Optional[Sequence[str]]" = None) -> argparse.Namespace:
    """Parse the command line arguments.

    Args:
        argv: The arguments to parse, or None to use ``sys.argv``.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyze the .out files produced by the LAPACK, BLAS "
        "and CBLAS test suites and print a summary of the test results.",
        epilog="By default all precisions and all test families are "
        "analyzed, each library is reported in its own section, and both "
        "the default API and extended API (_64) outputs are summarized "
        "when present.",
    )
    parser.add_argument(
        "-d",
        "--dir",
        default="TESTING",
        help="directory containing the LAPACK testing output (.out) files "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--blas-dir",
        default=str(Path("BLAS") / "TESTING"),
        help="directory containing the BLAS testing output (.out) files; "
        "skipped without warning if it does not exist (default: %(default)s)",
    )
    parser.add_argument(
        "--cblas-dir",
        default=str(Path("CBLAS") / "testing"),
        help="directory containing the CBLAS testing output (.out) files; "
        "skipped without warning if it does not exist (default: %(default)s)",
    )
    parser.add_argument(
        "-b",
        "--bin",
        default=None,
        help="directory containing the test drivers for --run; by default "
        "bin, bin/Release, bin/Debug, TESTING/LIN, TESTING/EIG, "
        "BLAS/TESTING and CBLAS/testing are probed",
    )
    parser.add_argument(
        "-r",
        "--run",
        action="store_true",
        help="run the LAPACK test drivers before analyzing their output "
        "(by default only existing .out files are analyzed)",
    )
    parser.add_argument(
        "-s",
        "--short",
        action="store_true",
        help="print only the summary table",
    )
    parser.add_argument(
        "-e",
        "--error",
        action="store_true",
        help="print only the error summary",
    )
    parser.add_argument(
        "-n",
        "--number",
        action="store_true",
        help="print only the numbers of failing tests (numerical failures "
        "and other errors, one per line)",
    )
    parser.add_argument(
        "-p",
        "--prec",
        choices=["s", "d", "sd", "c", "z", "cz", "x"],
        default="x",
        help="precisions to analyze: s=single, d=double, sd=single/double, "
        "c=complex, z=double complex, cz=complex/double complex, "
        "x=all (default)",
    )
    parser.add_argument(
        "-t",
        "--test",
        choices=list(ALL_FAMILIES) + ["lapack", "all"],
        default="all",
        help="test family to analyze: lin=linear equations, "
        "eig=eigenproblems (including balancing), mixed=mixed precision, "
        "rfp=RFP format, dmd=dynamic mode decomposition, blas=BLAS, "
        "cblas=CBLAS, lapack=all LAPACK families, all (default)",
    )
    parser.add_argument(
        "--suffix",
        action="append",
        choices=["none", "64"],
        default=None,
        help="API variant to analyze: 'none' for the default API, '64' for "
        "the index-64 extended API; may be given twice (default: analyze "
        "whichever variants have output files)",
    )
    parser.add_argument(
        "--merge-apis",
        action="store_true",
        help="when a library was analyzed for both the default and the "
        "extended API and both report the same errors, summarize them in a "
        "single section instead of one per API; affects only the summary "
        "table, not the detailed output or the exit status",
    )
    parser.add_argument(
        "--junit-xml",
        metavar="PATH",
        default=None,
        help="write a JUnit XML report of the analyzed output files to "
        "PATH (one testcase per output file), e.g. for GitLab CI test "
        "reports; written regardless of the display and --fail-* options",
    )
    parser.add_argument(
        "--fail-on-error",
        action="store_true",
        help="exit with a nonzero status if any test failure or error was "
        "found, a test driver could not be run, or unrecognized .out files "
        "were present in a testing directory",
    )
    parser.add_argument(
        "--fail-if-empty",
        action="store_true",
        help="exit with a nonzero status if no test results were analyzed",
    )
    parser.add_argument(
        "--fail-on-unrecognized",
        action="store_true",
        help="exit with a nonzero status if .out files not known to this "
        "script were present in a testing directory",
    )
    return parser.parse_args(argv)


def main(argv: "Optional[Sequence[str]]" = None) -> int:
    """Run the LAPACK, BLAS and CBLAS test summary tool.

    Args:
        argv: The command line arguments, or None to use ``sys.argv``.

    Returns:
        int: The process exit status. This is 2 for usage errors, 1 if a
        condition requested via ``--fail-on-error``, ``--fail-if-empty``
        or ``--fail-on-unrecognized`` occurred or a report requested via
        ``--junit-xml`` could not be written, and 0 otherwise.
    """
    args = parse_args(argv)
    short_summary: bool = args.short or args.number
    just_errors: bool = args.error

    test_dir = Path(args.dir)
    if not test_dir.is_dir():
        print(
            "lapack_testing.py: testing directory {} not found".format(test_dir),
            file=sys.stderr,
        )
        return 2

    # The BLAS tests are absent from builds that use an optimized BLAS,
    # and CBLAS is off by default, so a missing directory is normal and
    # is skipped silently.  An explicitly selected library that has no
    # directory is a usage error, though.
    directories: "Dict[str, Path]" = {LIBRARY_LAPACK: test_dir}
    for library, option in (
        (LIBRARY_BLAS, args.blas_dir),
        (LIBRARY_CBLAS, args.cblas_dir),
    ):
        directory = Path(option)
        if directory.is_dir():
            directories[library] = directory
        elif FAMILY_LIBRARY.get(args.test) == library:
            print(
                "lapack_testing.py: {} testing directory {} not found".format(
                    library, directory
                ),
                file=sys.stderr,
            )
            return 2

    letters = "sdcz" if args.prec == "x" else args.prec
    if args.test == "mixed":
        # The mixed-precision drivers exist only for d (ds) and z (zc);
        # like the old script, -t mixed analyzes both regardless of -p.
        if args.prec not in ("x", "dz"):
            print(
                "lapack_testing.py: -t mixed always analyzes the d and z "
                "mixed-precision tests; ignoring -p {}".format(args.prec),
                file=sys.stderr,
            )
        letters = "dz"
    if args.test == "all":
        families = list(ALL_FAMILIES)
    elif args.test == "lapack":
        families = list(LAPACK_FAMILIES)
    else:
        families = [args.test]
    # Drop the families of libraries that were not built.
    families = [family for family in families if FAMILY_LIBRARY[family] in directories]
    cases = build_test_cases(letters, families)
    if not cases:
        print(
            "lapack_testing.py: no test cases match -p {} -t {}".format(
                args.prec, args.test
            ),
            file=sys.stderr,
        )
        return 2

    forced_suffixes: "Optional[List[str]]" = None
    if args.suffix is not None:
        forced_suffixes = []
        for choice in args.suffix:
            suffix = "" if choice == "none" else "_64"
            if suffix not in forced_suffixes:
                forced_suffixes.append(suffix)
    elif args.run:
        forced_suffixes = [""]

    # One (library, API) pair per summary section, in reporting order.
    sections: "List[Tuple[str, str]]" = []
    for library in LIBRARIES:
        library_cases = [case for case in cases if case.library == library]
        if not library_cases:
            continue
        suffixes = forced_suffixes
        if suffixes is None:
            suffixes = discover_suffixes(library_cases, directories[library])
        for suffix in suffixes:
            sections.append((library, suffix))

    results_path = test_dir / RESULTS_FILENAME
    try:
        results_handle: "Optional[TextIO]" = open(
            str(results_path), "w", encoding="utf-8"
        )
    except OSError as error:
        print(
            "lapack_testing.py: cannot write {}: {}".format(results_path, error),
            file=sys.stderr,
        )
        results_handle = None
    log = SummaryLog(results_handle)

    if not short_summary:
        print(" ")
        print("-->  Testing LAPACK Routines  <--".center(80).rstrip())
        print(" ")
        print("Detailed results are stored in:".center(80).rstrip())
        print(str(results_path.resolve()).center(80).rstrip())

    summary = "\n" + "-->  LAPACK TESTING SUMMARY  <--".center(80).rstrip() + "\n"
    if not args.run:
        summary += "Processing LAPACK Testing output found in:".center(80).rstrip()
        summary += "\n" + str(test_dir.resolve()).center(80).rstrip() + "\n"

    grand_total = Counts()
    missing_files = 0
    run_failures = 0
    outcomes: "List[CaseOutcome]" = []

    # Counts are collected per (library, API) section first and rendered
    # afterwards, so that --merge-apis can decide to report two API
    # variants in one table.  The detailed output below stays per section.
    results: "Dict[Tuple[str, str], SectionResult]" = {}

    for library, suffix in sections:
        directory = directories[library]
        if len(sections) > 1 or suffix:
            if not short_summary:
                print(" ")
                print(section_heading(section_title(library, [suffix])))
        result = SectionResult()
        results[(library, suffix)] = result

        for letter, precision_name in PRECISIONS:
            precision_cases = [
                case
                for case in cases
                if case.precision == letter and case.library == library
            ]
            if not precision_cases:
                continue
            precision_total = Counts()

            for case in precision_cases:
                output_name = case.suffixed_output(suffix)
                run_error: "Optional[str]" = None
                duration: "Optional[float]" = None
                if not just_errors and not short_summary:
                    print(
                        "Testing {} '{}' ({})".format(
                            precision_name, case.description, output_name
                        ),
                        end=" ",
                    )
                if args.run:
                    start = time.monotonic()
                    run_error = run_test_case(case, suffix, directory, args.bin)
                    duration = time.monotonic() - start
                    if run_error is not None:
                        run_failures += 1
                        print(
                            "---- TESTING {}... FAILED({})!".format(
                                case.suffixed_executable(suffix), run_error
                            )
                        )
                lines = read_output_file(directory / output_name)
                if lines is None:
                    # A file that exists but cannot be read is a broken
                    # run, not a missing one; report it as an error.
                    if run_error is None and (directory / output_name).is_file():
                        run_error = (
                            "output file {} exists but could not be read".format(
                                output_name
                            )
                        )
                    outcomes.append(
                        CaseOutcome(case, suffix, run_error, None, duration)
                    )
                    missing_files += 1
                    if not short_summary:
                        print(
                            "---- WARNING: please check that you have the {} "
                            "output {}!".format(library, output_name)
                        )
                        print(
                            "---- WARNING: with the option -r, we can run the "
                            "testing for you"
                        )
                    continue
                log.record(
                    (
                        output_name
                        if library == LIBRARY_LAPACK
                        else str(directory / output_name)
                    ),
                    lines,
                )
                report = parse_lines(case.parser, lines)
                outcomes.append(CaseOutcome(case, suffix, run_error, report, duration))
                precision_total.add(report.counts)
                result.case_counts[case.output_name] = report.counts

                if not short_summary:
                    if not just_errors:
                        # Finish the "Testing ..." progress line.
                        if report.counts.runs > 0:
                            print("- passed: {}".format(report.counts.runs))
                        else:
                            print("")
                    for line in report.notable_lines:
                        print("-->  {}".format(line.strip()))
                    if report.counts.numerical > 0:
                        print(
                            "failing to pass the threshold: {}".format(
                                report.counts.numerical
                            )
                        )
                    if report.counts.illegal > 0:
                        print("Illegal Error: {}".format(report.counts.illegal))
                    if report.counts.info > 0:
                        print("Info Error: {}".format(report.counts.info))
                    if just_errors and report.counts.errors > 0:
                        print(
                            "ERROR IS LOCATED IN {} {} [ {} ]".format(
                                precision_name, case.description, output_name
                            )
                        )
                    if not just_errors or report.counts.errors > 0:
                        print("")
                sys.stdout.flush()

            result.precisions.append((precision_name, precision_total))
            result.total.add(precision_total)

        grand_total.add(result.total)

    log.close()

    # Group the sections for rendering, collapsing a library's API
    # variants into one table when --merge-apis is given and their errors
    # agree.  A missing output file leaves no entry to compare, so a
    # partial run never collapses.
    rendered: "List[Tuple[str, List[str], SectionResult]]" = []
    for library in LIBRARIES:
        library_suffixes = [
            suffix for section_library, suffix in sections if section_library == library
        ]
        if not library_suffixes:
            continue
        first = results[(library, library_suffixes[0])]
        if (
            args.merge_apis
            and len(library_suffixes) > 1
            and all(
                results[(library, suffix)].error_map() == first.error_map()
                for suffix in library_suffixes[1:]
            )
        ):
            rendered.append((library, library_suffixes, first))
            continue
        for suffix in library_suffixes:
            rendered.append((library, [suffix], results[(library, suffix)]))

    for library, shown_suffixes, result in rendered:
        if len(rendered) > 1 or shown_suffixes != [""]:
            summary += (
                "\n" + section_heading(section_title(library, shown_suffixes)) + "\n"
            )
        summary += SUMMARY_HEADER + "\n"
        summary += SUMMARY_RULE + "\n"
        for precision_name, precision_total in result.precisions:
            summary += format_summary_row(precision_name, precision_total) + "\n"
        if args.prec == "x":
            summary += (
                "\n" + format_summary_row("--> ALL PRECISIONS", result.total) + "\n"
            )
        if len(shown_suffixes) > 1 and any(
            results[(library, suffix)].case_counts[name].runs != counts.runs
            for suffix in shown_suffixes[1:]
            for name, counts in result.case_counts.items()
        ):
            # The errors agree, which is what the sections were collapsed
            # on, but the run counts do not; say which one is shown.
            summary += SUMMARY_INDENT + (
                "test counts differ between the APIs; those shown are the "
                "{}'s\n".format(api_name(shown_suffixes[0]))
            )

    if args.number:
        print(grand_total.numerical)
        print(grand_total.other)
    else:
        print(summary)
        if grand_total.runs == 0:
            print(
                "NO TESTS WERE ANALYZED, please use the -r option to run "
                "the LAPACK TESTING"
            )
    if missing_files > 0 and short_summary:
        print(
            "lapack_testing.py: {} expected output file(s) were missing "
            "(rerun without -s/-n for details)".format(missing_files),
            file=sys.stderr,
        )

    unrecognized = find_unrecognized_outputs(directories)
    if unrecognized:
        print(
            "lapack_testing.py: {} .out file(s) in {} are not known to this "
            "script and were NOT analyzed (new tests must be added to the "
            "test tables in this script):".format(
                len(unrecognized),
                ", ".join(str(directory) for directory in directories.values()),
            ),
            file=sys.stderr,
        )
        for name in unrecognized:
            print("    {}".format(name), file=sys.stderr)

    junit_error: "Optional[str]" = None
    if args.junit_xml is not None:
        junit_error = write_junit_xml(Path(args.junit_xml), outcomes, unrecognized)
        if junit_error is not None:
            print("lapack_testing.py: {}".format(junit_error), file=sys.stderr)

    if args.fail_if_empty and grand_total.runs == 0:
        return 1
    if args.fail_on_unrecognized and unrecognized:
        return 1
    if args.fail_on_error and (grand_total.errors > 0 or run_failures):
        return 1
    if junit_error is not None:
        return 1
    return 0


def _configure_output_streams() -> None:
    """Make stdout/stderr replace unencodable characters instead of dying.

    Fortran test output is not guaranteed to be encodable in the console
    encoding (notably on Windows with a redirected stdout).
    """
    for stream in (sys.stdout, sys.stderr):
        if isinstance(stream, io.TextIOWrapper):
            stream.reconfigure(errors="replace")


if __name__ == "__main__":
    _configure_output_streams()
    sys.exit(main())
