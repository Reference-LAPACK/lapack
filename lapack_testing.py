#!/usr/bin/env python3
"""Summarize (and optionally run) the LAPACK Fortran test suite.

This script analyzes the ``.out`` files written by the LAPACK testing
drivers (``xlintst*``, ``xeigtst*`` and ``xdmdeigtst*``) and prints a
summary table of the number of tests run and the number of failures per
precision (s/d/c/z).  With ``--run`` it executes the testing drivers
first and then analyzes their output.

When index-64 extended API outputs (``*_64.out``, produced by CMake
builds with ``BUILD_INDEX64_EXT_API=ON``) are present, they are analyzed
as well and reported in a separate "extended API" section so that the
default-API totals remain comparable across builds.

Examples:
    ./lapack_testing.py -n
        Print the numbers of failed tests by analyzing the LAPACK output.

    ./lapack_testing.py -n -r -p s
        Run the REAL precision tests, then print the numbers of failures.

    ./lapack_testing.py -n -p s -t eig
        Print the numbers of failures in REAL precision by analyzing only
        the eigenproblem test output.
"""

from __future__ import annotations

import argparse
import io
import math
import re
import subprocess
import sys
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

# All test families, in reporting order per precision.
ALL_FAMILIES: "Tuple[str, ...]" = ("eig",) + tuple(s[0] for s in LIN_SETS) + ("dmd",)

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

# Parser kinds, used by TestCase.parser.
PARSER_STANDARD = "standard"
PARSER_BALANCE = "balance"
PARSER_DMD = "dmd"


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


@dataclass(frozen=True)
class TestCase:
    """One LAPACK test driver invocation and its expected output file."""

    precision: str
    family: str
    description: str
    input_name: str
    output_name: str
    executable: str
    parser: str

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

    def suffixed_executable(self, suffix: str) -> str:
        """Return the test driver name for an API suffix.

        Args:
            suffix: The API suffix, either ``""`` or ``"_64"``.

        Returns:
            The executable name, e.g. ``xeigtsts_64``.
        """
        return self.executable + suffix


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


def parse_lines(parser: str, lines: "Sequence[str]") -> FileReport:
    """Parse test output lines with the parser kind of a test case.

    Args:
        parser: One of ``PARSER_STANDARD``, ``PARSER_BALANCE`` and
            ``PARSER_DMD``.
        lines: The lines of the output file.

    Returns:
        The counts and the notable (error) lines of the file.
    """
    if parser == PARSER_BALANCE:
        return parse_balance(lines)
    if parser == PARSER_DMD:
        return parse_dmd(lines)
    return parse_standard(lines)


def find_unrecognized_outputs(test_dir: Path) -> "List[str]":
    """Find ``.out`` files in the test directory this script cannot analyze.

    A file is unrecognized if its name matches no known test case in any
    precision, family or API variant — typically a test that was added to
    the harness without extending this script's test tables, or a renamed
    output such as those of ``make variants_testing``.  The current
    ``-p``/``-t`` selection is deliberately ignored: a deselected file is
    not an unrecognized one.

    Args:
        test_dir: The directory containing the ``.out`` files.

    Returns:
        The unrecognized file names, sorted alphabetically.
    """
    all_cases = build_test_cases("sdcz", ALL_FAMILIES)
    known = {
        case.suffixed_output(suffix) for case in all_cases for suffix in KNOWN_SUFFIXES
    }
    return sorted(
        path.name for path in test_dir.glob("*.out") if path.name not in known
    )


def discover_suffixes(test_dir: Path, cases: "Sequence[TestCase]") -> "List[str]":
    """Detect which API variants have output files in the test directory.

    Args:
        test_dir: The directory containing the ``.out`` files.
        cases: The selected test cases.

    Returns:
        The suffixes (out of ``""`` and ``"_64"``) for which at least one
        expected output file exists; ``[""]`` if none exist at all.
    """
    suffixes = [
        suffix
        for suffix in KNOWN_SUFFIXES
        if any((test_dir / case.suffixed_output(suffix)).is_file() for case in cases)
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
        ]
    for directory in directories:
        for filename in (name, name + ".exe"):
            candidate = directory / filename
            if candidate.is_file():
                return candidate.resolve()
    return None


def run_test_case(
    case: TestCase, suffix: str, test_dir: Path, bin_dir: "Optional[str]"
) -> "Optional[str]":
    """Run one test driver, redirecting its output to the ``.out`` file.

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
    input_path = test_dir / case.input_name
    if not input_path.is_file():
        # CMake build trees hold only the .out files; the .in files
        # live next to this script in the source tree.
        source_input = Path(__file__).resolve().parent / "TESTING" / case.input_name
        if source_input.is_file():
            input_path = source_input
        else:
            return "input file {} not found".format(input_path)
    output_path = test_dir / case.suffixed_output(suffix)
    # Write to a temporary file first so that a driver that cannot even
    # start does not clobber the results of an earlier run.
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    try:
        with open(str(input_path), "rb") as stdin, open(
            str(temporary_path), "wb"
        ) as stdout:
            process = subprocess.run(
                [str(executable)],
                stdin=stdin,
                stdout=stdout,
                stderr=subprocess.STDOUT,
                cwd=str(test_dir),
            )
    except OSError as error:
        try:
            temporary_path.unlink()
        except OSError:
            pass
        return "{} could not be run: {}".format(executable_name, error)
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


def section_title(suffix: str) -> str:
    """Return the human-readable name of an API section.

    Args:
        suffix: The API suffix, either ``""`` or ``"_64"``.

    Returns:
        The section name used in headings and messages.
    """
    if suffix:
        return "Extended API ({})".format(suffix)
    return "Default API"


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


def parse_args(argv: "Optional[Sequence[str]]" = None) -> argparse.Namespace:
    """Parse the command line arguments.

    Args:
        argv: The arguments to parse, or None to use ``sys.argv``.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyze the .out files produced by the LAPACK test "
        "suite and print a summary of the test results.",
        epilog="By default all precisions and all test families are "
        "analyzed, and both the default API and extended API (_64) "
        "outputs are summarized when present.",
    )
    parser.add_argument(
        "-d",
        "--dir",
        default="TESTING",
        help="directory containing the LAPACK testing output (.out) files "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "-b",
        "--bin",
        default=None,
        help="directory containing the LAPACK test drivers for --run; by "
        "default bin, bin/Release, bin/Debug, TESTING/LIN and TESTING/EIG "
        "are probed",
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
        choices=list(ALL_FAMILIES) + ["all"],
        default="all",
        help="test family to analyze: lin=linear equations, "
        "eig=eigenproblems (including balancing), mixed=mixed precision, "
        "rfp=RFP format, dmd=dynamic mode decomposition, all (default)",
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
        "--fail-on-error",
        action="store_true",
        help="exit with a nonzero status if any test failure or error was "
        "found, a test driver could not be run, or unrecognized .out files "
        "were present in the testing directory",
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
        "script were present in the testing directory",
    )
    return parser.parse_args(argv)


def main(argv: "Optional[Sequence[str]]" = None) -> int:
    """Run the LAPACK test summary tool.

    Args:
        argv: The command line arguments, or None to use ``sys.argv``.

    Returns:
        int: The process exit status. This is 2 for usage errors, 1 if a
        condition requested via ``--fail-on-error``, ``--fail-if-empty``
        or ``--fail-on-unrecognized`` occurred, and 0 otherwise.
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
    families = list(ALL_FAMILIES) if args.test == "all" else [args.test]
    cases = build_test_cases(letters, families)
    if not cases:
        print(
            "lapack_testing.py: no test cases match -p {} -t {}".format(
                args.prec, args.test
            ),
            file=sys.stderr,
        )
        return 2

    if args.suffix is not None:
        suffixes: "List[str]" = []
        for choice in args.suffix:
            suffix = "" if choice == "none" else "_64"
            if suffix not in suffixes:
                suffixes.append(suffix)
    elif args.run:
        suffixes = [""]
    else:
        suffixes = discover_suffixes(test_dir, cases)

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

    for suffix in suffixes:
        if len(suffixes) > 1 or suffix:
            summary += "\n" + section_heading(section_title(suffix)) + "\n"
            if not short_summary:
                print(" ")
                print(section_heading(section_title(suffix)))
        summary += SUMMARY_HEADER + "\n"
        summary += SUMMARY_RULE + "\n"
        section_total = Counts()

        for letter, precision_name in PRECISIONS:
            precision_cases = [case for case in cases if case.precision == letter]
            if not precision_cases:
                continue
            precision_total = Counts()

            for case in precision_cases:
                output_name = case.suffixed_output(suffix)
                if not just_errors and not short_summary:
                    print(
                        "Testing {} '{}' ({})".format(
                            precision_name, case.description, output_name
                        ),
                        end=" ",
                    )
                if args.run:
                    error_message = run_test_case(case, suffix, test_dir, args.bin)
                    if error_message is not None:
                        run_failures += 1
                        print(
                            "---- TESTING {}... FAILED({})!".format(
                                case.suffixed_executable(suffix), error_message
                            )
                        )
                lines = read_output_file(test_dir / output_name)
                if lines is None:
                    missing_files += 1
                    if not short_summary:
                        print(
                            "---- WARNING: please check that you have the LAPACK "
                            "output {}!".format(output_name)
                        )
                        print(
                            "---- WARNING: with the option -r, we can run the "
                            "LAPACK testing for you"
                        )
                    continue
                log.record(output_name, lines)
                report = parse_lines(case.parser, lines)
                precision_total.add(report.counts)

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

            summary += format_summary_row(precision_name, precision_total) + "\n"
            section_total.add(precision_total)

        if args.prec == "x":
            summary += (
                "\n" + format_summary_row("--> ALL PRECISIONS", section_total) + "\n"
            )
        grand_total.add(section_total)

    log.close()

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

    unrecognized = find_unrecognized_outputs(test_dir)
    if unrecognized:
        print(
            "lapack_testing.py: {} .out file(s) in {} are not known to this "
            "script and were NOT analyzed (new tests must be added to the "
            "test tables in this script):".format(len(unrecognized), test_dir),
            file=sys.stderr,
        )
        for name in unrecognized:
            print("    {}".format(name), file=sys.stderr)

    if args.fail_if_empty and grand_total.runs == 0:
        return 1
    if args.fail_on_unrecognized and unrecognized:
        return 1
    if args.fail_on_error and (grand_total.errors > 0 or run_failures):
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
