#!/usr/bin/env python3
"""Turn CTest's valgrind memory check logs into a JUnit XML report.

``ctest -T memcheck`` runs every test under valgrind and writes one log
per test to ``Testing/Temporary/MemoryChecker.<index>.log``, where
``<index>`` is the position CTest gave the test in its test list.  Those
logs hold the defects but not the names of the tests they belong to, and
CTest's own JUnit report (``--output-junit``) records whether the test
program passed, which under a memory check is not the question being
asked: valgrind passes the exit status of the program through, so a test
riddled with invalid reads still counts as passed.

This script joins the two.  It reads the test list of the same build
(``ctest --show-only=json-v1``), maps every log back to the test it
belongs to, and reports that test as failed when valgrind's ``ERROR
SUMMARY`` line counts errors -- the same verdict the memory check job of
``.github/workflows/special.yml`` gates the build on.  The result is a
JUnit XML report of the memory check itself, one test case per test, for
Codecov's test analytics.

Reporting is all this script does: a test with defects makes for a
failing test case in the report, not for a non-zero exit status here.
It fails only when it cannot produce a report at all, so that an empty
report is never uploaded in place of a real one.

Example:
    cd build
    ctest --show-only=json-v1 > ctest_tests.json
    ctest -T memcheck > memcheck.out 2>&1
    ../.github/scripts/memcheck_junit.py \
        --test-list ctest_tests.json \
        --log-dir Testing/Temporary \
        --ctest-output memcheck.out \
        --output memcheck_junit.xml
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Dict, List, Optional, Sequence, Tuple


# Characters that XML 1.0 does not allow in text.  A valgrind log quotes
# the strings of the program it checked, so it need not be clean.
RE_XML_FORBIDDEN = re.compile("[^\t\n\r\x20-\ud7ff\ue000-\ufffd\U00010000-\U0010ffff]")

# The name CTest gives the valgrind log of the test with that index.
RE_LOG_NAME = re.compile(r"^MemoryChecker\.(\d+)\.log$")

# valgrind's verdict, written once per checked process, as in
# "ERROR SUMMARY: 7 errors from 3 contexts (suppressed: 0 from 0)".
RE_ERROR_SUMMARY = re.compile(
    r"ERROR SUMMARY:\s+(\d+)\s+errors?\s+from\s+(\d+)\s+contexts?"
)

# The per-test line of the ctest console output, for the run times, as in
# "  1/122 Test   #7: BLAS-xblat1d ...........   Passed   43.21 sec".
# Only the index and the time are read; the name comes from the test list.
RE_CTEST_RESULT = re.compile(
    r"^\s*\d+/\d+\s+Test\s+#\s*(\d+):.*?([0-9]+\.[0-9]+)\s+sec\s*$"
)

# Cap on the text of one JUnit failure element, as in lapack_testing.py:
# a defect on a hot path is reported for every call and logs megabytes.
JUNIT_TEXT_LIMIT = 16 * 1024


@dataclass
class LogReport:
    """What one ``MemoryChecker.<index>.log`` says about its test.

    Attributes:
        path: The log file this was read from.
        errors: The number of errors valgrind reported, or None when the
            log carries no ``ERROR SUMMARY`` at all, which means
            valgrind did not get as far as writing its verdict.
        contexts: The number of distinct contexts those errors came
            from; zero when errors is None.
        text: The contents of the log, for the failure element.
        read_error: The reason the log could not be read, if it could
            not be; the other fields are then empty.
    """

    path: Path
    errors: "Optional[int]"
    contexts: int
    text: str
    read_error: "Optional[str]"


def sanitize_xml_text(text: str) -> str:
    """Replace characters that must not appear in XML 1.0 text.

    Args:
        text: The text to sanitize.

    Returns:
        The text with each forbidden character replaced by U+FFFD, the
        same replacement character used when decoding the logs.
    """
    return RE_XML_FORBIDDEN.sub("\ufffd", text)


def read_test_names(path: Path) -> "List[str]":
    """Read the test names of a build from a CTest test list.

    CTest identifies the log of a test by the index of the test in this
    list, counted from one, which is how the two are matched up here.

    Args:
        path: The file holding the output of
            ``ctest --show-only=json-v1``.

    Returns:
        The test names, in test list order.

    Raises:
        OSError: The file could not be read.
        ValueError: The file is not JSON, or not a CTest test list.
    """
    with path.open(encoding="utf-8") as stream:
        try:
            document = json.load(stream)
        except ValueError as error:
            raise ValueError("{}: not JSON: {}".format(path, error))
    tests = document.get("tests") if isinstance(document, dict) else None
    if not isinstance(tests, list):
        raise ValueError("{}: not the output of ctest --show-only=json-v1".format(path))
    names: "List[str]" = []
    for index, test in enumerate(tests, start=1):
        name = test.get("name") if isinstance(test, dict) else None
        if not isinstance(name, str) or not name:
            raise ValueError("{}: test #{} has no name".format(path, index))
        names.append(name)
    return names


def read_log(path: Path) -> LogReport:
    """Read the valgrind verdict out of one memory checker log.

    A log holds one checked process, because the test command execs the
    test program under a memory check rather than starting it as a child
    of its own (see CMAKE/RuntestCommand.cmake), and every process
    valgrind checks would otherwise write to this one file.  The counts
    of all verdicts in the log are nevertheless summed rather than the
    last one taken, so that an error reported by any process fails the
    test even if that ever changes.

    Args:
        path: The log file to read.

    Returns:
        What the log says about its test.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        return LogReport(
            path=path, errors=None, contexts=0, text="", read_error=str(error)
        )
    errors: "Optional[int]" = None
    contexts = 0
    for match in RE_ERROR_SUMMARY.finditer(text):
        errors = (errors or 0) + int(match.group(1))
        contexts += int(match.group(2))
    return LogReport(
        path=path, errors=errors, contexts=contexts, text=text, read_error=None
    )


def collect_logs(
    directory: Path, test_count: int
) -> "Tuple[Dict[int, LogReport], List[Path]]":
    """Read every memory checker log of a directory.

    Args:
        directory: The directory CTest wrote the logs to, normally
            ``Testing/Temporary`` of the build tree.
        test_count: The number of tests in the test list; a log whose
            index is beyond it belongs to no test of this build.

    Returns:
        The logs by test index, and the paths of the logs that could not
        be matched to a test.
    """
    logs: "Dict[int, LogReport]" = {}
    unmatched: "List[Path]" = []
    for path in sorted(directory.glob("MemoryChecker.*.log")):
        match = RE_LOG_NAME.match(path.name)
        index = int(match.group(1)) if match is not None else 0
        if not 1 <= index <= test_count:
            unmatched.append(path)
            continue
        logs[index] = read_log(path)
    return logs, unmatched


def read_durations(path: Path) -> "Dict[int, float]":
    """Read the run times of the tests from the ctest console output.

    Args:
        path: The file the ctest output of the memory check was saved to.

    Returns:
        The run time of each test that reported one, by test index.

    Raises:
        OSError: The file could not be read.
    """
    durations: "Dict[int, float]" = {}
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = RE_CTEST_RESULT.match(line)
            if match is not None:
                durations[int(match.group(1))] = float(match.group(2))
    return durations


def junit_testcase(
    name: str, log: "Optional[LogReport]", duration: "Optional[float]"
) -> "ET.Element":
    """Build the JUnit ``<testcase>`` element of one checked test.

    The element carries at most one status child: a ``<skipped>`` when
    the test has no log, because it did not run under the memory check;
    an ``<error>`` when its log could not be read or holds no verdict,
    because then nothing is known about the test; a ``<failure>`` when
    valgrind reported errors; and none when it reported none.  The text
    of the status child carries the message first and then the log, the
    way lapack_testing.py reports the output of a failing test, for the
    consumers that show the text and ignore the ``message`` attribute.

    Args:
        name: The name CTest knows the test by.
        log: What the log of the test says, or None when it has none.
        duration: The run time of the test in seconds, if known.

    Returns:
        The ``<testcase>`` element.
    """
    element = ET.Element("testcase", {"classname": "memcheck", "name": name})
    if duration is not None:
        element.set("time", "{:.3f}".format(duration))

    if log is None:
        status = "skipped"
        message = "no valgrind log; the test did not run under the memory check"
        details = ""
    elif log.read_error is not None:
        status = "error"
        message = "cannot read {}: {}".format(log.path, log.read_error)
        details = ""
    elif log.errors is None:
        status = "error"
        message = (
            "{} holds no valgrind ERROR SUMMARY; valgrind did not report a "
            "verdict on this test".format(log.path)
        )
        details = log.text
    elif log.errors > 0:
        status = "failure"
        message = "{} error(s) from {} context(s) reported by valgrind".format(
            log.errors, log.contexts
        )
        details = log.text
    else:
        return element

    child = ET.SubElement(element, status)
    child.set("message", sanitize_xml_text(message))
    text = message if not details else "{}\n{}".format(message, details)
    if len(text) > JUNIT_TEXT_LIMIT:
        text = text[:JUNIT_TEXT_LIMIT] + "\n... [log truncated]"
    child.text = sanitize_xml_text(text)
    return element


def build_junit_tree(
    names: "Sequence[str]",
    logs: "Dict[int, LogReport]",
    durations: "Dict[int, float]",
    unmatched: "Sequence[Path]",
) -> "ET.ElementTree":
    """Build the JUnit XML document of the memory check.

    The document holds one ``<testsuite>`` with one ``<testcase>`` per
    test of the test list, in test list order.  Logs that belong to no
    test of this build are reported as one extra failing test case, so
    that defects they report are not quietly dropped from the report.

    The suite and the document carry the totals of their test cases:
    ``tests``, ``failures``, ``errors``, ``skipped`` and their summed run
    time (``time``).  A test whose time is unknown contributes nothing to
    the sum rather than a zero.

    Args:
        names: The test names, in test list order.
        logs: The logs of the tests, by test index.
        durations: The run times of the tests, by test index.
        unmatched: The paths of the logs that belong to no test.

    Returns:
        The document; its root is a ``<testsuites>`` element.
    """
    root = ET.Element("testsuites", {"name": "lapack_memcheck"})
    suite = ET.SubElement(root, "testsuite", {"name": "memcheck"})

    failures = 0
    errors = 0
    skipped = 0
    total_time = 0.0
    timed = False
    for index, name in enumerate(names, start=1):
        duration = durations.get(index)
        element = junit_testcase(name, logs.get(index), duration)
        suite.append(element)
        if element.find("failure") is not None:
            failures += 1
        elif element.find("error") is not None:
            errors += 1
        elif element.find("skipped") is not None:
            skipped += 1
        if duration is not None:
            total_time += duration
            timed = True

    tests = len(names)
    if unmatched:
        message = (
            "{} memory checker log(s) belong to no test of this build and "
            "were not reported above".format(len(unmatched))
        )
        element = ET.SubElement(
            suite,
            "testcase",
            {"classname": "memcheck", "name": "unmatched memory checker logs"},
        )
        failure = ET.SubElement(element, "failure")
        failure.set("message", sanitize_xml_text(message))
        failure.text = sanitize_xml_text(
            "\n".join([message] + [str(path) for path in unmatched])
        )
        tests += 1
        failures += 1

    for element in (suite, root):
        element.set("tests", str(tests))
        element.set("failures", str(failures))
        element.set("errors", str(errors))
        element.set("skipped", str(skipped))
        if timed:
            element.set("time", "{:.3f}".format(total_time))
    return ET.ElementTree(root)


def write_junit_xml(path: Path, tree: "ET.ElementTree") -> "Optional[str]":
    """Write the JUnit XML report to its file.

    The file is written atomically, as in lapack_testing.py: first to a
    temporary file next to the target, which is renamed over it only when
    complete, so that an aborted run does not leave a truncated report
    behind for the upload to pick up.

    Args:
        path: The target path of the report; missing parent directories
            are created.
        tree: The document to write.

    Returns:
        An error message if the report could not be written, otherwise
        None.
    """
    # Path.with_name below would raise ValueError for such a path.
    if not path.name:
        return "cannot write {}: the path has no file name".format(path)
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


def main(argv: "Optional[Sequence[str]]" = None) -> int:
    """Write the JUnit XML report of a CTest memory check.

    Args:
        argv: The command line arguments, or None for sys.argv.

    Returns:
        0 when the report was written, 1 when it could not be: the
        defects it reports are the result, not the exit status.
    """
    parser = argparse.ArgumentParser(
        description="Turn CTest's valgrind memory check logs into a JUnit "
        "XML report.",
    )
    parser.add_argument(
        "--test-list",
        required=True,
        type=Path,
        help="the test list of the build that was checked, as written by "
        "'ctest --show-only=json-v1'; the names of the tests are taken "
        "from it",
    )
    parser.add_argument(
        "--log-dir",
        default=Path("Testing/Temporary"),
        type=Path,
        help="the directory holding the MemoryChecker.<index>.log files "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--ctest-output",
        type=Path,
        help="the saved console output of the 'ctest -T memcheck' run, for "
        "the run times of the tests; without it the report carries none",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="the path to write the JUnit XML report to",
    )
    args = parser.parse_args(argv)

    try:
        names = read_test_names(args.test_list)
    except (OSError, ValueError) as error:
        print("memcheck_junit.py: {}".format(error), file=sys.stderr)
        return 1
    if not names:
        print(
            "memcheck_junit.py: {} lists no tests".format(args.test_list),
            file=sys.stderr,
        )
        return 1

    logs, unmatched = collect_logs(args.log_dir, len(names))
    if not logs:
        print(
            "memcheck_junit.py: no MemoryChecker.<index>.log files in {}; "
            "valgrind did not run".format(args.log_dir),
            file=sys.stderr,
        )
        return 1
    for path in unmatched:
        print(
            "memcheck_junit.py: {} belongs to no test of this build".format(path),
            file=sys.stderr,
        )

    durations: "Dict[int, float]" = {}
    if args.ctest_output is not None:
        try:
            durations = read_durations(args.ctest_output)
        except OSError as error:
            print("memcheck_junit.py: {}".format(error), file=sys.stderr)
            return 1
        if not durations:
            print(
                "memcheck_junit.py: {} holds no test result lines; the report "
                "carries no run times".format(args.ctest_output),
                file=sys.stderr,
            )

    tree = build_junit_tree(names, logs, durations, unmatched)
    write_error = write_junit_xml(args.output, tree)
    if write_error is not None:
        print("memcheck_junit.py: {}".format(write_error), file=sys.stderr)
        return 1

    root = tree.getroot()
    print(
        "Wrote {}: {} test(s), {} with valgrind errors, {} error(s), "
        "{} not run".format(
            args.output,
            root.get("tests"),
            root.get("failures"),
            root.get("errors"),
            root.get("skipped"),
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
