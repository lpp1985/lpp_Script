
from __future__ import division, print_function
from collections import defaultdict
from xml.dom import minidom
import unittest
import sys


def iter_non_redundant_errors(errors):
    """
    Iterate over errors, skipping any that are judged to be duplicates
    based on their hashes.  Yields a tuple comprised of a 'unique' error
    and the number of occurences.
    """
    counts = defaultdict(int)
    for error in errors:
        counts[error] += 1
    displayed = set([])
    for error in errors:
        if (not error in displayed):
            displayed.add(error)
            yield error, counts[error]


def show_validation_errors(errors, out=sys.stdout, verbose=False,
                           use_termcolor=False):
    """
    Display a human-readable summary.  By default this will only print one
    error of each category; this can be overridden using the verbose flag.
    """
    try:
        import termcolor
    except ImportError as e:
        termcolor = None

    def _cprint(msg, c="yellow"):
        if not use_termcolor or termcolor is None:
            print(msg, file=out)
        else:
            termcolor.cprint(msg, c, attrs=['bold'], file=out)
    if (len(errors) == 0):
        if verbose:
            _cprint("(no spec violations detected)", c="green")
    else:
        _cprint("%d spec violation(s) detected:" % len(errors), c="red")
    if verbose:
        for error in errors:
            _cprint(str(error))
    else:
        for error, n_error in iter_non_redundant_errors(errors):
            _cprint("  " + str(error))
            if (n_error > 1):
                print("    [%d similar error(s) not displayed]" % (n_error - 1),
                      file=out)


def generate_multiple_file_junit_report(results, xml_out, skipped_files=()):
    """
    Produce JUnit-compatible output for multiple pbvalidate results (used in
    nightly builds).
    """
    doc = minidom.Document()
    n_failures = len([r for r in results if r.return_code != 0])
    root = doc.createElement("testsuite")
    doc.appendChild(root)
    root.setAttribute("errors", "0")
    root.setAttribute("failures", str(n_failures))
    root.setAttribute("name", "pbvalidate")
    root.setAttribute("skips", "0")
    root.setAttribute("tests", "2")
    root.setAttribute("time", str(sum([r.time for r in results])))
    for result in results:
        testcase = doc.createElement("testcase")
        testcase.setAttribute("classname", "pbvalidate.main.run_validator")
        testcase.setAttribute("name", result.file_name)
        testcase.setAttribute("time", str(result.time))
        if result.return_code != 0:
            failure = doc.createElement("failure")
            failure.setAttribute("message",
                                 "%d validation errors" % result.n_errors)
            failure.appendChild(doc.createCDATASection(result.error_string))
            testcase.appendChild(failure)
        root.appendChild(testcase)
    for file_name in skipped_files:
        testcase = doc.createElement("testcase")
        testcase.setAttribute("classname", "pbvalidate.main.run_validator")
        testcase.setAttribute("name", file_name)
        testcase.setAttribute("time", "0")
        testcase.appendChild(doc.createElement("skipped"))
        root.appendChild(testcase)
    xml_out.write(doc.toprettyxml(indent="  "))
