import logging
import argparse
import os
import pprint
import random
import sys
import time
import unittest

from pbcommand.cli import pacbio_args_runner, get_default_argparser

# the pbcommand version raise OSError for some reason
from pbcommand.cli.core import get_default_argparser_with_base_opts
from pbcommand.utils import (setup_console_and_file_logger, setup_logger,
                             Constants, compose)
from pbcommand.common_options import (add_log_debug_option,
                                      add_log_file_option,
                                      add_log_level_option)
from pbcommand.validators import validate_file

from pbsmrtpipe.engine import run_command_async
from pbsmrtpipe.cli import (LOG_LEVELS, resolve_dist_chunk_overrides)
from pbsmrtpipe.constants import SLOG_PREFIX

import pbsmrtpipe.testkit.butler as B
import pbsmrtpipe.testkit.loader as L
import pbsmrtpipe.testkit.xunit as X
import pbsmrtpipe.tools.utils as TU

from pbsmrtpipe.utils import setup_log, StdOutStatusLogFilter

log = logging.getLogger()
log.addHandler(logging.NullHandler())

slog = logging.getLogger(SLOG_PREFIX + __name__)
slog.addHandler(logging.NullHandler())


__version__ = '0.3.2'


def _patch_test_cases_with_job_dir(test_cases, job_dir):
    """This is dirrrttay. This must be called before the test cases are run"""
    for test_case in test_cases:
        test_case.__class__.job_dir = job_dir
        for t in test_case:
            t.__class__.job_dir = job_dir
            # log.debug("Setting job dir on {c}".format(c=t))


def _write_xunit_output(test_cases, result, output_xml, job_id):
    """Returns a XunitTestSuite instance"""
    xml = X.convert_suite_and_result_to_xunit(test_cases, result)

    log.debug("Writing Xunit XML output to {f}".format(f=output_xml))
    with open(output_xml, 'w+') as f:
        f.write(str(xml))

    xsuite = X.XunitTestSuite.from_xml(output_xml)
    # give the user some feedback that tests were run.
    # slog.info(str(xsuite))
    # slog.info(xsuite)
    nfailed_tests = xsuite.nfailure
    nerrors = xsuite.nerrors
    nskipped = xsuite.nskipped

    jenkins_xml = X.xunit_file_to_jenkins(output_xml, job_name=job_id)

    output_dir = os.path.dirname(output_xml)
    jenkins_name = "_".join(['jenkins', os.path.basename(output_xml)])
    jenkins_xml_file = os.path.join(output_dir, jenkins_name)

    log.info("Writing Jenkins XML output to {f}".format(f=jenkins_xml_file))
    with open(jenkins_xml_file, 'w') as f:
        f.write(str(jenkins_xml))
    
    log.info("Completed running {t} tests. {s} skipped, {e} errors, and {n} failed tests.".format(e=nerrors, n=nfailed_tests, s=nskipped, t=len(xsuite)))
    return xsuite


def run_butler_tests(test_cases, output_dir, output_xml, job_id):
    """

    :return:
    """

    # This is really hacky
    _patch_test_cases_with_job_dir(test_cases, output_dir)

    # This is the API to run directly from Unittest
    slog.info("Running test cases")
    log.debug(test_cases)
    result = unittest.TestResult()
    test_suite = unittest.TestSuite(test_cases)
    test_suite.run(result)
    #log.debug(result)

    xml = _write_xunit_output(test_cases, result, output_xml, job_id)
    log.info(str(xml))

    return 0 if result.wasSuccessful() else 1


def run_butler(butler, test_cases, output_xml,
               log_file,
               log_level=logging.DEBUG,
               force_distribute=None,
               force_chunk=None,
               ignore_test_failures=False):
    """
    Run a Butler instance.

    :param butler: Butler instance
    :return: exit code


    :rtype: int
    """
    started_at = time.time()

    if isinstance(force_distribute, bool):
        butler.force_distribute = force_distribute

    if isinstance(force_chunk, bool):
        butler.force_chunk = force_chunk

    cmd = butler.to_cmd()

    if not os.path.exists(butler.output_dir):
        os.mkdir(butler.output_dir)

    if isinstance(log_level, str):
        log_levels_d = {attr: getattr(logging, attr) for attr in LOG_LEVELS}
        console_level = log_levels_d[log_level]
    else:
        console_level = log_level

    # Always Write the log will full debug mode.
    setup_console_and_file_logger(console_level, Constants.LOG_FMT_SIMPLE, log_file, logging.DEBUG, Constants.LOG_FMT_STD)

    slog.info("completed setting up log to console and {f}".format(f=log_file))
    slog.info("Running butler with test id {i}".format(i=butler.job_id))
    slog.info("Running cmd '{c}'".format(c=cmd))

    def _to_p(file_name):
        return os.path.join(butler.output_dir, file_name)

    # pbsmrtpipe stdout and stderr
    stdout = _to_p('job.stdout')
    stderr = _to_p('job.stderr')

    with open(stdout, 'w') as stdout_fh:
        with open(stderr, 'w') as stderr_fh:

            rcode, stdout, stderr, run_time = run_command_async(cmd, stdout_fh, stderr_fh)

            rmessage = "was successful" if rcode == 0 else " failed"
            msg = "pbsmrtpipe command {m} ({s:.2f} sec) exit code {e}.'".format(e=rcode, s=run_time, m=rmessage)
            if rcode != 0:
                slog.error(msg)
                log.error(str(stderr))
            else:
                slog.info(msg)

        if test_cases is not None:
            slog.info("Running in test-only mode")
            trcode = run_butler_tests(test_cases, butler.output_dir, output_xml, butler.job_id)
        else:
            trcode = 0

    run_time = time.time() - started_at
    slog.info("Exiting testkit runner in {s:.2f} sec.".format(s=run_time))

    # Was butler successful + all the tests pass
    if ignore_test_failures:
        return rcode
    return rcode | trcode


def _args_run_butler(args):

    butler = B.config_parser_to_butler(args.testkit_cfg)

    test_cases = L.parse_cfg_file(args.testkit_cfg)

    if not os.path.exists(butler.output_dir):
        os.mkdir(butler.output_dir)

    output_xml = args.output_xml
    if output_xml is None:
        output_xml = os.path.join(butler.output_dir, 'testkit_xunit.xml')

    if args.log_file is None:
        log_file = os.path.join(butler.output_dir, 'testkit.log')
    else:
        log_file = args.log_file

    force_distribute, force_chunk = resolve_dist_chunk_overrides(args)

    log_level = args.log_level

    # Short hand for --log-level=DEBUG
    if args.debug:
        log_level = logging.DEBUG

    if log_level == logging.DEBUG:
        # The logger isn't setup yet
        print "Args", args

    if args.only_tests:
        # in test only mode, only emit to stdout (to avoid overwritten the
        # log file
        setup_logger(None, level=log_level)
        return run_butler_tests(test_cases, butler.output_dir, output_xml, butler.job_id)
    else:
        rcode = run_butler(butler, test_cases, output_xml, log_file,
                           log_level=log_level,
                           force_distribute=force_distribute,
                           force_chunk=force_chunk,
                           ignore_test_failures=args.ignore_test_failures)
        return rcode


def add_tests_only_option(p):
    p.add_argument('--only-tests', action='store_true', help="Only run the tests.")
    return p


def _add_config_file_option(p):
    p.add_argument('testkit_cfg', type=validate_file, help="Path to testkit.cfg file.")
    return p


def add_ignore_test_failures_option(p):
    p.add_argument("--ignore-test-failures", action="store_true",
                   help="Exit with code 0 if pbsmrtpipe ran successfully, even if " +
                        "some tests fail")
    return p


def add_output_xml_option(p):
    p.add_argument("--output-xml", action="store", dest="output_xml",
                   default=None, help="Path to output XUnit XML")
    return p


def get_parser():
    desc = "Testkit Tool to run pbsmrtpipe jobs."
    p = get_default_argparser_with_base_opts(__version__, desc)

    funcs = [TU.add_override_chunked_mode,
             TU.add_override_distribute_option,
             _add_config_file_option,
             add_tests_only_option,
             add_ignore_test_failures_option,
             add_output_xml_option]

    f = compose(*funcs)
    p = f(p)

    return p


def main(argv=sys.argv):
    parser = get_parser()
    # note, this logger will overwritten after the output dir is created
    return pacbio_args_runner(argv[1:], parser, _args_run_butler, log, setup_log_func=setup_log)
