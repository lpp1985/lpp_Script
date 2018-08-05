
"""
Parallel wrapper for pbtestkit-service-runner, similar to pbtestkit-multirunner
"""

import multiprocessing
import functools
import operator
import logging
import time
import sys
import os

from pbcommand.cli import (get_default_argparser_with_base_opts,
                           pacbio_args_runner)
from pbcommand.utils import setup_log

from pbsmrtpipe.testkit.multirunner import (validate_testkit_cfg_fofn,
    testkit_cfg_fofn_to_files)
from pbsmrtpipe.engine import backticks

__version__ = '0.1.0'
log = logging.getLogger(__name__)

class Constants(object):
    PORT = int(os.environ.get("PB_SERVICE_PORT", "8081"))
    HOST = os.environ.get("PB_SERVICE_HOST", "http://localhost")
    EXE = "pbtestkit-service-runner"
    NPROC = 1
    MISC_OPTS = ""
    SLEEP_TIME = 4 # seconds


# FIXME(nechols)(2016-01-22): this should use API calls, but I'd like to
# figure out how to deal sensibly with log output first
def _run_testkit_cfg(testkit_cfg, debug=False, misc_opts=Constants.MISC_OPTS,
                     sleep_time=0):
    time.sleep(sleep_time)
    os.chdir(os.path.dirname(testkit_cfg))
    cmd = "{e} --debug {m} {c}".format(c=testkit_cfg, e=Constants.EXE,
        m=misc_opts)
    rcode, stdout, stderr, run_time = backticks(cmd)
    if debug:
        log.debug(" ".join([str(i) for i in [cmd, rcode, stdout, stderr]]))
    # Returning the butler cfg is a bit odd, but this is necessary for the
    # parallelism to work consistently with the serial version
    return testkit_cfg, rcode, stdout, stderr, run_time


def run_services_testkit_jobs(host, port, testkit_cfg_fofn, nworkers=1,
                              ignore_test_failures=False, time_out=1800,
                              sleep_time=2, import_only=False):
    testkit_cfgs = testkit_cfg_fofn_to_files(testkit_cfg_fofn)
    nworkers = min(len(testkit_cfgs), nworkers)
    results = []
    started_at = time.time()
    misc_opts = [
        "--host", host,
        "--port", str(port),
        "--timeout", str(time_out),
        "--sleep", str(sleep_time)
    ]
    if ignore_test_failures:
        misc_opts.append("--ignore-test-failures")
    if import_only:
        misc_opts.append("--import-only")
    misc_opts = " ".join(misc_opts)
    if nworkers == 1:
        log.info("Running in serial mode.")
        for testkit_cfg in testkit_cfgs:
            bcfg, rcode, stdout, stderr, job_run_time = _run_testkit_cfg(
                testkit_cfg, misc_opts=misc_opts)
            d = dict(x=testkit_cfg, r=rcode, s=int(job_run_time),
                     m=job_run_time / 60.0)
            log.info("Completed running {x}. exit code {r} in {s} sec ({m:.2f} min).".format(**d))
            results.append((testkit_cfg, rcode, stdout, stderr, job_run_time))
    else:
        _results = []
        # XXX to avoid connection errors when submitting a job to services,
        # the calls to pbtestkit-service-runner run staggered with a sleep
        # time of Constants.SLEEP_TIME between each job.  This allows us to
        # start many more near-simultaneous pbsmrtpipe jobs than would
        # otherwise be the case.
        pool = multiprocessing.Pool(nworkers)
        for i_cfg, testkit_cfg in enumerate(testkit_cfgs):
            sleep_time = (i_cfg % nworkers) * Constants.SLEEP_TIME
            __run_testkit_cfg = functools.partial(_run_testkit_cfg,
                misc_opts=misc_opts, sleep_time=sleep_time)
            log.debug("Running {c} with sleep time {t}".format(
                      c=testkit_cfg, t=sleep_time))
            _results.append(pool.apply_async(__run_testkit_cfg, (testkit_cfg,)))
        pool.close()
        pool.join()
        results = [r.get() for r in _results]
            #results = pool.map_async(__run_testkit_cfg, testkit_cfgs).get(9999999)
    log.info("Results:")
    for bcfg, rcode, _, _, job_run_time in results:
        d = dict(r=rcode, j=bcfg, s=int(job_run_time), m=job_run_time / 60.0)
        log.info("exit code {r} in {s} sec ({m:.2f} min). job {j}".format(**d))
    njobs = len(results)
    rcodes = [operator.getitem(r, 1) for r in results]
    nfailed = len([r for r in rcodes if r != 0])
    run_time = time.time() - started_at
    d = dict(n=njobs, x=nfailed, s=int(run_time), m=run_time / 60.0)
    msg = "Completed {n} jobs in {s} sec ({m:.2f} min) {x} failed.".format(**d)
    log.info(msg)
    # should this propagate the rcodes from siv_butler calls?
    return 0 if nfailed == 0 else -1


def args_runner(args):
    return run_services_testkit_jobs(
        host=args.host,
        port=args.port,
        testkit_cfg_fofn=args.testkit_cfg_fofn,
        nworkers=args.nworkers,
        ignore_test_failures=args.ignore_test_failures,
        time_out=args.time_out,
        sleep_time=args.sleep,
        import_only=args.import_only)


def get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__)
    p.add_argument("testkit_cfg_fofn", type=validate_testkit_cfg_fofn,
                  help="Text file listing testkit.cfg files to run")
    p.add_argument("-u", "--host", dest="host", action="store",
                   default=Constants.HOST)
    p.add_argument("-p", "--port", dest="port", action="store",
                   default=Constants.PORT, help="Port number")
    p.add_argument("-n", "--nworkers", type=int, default=Constants.NPROC,
                   help="Number of jobs to concurrently run.")
    p.add_argument("-t", "--timeout", dest="time_out", type=int, default=1800,
                   help="Timeout for blocking after job submission")
    p.add_argument("-s", "--sleep", dest="sleep", type=int, default=2,
                   help="Sleep time after job submission")
    p.add_argument("--ignore-test-failures", dest="ignore_test_failures",
                   action="store_true",
                   help="Only exit with non-zero return code if the job "+
                        "itself failed, regardless of test outcome")
    p.add_argument("--import-only", dest="import_only", action="store_true",
                   help="Import datasets without running pipelines")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
