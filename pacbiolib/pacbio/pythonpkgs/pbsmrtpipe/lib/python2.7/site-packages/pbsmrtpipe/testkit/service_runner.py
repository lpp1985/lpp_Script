
"""
Utility for running a testkit job through services as an alternative to
pbtestkit-runner.
"""

import xml.etree.ElementTree as ET
import argparse
import logging
import os
import sys

from pbcommand.cli import (get_default_argparser_with_base_opts,
                           pacbio_args_runner)
from pbcommand.utils import setup_log
from pbcommand.services import ServiceAccessLayer, ServiceEntryPoint
from pbcommand.services.cli import run_analysis_job
import pbcommand.services

from pbsmrtpipe.pb_io import parse_pipeline_preset_xml, parse_pipeline_preset_json
from pbsmrtpipe.testkit.butler import config_parser_to_butler
from pbsmrtpipe.testkit.loader import (parse_cfg_file,
    dtype_and_uuid_from_dataset_xml)
from pbsmrtpipe.testkit.runner import run_butler_tests
from pbsmrtpipe.constants import to_opt_type_ns

log = logging.getLogger(__name__)


# FIXME(nechols)(2016-01-22): these utility functions should live elsewhere
def get_entrypoints(testkit_cfg):
    parsed_cfg = config_parser_to_butler(testkit_cfg)
    entrypoints = parsed_cfg.entry_points
    return entrypoints


def get_task_and_workflow_options(testkit_cfg):
    parsed_cfg = config_parser_to_butler(testkit_cfg)
    workflow_options, task_options = [], []
    def __get_option_type(val):
        option_type = to_opt_type_ns("string")
        if isinstance(val, bool):
            option_type = to_opt_type_ns("boolean")
        elif isinstance(val, int):
            option_type = to_opt_type_ns("integer")
        elif isinstance(val, float):
            option_type = to_opt_type_ns("float")
        elif val is None:
            val = ""
        return option_type, val
    if not parsed_cfg.preset_xml in [None, '']:
        if not parsed_cfg.preset_json in [None, '']:
            raise ValueError("Please use either preset_json or preset_xml, not both")
        presets = parse_pipeline_preset_xml(parsed_cfg.preset_xml)
        for option_id, option_value in presets.task_options:
            log.info("task_option: {i} = {v}".format(i=option_id,
                                                     v=option_value))
            option_type, option_value = __get_option_type(option_value)
            task_options.append(dict(
                optionId=option_id,
                value=option_value,
                optionTypeId=option_type))
        for option_id, option_value in presets.workflow_options:
            log.info("workflow_option: {i} = {v}".format(i=option_id,
                                                         v=option_value))
            workflow_options.append(dict(
                optionId=option_id,
                value=option_value,
                optionTypeId=__get_option_type(option_value)[0]))
    elif not parsed_cfg.preset_json in [None, '']:
        presets = parse_pipeline_preset_json(parsed_cfg.preset_json)
        for option_id, option_value in presets.task_options:
            log.info("task_option: {i} = {v}".format(i=option_id,
                                                     v=option_value))
            option_type, option_value = __get_option_type(option_value)
            task_options.append(dict(
                optionId=option_id,
                value=option_value,
                optionTypeId=option_type))
        for option_id, option_value in presets.workflow_options:
            log.info("workflow_option: {i} = {v}".format(i=option_id,
                                                         v=option_value))
            workflow_options.append(dict(
                optionId=option_id,
                value=option_value,
                optionTypeId=__get_option_type(option_value)[0]))
    return task_options, workflow_options


def entrypoints_dicts(entrypoints):
    """
    Extract dataset info from a list of entrypoints.
    """
    eps = []
    for entrypoint, dataset_xml in entrypoints.iteritems():
        dtype, unique_id = dtype_and_uuid_from_dataset_xml(dataset_xml)
        entry = {"_comment": "pbservice auto-job",
                 "datasetId": "{u}".format(u=unique_id),
                 "entryId": "{k}".format(k=entrypoint),
                 "fileTypeId": "{t}".format(t=dtype)}
        eps.append(entry)
    return eps


def pipeline_id_from_testkit_cfg(testkit_cfg):
    parsed_cfg = config_parser_to_butler(testkit_cfg)
    if parsed_cfg.workflow_xml is not None:
        tree = ET.parse(parsed_cfg.workflow_xml)
        root = tree.getroot()
        return root[0].attrib['id']
    else:
        return parsed_cfg.pipeline_id


def job_id_from_testkit_cfg(testkit_cfg):
    parsed_cfg = config_parser_to_butler(testkit_cfg)
    return parsed_cfg.job_id


# FIXME evil dwells here
def _patch_test_cases_with_service_access_layer(test_cases,
                                                service_access_layer, job_id):
    """This must be called before the test cases are run"""
    for test_case in test_cases:
        test_case.__class__.service_access_layer = service_access_layer
        test_case.__class__.job_id = job_id
        for t in test_case:
            t.__class__.service_access_layer = service_access_layer
            t.__class__.job_id = job_id


def run_butler_tests_from_cfg(testkit_cfg, output_dir, output_xml,
                              service_access_layer, services_job_id=None):
    job_id = job_id_from_testkit_cfg(testkit_cfg)
    butler = config_parser_to_butler(testkit_cfg)
    test_cases = parse_cfg_file(testkit_cfg)
    _patch_test_cases_with_service_access_layer(test_cases,
                                                service_access_layer,
                                                job_id=services_job_id)
    log.info("running tests...")
    exit_code = run_butler_tests(
        test_cases=test_cases,
        output_dir=output_dir,
        output_xml=output_xml,
        job_id=job_id)
    return exit_code


def run_services_testkit_job(host, port, testkit_cfg,
                             xml_out="test-output.xml",
                             ignore_test_failures=False,
                             time_out=1800, sleep_time=2,
                             import_only=False, test_job_id=None):
    """
    Given a testkit.cfg and host/port parameters:
        1. convert the .cfg to a JSON file
        2. connect to the SMRTLink services and start the job, then block
           until it finishes
        3. run the standard test suite on the job output
    """
    sal = ServiceAccessLayer(host, port, sleep_time=sleep_time)
    if test_job_id is not None:
        engine_job = sal.get_job_by_id(test_job_id)
        return run_butler_tests_from_cfg(
            testkit_cfg=testkit_cfg,
            output_dir=engine_job.path,
            output_xml=xml_out,
            service_access_layer=sal,
            services_job_id=test_job_id)
    entrypoints = get_entrypoints(testkit_cfg)
    pipeline_id = pipeline_id_from_testkit_cfg(testkit_cfg)
    job_id = job_id_from_testkit_cfg(testkit_cfg)
    log.info("job_id = {j}".format(j=job_id))
    log.info("pipeline_id = {p}".format(p=pipeline_id))
    log.info("url = {h}:{p}".format(h=host, p=port))
    task_options, workflow_options = get_task_and_workflow_options(testkit_cfg)
    service_entrypoints = [ServiceEntryPoint.from_d(x) for x in
                           entrypoints_dicts(entrypoints)]
    for ep, dataset_xml in entrypoints.iteritems():
        log.info("Importing {x}".format(x=dataset_xml))
        sal.run_import_local_dataset(dataset_xml)
    if import_only:
        log.info("Skipping job execution")
        return 0
    log.info("starting anaylsis job...")
    # XXX note that workflow options are currently ignored
    engine_job = run_analysis_job(sal, job_id, pipeline_id,
                                  service_entrypoints, block=True,
                                  time_out=time_out,
                                  task_options=task_options)
    exit_code = run_butler_tests_from_cfg(
        testkit_cfg=testkit_cfg,
        output_dir=engine_job.path,
        output_xml=xml_out,
        service_access_layer=sal,
        services_job_id=engine_job.id)
    if ignore_test_failures and engine_job.was_successful():
        return 0
    return exit_code


def args_runner(args):
    return run_services_testkit_job(
        host=args.host,
        port=args.port,
        testkit_cfg=args.testkit_cfg,
        xml_out=args.xml_out,
        ignore_test_failures=args.ignore_test_failures,
        time_out=args.time_out,
        sleep_time=args.sleep,
        import_only=args.import_only,
        test_job_id=args.test_job_id)


def get_parser():
    p = get_default_argparser_with_base_opts(
        version="0.1",
        description=__doc__)
    p.add_argument("testkit_cfg")
    p.add_argument("-u", "--host", dest="host", action="store",
                   default=os.environ.get("PB_SERVICE_HOST", "http://localhost"))
    p.add_argument("-p", "--port", dest="port", action="store", type=int,
                   default=int(os.environ.get("PB_SERVICE_PORT", "8081")),
                   help="Services port number")
    p.add_argument("-x", "--xunit", dest="xml_out", default="test-output.xml",
                   help="Output XUnit test results")
    p.add_argument("-t", "--timeout", dest="time_out", type=int, default=1800,
                   help="Timeout for blocking after job submission")
    p.add_argument("-s", "--sleep", dest="sleep", type=int, default=2,
                   help="Sleep time after job submission")
    p.add_argument("--ignore-test-failures", dest="ignore_test_failures",
                   action="store_true",
                   help="Only exit with non-zero return code if the job "+
                        "itself failed, regardless of test outcome")
    p.add_argument("--import-only", dest="import_only", action="store_true",
                   help="Import datasets without running pipeline")
    p.add_argument("--only-tests", dest="test_job_id", action="store",
                   type=int, default=None,
                   help="Run tests on an existing smrtlink job")
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
