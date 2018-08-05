"""CLI tool for helping getting metrics about the current configuration and setup


Runs the pbsmrtpipe.pipelines.dev_diagnostic pipelines.

ReferenceSet(eid_ref_dataset) -> Report


run-diagnostics preset.xml --simple --output-dir=./test-output # run the simplest hello-world job to queue
run-diagnostics preset.xml --output-dir=./test-output # run dev_diagnostic pipeline
"""
import os
import shutil
import sys
import logging
import argparse

from pbcommand.cli.utils import main_runner_default
from pbcommand.validators import validate_file
from pbcommand.cli import get_default_argparser

import pbsmrtpipe
from pbsmrtpipe.engine import run_command
from pbsmrtpipe.pb_io import parse_pipeline_preset_xml
from pbsmrtpipe.cluster import load_cluster_templates


log = logging.getLogger(__name__)


def _add_preset_xml_option(p):
    p.add_argument('preset_xml', type=validate_file, help="Path to Preset XML file.")
    return p


def _add_simple_mode_option(p):
    # Run the Full diagnostics suite
    p.add_argument('--simple', action='store_true', help="Perform full diagnostics tests (e.g., submit test job to cluster).")
    return p


def _test_can_find_temp_dir(path):
    if os.path.exists(path):
        log.info("Successfully found tmp-dir {t}".format(t=path))

    return True


def _test_can_load_cluster_templates(path):
    try:
        t = load_cluster_templates(path)
        return True
    except Exception as e:
        log.error("Failed to load cluster templates {t}".format(t=path))
        return False


def _to_echo_hello_world(output_file):
    return "echo 'Running hello world' > {o}".format(o=output_file)


def _write_echo_hello_world(output_file, run_sh):
    s = _to_echo_hello_world(output_file)
    with open(run_sh, 'w') as f:
        f.write(s)
    return s


def _to_path(root_dir):
    def _w(path):
        return os.path.join(root_dir, path)
    return _w


def run_simple_diagnostics(preset_xml, output_dir):
    """Setup simple job to run"""
    precord = parse_pipeline_preset_xml(preset_xml)

    wopts = precord.to_workflow_level_opt()
    to_p = _to_path(output_dir)

    ts = load_cluster_templates(wopts.cluster_manager_path)
    run_sh = to_p('run.sh')
    cluster_sh = to_p('cluster.sh')
    output_file = to_p('hello-world-output.txt')
    _write_echo_hello_world(output_file, run_sh)

    cluster_stderr = to_p("cluster.stderr")
    cluster_stdout = to_p("cluster.stdout")
    cluster_cmd = ts.render("start", run_sh, "job.dev-diagnostic-hello-world",
                            stdout=cluster_stdout, stderr=cluster_stderr)
    with open(cluster_sh, 'w') as f:
        f.write(cluster_cmd)

    print "Run.sh command {r}".format(r=run_sh)
    print "Exe'ing Cluster command {c}".format(c=cluster_cmd)
    rcode, stdout, stderr, run_time = run_command(cluster_cmd, sys.stdout, sys.stderr)

    if rcode == 0:
        print "Successfully submitted cluster job using templates {p}".format(p=wopts.cluster_manager_path)
    return rcode


def run_diagnostics(preset_xml, output_dir):
    """Run Hello World pipeline

    Submit to the cluster if configured
    """
    precord = parse_pipeline_preset_xml(preset_xml)
    wopts = precord.to_workflow_level_opt()

    to_p = _to_path(output_dir)

    input_txt = to_p("e-01_input.txt")
    with open(input_txt, 'w') as f:
        f.write("Mock data\n")

    job_preset_xml = to_p("preset.xml")
    shutil.copyfile(preset_xml, job_preset_xml)

    _d = dict(f=input_txt, p=job_preset_xml, d=output_dir)
    cmd_str = "pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.dev_dist -e \"e_01:{f}\" --preset-xml={p} --output-dir={d}"
    cmd = cmd_str.format(**_d)

    print "Running command {c}".format(c=cmd)
    rcode, stdout, stderr, run_time = run_command(cmd, sys.stdout, sys.stderr)

    if rcode == 0:
        print "Successfully submitted cluster job using templates {p}".format(p=wopts.cluster_manager_path)

    return rcode
