
"""
Alternate report gather for the classify task, which deals sensibly with
averaged values.
"""

import logging
import sys

from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.cli import pbparser_runner
from pbcommand.models.report import Report
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import get_main_runner

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbtranscript.tasks.gather_isoseq_classify_report"
    CHUNK_KEY = "$chunk.report_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbtranscript.tasks.gather_isoseq_classify_report --resolved-tool-contract "
    OPT_CHUNK_KEY = "pbsmrtpipe.task_options.gather_report_chunk_key"

def gather_report(json_files, output_file):
    """
    Combines statistics (usually raw counts) stored as JSON files.
    Data models: pbcommand.models.report
    """
    reports = [ load_report_from_json(fn) for fn in json_files ]
    merged = Report.merge(reports)
    total_num_flnc_bases = 0
    total_num_flnc = 0
    for r in reports:
        attrs = {a.id:a.value for a in r.attributes}
        num_flnc = attrs["num_flnc"]
        num_flnc_bases = attrs["num_flnc_bases"]
        total_num_flnc += num_flnc
        total_num_flnc_bases += num_flnc_bases
    if total_num_flnc > 0:
        for a in merged.attributes:
            if a.id == "avg_flnc_len":
                # mimicking pbtranscript.io.Summary
                a._value = int(total_num_flnc_bases / total_num_flnc)
                log.info("Setting avg_flnc_len = {v}".format(v=a.value))
    with open(output_file, "w") as writer:
        writer.write(merged.to_json())
    return output_file


run_main_gather_report = get_main_runner(gather_report)


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev JSON Gather",
                            "General Chunk JSON Statistics Gather",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with Json chunk key")
    p.add_output_file_type(FileTypes.JSON, "json_out",
                           "JSON",
                           "Gathered JSON", "gathered")
    # Only need to add to argparse layer for the commandline
    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         Constants.CHUNK_KEY,
                         "Chunk key",
                         "Chunk key to use (format $chunk.{chunk-key}")
    return p


def args_runner(args):
    return run_main_gather_report(args.cjson_in, args.json_out, args.chunk_key)


def rtc_runner(rtc):
    return run_main_gather_report(rtc.task.input_files[0], rtc.task.output_files[0], Constants.CHUNK_KEY)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
