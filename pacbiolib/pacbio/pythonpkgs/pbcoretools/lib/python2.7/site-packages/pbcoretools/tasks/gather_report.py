
import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import run_main_gather_report

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_report"
    CHUNK_KEY = "$chunk.report_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_report --resolved-tool-contract "
    OPT_CHUNK_KEY = 'pbcoretools.task_options.gather_report_chunk_key'
    REPORT_TYPE = FileTypes.REPORT


def get_parser(constants=Constants):
    p = get_gather_pbparser(constants.TOOL_ID,
                            constants.VERSION,
                            "Dev JSON Gather",
                            "General Chunk JSON Statistics Gather",
                            constants.DRIVER,
                            is_distributed=False)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with Json chunk key")

    p.add_output_file_type(constants.REPORT_TYPE, "json_out",
                           "JSON",
                           "Gathered JSON", "gathered")

    # Only need to add to argparse layer for the commandline
    p.arg_parser.add_str(constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         constants.CHUNK_KEY,
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
