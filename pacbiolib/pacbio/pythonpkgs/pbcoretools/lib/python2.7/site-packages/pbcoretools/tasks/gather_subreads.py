import logging
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
import sys
from pbcommand.utils import setup_log
from pbcoretools.chunking.gather import run_main_gather_subreadset

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_subreadset"
    CHUNK_KEY = "$chunk.subreadset_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_subreads --resolved-tool-contract "


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev SubreadSet Gather",
                            "General Chunk SubreadSet Gather",
                            Constants.DRIVER,
                            is_distributed=True)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with SubreadSet chunk key")

    p.add_output_file_type(FileTypes.DS_SUBREADS, "ds_out", "SubreadSet",
                           "Gathered SubreadSet", default_name="gathered")
    return p


def args_runner(args):
    return run_main_gather_subreadset(args.cjson_in, args.ds_out,
                                      Constants.CHUNK_KEY)

def rtc_runner(rtc):
    return run_main_gather_subreadset(rtc.task.input_files[0],
        rtc.task.output_files[0], rtc.task.chunk_key)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
