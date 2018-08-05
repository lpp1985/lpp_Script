import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import run_main_gather_txt

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_txt"
    CHUNK_KEY = "$chunk.txt_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_txt --resolved-tool-contract "
    OPT_CHUNK_KEY = 'pbcoretools.task_options.gather_txt_chunk_key'


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev text Gather",
                            "General Chunk text Gather",
                            Constants.DRIVER,
                            is_distributed=False)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with text chunk key")

    p.add_output_file_type(FileTypes.TXT, "txt_out",
                           "Text",
                           "Gathered text", "gathered")

    # Only need to add to argparse layer for the commandline
    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         "$chunk.txt_id",
                         "Chunk key",
                         "Chunk key to use (format $chunk.{chunk-key}")

    return p


def args_runner(args):
    return run_main_gather_txt(args.cjson_in, args.txt_out, args.chunk_key)


def rtc_runner(rtc):
    return run_main_gather_txt(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
