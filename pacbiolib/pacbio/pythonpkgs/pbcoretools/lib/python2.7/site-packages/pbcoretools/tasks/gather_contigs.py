import sys
import logging

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import run_main_gather_contigset

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_contigset"
    CHUNK_KEY = '$chunk.contigset_id'
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_contigs --resolved-tool-contract "
    OPT_CHUNK_KEY = '$chunk.contigset_id'


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev ContigSet Gather",
                            "General Chunk ContigSet Gather",
                            Constants.DRIVER,
                            is_distributed=True)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with ContigSet chunk key")

    p.add_output_file_type(FileTypes.DS_CONTIG,
                           "contigset",
                           "ContigSet",
                           "Gathered ContigSet",
                           "gathered")

    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         Constants.CHUNK_KEY,
                         "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def args_runner(args):
    return run_main_gather_contigset(args.cjson_in,
                                     args.contigset,
                                     args.chunk_key)


def rtc_runner(rtc):
    return run_main_gather_contigset(rtc.task.input_files[0],
                                     rtc.task.output_files[0],
                                     rtc.task.chunk_key)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
