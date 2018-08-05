import logging
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import FileTypes, get_gather_pbparser

from pbcoretools.chunking.gather import run_main_gather_fastq_contigset

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_fastq"
    CHUNK_KEY = "$chunk.fastq_id"
    CHUNK_KEY_ID = "pbcoretools.task_options.dev_scatter_chunk_key"


def get_contract_parser():
    driver = "python -m pbcoretools.tasks.gather_fastq --resolved-tool-contract "

    p = get_gather_pbparser(Constants.TOOL_ID, "0.1.3", "Gather Fastq",
                            "Gather Fastq", driver, is_distributed=True)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "Gather ChunkJson",
                          "Fastq Gather Chunk JSON")

    p.add_output_file_type(FileTypes.FASTQ, "fastq_out", "Fastq Gathered",
                           "Fastq Gathered",
                           default_name="gathered")

    p.add_str(Constants.CHUNK_KEY_ID, "chunk_key",
              "$chunk:fastq_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def args_runner(args):
    return run_main_gather_fastq_contigset(args.cjson_in, args.fastq_out,
        args.chunk_key)


def rtc_runner(rtc):
    """
    :type rtc: pbcommand.models.ResolvedToolContract
    :return:
    """
    # the input file is just a sentinel file
    return run_main_gather_fastq_contigset(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp, args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
