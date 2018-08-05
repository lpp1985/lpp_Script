import logging
import os
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.scatter_contigset"


class Constants(object):
    DEFAULT_NCHUNKS = 24
    CHUNK_KEY = '$chunk.contigset_id'


def get_contract_parser():
    driver = "python -m pbcoretools.tasks.scatter_contigset --resolved-tool-contract "

    chunk_keys = ("$chunk.contigset_id", )
    p = get_scatter_pbparser(TOOL_ID, "0.1.3", "Scatter ContigSet",
                             "Scatter ContigSet", driver, chunk_keys,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_CONTIG, "dataset_in", "ContigSet In",
                          "PacBio ContigSet")
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out", "Chunk JSON Filtered Fasta",
                           "Chunked JSON ContigSet",
                           "fasta.chunked")
    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.dev_scatter_max_nchunks", "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    p.add_str("pbcoretools.task_options.dev_scatter_chunk_key", "chunk_key",
              Constants.CHUNK_KEY, "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(dataset_file, output_json, max_nchunks, chunk_key):
    log.info("Running {f} into {n} chunks".format(f=dataset_file, n=max_nchunks))
    output_dir = os.path.dirname(output_json)
    CU.write_contigset_chunks_to_file(output_json, dataset_file, max_nchunks,
                                      output_dir, "scattered-contigset",
                                      FileTypes.DS_CONTIG.ext)
    return 0


def _args_run(args):
    return run_main(args.contigset_in, args.cjson_out, args.max_nchunks, args.chunk_key)


def _rtc_runner(rtc):
    # the chunk key isn't really something that can be tweaked here.
    return run_main(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.max_nchunks, Constants.CHUNK_KEY)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_run,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
