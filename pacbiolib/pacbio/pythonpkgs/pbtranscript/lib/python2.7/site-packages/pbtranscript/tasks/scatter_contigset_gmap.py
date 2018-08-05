
"""
Specialized ContigSet scatter, incorporating unchunked reference required by
GMAP.
"""

import logging
import os
import sys

import pbcoretools.chunking.chunk_utils as CU
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbtranscript.tasks.scatter_contigset_gmap"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m pbtranscript.tasks.scatter_contigset_gmap --resolved-tool-contract "
    CHUNK_KEYS = ('$chunk.contigset_id', '$chunk.reference_id',)


def get_contract_parser():

    p = get_scatter_pbparser(Constants.TOOL_ID, Constants.VERSION,
                             "Scatter ContigSet for GMAP",
                             __doc__, Constants.DRIVER_EXE,
                             chunk_keys=Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_CONTIG, "dataset_in", "ContigSet In",
                          "PacBio ContigSet")
    p.add_input_file_type(FileTypes.DS_REF, "ref_in", "ReferenceSet",
                          "PacBio ReferenceSet")
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk JSON Filtered Fasta",
                           "Chunked JSON ContigSet",
                           "fasta.chunked")
    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.dev_scatter_max_nchunks", "max_nchunks",
              Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(dataset_file, ref_dataset_file, output_json,
             max_nchunks):
    log.info("Running {f} into {n} chunks".format(f=dataset_file,
                                                  n=max_nchunks))
    output_dir = os.path.dirname(output_json)
    CU.write_contigset_chunks_to_file(
        chunk_file=output_json,
        dataset_path=dataset_file,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="scattered-contigset",
        chunk_ext="contigset.xml",
        extra_chunk_keys={"$chunk.reference_id": ref_dataset_file})
    return 0


def _args_run(args):
    raise NotImplementedError()


def _rtc_runner(rtc):
    # the chunk key isn't really something that can be tweaked here.
    return run_main(
        dataset_file=rtc.task.input_files[0],
        ref_dataset_file=rtc.task.input_files[1],
        output_json=rtc.task.output_files[0],
        max_nchunks=rtc.task.max_nchunks)


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
