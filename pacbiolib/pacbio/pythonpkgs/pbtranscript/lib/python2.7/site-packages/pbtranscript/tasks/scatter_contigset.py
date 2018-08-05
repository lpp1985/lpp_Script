
"""
Specialized ContigSet scatter, incorporating unchunked files required by
ice_partial
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
    TOOL_ID = "pbtranscript.tasks.scatter_contigset"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m pbtranscript.tasks.scatter_contigset --resolved-tool-contract "
    CHUNK_KEYS = ('$chunk.contigset_id', '$chunk.ref_contigset_id',
                  '$chunk.ccsset_id')


def get_contract_parser():

    p = get_scatter_pbparser(Constants.TOOL_ID, Constants.VERSION,
                             "Scatter ContigSet",
                             __doc__, Constants.DRIVER_EXE,
                             chunk_keys=Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_CONTIG, "dataset_in", "ContigSet In",
                          "PacBio ContigSet")
    p.add_input_file_type(FileTypes.DS_CONTIG, "ref_in", "Reference ContigSet",
                          "PacBio ContigSet")
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_in", "ConsensusReadSet",
                          "PacBio ConsensusRead DataSet")
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk JSON Filtered Fasta",
                           "Chunked JSON ContigSet",
                           "fasta.chunked")
    # max nchunks for this specific task
    p.add_int("pbsmrtpipe.task_options.dev_scatter_max_nchunks", "max_nchunks",
              Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(dataset_file, ref_dataset_file, ccs_file, output_json,
             max_nchunks):
    # FIXME this makes unit testing difficult
    log.info("Running {f} into {n} chunks".format(f=dataset_file,
                                                  n=max_nchunks))
    output_dir = os.path.dirname(output_json)
    CU.write_contigset_chunks_to_file(output_json, dataset_file, max_nchunks,
                                      output_dir, "scattered-contigset", "contigset.xml",
                                      extra_chunk_keys={
                                          '$chunk.ref_contigset_id': ref_dataset_file,
                                          '$chunk.ccsset_id': ccs_file,
                                      })
    return 0


def _args_run(args):
    raise NotImplementedError()


def _rtc_runner(rtc):
    # the chunk key isn't really something that can be tweaked here.
    return run_main(rtc.task.input_files[0], rtc.task.input_files[1], rtc.task.input_files[2], rtc.task.output_files[0], rtc.task.max_nchunks)


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
