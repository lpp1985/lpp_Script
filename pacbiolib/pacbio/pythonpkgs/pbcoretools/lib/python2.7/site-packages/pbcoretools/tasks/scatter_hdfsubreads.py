import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes, SymbolTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.h5_subreadset_scatter"


class Constants(object):
    DRIVER = "python -m pbcoretools.tasks.scatter_hdfsubreads --resolved-tool-contract "
    DEFAULT_NCHUNKS = 5
    CHUNK_KEYS = ('$chunk.hdf5subreadset_id', )


def get_contract_parser():
    p = get_scatter_pbparser(TOOL_ID, "0.1.3", "H5 SubreadSet scatter",
                             "Scatter Hdf5 Subread DataSet",
                             Constants.DRIVER,
                             Constants.CHUNK_KEYS,
                             is_distributed=True,
                             nchunks=SymbolTypes.MAX_NCHUNKS)

    p.add_input_file_type(FileTypes.DS_SUBREADS_H5, "h5_subreads", "HdfSubreadSet",
                          "Pac Bio Fasta format")

    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk HdfSubreadSet",
                           "PacBio Chunked JSON HdfSubread Set",
                           "hdfsubreadset_chunked")

    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_hdfsubread_max_nchunks", "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")

    p.add_str("pbcoretools.task_options.dev_scatter_chunk_key", "chunk_key",
              "$chunk:fasta_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(chunk_output_json, hdfsubread_xml, max_nchunks, output_dir):
    return CU.write_hdfsubreadset_chunks_to_file(chunk_output_json, hdfsubread_xml,
                                                 max_nchunks,
                                                 output_dir,
                                                 "chunk_hdfsubreadset",
                                                 FileTypes.DS_SUBREADS_H5.ext)


def _args_run_to_random_fasta_file(args):
    return run_main(args.chunk_report_json, args.hdfsubreadset,
                    args.max_total_chunks, args.output_dir)


def _rtc_runner(rtc):
    #max_nchunks = rtc.task.options['pbcoretools.task_options.dev_scatter_max_nchunks']
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0], max_nchunks, output_dir)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_run_to_random_fasta_file,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
