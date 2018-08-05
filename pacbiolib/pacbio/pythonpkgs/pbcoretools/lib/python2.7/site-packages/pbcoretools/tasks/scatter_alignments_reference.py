"""
Scatter AlignmentSet, ReferenceSet -> Chunk.JSON
"""
import logging
import os
import sys

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.alignment_contig_scatter"
MODULE_NAME = "pbcoretools.tasks.scatter_alignments_reference"


class Constants(object):
    DEFAULT_NCHUNKS = 12
    CHUNK_KEYS = ('$chunk.alignmentset_id', "$chunk.reference_id")
    DRIVER_BASE = "python -m {module} --resolved-tool-contract "
    OPT_CHUNK_KEY = "pbcoretools.task_options.dev_scatter_chunk_key"
    OPT_MAX_NCHUNKS = 'pbcoretools.task_options.scatter_alignments_reference_max_nchunks'


def get_contract_parser(tool_id=TOOL_ID, module_name=MODULE_NAME):
    p = get_scatter_pbparser(tool_id, "0.1.3",
                             "Scatter AlignmentSet",
                             "Pacbio DataSet AlignmentSet",
                             Constants.DRIVER_BASE.format(module=module_name),
                             Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_ALIGN,
                          "alignment_ds",
                          "AlignmentSet",
                          "Pacbio DataSet AlignmentSet")

    p.add_input_file_type(FileTypes.DS_REF,
                          "ds_reference",
                          "ReferenceSet",
                          "Pac Bio Fasta format")

    p.add_output_file_type(FileTypes.CHUNK,
                           "cjson_out",
                           "Chunk JSON Filtered Fasta",
                           "Chunked JSON Filtered Fasta",
                           "alignments_reference.chunked")

    # max nchunks for this specific task
    p.add_int(Constants.OPT_MAX_NCHUNKS,
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(ds_xml, reference_set_xml, output_json, max_nchunks, output_dir):
    CU.write_alignmentset_chunks_to_file(output_json,
                                         ds_xml,
                                         reference_set_xml,
                                         max_nchunks,
                                         output_dir,
                                         "chunk_alignmentset",
                                         FileTypes.DS_ALIGN.ext)
    return 0


def args_runner(args):
    # FIXME. The chunk needs to be passed directly the func
    chunk_key = args.chunk_key
    #chunk_key = "alignmentset_id"
    output_dir = os.path.dirname(args.cjson_out)
    return run_main(args.alignment_ds, args.ds_reference, args.cjson_out, args.max_nchunks, output_dir)


def rtc_runner(rtc):
    return run_main(rtc.task.input_files[0],
                    rtc.task.input_files[1],
                    rtc.task.output_files[0],
                    rtc.task.max_nchunks,
                    os.path.dirname(rtc.task.output_files[0]))


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
