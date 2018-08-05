
"""
Scatter AlignmentSet, ReferenceSet -> Chunk.JSON

This is effectively an 
"""
import logging
import os
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.tasks import scatter_alignments_reference

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.alignment_contig_scatter_basemods"
MODULE_NAME = "pbcoretools.tasks.scatter_alignments_reference_basemods"

def main(argv=sys.argv):
    mp = scatter_alignments_reference.get_contract_parser(
        tool_id=TOOL_ID,
        module_name=MODULE_NAME)
    return pbparser_runner(argv[1:],
                           mp,
                           scatter_alignments_reference.args_runner,
                           scatter_alignments_reference.rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
