"""
Combine consensus isoforms of all cluster bins.
"""

import logging
import os.path as op
import sys

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key

from pbtranscript.__init__ import get_version
from pbtranscript.Utils import as_contigset
from pbtranscript.CombineUtils import CombinedFiles, combine_consensus_isoforms

log = logging.getLogger(__name__)


class Constants(object):
    """Constants used in pbtranscript.tasks.gather_consensus_isoforms_of_cluster_bins"""
    TOOL_ID = "pbtranscript.tasks.gather_consensus_isoforms_of_cluster_bins"
    VERSION = get_version()
    DRIVER = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__
    CHUNK_KEY = "$chunk.contigset_id"


def get_contract_parser():
    """Return contract parser.
    input idx 0: chunk json
    output idx 0: combined consensus isoforms contigset from all cluster bins
    """
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Gather consensus isoforms ContigSets",
                            "Chunk ContigSet Gather for Consensus Isoforms",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "Gather CHUNK Json",
                          "Gathered CHUNK Json with ContigSet chunk key")
    p.add_output_file_type(FileTypes.DS_CONTIG, "consensus_isoforms_out", "ContigSet",
                           "Combined ContigSet", "Combined ContigSet")
    return p


def run_main(chunk_json, contigset_output, chunk_key):
    """run main"""
    chunks = load_pipeline_chunks_from_json(chunk_json)

    # Allow looseness
    if not chunk_key.startswith('$chunk.'):
        chunk_key = '$chunk.' + chunk_key
        log.warn("Prepending chunk key with '$chunk.' to '%s'", str(chunk_key))

    fasta_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    log.debug("Chunked consensus isoforms files are %s.", (', '.join(fasta_files)))

    out_fa = CombinedFiles(combined_dir=op.dirname(contigset_output)).all_consensus_isoforms_fa
    combine_consensus_isoforms(split_indices=range(0, len(fasta_files)),
                               split_files=fasta_files,
                               combined_consensus_isoforms_fa=out_fa)
    log.info("Combining files to %s.", out_fa)

    log.info("Writing contigset %s", contigset_output)
    assert contigset_output.endswith('xml')
    as_contigset(out_fa, contigset_output)

    #cs = ContigSet(*fasta_files)
    #cs.newUuid()
    #cs.write(contigset_output)
    return 0


def args_runner(args):
    """Args runner"""
    raise NotImplementedError()


def rtc_runner(rtc):
    """
    :type rtc: pbcommand.models.ResolvedToolContract
    :return:
    """
    return run_main(chunk_json=rtc.task.input_files[0],
                    contigset_output=rtc.task.output_files[0],
                    chunk_key=Constants.CHUNK_KEY)#rtc.task.chunk_key


def main():
    """main"""
    mp = get_contract_parser()
    return pbparser_runner(sys.argv[1:],
                           mp, args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
