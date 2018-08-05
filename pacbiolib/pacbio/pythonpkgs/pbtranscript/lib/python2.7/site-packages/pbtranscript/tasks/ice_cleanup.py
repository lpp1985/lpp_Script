# clean up intermediate files after ice_quiver_postprocess.

"""
Clean up intermediate files after ice_polish_cluster_bins complete.
Inputs: separate_flnc.pickle
Outputs: None
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.Utils import execute
from pbtranscript.ice.IceFiles import IceFiles
from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.tasks.TPickles import ChunkTasksPickle, ClusterChunkTask

log = logging.getLogger(__name__)


class Constants(object):
    """Constants used TOOL_ID"""
    TOOL_ID = "pbtranscript.tasks.ice_cleanup"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__


def args_runner(args):
    """args runner."""
    raise NotImplementedError()


def get_contract_parser():
    """ Return resolved tool contract.
    input files:
        idx 0 - ChunkTasksPickle of ClusterChunkTask objects
        idx 1 - ice_polish_cluster_bins sentinel txt
    output files:
        idx 0 - sentinal
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    p.add_input_file_type(FileTypes.PICKLE, "cluster_chunks_pickle", "Pickle In",
                          "Cluster chunks pickle file") # input idx 0
    p.add_input_file_type(FileTypes.TXT, "ice_polish_cluster_bins_sentinel", "Sentinel In",
                          "Setinel file") # input idx 1

    # output idx 0, consensus_isoforms.contigset,
    p.add_output_file_type(FileTypes.TXT, "ice_cleanup_sentinel",
                           name="Sentinel file",
                           description="Output sentinel file",
                           default_name="ice_cleanup_done")
    return p


def resolved_tool_contract_runner(rtc):
    """
    For each cluster bin, clean up intermediate files under tmp.
    """
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    assert all([isinstance(task, ClusterChunkTask) for task in p])

    cluster_bin_indices = [task.cluster_bin_index for task in p]
    # sanity check that Cluster indices are unique!
    assert len(set(cluster_bin_indices)) == len(cluster_bin_indices)

    sentinel_out = rtc.task.output_files[0]
    with open(sentinel_out, 'w') as writer:
        for task in p:
            icef = IceFiles(prog_name="ice_cleanup",
                            root_dir=task.cluster_out_dir)
            tmp_dir = icef.tmp_dir
            log.info("Cleaning up, removing %s", tmp_dir)
            writer.write("removing %s\n" % tmp_dir)
            execute("rm -rf %s" % tmp_dir)

            quivered_dir = icef.quivered_dir
            log.info("Cleaning up, removing %s", quivered_dir)
            writer.write("removing %s\n" % quivered_dir)
            execute("rm -rf %s" % quivered_dir)


def main():
    """Main"""
    mp = get_contract_parser()
    return pbparser_runner(argv=sys.argv[1:],
                           parser=mp,
                           args_runner_func=args_runner,
                           contract_runner_func=resolved_tool_contract_runner,
                           alog=log,
                           setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
