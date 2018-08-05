
"""
Calls the ICE algorithm to identify de novo consensus isoforms
of cluster bins.
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import BaseConstants, \
        get_base_contract_parser, get_argument_parser
from pbtranscript.tasks.TPickles import ClusterChunkTask, \
        ChunkTasksPickle

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    """Constants used in pbtranscript.tasks.cluster_bins"""
    TOOL_ID = "pbtranscript.tasks.cluster_bins"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__


def get_contract_parser():
    """
    Return tool contract parser for running ICE clustering
    no quiver, no sge setting needed.
    The tool contract parser has 2 inputs:
        idx 0 - cluster_chunks.pickle (original or spawned)
        idx 1 - ccs,
    and has one outputs:
        idx 0 - cluster_chunk_done.txt, a sential file which does nothing
                but connect ICE cluster and the subsequence ice_partial
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    p.add_input_file_type(FileTypes.PICKLE, "cluster_chunks_pickle", "Pickle In",
                          "Cluster chunks pickle file") # input 0
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_in", "ConsensusReadSet In",
                          "PacBio ConsensusReadSet") # input 1
    p.add_output_file_type(FileTypes.TXT, "cluster done txt",
                           name="Cluster Done Txt file",
                           description="Cluster Done Txt file.",
                           default_name="cluster_chunks_done")
    return p


def args_runner(args):
    """args runner"""
    raise NotImplementedError()


def task_to_args(task, ccs_file, nproc, use_finer_qv=False):
    """Convert a ClusterChunkTask to 'pbtranscript cluster' args."""
    args = ["--verbose",
            "cluster",
            "--ccs_fofn", ccs_file,
            "--blasr_nproc", str(nproc)]
    if use_finer_qv:
        args.append("--use_finer_qv")

    assert op.exists(task.flnc_file)

    args.extend(["-d", task.cluster_out_dir,
                 task.flnc_file,
                 task.consensus_isoforms_file])
    return get_argument_parser().parse_args(args)


def resolved_tool_contract_runner(rtc):
    """run all tasks in cluster_chunks.pickle given rtc"""
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    assert all([isinstance(task, ClusterChunkTask) for task in p])

    ccs_file = rtc.task.input_files[1]
    assert op.exists(ccs_file)

    nproc = rtc.task.nproc
    use_finer_qv = False
    #if rtc.task.options.get(Constants.USE_FINER_QV_ID, False):
    #    use_finer_qv = True

    with open(rtc.task.output_files[0], 'w') as writer:
        for i, task in enumerate(p):
            args = task_to_args(task=task, ccs_file=ccs_file,
                                nproc=nproc, use_finer_qv=use_finer_qv)
            log.info("ARGUMENTS of Task %s/%s:\n%s", str(i), str(len(p)), str(args))
            log.info("Running ICE on cluster bin %s", task.cluster_bin_index)
            PBTranscript(args, subCommand="cluster").start()
            writer.write("ICE of cluster bin %s in %s is DONE: %s\n" %
                         (task.cluster_bin_index, task.cluster_out_dir,
                          task.consensus_isoforms_file))


def main():
    """main"""
    mp = get_contract_parser()
    return pbparser_runner(argv=sys.argv[1:],
                           parser=mp,
                           args_runner_func=args_runner,
                           contract_runner_func=resolved_tool_contract_runner,
                           alog=log,
                           setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
