"""
Given a ChunkTasksPickle of a list of PolishChunkTask objects,
gather polished HQ|LQ isoforms in each bin.
"""

import logging
import sys
import os
import os.path as op
#from itertools import groupby

from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptOptions import BaseConstants
from pbtranscript.ClusterOptions import  IceQuiverHQLQOptions
from pbtranscript.ice.IceQuiverPostprocess import IceQuiverPostprocess
from pbtranscript.tasks.TPickles import ChunkTasksPickle, PolishChunkTask


log = logging.getLogger(__name__)


class Constants(BaseConstants):
    """Constants used for TOOL_ID"""
    TOOL_ID = "pbtranscript.tasks.gather_polished_isoforms_in_each_bin"
    DRIVER = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__
    VERSION = "0.1.0"


def get_contract_parser():
    """Tool contract should have the following inputs and outputs.
    Input:
        idx 0 - ChunkTasksPickle of PolishChunkTask objects.
        idx 1 - sentinel txt
    Output:
        idx 0 - gather_polished_isoforms_in_each_bin_done.txt,
                a sential file which does nothing
                but connect this task and the subsequence task
    """
    p = get_pbparser(tool_id=Constants.TOOL_ID,
                     version=Constants.VERSION,
                     name=Constants.TOOL_ID,
                     description=__doc__,
                     driver_exe=Constants.DRIVER)
    p.add_input_file_type(FileTypes.PICKLE, "polish_chunks_pickle", "Pickle In",
                          "Polish chunks pickle file") # input 0
    p.add_input_file_type(FileTypes.TXT, "polish_sentinel_in", "Sentinel In",
                          "Setinel file") # input idx 1
    p.add_output_file_type(FileTypes.TXT,
                           "gather polished isoforms in each bin done txt",
                           name="Gather polished isoforms in each bin Done Txt file",
                           description="Gather Done Txt file.",
                           default_name="gather_polished_isoforms_in_each_bin_done")

    p.add_float(BaseConstants.HQ_QUIVER_MIN_ACCURACY_ID, "hq_quiver_min_accuracy",
                default=BaseConstants.HQ_QUIVER_MIN_ACCURACY_DEFAULT,
                name="Minimum Accuracy of polished isoforms",
                description="Minimum Acuuracy of polished isoforms")
    p.add_int(BaseConstants.QV_TRIM_FIVEPRIME_ID, "qv_trim_5",
              default=BaseConstants.QV_TRIM_FIVEPRIME_DEFAULT,
              name="Trim QVs 5'", description="Ignore QV of n bases in the 5' end.")
    p.add_int(BaseConstants.QV_TRIM_THREEPRIME_ID, "qv_trim_3",
              default=BaseConstants.QV_TRIM_THREEPRIME_DEFAULT,
              name="Trim QVs 3'", description="Ignore QV of n bases in the 3' end.")
    return p


def args_runner(args):
    """args runner"""
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    """
    For each cluster bin, create summary.json, cluster_report.csv,
    hq_isoforms.fa|fq, lq_isoforms.fa|fq
    """
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    assert all([isinstance(task, PolishChunkTask) for task in p])
    p.sorted_by_attr(attr='cluster_bin_index')

    opts = rtc.task.options
    ipq_opts = IceQuiverHQLQOptions(qv_trim_5=opts[Constants.QV_TRIM_FIVEPRIME_ID],
                                    qv_trim_3=opts[Constants.QV_TRIM_THREEPRIME_ID],
                                    hq_quiver_min_accuracy=opts[Constants.HQ_QUIVER_MIN_ACCURACY_ID])

    with open(rtc.task.output_files[0], 'w') as writer:
        for cluster_bin_index, cluster_out_dir in p.sorted_no_redundant_cluster_bins():
            log.info("ice_quiver_postprocess of cluster bin index %s in %s.",
                     str(cluster_bin_index), str(cluster_out_dir))
            good_hq, bad_hq = \
            ice_quiver_postprocess_a_cluster_bin(cluster_out_dir=cluster_out_dir,
                                                 ipq_opts=ipq_opts)
            writer.write("ice_quiver_postprocess of cluster bin index %s in %s DONE:\n%s\n%s\n" %
                         (cluster_bin_index, cluster_out_dir, good_hq, bad_hq))


def ice_quiver_postprocess_a_cluster_bin(cluster_out_dir, ipq_opts):
    """
    ice_quiver_postprocess a cluster bin, create summary.json,
    cluster_report.csv, hq|lq_isoforms.fa|fq.
    return hq.fq and lq.fq
    Parameters:
      cluster_out_dir - root dir running ICE, ice_partial and ice_quiver
                        for this cluster bin.
    """
    _jobs_log = op.join(cluster_out_dir, "log", "submitted_quiver_jobs.txt")
    shell_scripts = []
    for file_name in os.listdir(op.join(cluster_out_dir, "quivered")):
        if file_name.endswith(".sh"):
            shell_scripts.append(file_name)
    with open(_jobs_log, 'w') as f:
        f.write("\n".join(["\t".join(["local", s]) for s in shell_scripts]))

    icep = IceQuiverPostprocess(root_dir=cluster_out_dir,
                                ipq_opts=ipq_opts)
    for file_name in os.listdir(op.join(cluster_out_dir, "quivered")):
        if file_name.endswith(".sh"):
            shell_scripts.append(file_name)

    icep.run()
    return (icep.quivered_good_fq, icep.quivered_bad_fq)


def main():
    """main"""
    mp = get_contract_parser()
    return pbparser_runner(
        argv=sys.argv[1:],
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
