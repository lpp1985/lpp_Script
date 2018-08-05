"""
Combine results from all cluster bins, including consensus isoforms,
polished HQ|LQ isoforms, cluster_summary.csv and cluster_report.json
"""

import logging
import os.path as op
import sys
import shutil

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.Utils import as_contigset, mkdir, get_sample_name, ln
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser,
                                              add_cluster_summary_report_arguments,
                                              add_ice_post_quiver_hq_lq_arguments)
from pbtranscript.ClusterOptions import  IceQuiverHQLQOptions
from pbtranscript.ice.IceQuiverPostprocess import IceQuiverPostprocess
from pbtranscript.io.Summary import write_cluster_summary
from pbtranscript.CombineUtils import CombinedFiles, combine_consensus_isoforms, \
        combine_polished_isoforms, write_combined_cluster_report
from pbtranscript.tasks.TPickles import ChunkTasksPickle, ClusterChunkTask


log = logging.getLogger(__name__)


class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.combine_cluster_bins"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__


def args_runner(args):
    """args runner."""
    raise NotImplementedError()


def get_contract_parser():
    """ Return resolved tool contract.
    input files:
        idx 0 - ChunkTasksPickle of ClusterChunkTask objects
        idx 1 - sentinel txt
    output files:
        idx 0 - consensus_isoforms.fa
        idx 1 - summary.json
        idx 2 - cluster_report.csv
        idx 3 - hq_isoforms.contigset.xml
        idx 4 - hq_isoforms.fq
        idx 5 - lq_isoforms.contigset.xml
        idx 6 - lq_isoforms.fq
        idx 7 - hq_lq_prefix_dict.pickle
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    p.add_input_file_type(FileTypes.PICKLE, "cluster_chunks_pickle", "Pickle In",
                          "Cluster chunks pickle file") # input 0
    p.add_input_file_type(FileTypes.TXT, "cluster_sentinel_in", "Sentinel In",
                          "Setinel file") # input idx 1

    # output idx 0, consensus_isoforms.contigset,
    p.add_output_file_type(FileTypes.DS_CONTIG, "consensus_isoforms",
                           name="Unpolished Consensus Isoforms",
                           description="Output unpolished consensus isoforms",
                           default_name="consensus_isoforms")

    # output idx 1, summary.json,
    # output idx 2, cluster_report.csv
    add_cluster_summary_report_arguments(p)

    # output idx 3, hq_isoforms_fa, idx 4, hq_isoforms_fq,
    # idx 5 lq_isoforms_fa, idx 6 lq_isoforms_fq
    add_ice_post_quiver_hq_lq_arguments(p)

    # output idx 7, hq_lq_prefix_dict.pickle
    p.add_output_file_type(FileTypes.PICKLE, "hq_lq_prefix_dict",
                           name="HQ LQ Sample Prefix Dict",
                           description="Pickle mapping HQ (LQ) sample prefixes to ICE dir",
                           default_name="hq_lq_prefix_dict")

    # user specified sample name.
    p.add_str(option_id=Constants.SAMPLE_NAME_ID, option_str="sample_name",
              default=Constants.SAMPLE_NAME_DEFAULT, name="sample Name",
              description="Sample Name")
    return p


def resolved_tool_contract_runner(rtc):
    """
    For each cluster bin, create summary.json, cluster_report.csv,
    hq_isoforms.fa|fq, lq_isoforms.fa|fq
    Finally, merge all cluster bins and save all outputs to 'combined'.
    """
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    assert all([isinstance(task, ClusterChunkTask) for task in p])
    p.sorted_by_attr(attr='cluster_bin_index')

    opts = rtc.task.options
    ipq_opts = IceQuiverHQLQOptions(qv_trim_5=opts[Constants.QV_TRIM_FIVEPRIME_ID],
                                    qv_trim_3=opts[Constants.QV_TRIM_THREEPRIME_ID],
                                    hq_quiver_min_accuracy=opts[Constants.HQ_QUIVER_MIN_ACCURACY_ID])
    sample_name = get_sample_name(input_sample_name=opts[Constants.SAMPLE_NAME_ID])

    out_consensus_isoforms_cs = rtc.task.output_files[0]
    out_summary = rtc.task.output_files[1]
    out_report = rtc.task.output_files[2]
    out_hq_cs = rtc.task.output_files[3]
    out_hq_fq = rtc.task.output_files[4]
    out_lq_cs = rtc.task.output_files[5]
    out_lq_fq = rtc.task.output_files[6]
    out_hq_lq_prefix_dict_pickle = rtc.task.output_files[7]

    assert out_consensus_isoforms_cs.endswith(".contigset.xml")
    assert out_hq_cs.endswith(".contigset.xml")
    assert out_lq_cs.endswith(".contigset.xml")
    out_consensus_isoforms_fa = out_consensus_isoforms_cs.replace(".contigset.xml", ".fasta")
    out_hq_fa = out_hq_cs.replace('.contigset.xml', '.fasta')
    out_lq_fa = out_lq_cs.replace('.contigset.xml', '.fasta')

    hq_fq_fns, lq_fq_fns = [], []
    split_uc_pickles, split_partial_uc_pickles = [], []
    split_consensus_isoforms = []

    cluster_bin_indices = [task.cluster_bin_index for task in p]
    cluster_out_dirs = [task.cluster_out_dir for task in p]
    # sanity check that Cluster indices are unique!
    assert len(set(cluster_bin_indices)) == len(cluster_bin_indices)

    for task in p:
        ice_pq = IceQuiverPostprocess(root_dir=task.cluster_out_dir,
                                      ipq_opts=ipq_opts)
        hq_fq_fns.append(ice_pq.quivered_good_fq)
        lq_fq_fns.append(ice_pq.quivered_bad_fq)
        split_uc_pickles.append(ice_pq.final_pickle_fn)
        split_partial_uc_pickles.append(ice_pq.nfl_all_pickle_fn)
        split_consensus_isoforms.append(ice_pq.final_consensus_fa)

    combined_dir = op.join(op.dirname(op.dirname(cluster_out_dirs[0])), "combined")
    mkdir(combined_dir)
    combined_files = CombinedFiles(combined_dir)
    log.info("Combining results of all cluster bins to %s.", combined_dir)
    log.info("Merging HQ|LQ isoforms from all cluster bins.")
    log.info("HQ isoforms are: %s.", ",".join(hq_fq_fns))
    log.info("LQ isoforms are: %s.", ",".join(lq_fq_fns))
    combine_polished_isoforms(split_indices=cluster_bin_indices,
                              split_hq_fns=hq_fq_fns,
                              split_lq_fns=lq_fq_fns,
                              combined_hq_fa=combined_files.all_hq_fa,
                              combined_hq_fq=combined_files.all_hq_fq,
                              combined_lq_fa=combined_files.all_lq_fa,
                              combined_lq_fq=combined_files.all_lq_fq,
                              hq_lq_prefix_dict_pickle=combined_files.hq_lq_prefix_dict_pickle,
                              sample_name=sample_name)

    ln(combined_files.all_hq_fa, out_hq_fa) #'HQ isoforms'
    ln(combined_files.all_hq_fq, out_hq_fq) #'HQ isoforms'
    ln(combined_files.all_lq_fa, out_lq_fa) #'LQ isoforms'
    ln(combined_files.all_lq_fq, out_lq_fq) #'LQ isoforms'
    ln(combined_files.hq_lq_prefix_dict_pickle, out_hq_lq_prefix_dict_pickle)

    as_contigset(out_hq_fa, out_hq_cs)
    as_contigset(out_lq_fa, out_lq_cs)

    log.info("Merging consensus isoforms from all cluster bins.")
    combine_consensus_isoforms(split_indices=cluster_bin_indices,
                               split_files=split_consensus_isoforms,
                               combined_consensus_isoforms_fa=combined_files.all_consensus_isoforms_fa,
                               sample_name=sample_name)
    ln(combined_files.all_consensus_isoforms_fa, out_consensus_isoforms_fa)
    #consensus isoforms
    as_contigset(out_consensus_isoforms_fa, out_consensus_isoforms_cs)

    log.info("Writing cluster summary to %s", combined_files.all_cluster_summary_fn)
    write_cluster_summary(summary_fn=combined_files.all_cluster_summary_fn,
                          isoforms_fa=out_consensus_isoforms_cs,
                          hq_fa=out_hq_fa,
                          lq_fa=out_lq_fa)
    ln(combined_files.all_cluster_summary_fn, out_summary) # "cluster summary"

    log.info("Writing cluster report to %s", combined_files.all_cluster_report_fn)
    write_combined_cluster_report(split_indices=cluster_bin_indices,
                                  split_uc_pickles=split_uc_pickles,
                                  split_partial_uc_pickles=split_partial_uc_pickles,
                                  report_fn=combined_files.all_cluster_report_fn,
                                  sample_name=sample_name)
    ln(combined_files.all_cluster_report_fn, out_report) # "cluster report"


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
