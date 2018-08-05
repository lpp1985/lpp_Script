#! python
"""
Class DalignerRunner which aligns a query FASTA file and
a target FASTA file using daligner and interpret las output
using LA4Ice.

"""

import os.path as op
import sys
import logging
import argparse

from pbcore.util.ToolRunner import PBToolRunner
from pbtranscript.__init__ import get_version
from pbtranscript.ClusterOptions import IceQuiverHQLQOptions
from pbtranscript.Utils import mkdir, get_sample_name, ln
from pbtranscript.PBTranscriptOptions import BaseConstants
from pbtranscript.ice.IceQuiverPostprocess import IceQuiverPostprocess
from pbtranscript.separate_flnc import convert_pickle_to_sorted_flnc_files
from pbtranscript.io.Summary import write_cluster_summary
from pbtranscript.CombineUtils import CombinedFiles, combine_consensus_isoforms, \
        combine_polished_isoforms, write_combined_cluster_report


log = logging.getLogger(__name__)


__author__ = "yli@pacificbiosciences.com"


def add_ice_combine_cluster_bins_arguments(parser):
    """Set up argument parser."""
    helpstr = "A pickle file generated by separate_flnc, from which " + \
              "cluster bin output directory of each separated flnc file " + \
              "will be inferred."
    parser.add_argument("separate_flnc_pickle", type=str, help=helpstr)

    helpstr = "Output directory to save all combined files"
    parser.add_argument("combined_dir", type=str, help=helpstr)

    helpstr = "A list of comma delimited cluster bin directories. " + \
              "If set, will overwrite the default cluster bin directories " + \
              "inferred from separated flnc fasta files " + \
              "(e.g., 3to4kb_part0/cluster_out,4to5kb_part0/cluster_out)"
    parser.add_argument("--cluster_bin_dirs", type=str, default=None, help=helpstr)

    # user specified sample name.
    parser.add_argument("--sample_name", help="Sample Name", type=str, default=None)

    # combined output files
    ogroup = parser.add_argument_group("Combined output files")
    helpstr = "Combined consensus isoforms fasta, default: ${combined_dir}/all.consensus_isoforms.fasta"
    ogroup.add_argument("--consensus_isoforms_fa", type=str, default=None, help=helpstr)

    helpstr = "Combined cluster report csv, default: ${combined_dir}/all.cluster_report.csv"
    ogroup.add_argument("--report", default=None, type=str, dest="report_fn", help=helpstr)

    helpstr = "Combined cluster summary json, default: ${combined_dir}/cluster_summary.json"
    ogroup.add_argument("--summary", default=None, type=str, dest="summary_fn", help=helpstr)

    # hq_isoforms_fa, hq_isoforms_fq, lq_isoforms_fa, lq_isoforms_fq
    helpstr = "Combined HQ isoforms fasta, default: ${combined_dir}/all.polished_hq.fasta"
    ogroup.add_argument("--hq_isoforms_fa", default=None, type=str, help=helpstr)
    helpstr = "Combined HQ isoforms fastq, default: ${combined_dir}/all.polished_hq.fastq"
    ogroup.add_argument("--hq_isoforms_fq", default=None, type=str, help=helpstr)

    helpstr = "Combined LQ isoforms fasta, default: ${combined_dir}/all.polished_lq.fasta"
    ogroup.add_argument("--lq_isoforms_fa", default=None, type=str, help=helpstr)
    helpstr = "Combined LQ isoforms fastq, default: ${combined_dir}/all.polished_lq.fastq"
    ogroup.add_argument("--lq_isoforms_fq", default=None, type=str, help=helpstr)

    icq_gp = parser.add_argument_group("HQ LQ QV arguments used in cluster bins", argparse.SUPPRESS)
    icq_gp.add_argument("--hq_quiver_min_accuracy", type=float,
                        default=BaseConstants.HQ_QUIVER_MIN_ACCURACY_DEFAULT,
                        dest="hq_quiver_min_accuracy", help=argparse.SUPPRESS)
    icq_gp.add_argument("--qv_trim_5", type=int,
                        default=BaseConstants.QV_TRIM_FIVEPRIME_DEFAULT,
                        dest="qv_trim_5", help=argparse.SUPPRESS)
    icq_gp.add_argument("--qv_trim_3", type=int,
                        default=BaseConstants.QV_TRIM_THREEPRIME_DEFAULT,
                        dest="qv_trim_3", help=argparse.SUPPRESS)
    return parser


class CombineClusterBinsRunner(PBToolRunner):
    """
    Combine cluster bin output from multiple separated flnc reads.
    """
    def __init__(self):
        desc = __doc__
        PBToolRunner.__init__(self, desc)
        add_ice_combine_cluster_bins_arguments(self.parser)

    def getVersion(self):
        """Get version string."""
        return get_version()

    def get_cluster_bin_dirs(self, separate_flnc_pickle, cluster_bin_dirs):
        """
        Get a list of cluster bin directories.

        If cluster_bin_dirs is not None, it is a comma delimited list of cluster bin dirs.
        Otherwise, cluster bin dirs can be inferred from separated flnc fasta files, while
        each sepated flnc fasta file must correpsond to a cluster bin.

            3to4kb_part0
            3to4kb_part0/isoseq_flnc.fasta # sepearted flnc fasta
            3to4kb_part0/cluster_out       # cluster bin dir
            4to5kb_part0
            4to5kb_part0/isoseq_flnc.fasta # sepearted flnc fasta
            4to5kb_part0/cluster_out       # cluster bin dir
        """
        ret = []
        if cluster_bin_dirs is None:
            # Get all separated flnc fasta file.
            flnc_fas = convert_pickle_to_sorted_flnc_files(separate_flnc_pickle)
            ret = [op.join(op.dirname(flnc_fa), "cluster_out") for flnc_fa in flnc_fas]
        else:
            ret = cluster_bin_dirs.split(',')

        for d in ret:
            if not op.exists(d):
                raise IOError("Could not find cluster bin directory %s" % d)
        return ret

    def run(self):
        """
        For each cluster bin, create summary.json, cluster_report.csv,
        hq_isoforms.fa|fq, lq_isoforms.fa|fq
        Finally, merge all cluster bins and save all outputs to 'combined'.
        """
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args

        # Get cluster bins directories as input
        cluster_bin_dirs = self.get_cluster_bin_dirs(separate_flnc_pickle=args.separate_flnc_pickle,
                                                     cluster_bin_dirs=args.cluster_bin_dirs)
        cluster_bin_indices = range(0, len(cluster_bin_dirs))

        # Create output dir
        combined_dir = args.combined_dir
        mkdir(combined_dir)

        # Get combined output filenames
        def f(input_fn, default_fn):
            if input_fn is None:
                return op.join(combined_dir, default_fn)

        out_consensus_isoforms_fa = f(args.consensus_isoforms_fa, "all.consensus_isoforms.fasta")
        out_summary = f(args.summary_fn, "all.cluster_summary.json")
        out_report = f(args.report_fn, "all.cluster_report.csv")
        out_hq_fa = f(args.hq_isoforms_fa, "all.polished_hq.fasta")
        out_lq_fa = f(args.lq_isoforms_fa, "all.polished_lq.fasta")
        out_hq_fq = f(args.hq_isoforms_fq, "all.polished_hq.fastq")
        out_lq_fq = f(args.lq_isoforms_fq, "all.polished_lq.fastq")

        ipq_opts = IceQuiverHQLQOptions(qv_trim_5=args.qv_trim_5,
                                        qv_trim_3=args.qv_trim_3,
                                        hq_quiver_min_accuracy=args.hq_quiver_min_accuracy)
        sample_name = get_sample_name(input_sample_name=args.sample_name)


        hq_fq_fns, lq_fq_fns = [], []
        split_uc_pickles, split_partial_uc_pickles = [], []
        split_consensus_isoforms = []

        for cluster_bin_dir in cluster_bin_dirs:
            ice_pq = IceQuiverPostprocess(root_dir=cluster_bin_dir, ipq_opts=ipq_opts)
            hq_fq_fns.append(ice_pq.quivered_good_fq)
            lq_fq_fns.append(ice_pq.quivered_bad_fq)
            split_uc_pickles.append(ice_pq.final_pickle_fn)
            split_partial_uc_pickles.append(ice_pq.nfl_all_pickle_fn)
            split_consensus_isoforms.append(ice_pq.final_consensus_fa)

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

        log.info("Merging consensus isoforms from all cluster bins.")
        combine_consensus_isoforms(split_indices=cluster_bin_indices,
                                   split_files=split_consensus_isoforms,
                                   combined_consensus_isoforms_fa=combined_files.all_consensus_isoforms_fa,
                                   sample_name=sample_name)
        ln(combined_files.all_consensus_isoforms_fa, out_consensus_isoforms_fa)

        log.info("Writing cluster summary to %s", combined_files.all_cluster_summary_fn)
        write_cluster_summary(summary_fn=combined_files.all_cluster_summary_fn,
                              isoforms_fa=out_consensus_isoforms_fa,
                              hq_fa=out_hq_fa, lq_fa=out_lq_fa)
        ln(combined_files.all_cluster_summary_fn, out_summary) # "cluster summary"

        log.info("Writing cluster report to %s", combined_files.all_cluster_report_fn)
        write_combined_cluster_report(split_indices=cluster_bin_indices,
                                      split_uc_pickles=split_uc_pickles,
                                      split_partial_uc_pickles=split_partial_uc_pickles,
                                      report_fn=combined_files.all_cluster_report_fn,
                                      sample_name=sample_name)
        ln(combined_files.all_cluster_report_fn, out_report) # "cluster report"


def main():
    """Main function to call CombineClusterBinsRunner"""
    runner = CombineClusterBinsRunner()
    return runner.start()


if __name__ == "__main__":
    sys.exit(main())