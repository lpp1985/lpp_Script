
"""
Generates a report of statistics for CCS reads mapped to a reference genome
with Blasr/pbalign.
"""

from collections import OrderedDict
import logging
import sys
import os
import os.path as op
import numpy as np

from pbcommand.models import get_pbparser, FileTypes

from pbreports.report.streaming_utils import PlotViewProperties
from pbreports.plot.helper import (get_fig_axes_lpr, apply_line_data, Line)
from pbreports.report.mapping_stats import *
from pbreports.report.mapping_stats import Constants as BaseConstants
from pbreports.report.ccs import create_plot
from pbreports.io.specs import *


class Constants(BaseConstants):
    TOOL_ID = "pbreports.tasks.mapping_stats_ccs"
    DRIVER_EXE = "python -m pbreports.report.mapping_stats_ccs --resolved-tool-contract"
    R_ID = "mapping_stats_ccs"

    T_STATS = "mapping_stats_table"
    C_MOVIE = "movie"
    C_READS = "mapped_reads"
    C_NBASES = "mapped_bases"
    C_READ_LENGTH = "mapped_read_length_mean"
    C_READ_LENGTH_N50 = "mapped_read_length_n50"
    C_SUBREAD_CONCORDANCE = "mapped_read_concordance_mean"

    PG_READ_CONCORDANCE = "ccs_read_concordance_group"
    P_READ_CONCORDANCE = "concordance_plot"
    PG_READLENGTH = "ccs_readlength_group"
    PG_QV_CALIBRATION = "qv_calibration_group"
    P_QV_CALIBRATION = "qv_calibration_plot"

    A_NREADS = "mapped_reads_n"
    A_READLENGTH = "mapped_readlength_mean"
    A_READLENGTH_N50 = "mapped_readlength_n50"
    A_NBASES = "mapped_bases_n"
    A_READ_CONCORDANCE = "mapped_read_concordance_mean"

__version__ = "0.1"
log = logging.getLogger(__name__)
spec = load_spec(Constants.R_ID)


def scatter_plot_accuracy_vs_concordance(
        data,
        axis_labels=(
            get_plot_xlabel(spec, Constants.PG_QV_CALIBRATION,
                            Constants.P_QV_CALIBRATION),
            get_plot_ylabel(spec, Constants.PG_QV_CALIBRATION,
                            Constants.P_QV_CALIBRATION)),
        nbins=None,
        barcolor=None):
    accuracy, concordance = data
    fig, ax = get_fig_axes_lpr()
    data = [Line(xData=accuracy,
                 yData=concordance,
                 style='+')]
    apply_line_data(
        ax=ax,
        line_models=data,
        axis_labels=axis_labels,
        only_whole_ticks=False)
    xlim = ax.get_xlim()
    xy = np.linspace(xlim[0], xlim[1])
    ax.plot(xy, xy, '-', color='r')
    return fig, ax


class CCSMappingStatsCollector(MappingStatsCollector):
    COLUMN_ATTR = [
        Constants.A_NREADS, Constants.A_READLENGTH, Constants.A_READLENGTH_N50,
        Constants.A_NBASES, Constants.A_READ_CONCORDANCE
    ]
    HISTOGRAM_IDS = {
        Constants.P_READ_CONCORDANCE: "read_concordance_histogram",
        Constants.P_READLENGTH: "readlength_histogram",
    }

    COL_IDS = [
        Constants.C_MOVIE,
        Constants.C_READS,
        Constants.C_READ_LENGTH,
        Constants.C_READ_LENGTH_N50,
        Constants.C_NBASES,
        Constants.C_READ_CONCORDANCE
    ]

    COLUMN_AGGREGATOR_CLASSES = [
        ReadCounterAggregator,
        MeanReadLengthAggregator,
        N50Aggreggator,
        NumberSubreadBasesAggregator,
        MeanSubreadConcordanceAggregator
    ]

    def _get_plot_view_configs(self):
        """
        There are just two histogram plots in this report:

        1. Mapped concordance
        2. Read length
        """
        _p = [
            PlotViewProperties(
                Constants.P_READ_CONCORDANCE,
                Constants.PG_READ_CONCORDANCE,
                generate_plot,
                'mapped_read_concordance_histogram.png',
                xlabel=get_plot_xlabel(spec, Constants.PG_READ_CONCORDANCE,
                                       Constants.P_READ_CONCORDANCE),
                ylabel=get_plot_ylabel(spec, Constants.PG_READ_CONCORDANCE,
                                       Constants.P_READ_CONCORDANCE),
                color=get_green(3),
                edgecolor=get_green(2),
                use_group_thumb=True,
                plot_group_title=""),
            PlotViewProperties(
                Constants.P_READLENGTH,
                Constants.PG_READLENGTH,
                generate_plot,
                'mapped_readlength_histogram.png',
                xlabel=get_plot_xlabel(spec, Constants.PG_READLENGTH,
                                       Constants.P_READLENGTH),
                ylabel=get_plot_ylabel(spec, Constants.PG_READLENGTH,
                                       Constants.P_READLENGTH),
                color=get_blue(3),
                edgecolor=get_blue(2),
                use_group_thumb=True,
                plot_group_title="")
        ]
        return {v.plot_id: v for v in _p}

    def _get_total_aggregators(self):
        return OrderedDict([
            (Constants.A_READ_CONCORDANCE, MeanSubreadConcordanceAggregator()),
            (Constants.A_NREADS, ReadCounterAggregator()),
            (Constants.A_NBASES, NumberBasesAggregator()),
            (Constants.A_READLENGTH, MeanReadLengthAggregator()),
            (Constants.A_READLENGTH_N50, N50Aggreggator()),
            (Constants.A_READLENGTH_Q95, MappedReadLengthQ95(dx=10,
                                                             nbins=10000)),
            (Constants.A_READLENGTH_MAX, MaxReadLengthAggregator()),
            (self.HISTOGRAM_IDS[Constants.P_READLENGTH],
             ReadLengthHistogram()),
            (self.HISTOGRAM_IDS[Constants.P_READ_CONCORDANCE],
             SubReadConcordanceHistogram(dx=0.005, nbins=1001))
        ])

    def add_more_plots(self, plot_groups, output_dir):
        try:
            ds = ConsensusAlignmentSet(self.alignment_file)
            accuracy, concordance = [], []
            for bam in ds.resourceReaders():
                accuracy.extend(list(bam.pbi.readQual))
                concordance.extend(list(bam.identity))
            qv_validation_plot = create_plot(
                _make_plot_func=scatter_plot_accuracy_vs_concordance,
                plot_id=Constants.P_QV_CALIBRATION,
                axis_labels=(
                    get_plot_xlabel(spec, Constants.PG_QV_CALIBRATION,
                                    Constants.P_QV_CALIBRATION),
                    get_plot_ylabel(spec, Constants.PG_QV_CALIBRATION,
                                    Constants.P_QV_CALIBRATION)),
                nbins=None,
                plot_name="mapped_qv_calibration.png",
                barcolor="#A0A0FF",
                data=(accuracy, concordance),
                output_dir=output_dir)
            pg = PlotGroup(Constants.PG_QV_CALIBRATION,
                           plots=[qv_validation_plot],
                           thumbnail=qv_validation_plot.thumbnail)
            plot_groups.append(pg)
        except Exception as e:
            log.exception(e)
        return plot_groups


def to_report(alignment_file, output_dir, subreads_file=None):
    return spec.apply_view(CCSMappingStatsCollector(alignment_file).to_report(output_dir, Constants.R_ID))


def _args_runner(args):
    return run_and_write_report(args.alignment_file, args.report_json,
                                report_func=to_report)


def _resolved_tool_contract_runner(resolved_contract):
    """
    Run the mapping report from a resolved tool contract.

    :param resolved_contract:
    :type resolved_contract: ResolvedToolContract
    :return: Exit code
    """

    alignment_path = resolved_contract.task.input_files[0]

    # This is a bit different than the output dir oddness that Johann added.
    # the reports dir should be determined from the dir of the abspath of the report.json file
    # resolved values will always be absolute paths.
    output_report = resolved_contract.task.output_files[0]

    return run_and_write_report(alignment_path, output_report,
                                report_func=to_report)


def _get_parser():
    parser = get_pbparser(Constants.TOOL_ID, __version__,
                          spec.title, __doc__,
                          Constants.DRIVER_EXE)
    parser.add_input_file_type(FileTypes.DS_ALIGN_CCS, "alignment_file",
                               "ConsensusAlignment XML DataSet",
                               "BAM, SAM or ConsensusAlignment DataSet")
    parser.add_output_file_type(FileTypes.REPORT, "report_json",
                                "Mapping Statistics Report",
                                "Summary of alignment results",
                                default_name=Constants.R_ID)
    return parser


if __name__ == '__main__':
    sys.exit(main(get_parser_func=_get_parser,
                  args_runner_func=_args_runner,
                  rtc_runner_func=_resolved_tool_contract_runner))
