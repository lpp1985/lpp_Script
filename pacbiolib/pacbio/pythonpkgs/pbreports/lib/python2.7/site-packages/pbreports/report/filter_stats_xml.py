
# TODO(nechols)(2016-04-22) rename this to raw_data_report or similar

"""
Generates:
 - Filter ReadLength histograms with SDF (with thumbnails)
 - Filter ReadScore Histogram with SDF (with thumbnails)
 - Table of total bases, # of reads, mean readlengh, mean
   readscore
"""

from collections import OrderedDict
import functools
import logging
import os
import os.path as op
import sys

import numpy as np

from pbcommand.utils import setup_log
from pbcommand.models.report import (Report, Table, Column, Plot, Attribute,
                                     PlotGroup)
from pbcommand.models import TaskTypes, FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcore.io import SubreadSet

from pbreports.plot.helper import (get_fig_axes_lpr,
                                   save_figure_with_thumbnail, get_green)
from pbreports.util import compute_n50, continuous_dist_shaper
from pbreports.io.specs import *

__version__ = '0.1.0'


class Constants(object):
    TOOL_ID = "pbreports.tasks.filter_stats_report_xml"
    DRIVER_EXE = ("python -m pbreports.report.filter_stats_xml "
                  "--resolved-tool-contract ")
    R_ID = "raw_data_report"

    A_NBASES = "nbases"
    A_NREADS = "nreads"
    A_READ_N50 = "read_n50"
    A_READ_LENGTH = "read_length"
#    A_READ_QUALITY = "read_quality"
    A_INSERT_LENGTH = "insert_length"
#    A_INSERT_QUALITY = "insert_quality"

    READ_ATTR = [A_NBASES, A_NREADS, A_READ_LENGTH, A_READ_N50]
    INSERT_ATTR = [A_INSERT_LENGTH]
    PG_RL = "read_length_plot_group"
    P_RL = "read_length_plot"
    PG_IL = "insert_length_plot_group"
    P_IL = "insert_length_plot"


class ReadStatsPlots(object):
    P_LENGTH = "read_length_plot"
    P_LENGTH_PREFIX = "readLenDist"
    P_QUAL = "read_quality_plot"
    P_QUAL_PREFIX = "readQualDist"
    PG_LENGTH = "read_length_plot_group"
    PG_QUAL = "read_quality_plot_group"


class InsertStatsPlots(object):
    P_LENGTH_PREFIX = "insertLenDist"
    P_LENGTH = "insert_length_plot"
    P_QUAL = "insert_quality_plot"
    P_QUAL_PREFIX = "insertQualDist"
    PG_LENGTH = "insert_length_plot_group"
    PG_QUAL = "insert_quality_plot_group"


log = logging.getLogger(__name__)
spec = load_spec(Constants.R_ID)


def _total_from_bins(bins, min_val, bin_width):
    _min = min_val
    _wid = bin_width
    bin_means = [_min + (_wid * i) + _wid / 2 for i in range(len(bins))]
    bin_totals = [count * mean for count, mean in zip(bins, bin_means)]
    return sum(bin_totals)


def _to_read_stats_attributes(readLenDists, readQualDists):
    # Build the stats table:
    nbases = 0
    nreads = 0
    n50 = 0
    readscoretotal = 0
    readscorenumber = 0
    approx_read_lens = []

    # if a merge failed there may be more than one dist:
    for rlendist in readLenDists:
        nbases += rlendist.sampleMean * rlendist.sampleSize
        nreads += rlendist.sampleSize

        # N50:
        for i, lbin in enumerate(rlendist.bins):
            # use the average, except for the last bin
            if i != len(rlendist.bins) - 1:
                value = ((i * rlendist.binWidth) + rlendist.minBinValue +
                         rlendist.binWidth / 2)
            # for the last bin, just use the value
            else:
                value = (i * rlendist.binWidth) + rlendist.minBinValue
            approx_read_lens.extend([value] * lbin)
            # TODO(mdsmith)(2016-02-09) make sure maxOutlierValue is updated
            # during a merge /todo
            # but pop off that last value and replace it with the
            # maxOutlierValue:
            # approx_read_lens.pop()
            # approx_read_lens.append(rlendist.maxBinValue)
    n50 = int(np.round(compute_n50(approx_read_lens), decimals=0))
    for rqualdist in readQualDists:
        readscoretotal += _total_from_bins(rqualdist.bins,
                                           rqualdist.minBinValue,
                                           rqualdist.binWidth)
        readscorenumber += sum(rqualdist.bins)

    readlen = 0
    if nreads != 0:
        readlen = nbases / nreads
    readlen = int(np.round(readlen, decimals=0))
    readQuality = 0
    if readscorenumber != 0:
        readQuality = readscoretotal / readscorenumber
    return [int(np.round(nbases, decimals=0)),
            nreads,
            readlen,
            n50,
            ]  # readQuality]


def _make_attributes(read_attr, attr_values):
    assert len(read_attr) == len(attr_values)
    return [Attribute(attr_id, value=value)
            for attr_id, value in zip(read_attr, attr_values)]


def to_read_stats_attributes(readLenDists, readQualDists):
    return _make_attributes(
        read_attr=Constants.READ_ATTR,
        attr_values=_to_read_stats_attributes(readLenDists, readQualDists))


# XXX just Insert Length!
def to_insert_stats_attributes(readLenDists, readQualDists):
    return _make_attributes(
        read_attr=Constants.INSERT_ATTR,
        attr_values=[_to_read_stats_attributes(readLenDists, readQualDists)[2]])


def _to_read_stats_plots(PlotConstants, title, readLenDists, readQualDists,
                         output_dir, dpi=72, lenDistShaper=None):
    length_plots = []
    # ReadLen distribution to barplot:
    if lenDistShaper is None:
        lenDistShaper = continuous_dist_shaper(readLenDists, trim_excess=True)
    for i, orig_rlendist in enumerate(readLenDists):
        rlendist = lenDistShaper(orig_rlendist)
        assert sum(orig_rlendist.bins) == sum(rlendist.bins)
        len_fig, len_axes = get_fig_axes_lpr()
        len_axes.bar(rlendist.labels, rlendist.bins,
                     color=get_green(0), edgecolor=get_green(0),
                     width=(rlendist.binWidth * 0.75))
        len_axes.set_xlabel(get_plot_xlabel(spec, PlotConstants.PG_LENGTH,
                                            PlotConstants.P_LENGTH))
        len_axes.set_ylabel(get_plot_ylabel(spec, PlotConstants.PG_LENGTH,
                                            PlotConstants.P_LENGTH))
        png_fn = os.path.join(output_dir, "{p}{i}.png".format(i=i,
                                                              p=PlotConstants.P_LENGTH_PREFIX))
        png_base, thumbnail_base = save_figure_with_thumbnail(len_fig, png_fn,
                                                              dpi=dpi)
        length_plots.append(
            Plot("{p}_{i}".format(i=i, p=PlotConstants.P_LENGTH),
                 os.path.relpath(png_base, output_dir),
                 title=title, caption=title,
                 thumbnail=os.path.relpath(thumbnail_base, output_dir)))

    plot_group_thumb = length_plots[-1].thumbnail if length_plots else None

    plot_groups = [
        PlotGroup(PlotConstants.PG_LENGTH,
                  title=title,
                  plots=length_plots,
                  thumbnail=plot_group_thumb)
    ]
    return plot_groups

    # FIXME these aren't useful yet
    qual_plots = []
    # ReadQual distribution to barplot:
    shaper = continuous_dist_shaper(readQualDists, trim_excess=True)
    for i, orig_rqualdist in enumerate(readQualDists):
        rqualdist = shaper(orig_rqualdist)
        qual_fig, qual_axes = get_fig_axes_lpr()
        qual_axes.bar(rqualdist.labels, rqualdist.bins,
                      color=get_green(0), edgecolor=get_green(0),
                      width=(rqualdist.binWidth * 0.75))
        qual_axes.set_xlabel(get_plot_xlabel(spec, PlotConstants.PG_QUAL,
                                             PlotConstants.P_QUAL))
        qual_axes.set_ylabel(get_plot_ylabel(spec, PlotConstants.PG_QUAL,
                                             PlotConstants.P_QUAL))
        png_fn = os.path.join(output_dir, "{p}{i}.png".format(i=i,
                                                              p=PlotConstants.P_QUAL_PREFIX))
        png_base, thumbnail_base = save_figure_with_thumbnail(qual_fig, png_fn,
                                                              dpi=dpi)
        qual_plots.append(
            Plot("{p}_{i}".format(i=i, p=PlotConstants.P_QUAL),
                 os.path.relpath(png_base, output_dir),
                 thumbnail=os.path.relpath(thumbnail_base, output_dir)))
    plot_groups.append(
        PlotGroup(PlotConstants.PG_QUAL,
                  plots=qual_plots))
    return plot_groups

to_read_stats_plots = functools.partial(_to_read_stats_plots, ReadStatsPlots,
                                        get_plotgroup_title(spec, Constants.PG_RL))
to_insert_stats_plots = functools.partial(_to_read_stats_plots, InsertStatsPlots,
                                          get_plotgroup_title(spec, Constants.PG_IL))


def to_report(stats_xml, output_dir, dpi=72):
    """Main point of entry

    :type stats_xml: str
    :type output_dir: str
    :type dpi: int

    :rtype: Report
    """
    log.info("Analyzing XML {f}".format(f=stats_xml))
    # stats_xml should be a dataset:
    dset = SubreadSet(stats_xml)

    dataset_uuids = [dset.uuid]
    # but if it isn't, no problem:
    if not dset.metadata.summaryStats:
        dset.loadStats(stats_xml)
        # an sts file was provided which will generate a new random uuid
        dataset_uuids = []
    if not dset.metadata.summaryStats.readLenDists:
        raise IOError("Pipeline Summary Stats (sts.xml) not found or missing "
                      "key distributions")

    # we want all of the length distributions in this report to look the same,
    # so we make the shaper here and pass it around:
    alldists = (dset.metadata.summaryStats.readLenDists[:] +
                dset.metadata.summaryStats.insertReadLenDists[:])
    len_dist_shaper = continuous_dist_shaper(alldists, trim_excess=True)

    attr = to_read_stats_attributes(
        readLenDists=dset.metadata.summaryStats.readLenDists,
        readQualDists=dset.metadata.summaryStats.readQualDists)
    attr.extend(to_insert_stats_attributes(
        readLenDists=dset.metadata.summaryStats.insertReadLenDists,
        readQualDists=dset.metadata.summaryStats.insertReadQualDists))

    plot_groups = to_read_stats_plots(
        readLenDists=dset.metadata.summaryStats.readLenDists,
        readQualDists=dset.metadata.summaryStats.readQualDists,
        output_dir=output_dir,
        lenDistShaper=len_dist_shaper)
    plot_groups.extend(to_insert_stats_plots(
        readLenDists=dset.metadata.summaryStats.insertReadLenDists,
        readQualDists=dset.metadata.summaryStats.insertReadQualDists,
        output_dir=output_dir,
        lenDistShaper=len_dist_shaper))

    # build the report:
    report = Report(Constants.R_ID,
                    attributes=attr,
                    plotgroups=plot_groups,
                    dataset_uuids=dataset_uuids)

    return spec.apply_view(report)


def args_runner(args):
    log.info("Starting {f} v{v}".format(f=os.path.basename(__file__),
                                        v=__version__))
    output_dir = os.path.dirname(args.report)
    report = to_report(args.subread_set, output_dir)
    report.write_json(args.report)
    return 0


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    log.info("Starting {f} v{v}".format(f=os.path.basename(__file__),
                                        v=__version__))
    output_dir = os.path.dirname(rtc.task.output_files[0])
    report = to_report(rtc.task.input_files[0], output_dir)
    report.write_json(rtc.task.output_files[0])
    return 0


def _add_options_to_parser(p):
    p.add_input_file_type(
        FileTypes.DS_SUBREADS,
        file_id="subread_set",
        name="SubreadSet",
        description="SubreadSet")
    p.add_output_file_type(
        FileTypes.REPORT, "report", spec.title,
        description=("Filename of JSON output report. Should be name only, "
                     "and will be written to output dir"),
        default_name="report")


def add_options_to_parser(p):
    """
    API function for extending main pbreport arg parser (independently of
    tool contract interface).
    """
    p_wrap = _get_parser_core()
    p_wrap.arg_parser.parser = p
    p.description = __doc__
    add_debug_option(p)
    _add_options_to_parser(p_wrap)
    p.set_defaults(func=args_runner)
    return p


def _get_parser_core():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        "Raw Dataset Statistics XML Report",
        __doc__,
        Constants.DRIVER_EXE,
        is_distributed=True)
    return p


def get_parser():
    p = _get_parser_core()
    _add_options_to_parser(p)
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)


# for 'python -m pbreports.report.sat ...'
if __name__ == "__main__":
    sys.exit(main())
