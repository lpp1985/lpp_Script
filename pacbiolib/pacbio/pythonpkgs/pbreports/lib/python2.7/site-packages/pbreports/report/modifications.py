
# FIXME we might want to move this to kineticsTools someday, just to keep the
# HDF5 dependency contained (but as long as it's in pbcore it doesn't matter)

"""
Generates plots showing the distribution of kinetics across all bases, taken
from ipdSummary output.
"""

import collections
import argparse
import logging
import gzip
import os
import os.path as op
import sys

from pylab import legend, arange
import numpy as np
import h5py

from pbcommand.models.report import Report, PlotGroup, Plot
from pbcommand.models import TaskTypes, FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcommand.utils import setup_log

import pbreports.plot.helper as PH
from pbreports.util import (add_base_and_plot_options,
                            add_base_options_pbcommand)
from pbreports.util import Constants as BaseConstants
from pbreports.io.specs import *

log = logging.getLogger(__name__)

__version__ = '2.1'

spec = load_spec("modifications")


class Constants(BaseConstants):
    TOOL_ID = "pbreports.tasks.modifications_report"
    DRIVER_EXE = "python -m pbreports.report.modifications --resolved-tool-contract"
    PG_KIN = "kinetic_detections"
    P_SCAT = "kinetic_scatter"
    P_HIST = "kinetic_histogram"


def _create_fig_template(dims=(8, 6), facecolor='#ffffff', gridcolor='#e0e0e0'):
    fig, ax = PH.get_fig_axes_lpr(dims=dims)
    ax = fig.add_subplot(111)

    ax.axesPatch.set_facecolor(facecolor)
    ax.grid(color=gridcolor, linewidth=0.5, linestyle='-')
    ax.set_axisbelow(True)
    PH.set_tick_label_font_size(ax, 12, 12)
    PH.set_axis_label_font_size(ax, 16)
    return fig, ax


def _get_datasets(basemods_h5, ds_id):
    log.debug("extracting {d}...".format(d=ds_id))
    data = []
    for ref_name in basemods_h5.keys():
        data.append(basemods_h5[ref_name][ds_id])
    log.debug("  concatenating arrays...")
    return np.concatenate(data)


def plot_kinetics_scatter(basemods_h5, ax):

    handles = []
    colors = ['red', 'green', 'blue', 'magenta']
    base_ids = ['A', 'C', 'G', 'T']

    coverage = _get_datasets(basemods_h5, 'coverage')
    score = _get_datasets(basemods_h5, 'score')
    bases = _get_datasets(basemods_h5, 'base').__array__()

    for base, color in zip(base_ids, colors):
        baseHits = bases == base
        n_bases = np.count_nonzero(baseHits)

        if n_bases > 0:
            # Add a bit of scatter to avoid ugly aliasing in plot due to
            # integer quantization
            coverage_ = coverage[baseHits] + 0.25 * np.random.randn(n_bases)
            score_ = score[baseHits] + 0.25 * np.random.randn(n_bases)

            pl = ax.scatter(coverage_, score_, c=color, label=base,
                            lw=0, alpha=0.3, s=12)
        else:
            log.warn("Base {b} not found".format(b=base))

    ax.set_xlabel(get_plot_xlabel(spec, Constants.PG_KIN, Constants.P_SCAT))
    ax.set_ylabel(get_plot_ylabel(spec, Constants.PG_KIN, Constants.P_SCAT))
    ax.legend(loc='upper left')

    if len(coverage) > 0:
        ax.set_xlim(0, np.percentile(coverage, 95.0) * 1.4)
        ax.set_ylim(0, np.percentile(score, 99.9) * 1.3)


def plot_kinetics_hist(basemods_h5, ax):

    colors = ['red', 'green', 'blue', 'magenta']
    base_ids = ['A', 'C', 'G', 'T']

    # Check for empty or peculiar modifications report:
    bases = _get_datasets(basemods_h5, 'base')
    scores = _get_datasets(basemods_h5, 'score')
    if len(scores) == 0:
        binLim = 1.0
    elif np.isnan(np.sum(scores)):
        binLim = np.nanmax(scores)
    else:
        binLim = np.percentile(scores, 99.9) * 1.2
    log.debug("binLim = {l}".format(l=binLim))
    ax.set_xlim(0, binLim)
    bins = arange(0, binLim, step=binLim / 75)

    for base, color in zip(base_ids, colors):
        baseHits = bases == base
        if np.count_nonzero(baseHits) > 0:
            pl = ax.hist(scores[baseHits], color=color, label=base,
                         bins=bins, histtype="step", log=True)
        else:
            log.warn("Base {b} not found".format(b=base))

    ax.set_xlabel(get_plot_xlabel(spec, Constants.PG_KIN, Constants.P_HIST))
    ax.set_ylabel(get_plot_ylabel(spec, Constants.PG_KIN, Constants.P_HIST))

    if len(scores) > 0:
        ax.legend(loc='upper right')


def get_qmod_plot(basemods_h5, output_dir, dpi):
    """
    Return a plot object
    """
    fig, ax = _create_fig_template()

    plot_kinetics_scatter(basemods_h5, ax)

    png_path = os.path.join(output_dir, "kinetic_detections.png")
    png, thumbpng = PH.save_figure_with_thumbnail(fig, png_path, dpi=dpi)

    return Plot(Constants.P_SCAT, os.path.basename(png),
                thumbnail=os.path.basename(thumbpng))


def get_qmod_hist(basemods_h5, output_dir, dpi):
    """
    Return a plot object
    """
    fig, ax = _create_fig_template()

    plot_kinetics_hist(basemods_h5, ax)

    png_path = os.path.join(output_dir, "kinetic_histogram.png")
    png, thumbpng = PH.save_figure_with_thumbnail(fig, png_path, dpi=dpi)

    return Plot(Constants.P_HIST, os.path.basename(png),
                thumbnail=os.path.basename(thumbpng))


def make_modifications_report(modifications_h5, report, output_dir, dpi=72):
    """
    Entry point to report generation.
    """
    basemods_h5 = h5py.File(modifications_h5)
    scatter = get_qmod_plot(basemods_h5, output_dir, dpi)
    hist = get_qmod_hist(basemods_h5, output_dir, dpi)
    pg = PlotGroup(Constants.PG_KIN,
                   title=get_plotgroup_title(spec, Constants.PG_KIN),
                   thumbnail=scatter.thumbnail,
                   plots=[scatter, hist])
    rpt = Report(spec.id, plotgroups=[pg])
    rpt = spec.apply_view(rpt)
    rpt.write_json(os.path.join(output_dir, report))
    return 0


def add_options_to_parser(p):
    from pbreports.io.validators import validate_file
    p.description = __doc__  # FIXME which is probably wrong
    p.version = __version__
    p = add_base_and_plot_options(p)
    p.add_argument("basemods_h5", help="modifications.h5", type=validate_file)
    p.set_defaults(func=_args_runner)
    return p


def args_runner(args):
    return make_modifications_report(
        modifications_h5=args.basemods_h5,
        report=os.path.basename(args.report),
        output_dir=os.path.dirname(args.report))


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    return make_modifications_report(
        modifications_h5=rtc.task.input_files[0],
        report=os.path.basename(rtc.task.output_files[0]),
        output_dir=os.path.dirname(rtc.task.output_files[0]))


def get_parser():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        "Modifications Report",
        __doc__,
        Constants.DRIVER_EXE,
        is_distributed=True)
    p.add_input_file_type(FileTypes.H5, "basemods_h5", "HDF5 file",
                          "HDF5 file of base modifications from ipdSummary")
    p.add_output_file_type(FileTypes.REPORT, "report", "Basemods report",
                           description="Summary of basemod results",
                           default_name="report")
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)


if __name__ == "__main__":
    sys.exit(main())
