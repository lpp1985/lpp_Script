
"""
Generate a report for a Iso-Seq Classify run, given both primer-trimmed,
non-chimeric, full-length reads and a classify summary.
"""
import functools
import os
import os.path as op
import sys
import logging
import argparse
from pprint import pformat

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import (Report, Table, Column, Attribute, Plot,
                                     PlotGroup)
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import ContigSet

from pbreports.plot.helper import (get_fig_axes_lpr, apply_histogram_data,
                                   get_blue)
from pbreports.util import validate_file
from pbreports.io.specs import *

log = logging.getLogger(__name__)

__version__ = '0.1.0.132615'  # The last 6 digits is changelist


class Constants(object):
    TOOL_ID = "pbreports.tasks.isoseq_classify"
    DRIVER_EXE = "python -m pbreports.report.isoseq_classify --resolved-tool-contract "
    R_ID = "isoseq_classify"

    """Names used within plot groups."""

    # PlotGroup
    PG_READLENGTH = "fulllength_nonchimeric_readlength_group"

    # Plots
    P_READLENGTH = "fulllength_nonchimeric_readlength_hist"

    # Table
    T_ATTR = "isoseq_classify_table"

spec = load_spec(Constants.R_ID)


def _report_to_attributes(summary_json):
    report = load_report_from_json(summary_json)
    return [Attribute(attr.id, attr.value) for attr in report.attributes]


def _attributes_to_table(attributes):
    """Build a report table from Iso-Seq Classify attributes.

    """
    columns = [Column(x.id) for x in attributes]

    table = Table(Constants.T_ATTR,
                  columns=columns)

    for x in attributes:
        table.add_data_by_column_id(x.id, x.value)

    return table


def _make_histogram(datum, axis_labels, nbins, barcolor):
    """Create a fig, ax instance and generate a histogram.

    :param datum: np.array
    :param axis_labels: (tuple of str) (axis label, y axis label)
    :return: matplotlib fig, ax
    """
    # axis_labels = ('Median Distance Between Adapters', 'Pre-Filter Reads')
    fig, ax = get_fig_axes_lpr()
    apply_histogram_data(ax, datum, nbins, axis_labels=axis_labels,
                         barcolor=barcolor)
    return fig, ax


def _make_histogram_with_cdf(datum, axis_labels, nbins, barcolor):
    """
    Make a histogram png file with cdf.
    """
    fig, ax = _make_histogram(datum, axis_labels, nbins, barcolor)

    bins, bin_edges = np.histogram(datum, bins=nbins)
    bin_edges = np.array(bin_edges)

    rax = ax.twinx()

    def to_cdf(points):
        """Given a list of points, return its cdf."""
        _total = 0
        datum = []
        for x, y in points:
            _total += int(x * y)
            datum.append(_total)
        return datum

    log.debug("Min edges {e} bins {b}".format(e=len(bin_edges), b=len(bins)))

    cdf = to_cdf(zip(bin_edges[:-1], bins))
    max_cdf = max(cdf)
    sdf = [max_cdf - i for i in cdf]

    log.debug((len(bin_edges), len(sdf)))

    # Plot the data
    rax.plot(bin_edges[:-1], sdf, 'k')
    rax.set_xlim(bin_edges.min(), bin_edges.max())

    if len(axis_labels) == 3:
        rax.set_ylabel(axis_labels[2])

    return fig, ax


def __create_plot(_make_plot_func, plot_id, axis_labels, nbins,
                  plot_name, barcolor, datum, output_dir, dpi=72):
    """Internal function used to create Plot instances.

    This should probably have a special container class to capture all the
    plot config options.
    """

    fig, _ax = _make_plot_func(datum, axis_labels, nbins, barcolor)
    path = os.path.join(output_dir, plot_name)
    try:
        fig.tight_layout()
    except AttributeError as e:  # FIXME bug 25872
        log.warn("figure.tight_layout() not available")
        log.warn(str(e))
    except ValueError as e:
        log.error(str(e))
    fig.savefig(path, dpi=dpi)
    log.debug("Saved plot with id {i} to {p}".format(p=path, i=plot_id))
    thumbnail = plot_name.replace(".png", "_thumb.png")

    fig.savefig(os.path.join(output_dir, thumbnail), dpi=20)
    plt.close(fig)
    log.debug("Saved plot to {p}".format(p=thumbnail))
    plot = Plot(plot_id, os.path.basename(plot_name),
                thumbnail=os.path.basename(thumbnail))

    return plot

create_readlength_plot = functools.partial(
    __create_plot, _make_histogram_with_cdf, Constants.P_READLENGTH,
    ("Read Length", "Reads", "Reads > Read Length"),
    80, "fulllength_nonchimeric_readlength_hist.png", get_blue(3))


def make_report(contig_set, summary_txt, output_dir):
    """
    Generate a report with ID, tables, attributes and plot groups.

    :param contig_set: an input FASTA file which has all full-length,
    non-chimeric reads produced by pbtranscript.py classify.

    This file is required to plot a read length histogram as part of
    the report:
         fulllength_nonchimeric_readlength_hist.png

    :param summary_txt: a summary TXT file with classify attributes,
    including 6 attributes,
        number of consensus reads
        number of five prime reads,
        number of three prime reads,
        number of poly-A reads,
        number of full-length non-chimeric reads,
        average full-length n on-chimeric read length

    Attributes of the report are extracted from this file.

    :type contig_set: str
    :type summary_txt: str
    :type output_dir: str

    :rtype: Report
    """
    log.info("Plotting read length histogram from file: {f}".
             format(f=contig_set))

    # Collect read lengths of
    def _get_reads():
        with ContigSet(contig_set) as f:
            for record in f:
                yield len(record.sequence)

    readlengths = np.fromiter(_get_reads(), dtype=np.int64, count=-1).astype(float)

    # Plot read length histogram
    readlength_plot = create_readlength_plot(readlengths, output_dir)
    readlength_group = PlotGroup(Constants.PG_READLENGTH,
                                 plots=[readlength_plot],
                                 thumbnail=readlength_plot.thumbnail)

    log.info("Plotting summary attributes from file: {f}".
             format(f=summary_txt))
    dataset_uuids = [ContigSet(contig_set).uuid]
    # Produce attributes based on summary.
    attributes = _report_to_attributes(summary_txt)
    r = load_report_from_json(summary_txt)
    # FIXME(nechols)(2016-03-22) not using the dataset UUIDs from these
    # reports; should we be?

    table = _attributes_to_table(attributes)
    log.info(str(table))

    # A report is consist of ID, tables, attributes, and plotgroups.
    report = Report(Constants.R_ID,
                    tables=[table],
                    attributes=attributes,
                    plotgroups=[readlength_group],
                    dataset_uuids=dataset_uuids)

    return spec.apply_view(report)


def _run(contig_set, summary_txt, output_dir, json_report):
    if output_dir in ["", None]:
        output_dir = os.getcwd()
    report = make_report(contig_set, summary_txt, output_dir)
    log.info(pformat(report.to_dict()))
    report.write_json(json_report)
    return 0


def args_runner(args):
    return _run(
        contig_set=args.inReadsFN,
        summary_txt=args.inSummaryFN,
        json_report=args.outJson,
        output_dir=os.path.dirname(args.outJson))


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    return _run(
        contig_set=rtc.task.input_files[0],
        summary_txt=rtc.task.input_files[1],
        json_report=rtc.task.output_files[0],
        output_dir=os.path.dirname(rtc.task.output_files[0]))

# XXX this module is *not* used by the main 'pbreports' app, so we can skip
# the separate argparse-only parser


def get_contract_parser():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        spec.title,
        __doc__,
        Constants.DRIVER_EXE)

    p.add_input_file_type(FileTypes.DS_CONTIG, "inReadsFN", "Fasta reads",
                          description="Reads in FASTA format, usually are full-length, " +
                          "non-chimeric, primer-trimmed reads produced by Iso-Seq classify.")
    p.add_input_file_type(FileTypes.JSON, "inSummaryFN", "Summary file",
                          description="A summary produced by Iso-Seq Classify, e.g. " +
                          "classify_summary.json")
    p.add_output_file_type(FileTypes.REPORT, "outJson", "Iso-Seq Classification Report",
                           description="Summary of classification results from pbtranscript",
                           default_name="isoseq_classify_report")

    return p


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
