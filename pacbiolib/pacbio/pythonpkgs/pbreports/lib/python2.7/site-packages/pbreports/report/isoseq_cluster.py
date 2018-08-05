"""
Generate a report for a Iso-Seq Cluster run, given both
consensus isoforms reads in Fasta file and a
cluster summary.
"""

from pprint import pformat
import functools
import os
import sys
import logging

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import (Report, Table, Column, Plot,
                                     PlotGroup)
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import ContigSet

from pbreports.plot.helper import (get_fig_axes_lpr, apply_histogram_data,
                                   get_blue)
from pbreports.io.specs import *

log = logging.getLogger(__name__)

__version__ = '0.1.0.132615'  # The last 6 digits is changelist


class Constants(object):
    TOOL_ID = "pbreports.tasks.isoseq_cluster"
    DRIVER_EXE = "python -m pbreports.report.isoseq_cluster --resolved-tool-contract"
    R_ID = "isoseq_cluster"

    # Attributes
    A_LENGTH = "avg_consensus_isoform_length"
    A_CONSENSUS = "num_consensus_isoforms"
    A_BASES = "num_total_bases"

    # PlotGroup
    PG_READLENGTH = "consensus_isoforms_readlength_group"
    PG_AVGQV = "hq_lq_isoforms_avgqv_group"

    # Plots
    P_READLENGTH = "consensus_isoforms_readlength_hist"
    P_AVGQV = "hq_lq_isoforms_avgqv_hist"

    # Table
    T_ATTR = "isoseq_classify_table"

spec = load_spec(Constants.R_ID)


def _report_to_attributes(report_file):
    report = load_report_from_json(report_file)
    return report.attributes


def attributesToTable(attributes):
    """Build a report table from Iso-Seq cluster attributes."""
    columns = [Column(x.id, header="") for x in attributes]

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
    ("Read Length", "Reads", "Reads > Read Length"), 80,
    "consensus_isoforms_readlength_hist.png", get_blue(3))

create_avgqv_plot = functools.partial(
    __create_plot, _make_histogram_with_cdf, Constants.P_AVGQV,
    ("HQ LQ Isoform Average QV", "Isoforms", "Isoforms > Average QV"), 80,
    "hq_lq_isoforms_avgqv_hist.png", get_blue(3))


def makeReport(inReadsFN, hq_isoforms_fq, lq_isoforms_fq, inSummaryFN, outDir):
    """
    Generate a report with ID, tables, attributes and plot groups.

    inReadsFN --- an input FASTA file which has all consensus
    isoforms produced by pbtranscript.py cluster.
    This file is required to plot a read length histogram as part of
    the report:
         consensus_isoforms_readlength_hist.png

    hq_isoforms_fq/lq_isoforms_lq --- input FASTQ files which has
    all HQ/LQ isoforms produced by pbtranscript.py cluster.
    These two files will be required to plot the average QV histograms:
         hq_lq_isoforms_avgqv_hist.png

    inSummaryFN --- a summary TXT file with cluster attributes,
    including two attributes:
         number of consensus isoforms
         average length of consensus isoforms
    Attributes of the report are extracted from this file.

    """
    log.info("Plotting read length histogram from file: {f}".
             format(f=inReadsFN))

    # Collect read lengths of
    reader = ContigSet(inReadsFN)
    rs = [len(r.sequence) for r in reader]
    reader.close()
    readlengths = np.array(rs).astype(float)

    # Plot read length histogram
    readlength_plot = create_readlength_plot(readlengths, outDir)
    readlength_group = PlotGroup(Constants.PG_READLENGTH,
                                 plots=[readlength_plot],
                                 thumbnail=readlength_plot.thumbnail)

    # Collect average qvs
    hq_qvs = [np.mean(r.quality) for r in ContigSet(hq_isoforms_fq)]
    lq_qvs = [np.mean(r.quality) for r in ContigSet(lq_isoforms_fq)]
    avgqvs = np.array(hq_qvs + lq_qvs)

    # Plot average qv histogram
    avgqv_plot = create_avgqv_plot(avgqvs, outDir)
    avgqv_group = PlotGroup(Constants.PG_AVGQV,
                            plots=[avgqv_plot],
                            thumbnail=avgqv_plot.thumbnail)

    log.info("Plotting summary attributes from file: {f}".
             format(f=inSummaryFN))
    # Produce attributes based on summary.
    dataset_uuids = [ContigSet(inReadsFN).uuid]
    attributes = _report_to_attributes(inSummaryFN)
    r = load_report_from_json(inSummaryFN)
    # FIXME(nechols)(2016-03-22) not using the dataset UUIDs from these
    # reports; should we be?

    table = attributesToTable(attributes)
    log.info(str(table))

    # A report is consist of ID, tables, attributes, and plotgroups.
    report = Report(Constants.R_ID,
                    attributes=attributes,
                    plotgroups=[readlength_group, avgqv_group],
                    dataset_uuids=dataset_uuids)

    return spec.apply_view(report)


def _run(fasta_file, hq_isoforms_fq, lq_isoforms_fq, summary_txt, json_report, outDir):
    if outDir in ["", None]:
        outDir = os.getcwd()
    report = makeReport(
        inReadsFN=fasta_file,
        hq_isoforms_fq=hq_isoforms_fq,
        lq_isoforms_fq=lq_isoforms_fq,
        inSummaryFN=summary_txt,
        outDir=outDir)
    log.info(pformat(report.to_dict()))
    report.write_json(json_report)
    return 0


def args_runner(args):
    return _run(
        fasta_file=args.inReadsFN,
        hq_isoforms_fq=args.hq_isoforms_fq,
        lq_isoforms_fq=args.lq_isoforms_fq,
        summary_txt=args.inSummaryFN,
        json_report=args.outJson,
        outDir=os.path.dirname(args.outJson))


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    return _run(
        fasta_file=rtc.task.input_files[0],
        hq_isoforms_fq=rtc.task.input_files[2],
        lq_isoforms_fq=rtc.task.input_files[3],
        summary_txt=rtc.task.input_files[1],
        json_report=rtc.task.output_files[0],
        outDir=os.path.dirname(rtc.task.output_files[0]))


def get_contract_parser():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        "Iso-Seq Cluster Report",
        __doc__,
        Constants.DRIVER_EXE,
        is_distributed=True)

    p.add_input_file_type(FileTypes.DS_CONTIG, "inReadsFN", "Fasta reads",
                          description="Reads in FASTA format, usually are consensus, " +
                          "isoforms produced by Iso-Seq Cluster.")

    p.add_input_file_type(FileTypes.DS_CONTIG, "hq_isoforms_fq", "HQ isoforms in Fastq",
                          description="HQ isoforms in FASTQ format produced by Iso-Seq Cluster.")

    p.add_input_file_type(FileTypes.DS_CONTIG, "lq_isoforms_fq", "LQ isoforms in Fastq",
                          description="LQ isoforms in FASTQ format produced by Iso-Seq Cluster.")

    p.add_input_file_type(FileTypes.JSON, "inSummaryFN", "Summary text",
                          description="A summary produced by Iso-Seq Cluster, e.g. " +
                          "cluster_summary.txt")
    p.add_output_file_type(FileTypes.REPORT, "outJson", "Transcript Clustering Report",
                           description="Summary of results from pbtranscript",
                           default_name="isoseq_cluster_report")
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
