
# FIXME ideally we would use the .pbi as much as possible here.  however, among
# other problems, the pbi does not currently give us any way to retrieve CCS
# read length, since qStart/qEnd are always -1.

"""
Generate a report summarizing Circular Consensus Read (CCS) results.
"""

from collections import OrderedDict, Counter, defaultdict, namedtuple
import functools
import os
import os.path as op
import sys
import logging
import argparse
import time
from pprint import pformat

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import (Report, Table, Column, Attribute, Plot,
                                     PlotGroup)

from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import ConsensusReadSet, BarcodeSet

from pbreports.plot.helper import (get_fig_axes_lpr, apply_histogram_data,
                                   get_blue, get_green, Line, apply_line_data)
from pbreports.util import accuracy_as_phred_qv
from pbreports.io.specs import *

log = logging.getLogger(__name__)
__version__ = '0.44'

spec = load_spec("ccs")


class Constants(object):

    """ ids for the core Report objects (e.g., Plot, PlotGroup, etc...)"""
    TOOL_ID = "pbreports.tasks.ccs_report"
    TOOL_NAME = "ccs_report"
    DRIVER_EXE = "python -m pbreports.report.ccs --resolved-tool-contract"

    NO_BC_LABEL = "none"

    R_ID = "ccs"

    # PlotGroup
    PG_READLENGTH = 'readlength_group'
    PG_ACCURACY = "accuracy_group"
    PG_NPASSES = "npasses_hist"
    PG_SCATTER = "npasses_vs_accuracy"

    # Plots
    P_READLENGTH = "readlength_hist"
    P_ACCURACY = "accuracy_hist"
    P_NPASSES = "npasses_hist"
    P_SCATTER = "npasses_vs_accuracy"

    I_CCS_READ_LENGTH_HIST = "ccs_readlength_hist.png"
    I_CCS_READ_ACCURACY_HIST = "ccs_accuracy_hist.png"
    I_CCS_NUM_PASSES_HIST = "ccs_npasses_hist.png"
    I_CCS_SCATTER_PLOT = "ccs_npasses_vs_accuracy.png"

    # Attributes
    A_NREADS = 'number_of_ccs_reads'
    A_TOTAL_BASES = 'total_number_of_ccs_bases'
    A_MEAN_READLENGTH = 'mean_ccs_readlength'
    A_MEAN_ACCURACY = 'mean_accuracy'
    #A_MEAN_QV = 'mean_qv'
    A_MEAN_NPASSES = 'mean_ccs_num_passes'

    # Table
    T_ID = 'ccs_table'
    T_BARCODES = "ccs_barcodes"

    # Columns
    C_MOVIE_NAME = 'movie_name'
    C_NREADS = "number_of_ccs_reads"
    C_TOTAL_BASES = 'total_number_of_ccs_bases'
    C_MEAN_READLENGTH = 'ave_ccs_readlength'
    C_MEAN_ACCURACY = 'ave_ccs_accuracy'
    #C_MEAN_QV = 'mean_ccs_qv'
    C_MEAN_NPASSES = 'mean_ccs_num_passes'

    C_BARCODE_ID = "barcode_name"
    C_BARCODE_COUNTS = "number_of_ccs_reads"
    C_BARCODE_NBASES = "total_number_of_ccs_bases"
    C_BARCODE_QUALITY = "mean_ccs_accuracy"
    C_BARCODE_NPASSES = "mean_ccs_num_passes"
    C_BARCODE_READLENGTH = "mean_ccs_readlength"

MovieResult = namedtuple("MovieResult", ["movie_name", "read_lengths",
                                         "accuracies", "num_passes"])
BamStats = namedtuple("BamStats", ["qLen", "numPasses", "readScore",
                                   "movieName", "bcForward", "bcReverse"])


def _stats_from_dataset(ccs_set):
    results = []
    movie_names = set()
    for bam in ccs_set.resourceReaders():
        for rg in bam.readGroupTable:
            assert rg["ReadType"] == "CCS"
            movie_names.add(rg["MovieName"])
        is_barcoded = hasattr(bam.pbi, "bcForward")
        for k, r in enumerate(bam):
            bcForward = bam.pbi.bcForward[k] if is_barcoded else None
            bcReverse = bam.pbi.bcReverse[k] if is_barcoded else None
            results.append(BamStats(r.qLen, r.numPasses, r.readScore,
                                    r.movieName, bcForward, bcReverse))
    return results, movie_names


def _stats_to_movie_results(bam_stats, movie_names):
    """
    Separate out per-movie results from process stats.
    """
    results = []
    movies = sorted(list(movie_names))
    for movie_name in movies:
        def _base_calls():
            for r in bam_stats:
                if r.movieName == movie_name:
                    yield r.qLen

        def _num_passes():
            for r in bam_stats:
                if r.movieName == movie_name:
                    yield r.numPasses

        def _accuracy():
            for r in bam_stats:
                if r.movieName == movie_name:
                    yield r.readScore

        read_lengths = np.fromiter(_base_calls(), dtype=np.int64, count=-1)
        num_passes = np.fromiter(_num_passes(), dtype=np.int64, count=-1)
        accuracy = np.fromiter(_accuracy(), dtype=np.float, count=-1)

        results.append(MovieResult(
            movie_name, read_lengths, accuracy, num_passes))
    return results


def _movie_results_to_attributes(movie_results):
    """Create the necessary attributes for the CCS report"""
    rs = [m.read_lengths for m in movie_results]
    read_lengths = np.concatenate(rs)
    ac = [m.accuracies for m in movie_results]
    accuracies = np.concatenate(ac)
    npass = [m.num_passes for m in movie_results]
    num_passes = np.concatenate(npass)

    m_readlength = int(read_lengths.mean()) if read_lengths.size > 0 else 0
    m_accuracy = accuracies.mean() if accuracies.size > 0 else 0.0
    m_npasses = int(np.round(num_passes.mean(), decimals=0)
                    ) if num_passes.size > 0 else 0
    #m_qv = int(round(accuracy_as_phred_qv(float(m_accuracy))))

    attributes = []
    attributes.append(Attribute(Constants.A_NREADS, read_lengths.shape[0]))
    attributes.append(
        Attribute(Constants.A_TOTAL_BASES, int(read_lengths.sum())))
    attributes.append(Attribute(Constants.A_MEAN_READLENGTH, m_readlength))
    attributes.append(Attribute(Constants.A_MEAN_ACCURACY, float(m_accuracy)))
    attributes.append(Attribute(Constants.A_MEAN_NPASSES, m_npasses))

    return attributes


def _movie_results_to_table(movie_results):
    """Group movie results by movie name and build a report table.

    Table has movie name, # of CCS bases, Total CCS bases,
    mean CCS readlength and mean CCS accuracy.
    """

    columns = []
    columns.append(Column(Constants.C_MOVIE_NAME, values=[]))
    columns.append(Column(Constants.C_NREADS, values=[]))
    columns.append(Column(Constants.C_TOTAL_BASES, values=[]))
    columns.append(Column(Constants.C_MEAN_READLENGTH, values=[]))
    columns.append(Column(Constants.C_MEAN_ACCURACY, values=[]))
    columns.append(Column(Constants.C_MEAN_NPASSES, values=[]))
    table = Table(Constants.T_ID, columns=columns)

    movie_names = {m.movie_name for m in movie_results}

    for movie_name in movie_names:
        rs = [
            m.read_lengths for m in movie_results if m.movie_name == movie_name]
        read_lengths = np.concatenate(rs)
        ac = [
            m.accuracies for m in movie_results if m.movie_name == movie_name]
        accuracies = np.concatenate(ac)
        npass = [
            m.num_passes for m in movie_results if m.movie_name == movie_name]
        num_passes = np.concatenate(npass)

        m_readlength = int(
            read_lengths.mean()) if read_lengths.size > 0 else 0
        m_accuracy = accuracies.mean() if accuracies.size > 0 else 0.0
        m_npasses = int(np.round(num_passes.mean(), decimals=0)
                        ) if num_passes.size > 0 else 0
        #m_qv = int(round(accuracy_as_phred_qv(float(accuracies.mean()))))

        table.add_data_by_column_id(Constants.C_MOVIE_NAME, movie_name)
        table.add_data_by_column_id(Constants.C_NREADS, read_lengths.shape[0])
        table.add_data_by_column_id(
            Constants.C_TOTAL_BASES, int(read_lengths.sum()))
        table.add_data_by_column_id(Constants.C_MEAN_READLENGTH, m_readlength)
        table.add_data_by_column_id(Constants.C_MEAN_ACCURACY, m_accuracy)
        #table.add_data_by_column_id(Constants.A_MEAN_QV, m_qv)
        table.add_data_by_column_id(Constants.C_MEAN_NPASSES, m_npasses)

    return table


def _make_barcode_table(bam_stats, ccs_set):
    """
    Generate a table of per-barcode results
    """
    barcode_counts = defaultdict(int)
    barcode_nbases = defaultdict(int)
    barcode_npasses = defaultdict(list)
    barcode_readscores = defaultdict(list)
    is_symmetric = all([r.bcForward == r.bcReverse for r in bam_stats])
    for r in bam_stats:
        key = r.bcForward
        if not is_symmetric:
            key = (r.bcForward, r.bcReverse)
        barcode_counts[key] += 1
        barcode_nbases[key] += r.qLen
        barcode_npasses[key].append(r.numPasses)
        barcode_readscores[key].append(r.readScore)
    barcode_labels = {}
    for er in ccs_set.externalResources:
        bcs = er.barcodes
        if bcs is not None:
            with BarcodeSet(bcs) as bc_set:
                for i_bc, rec in enumerate(bc_set):
                    if i_bc in barcode_labels:
                        assert barcode_labels[i_bc] == rec.id, "Barcode ID mismatch: {l} versus {r}".format(
                            l=barcode_labels[i_bc], r=rec.id)
                    else:
                        barcode_labels[i_bc] = rec.id
    barcode_ids = sorted(barcode_counts.keys())
    counts = [barcode_counts[i_bc] for i_bc in barcode_ids]
    nbases = [barcode_nbases[i_bc] for i_bc in barcode_ids]
    mean_length = [int(float(n) / c) for (c, n) in zip(counts, nbases)]
    labels = []
    for i_bc in barcode_ids:
        if is_symmetric:
            labels.append(barcode_labels.get(i_bc, Constants.NO_BC_LABEL))
        else:
            labels.append("{f}, {r}".format(
                          f=barcode_labels.get(i_bc[0], Constants.NO_BC_LABEL),
                          r=barcode_labels.get(i_bc[1], Constants.NO_BC_LABEL)))
    npasses = [sum(barcode_npasses[i_bc]) / len(barcode_npasses[i_bc])
               for i_bc in barcode_ids]
    readquals = [sum(barcode_readscores[i_bc]) / len(barcode_readscores[i_bc])
                 for i_bc in barcode_ids]
    assert len(labels) == len(counts) == len(nbases)
    columns = [
        Column(Constants.C_BARCODE_ID, values=labels),
        Column(Constants.C_BARCODE_COUNTS, values=counts),
        Column(Constants.C_BARCODE_NBASES, values=nbases),
        Column(Constants.C_BARCODE_READLENGTH, values=mean_length),
        Column(Constants.C_BARCODE_QUALITY, values=readquals),
        Column(Constants.C_BARCODE_NPASSES, values=npasses)
    ]
    return Table(Constants.T_BARCODES, columns=columns)


def _make_histogram(data, axis_labels, nbins, barcolor):
    """Create a fig, ax instance and generate a histogram.

    :param data: np.array
    :param axis_labels: (tuple of str) (axis label, y axis label)
    :return: matplotlib fig, ax
    """
    # axis_labels = ('Median Distance Between Adapters', 'Pre-Filter Reads')
    fig, ax = get_fig_axes_lpr()
    apply_histogram_data(
        ax, data, nbins, axis_labels=axis_labels, barcolor=barcolor)
    return fig, ax


def to_cdf(points):
    _total = 0
    data = []
    for x, y in points:
        _total += int(x * y)
        data.append(_total)
    return data


def _make_histogram_with_cdf(data, axis_labels, nbins, barcolor):
    """

    """
    fig, ax = _make_histogram(data, axis_labels, nbins, barcolor)

    bins, bin_edges = np.histogram(data, bins=nbins)

    rax = ax.twinx()

    log.debug(
        "Min edges {e} bins {b}".format(e=len(bin_edges), b=len(bins)))

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


def _custom_histogram_with_cdf(new_rlabel, threshold, data, axis_labels, nbins, barcolor):
    fig, ax = _make_histogram(data, axis_labels, nbins, barcolor)

    bins, bin_edges = np.histogram(data, bins=nbins)

    rax = ax.twinx()

    log.debug(
        "Min edges {e} bins {b}".format(e=len(bin_edges), b=len(bins)))

    cdf = to_cdf(zip(bin_edges[:-1], bins))
    max_cdf = max(cdf)

    exceeded_threshold = False
    if max_cdf > threshold:
        exceeded_threshold = True
        tmp_cdf = [x / float(threshold) for x in cdf]
        cdf = tmp_cdf
        max_cdf = max(cdf)

    sdf = [max_cdf - i for i in cdf]

    log.debug((len(bin_edges), len(sdf)))

    # Plot the data
    rax.plot(bin_edges[:-1], sdf, 'k')
    rax.set_xlim(bin_edges.min(), bin_edges.max())

    if len(axis_labels) == 3:
        if exceeded_threshold:
            rax.set_ylabel(new_rlabel)
        else:
            # use the default rlabel given
            rax.set_ylabel(axis_labels[2])

    return fig, ax


def scatter_plot_accuracy_vs_numpasses(
        data,
        axis_labels=(
            get_plot_xlabel(spec, Constants.PG_SCATTER, Constants.P_SCATTER),
            get_plot_ylabel(spec, Constants.PG_SCATTER, Constants.P_SCATTER)),
        nbins=None, barcolor=None):
    """
    """
    npasses, accuracy = data
    qvs = accuracy_as_phred_qv(accuracy)
    fig, ax = get_fig_axes_lpr()
    data = [Line(xData=npasses,
                 yData=qvs,
                 style='o')]
    apply_line_data(
        ax=ax,
        line_models=data,
        axis_labels=axis_labels,
        only_whole_ticks=False)
    return fig, ax


def create_plot(_make_plot_func, plot_id, axis_labels, nbins, plot_name, barcolor, data, output_dir, dpi=72):
    """Internal function used to create Plot instances.

    This should probably have a special container class to capture all the
    plot config options.
    """

    fig, ax = _make_plot_func(data, axis_labels, nbins, barcolor)
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

    to_b = lambda x: os.path.basename(x)
    fig.savefig(os.path.join(output_dir, thumbnail), dpi=20)
    plt.close(fig)
    log.debug("Saved plot to {p}".format(p=thumbnail))
    plot = Plot(plot_id, to_b(plot_name), thumbnail=to_b(thumbnail))

    return plot

# These functions create signatures (data, axis_labels, nbins, barcolor
_custom_read_length_histogram = functools.partial(
    _custom_histogram_with_cdf, "Mb > Read Length", 50)
_custom_read_accuracy_histogram = functools.partial(
    _custom_histogram_with_cdf, "Mb > Read Score", 100)


# These functions need to generate a function with signature (data,
# output_dir, dpi=)
create_readlength_plot = functools.partial(
    create_plot, _custom_read_length_histogram, Constants.P_READLENGTH,
    (get_plot_xlabel(spec, Constants.PG_READLENGTH, Constants.P_READLENGTH),
     "Reads", "bp > Read Length"),
    80, Constants.I_CCS_READ_LENGTH_HIST, get_blue(3))

create_accuracy_plot = functools.partial(
    create_plot, _custom_read_accuracy_histogram, Constants.P_ACCURACY,
    (get_plot_xlabel(spec, Constants.PG_ACCURACY, Constants.P_ACCURACY),
     "Reads", "bp > Read Score"),
    80, Constants.I_CCS_READ_ACCURACY_HIST, get_green(3))

create_npasses_plot = functools.partial(
    create_plot, _make_histogram, Constants.P_NPASSES,
    (get_plot_xlabel(spec, Constants.PG_NPASSES, Constants.P_NPASSES),
     get_plot_ylabel(spec, Constants.PG_NPASSES, Constants.P_NPASSES)),
    80, Constants.I_CCS_NUM_PASSES_HIST, "#F18B17")

create_scatter_plot = functools.partial(
    create_plot, scatter_plot_accuracy_vs_numpasses, Constants.P_SCATTER,
    (get_plot_xlabel(spec, Constants.PG_SCATTER, Constants.P_SCATTER),
     get_plot_ylabel(spec, Constants.PG_SCATTER, Constants.P_SCATTER)),
    None, Constants.I_CCS_SCATTER_PLOT, get_blue(3))


def to_report(ccs_set, output_dir):
    bam_files = list(ccs_set.toExternalFiles())
    log.info("Generating report from files: {f}".format(f=bam_files))
    bam_stats, movie_names = _stats_from_dataset(ccs_set)
    movie_results = _stats_to_movie_results(bam_stats, movie_names)
    log.debug("\n" + pformat(movie_results))

    rs = [m.read_lengths for m in movie_results]
    readlengths = np.concatenate(rs).astype(float)
    ac = [m.accuracies for m in movie_results]
    accuracies = np.concatenate(ac)
    ps = [m.num_passes for m in movie_results]
    num_passes = np.concatenate(ps)

    readlength_plot = create_readlength_plot(readlengths, output_dir)
    accuracy_plot = create_accuracy_plot(accuracies, output_dir)
    npasses_plot = create_npasses_plot(num_passes, output_dir)
    scatter_plot = create_scatter_plot((num_passes, accuracies), output_dir)

    readlength_group = PlotGroup(Constants.PG_READLENGTH,
                                 plots=[readlength_plot],
                                 thumbnail=readlength_plot.thumbnail)
    accuracy_group = PlotGroup(Constants.PG_ACCURACY, plots=[accuracy_plot],
                               thumbnail=accuracy_plot.thumbnail)

    npasses_group = PlotGroup(Constants.PG_NPASSES, plots=[npasses_plot],
                              thumbnail=npasses_plot.thumbnail)

    scatter_group = PlotGroup(Constants.PG_SCATTER, plots=[scatter_plot],
                              thumbnail=scatter_plot.thumbnail)

    movie_table = _movie_results_to_table(movie_results)
    log.debug(str(movie_table))
    tables = [movie_table]
    if ccs_set.isBarcoded:
        tables.append(_make_barcode_table(bam_stats, ccs_set))

    attributes = _movie_results_to_attributes(movie_results)

    report = Report(Constants.R_ID,
                    tables=tables, attributes=attributes,
                    plotgroups=[readlength_group, accuracy_group,
                                npasses_group, scatter_group],
                    dataset_uuids=(ccs_set.uuid,))

    return spec.apply_view(report)


def run_report(
        input_file,
        report_json,
        output_dir):
    log.info("Running {f} v{v}.".format(
        f=os.path.basename(__file__), v=__version__))
    report = None
    ds = ConsensusReadSet(input_file)
    report = to_report(ds, output_dir)
    log.info(pformat(report.to_dict()))
    report.write_json(report_json)
    return 0


def args_runner(args):
    return run_report(
        input_file=args.ccs_in,
        report_json=args.report_json,
        output_dir=args.output_dir)


def resolved_tool_contract_runner(rtc):
    return run_report(
        input_file=rtc.task.input_files[0],
        report_json=rtc.task.output_files[0],
        output_dir=os.path.dirname(rtc.task.output_files[0]))


def get_parser():
    p = get_pbparser(
        tool_id=Constants.TOOL_ID,
        version=__version__,
        name=Constants.TOOL_NAME,
        description=__doc__,
        driver_exe=Constants.DRIVER_EXE)
    ap = p.arg_parser.parser
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_in",
                          name="ConsensusReadSet",
                          description="ConsensusRead DataSet file")
    p.add_output_file_type(FileTypes.REPORT, "report_json",
                           name="CCS Report",
                           description="Summary of results from CCS2",
                           default_name="ccs_report")
    ap.add_argument('-o', '--output-dir', dest='output_dir',
                    default=os.getcwd(),
                    help="Path to write histogram images to.")
    # ap.add_argument('--debug', action='store_true',
    #               help='Flag to debug to stdout.')
    return p


def main(argv=sys.argv):
    """Main point of Entry"""
    return pbparser_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
