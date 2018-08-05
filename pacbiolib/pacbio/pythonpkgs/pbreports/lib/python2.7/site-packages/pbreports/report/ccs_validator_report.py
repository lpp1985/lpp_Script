
# TODO modernize or delete

from collections import OrderedDict
import functools
import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import Table, Column, Report, Plot, PlotGroup
from pbcommand.validators import validate_dir, validate_file
from pbcommand.cli.core import pacbio_args_runner
from pbcommand.utils import setup_log
from pbcore.io.FastqIO import FastqReader

from pbreports.plot.helper import get_fig_axes_lpr
from pbreports.io.specs import *

log = logging.getLogger(__name__)
__version__ = '1.2'


class Constants(object):
    R_ID = "ccs_validator"
    PG_CCS = "ccs_validator_group"
    P_RL = "readlength_hist"
    P_QV = "qv_hist"
    T_FASTQ = "fastq_table"
    C_FN = 'file_name'
    C_NREADS = 'n_reads'
    C_TOT_BASES = 'total_bases'
    C_READLENGTH = 'mean_readlength'
    C_QV = 'mean_qv'

spec = load_spec(Constants.R_ID)


class FastqStats(object):

    def __init__(self, reads, qvs, file_name):
        """Simple container class"""
        self.qvs = qvs
        # these are read lengths
        self.reads = reads
        self.file_name = file_name

    @staticmethod
    def from_file(file_name):
        qvs, reads = _get_stats(file_name)
        return FastqStats(reads, qvs, file_name)

    def __str__(self):
        outs = list()
        outs.append("Reads           :{n}".format(n=self.reads.shape[0]))
        outs.append("Mean readlength :{m}".format(int(np.sum(self.reads))))
        outs.append("Total bases     :{m}".format(m=int(np.sum(self.reads))))
        outs.append("Mean qv         :{m:.2f}".format(m=self.qvs.mean()))
        return "\n".join(outs)


def _get_stats(fastq_file_name):
    raw_qvs = np.array([r.quality for r in FastqReader(fastq_file_name)])
    qvs = np.hstack(raw_qvs)
    reads = np.array([len(r.sequence) for r in FastqReader(fastq_file_name)])
    return qvs, reads


def _generate_histogram(datum, title, xlabel, ylabel=None):
    fig, ax = get_fig_axes_lpr()
    fig.suptitle(title)
    ax.hist(datum)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    return fig, ax


def __generate_histogram_comparison(method_name, title, xlabel, list_fastq_stats):
    fig, ax = get_fig_axes_lpr()
    fig.suptitle(title)

    alpha = 0.3
    hs = OrderedDict()
    for fastq_stat in list_fastq_stats:
        label = os.path.basename(fastq_stat.file_name)
        h = ax.hist(getattr(fastq_stat, method_name),
                    alpha=alpha, bins=85, label=label)
        hs[label] = h

    ax.set_xlabel(xlabel)
    ax.legend(loc="best")
    return fig, ax

to_qv_histogram = functools.partial(
    __generate_histogram_comparison, 'qvs',
    get_plot_title(spec, Constants.PG_CCS, Constants.P_QV),
    get_plot_xlabel(spec, Constants.PG_CCS, Constants.P_QV))
to_read_length_histogram = functools.partial(
    __generate_histogram_comparison, 'reads',
    get_plot_title(spec, Constants.PG_CCS, Constants.P_RL),
    get_plot_xlabel(spec, Constants.PG_CCS, Constants.P_RL))


def _generate_table(list_fastq_stats):
    columns = [Column(Constants.C_FN),
               Column(Constants.C_NREADS),
               Column(Constants.C_TOT_BASES),
               Column(Constants.C_READLENGTH),
               Column(Constants.C_QV)]

    table = Table(Constants.T_FASTQ, columns=columns)

    for fastq_stat in list_fastq_stats:
        table.add_data_by_column_id(
            Constants.C_FN, os.path.basename(fastq_stat.file_name))
        table.add_data_by_column_id(
            Constants.C_NREADS, fastq_stat.reads.shape[0])
        table.add_data_by_column_id(
            Constants.C_TOT_BASES, int(np.sum(fastq_stat.reads)))
        table.add_data_by_column_id(
            Constants.C_READLENGTH, int(fastq_stat.reads.mean()))
        table.add_data_by_column_id(Constants.C_QV, fastq_stat.qvs.mean())

    return table


def fastq_files_to_stats(fastq_files):
    fastq_stats = {file_name: FastqStats.from_file(
        file_name) for file_name in fastq_files}
    return fastq_stats


def to_report(fastq_files, qv_hist=None, readlength_hist=None):
    """Generate a histogram of read lengths and quality values"""
    fastq_stats = fastq_files_to_stats(fastq_files)

    table = _generate_table(fastq_stats.values())
    log.debug(str(table))

    if qv_hist is not None:
        fig, ax = to_qv_histogram(fastq_stats.values())
        fig.savefig(qv_hist)
    if readlength_hist is not None:
        fig, ax = to_read_length_histogram(fastq_stats.values())
        fig.savefig(readlength_hist)
    plt.close(fig)
    readlength_hist_plot = Plot(Constants.P_RL, readlength_hist)
    plotgroup = PlotGroup(Constants.PG_RL, plots=[
                          readlength_hist_plot])
    report = Report(Constats.R_ID, tables=[table], plotgroups=[plotgroup])
    return spec.apply_view(report)


def args_runner(args):
    output_dir = args.output_dir
    json_report_name = args.report

    to_p = lambda x: os.path.join(output_dir, x)
    json_report = to_p(json_report_name)
    readlength_hist = to_p(
        'ccs_validation_readlength_histogram.png') if output_dir else None
    qv_hist = to_p('ccs_validation_qv_histogram.png') if output_dir else None

    log.info("Starting v{v} of {f}".format(v=__version__,
                                           f=os.path.basename(__file__)))
    fastq_files = [args.fastq_1, args.fastq_2]

    # weak attempt to make the plots labels show up consistently
    fastq_files.sort()
    report = to_report(fastq_files, qv_hist=qv_hist,
                       readlength_hist=readlength_hist)

    log.info("writing report to {j}".format(j=json_report))
    report.write_json(json_report)

    return 0


def get_parser():
    p = argparse.ArgumentParser(version=__version__)
    p.add_argument('fastq_1', type=validate_file)
    p.add_argument('fastq_2', type=validate_file)
    p.add_argument('--output-dir', required=True, type=validate_dir,
                   dest='output_dir',
                   help="Directory to write Read length and Quality Value histograms.")
    p.add_argument('-r', '--report', type=str,
                   default="ccs_validator_report.json",
                   help=spec.title)
    p.add_argument('--debug', action='store_true', help='Debug to stdout.')
    return p


def main(argv=sys.argv):
    """Main point of Entry"""
    log.info("Starting {f} version {v} report generation".format(
        f=__file__, v=__version__))
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == '__main__':
    sys.exit(main())
