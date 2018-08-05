
# TODO modernize this

import argparse
import logging
import os
import os.path as op
import math
import array
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mp

from pbcommand.models.report import (Attribute, Report, PlotGroup, Plot,
                                     PbReportError)
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.io import CmpH5Reader

from pbreports.io.validators import validate_dir
from pbreports.io.validators import validate_file
from pbreports.io.filtered_summary_reader import FilteredSummaryReader
from pbreports.plot.helper import (get_fig_axes_lpr,
                                   save_figure_with_thumbnail,
                                   set_tick_label_font_size,
                                   set_axis_label_font_size)
from pbreports.util import add_base_and_plot_options
from pbreports.io.specs import *


log = logging.getLogger(__name__)

__version__ = '0.1'

CSV_COLUMN_MAP = {"ReadId": ("|S128", str),
                  "Readlength": (int, int),
                  "ReadScore": (float, float),
                  "PassedFilter": (int, int)}


class Constants(object):
    R_ID = "control"
    A_CONTROL_SEQ = "control_sequence"
    A_NCONTROL = "n_control_reads"
    A_FRAC_CONTROL = "frac_control_reads"
    A_ACCURACY = "control_subread_accuracy"
    A_N50 = "control_n50"
    A_PCT95 = "control_95_percentile_readlength"
    A_LENGTH = "control_mean_readlength"

    PG_QUAL = "polymerase_read_quality"
    P_QUAL = "control_noncontrol_readquality"
    PG_LENGTH = "polymerase_read_length"
    P_LENGTH = "control_noncontrol_readlength"

spec = load_spec(Constants.R_ID)


def make_control_report(control_cmph5, filtered_subreads_csv, report,
                        output_dir, dpi, dumpdata):
    """
    Entry to report.
    :param control_cmph5: (str) path to control_reads.cmp.h5
    :param filtered_subreads_csv: (str) path to filtered_subread_summary.csv
    """
    _validate_inputs(control_cmph5, filtered_subreads_csv)
    name, control_reads = _get_control_reads(control_cmph5)
    filtered_reads = _get_filtered_reads(filtered_subreads_csv)
    control_data, sample_data = _process_reads(control_reads, filtered_reads)
    nr = _get_num_control_reads(control_data)
    if nr == 0:
        # Not sure this ever happens, but logic exists in makeControlReport.py
        r = _get_error_report()
        r.write_json(os.path.join(output_dir, report))
        return
    atts = _get_attributes(name, control_data, sample_data)
    pgs = [_get_plot_group_score(control_data,
                                 sample_data, output_dir),
           _get_plot_group_length(control_data,
                                  sample_data, output_dir)]
    r = Report(Constants.R_ID, attributes=atts, plotgroups=pgs)
    r = spec.apply_view(r)
    r.write_json(os.path.join(output_dir, report))


def _get_plot_group_length(control_data, sample_data, output_dir):
    """
    Create the quality plot group and return it.
    """
    fig = _create_length_figure(control_data, sample_data)
    fname = 'control_non-control_readlength.png'
    thumb = save_figure_with_thumbnail(fig, os.path.join(output_dir, fname))[1]
    plots = [Plot(Constants.P_LENGTH, fname)]
    pg = PlotGroup(Constants.PG_LENGTH,
                   thumbnail=os.path.basename(thumb), plots=plots)
    return pg


def _get_plot_group_score(control_data, sample_data, output_dir):
    """
    Create the length plot group and return it.
    """
    fig = _create_score_figure(control_data, sample_data)
    fname = 'control_non-control_readquality.png'
    thumb = save_figure_with_thumbnail(fig, os.path.join(output_dir, fname))[1]
    plots = [Plot(Constants.P_QUAL, fname)]
    pg = PlotGroup(Constants.PG_QUAL,
                   thumbnail=os.path.basename(thumb), plots=plots)
    return pg


def _get_attributes(name, control_data, sample_data):
    """
    Return a list of report attributes
    :param name: (:string) control contig name
    :param control_data: numpy array
    :param sample_data: numpy array
    """
    nc = _get_num_control_reads(control_data)
    ns = sample_data.shape[1]
    l = []
    l.append(Attribute(Constants.A_CONTROL_SEQ, name))
    l.append(_get_attr_n_control_reads(nc))
    l.append(_get_attr_fraction_control_reads(nc, ns))
    l.append(_get_attr_control_subread_acc(control_data))
    l.append(_get_attr_n50(control_data))
    l.append(_get_attr_control_95_readlength(control_data))
    l.append(_get_attr_control_mean_readlength(control_data))
    return l


def _get_num_control_reads(control_data):
    """
    Return the number of control reads
    """
    return control_data.shape[1]


def _validate_inputs(control_cmph5, filtered_subreads_csv):
    """
    Raise an Error if a required file is null or non-existent
    :param control_cmph5: (str) path to control_reads.cmp.h5
    :param filtered_subreads_csv: (str) path to filtered_subread_summary.csv
    """
    if control_cmph5 is None:
        raise PbReportError('control_cmph5 cannot be None')
    if not os.path.exists(control_cmph5):
        raise IOError(
            'control_cmph5 {g} does not exist: '.format(g=control_cmph5))
    if filtered_subreads_csv is None:
        raise PbReportError('filtered_subreads_csv cannot be None')
    if not os.path.exists(filtered_subreads_csv):
        raise IOError(
            'filtered_subreads_csv {g} does not exist: '.format(g=filtered_subreads_csv))


def _get_control_reads(control_cmph5):
    """
    Return a tuple of len == 2:
    Position 0: (string) control name 
    Position 1: (dict) dict of string to tuple (int,float) . The key is control readId,  
    position 0 of the tuple is accuracy, position 1 is length.
    :param control_cmph5: (str) path to control_reads.cmp.h5
    """
    control_reads = {}
    c = CmpH5Reader(control_cmph5)
    for ca in c:
        read_id = '%s/%d' % (ca.movieInfo.Name, ca.HoleNumber)
        if read_id in control_reads:
            log.warn(
                'read {i} is control read and has subreads?'.format(i=read_id))
        control_reads[read_id] = (ca.accuracy, ca.readLength)
    name = c.referenceInfo('ref000001').FullName
    return name, control_reads


def _get_filtered_reads(filtered_subreads_csv):
    """
    Return a numpy array of csv data filtered by
        PassedFilter > 0

    :param filtered_subreads_csv: path to filtered_summary.csv f
    """
    reader = FilteredSummaryReader(filtered_subreads_csv, CSV_COLUMN_MAP)
    reader.load()
    data = reader.data_as_numpy_array()
    data = data[data["PassedFilter"] > 0]
    log.info('Total # reads in {f}: {i}'.format(
        f=filtered_subreads_csv, i=reader.num_records))
    log.info('# reads that passed filter: {i}'.format(i=len(data)))
    return data


def _process_reads(control_reads, filtered_reads):
    """
    Return 2 numpy arrays, control_data and sample_data. Each array
    contains accuracy and read length metrics.
    :param control_reads: (dict) key = readId, value = accuracy & length tuple
    :param filtered_reads: (numpy array) see _get_filtered_reads
    """
    control_scores = array.array('f')
    sample_scores = array.array('f')
    control_pred_lens = array.array('i')
    control_obs_lens = array.array('i')
    sample_lens = array.array('i')
    control_acc = array.array('f')

    for row in filtered_reads:
        id, score, length = row['ReadId'],\
            row['ReadScore'],\
            row['Readlength']

        if id in control_reads:

            control_scores.append(score)
            control_pred_lens.append(length)
            control_acc.append(control_reads[id][0])
            control_obs_lens.append(control_reads[id][1])
        else:
            sample_scores.append(score)
            sample_lens.append(length)

    control_data = np.array([control_scores, control_pred_lens,
                             control_acc, control_obs_lens])
    sample_data = np.array([sample_scores, sample_lens])
    return control_data, sample_data


def _get_error_report():
    """
    Convenience function to return a report object. If
    num_control_reads is 0, returns a special report with
    a single "warning" attribute.
    """
    log.warn('Returning a report with a warning that 0 controls reads have '
             'been found.')
    a = Attribute('warning', 'No control reads found', 'Warning')
    return Report(Constants.R_ID, title=spec.title, attributes=[a])


def _create_score_figure(control_data, sample_data):
    """
    Returns a matplotlib fig object, with score histogram data applied
    """
    min_score = min(min(control_data[0, :]), min(sample_data[0, :]))
    if min_score > 0.6:
        min_score = 0.6
    x_data = np.arange(min_score, 1.0, 0.02)
    y1_data = control_data[0, :]
    y2_data = sample_data[0, :]
    xlabel = get_plot_xlabel(spec, Constants.PG_QUAL, Constants.P_QUAL)
    ylabel = get_plot_ylabel(spec, Constants.PG_QUAL, Constants.P_QUAL)
    return _apply_plot_data(x_data, y1_data, y2_data, (xlabel, ylabel), legend_loc='upper left')


def _create_length_figure(control_data, sample_data):
    """
    Returns a matplotlib fig object, with length histogram data applied
    """
    len_max = max(max(control_data[1, :]), max(sample_data[1, :]))
    num_len_bins = 30.0
    len_unit = _get_pretty_value(math.ceil(len_max / num_len_bins))
    x_data = np.arange(0, len_unit * num_len_bins, len_unit)
    y1_data = control_data[1, :]
    y2_data = sample_data[1, :]
    labels = (
        get_plot_xlabel(spec, Constants.PG_LENGTH, Constants.P_LENGTH),
        get_plot_ylabel(spec, Constants.PG_LENGTH, Constants.P_LENGTH))
    return _apply_plot_data(x_data, y1_data, y2_data, labels, legend_loc='upper right')


def _apply_plot_data(x_data, y1_data, y2_data, labels, legend_loc=None):
    """Default labels assume y1_data==control, y2_data==sample"""

    h1_color = '#5050f0'
    h2_color = '#f05050'

    # log option isn't really working yet ... funky with polygons
    # these are unintuitively inverted b/c we require the sample
    # on the left side of the doubleY axis
    h2, t = np.histogram(y1_data, bins=x_data)
    h1, t = np.histogram(y2_data, bins=x_data)

    fig, ax = get_fig_axes_lpr()
    x_data = x_data[:-1]
    y0 = np.zeros(len(x_data))
#    if log:
#        h1 = np.log10(h1)
#        h2 = np.log10(h2)

    ax.fill_between(x_data, y0, h1.T, alpha=0.6,
                    edgecolor=h1_color, facecolor=h1_color)

    fake_h2 = mp.Rectangle((0.1, 0.1), 0.1, 0.1, facecolor=h2_color,
                           edgecolor=h2_color, alpha=0.6)
    fake_h1 = mp.Rectangle((0.1, 0.1), 0.1, 0.1, facecolor=h1_color,
                           edgecolor=h1_color, alpha=0.6)

    ax.set_xlabel(labels[0])
    ax.set_ylabel('%s (Sample)' % labels[1])
    ax.legend([fake_h2, fake_h1],
              ['Control', 'Sample'], loc=legend_loc)
    # gray border around legend
    ax.get_legend().get_frame().set_edgecolor('#a0a0a0')
    set_tick_label_font_size(ax, 12, 12)
    set_axis_label_font_size(ax, 16)

    ax2 = ax.twinx()
    ax2.fill_between(x_data, y0, h2.T, alpha=0.6,
                     edgecolor=h2_color, facecolor=h2_color)
    ax2.set_ylabel('%s (Control)' % labels[1])
    set_tick_label_font_size(ax2, 12, 12)
    set_axis_label_font_size(ax2, 16)
    return fig


def _get_attr_control_95_readlength(control_data):
    """
    Return the 95% readlength
    :param control_data: numpy array
    """
    sorted = np.sort(control_data[3, :])
    index = np.ceil((95 / 100.0) * len(sorted))
    val = None
    if index >= len(sorted):
        val = sorted[-1]
    else:
        val = sorted[index]
    return Attribute(Constants.A_PCT95, int(val))


def _get_attr_control_mean_readlength(control_data):
    """
    :param control_data: numpy array
    """
    rl = np.mean(control_data[3, :])
    return Attribute(Constants.A_LENGTH, int(rl))


def _get_attr_control_subread_acc(control_data):
    """
    :param control_data: numpy array
    """
    acc = np.mean(control_data[2, :])
    return Attribute(Constants.A_ACCURACY, acc)


def _get_attr_n_control_reads(num_control):
    return Attribute(Constants.A_NCONTROL, num_control)


def _get_attr_fraction_control_reads(num_control, num_sample):
    return Attribute(Constants.A_FRAC_CONTROL, num_control / float(num_sample + num_control))


def _get_attr_n50(control_data):

    readlengths = [x for x in control_data[3]]
    r = sorted(readlengths, reverse=True)
    testn50 = sum(r) / 2.0
    testSum = 0
    n50 = 0
    for c in r:
        testSum += c
        if testn50 < testSum:
            n50 = c
            break

    return Attribute(Constants.A_N50, n50)


def _get_pretty_value(raw_value):
    """Returns a number like 1, 2, 5, 20, 50, 100, 200, 500, 1000 which is 
    closest to the input value (in log space)"""
    if raw_value <= 0.0:
        return 0
    v_log10 = math.log(raw_value) / math.log(10.0)
    iv = math.floor(v_log10)
    mantissa = v_log10 - iv
    targets = [1.0, 2.0, 5.0, 10.0]
    best_dist = 100.0
    best_target = 0
    for t in map(math.log10, targets):
        d = abs(mantissa - t)
        if d < best_dist:
            best_dist = d
            best_target = t
    return int(round(math.pow(10.0, iv + best_target)))


def args_runner(args):
    return make_control_report(args.cntrCmpH5, args.csv,
                               args.report, args.output, args.dpi, args.dumpdata)


def add_options_to_parser(p):
    desc = spec.description
    p.description = desc
    p.version = __version__
    p = add_base_and_plot_options(p)
    p.add_argument("cntrCmpH5", help="control_reads.cmp.h5",
                   type=validate_file)
    p.add_argument("csv", help="filtered_summary.csv", type=validate_file)

    p.set_defaults(func=args_runner)
    return p


def get_parser():
    # for standalone exe
    p = argparse.ArgumentParser(version=__version__)
    p = add_options_to_parser(p)
    return p


def main(argv=sys.argv[1:]):
    """Main point of Entry"""
    return pacbio_args_runner(
        argv=argv,
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == '__main__':
    sys.exit(main())
