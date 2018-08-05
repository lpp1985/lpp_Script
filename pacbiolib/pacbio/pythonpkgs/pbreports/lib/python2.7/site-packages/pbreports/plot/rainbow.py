
"""
Convert an input bam or DataSet XML file (used to be a cmp.h5 file, thus the
name) to a figure of Concordance vs. Subread Length. Each point on the graph
represents the concordance and length of a single subread as measured by a local
alignment to the reference. The points are colored by qv-score.
"""

import argparse
import logging
import time
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcore.io import openDataFile, DataSet, CmpH5Reader
from pbcommand.models.report import Report, PlotGroup, Plot

from pbreports.plot.helper import save_figure_with_thumbnail
from pbreports.io.validators import (validate_file,
                                     validate_output_dir)

log = logging.getLogger(__name__)


class Constants(object):
    PLOT_GROUP_ID = "alignment_to_png_plot_group"
    CAPTION_STRING = """\
Each point on the graph represents the concordance and length of a single subread
as measured by a local alignment to the reference.  The points are colored by
Z-score, a measure of the significance of each alignment."""


def _read_in_file(in_fn, reference=None):
    """ Read in a file, compute statistics per reference or for a particular
    reference. MapQV coloring used to be z_score coloring.

    Args:
        in_fn: The name of an indexed BAM file, DataSet XML or cmp.h5 file
        reference: the reference to compute statistics for (all references if
                   None)

    Returns:
        A 2D array of lengths, percent concordance and color by MapQV
    """
    if in_fn.endswith(".xml"):
        return _read_in_indexed_alignmentset(in_fn, reference)

    def _openAlignments():
        if in_fn.endswith(".cmp.h5"):
            return CmpH5Reader(in_fn)
        else:
            return openDataFile(in_fn)
    lengths, percent_accs, map_qvs = [], [], []
    with _openAlignments() as alignments:
        for row in alignments:
            if reference == None or row.referenceName == reference:
                try:
                    length = row.aEnd - row.aStart
                except (AttributeError, IndexError):
                    length = row.rEnd - row.rStart
                lengths.append(length)
                # if bam, breaks if cmp.h5:
                try:
                    n_ins = row.aEnd - row.aStart - row.nM - row.nMM
                    n_del = row.tEnd - row.tStart - row.nM - row.nMM
                    percent_accs.append(1.0 - (row.nMM + n_ins + n_del) /
                                        float(length))
                except (AttributeError, IndexError):
                    percent_accs.append(1.0 - (row.nMM + row.nIns + row.nDel) /
                                        float(length))
                map_qvs.append(float(row.MapQV))

    data = np.array([lengths, percent_accs, map_qvs])
    data = data.transpose()
    return data


def _read_in_indexed_alignmentset(in_fn, reference=None):
    """
    Extract data from the .pbi files in an AlignmentSet using numpy array
    operations.
    """
    lengths, percent_accs, map_qvs = [], [], []
    with openDataFile(in_fn) as ds:
        for bam in ds.resourceReaders():
            if len(bam) == 0:
                continue
            identities = bam.identity
            ref_name_to_id = {r.Name: r.ID for r in bam.referenceInfoTable}
            sel = np.full(len(identities), True, dtype=bool)
            bam_lengths = bam.pbi.aEnd - bam.pbi.aStart
            if reference is not None:
                ref_id = None
                # FIXME there must be a cleaner way to do this...
                for ref_info in bam.referenceInfoTable:
                    if ref_info.Name == reference:
                        ref_id = ref_info.ID
                        break
                sel = bam.pbi.tId == ref_id
            lengths.extend(bam_lengths[sel])
            percent_accs.extend(identities[sel])
            map_qvs.extend(bam.pbi.mapQV[sel])
    data = np.array([lengths, percent_accs, map_qvs])
    data = data.transpose()
    return data


def _make_plot(data, png_fn, bounds=None, dpi=60, nolegend=False):
    """Make a scatterplot of read length and concordance"""
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.axesPatch.set_facecolor('#ffffff')
    axes.grid(color='#e0e0e0', linewidth=0.5, linestyle='-')
    axes.set_axisbelow(True)

    # from color brewer
    # qv_colors = ['#a6cee3', '#1f77b4', '#b2df8a', '#33a02c', '#fb9a99',
    #'#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
    #'#ffff99']
    qv_colors = ['#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d']
    # qv_colors.extend(qv_colors)
    # qv_colors.extend(qv_colors)
    # plot by z-values
    qv_min = 1.0
    #qv_delta = 3.0
    handles = []
    labels = []
    # Make sure the max actually gets in a bin
    qv_max = max(data[:, 2]) + 1
    qv_delta = (qv_max - qv_min) / len(qv_colors)
    for qv_bin, color in zip(
            #np.arange(qv_min, qv_min + qv_delta * len(qv_colors), qv_delta),
            np.arange(qv_min, qv_max, qv_delta),
            qv_colors):
        if qv_bin > qv_max:
            break
        qv_bin_max = qv_bin + qv_delta
        points = data[(data[:, 2] >= qv_bin) * (data[:, 2] < qv_bin_max), :]
        if len(points[:, 0]) > 0:
            l, = axes.plot(points[:, 0], points[:, 1], 'o', c=color, mec=color,
                           alpha=0.1, ms=2.0)
            handles.append(l)
            labels.append('QV >= %d' % qv_bin)
    if not nolegend:
        axes.legend(handles, labels, loc='lower right', numpoints=1,
                    borderpad=0.3, markerscale=2.0, handletextpad=0.3,
                    labelspacing=0.3, handlelength=0.5)
        axes.get_legend().get_frame().set_edgecolor('#a0a0a0')

    if bounds:
        intbounds = map(int, bounds.split(":"))
        axes.set_xlim(xmin=intbounds[0], xmax=intbounds[1])
        axes.set_ylim(ymin=intbounds[2], ymax=intbounds[3])
    axes.set_xlabel('Subread Length / bp')
    axes.set_ylabel('% Concordance')
    save_figure_with_thumbnail(fig, png_fn, dpi=int(dpi))
    plt.close(fig)


def make_report(in_fn, out_dir='.', bounds=None, nolegend=False,
                reference=None, dpi=60, name=None):
    """AlignmentToPng Report

    Convert an input bam or DataSet XML file to a figure of Concordance vs.
    Subread Length.

    Args:
        in_fn: the bam, DataSet XML or cmp.h5 file to turn into a length vs
               concordance plot
        out_dir: the output directory to be used with the file name or default
        name: the file name to be used with the outdir or default (no full
              path filenames!)
        bounds: the figure limits (in xmin:xmax:ymin:ymax)
        nolegend: exclude the figure legend
        reference: the reference to use in the figure. Default of all
                   references
        dpi: the dots per inch (resolution) of the figure
    """

    data = _read_in_file(in_fn, reference)
    report = Report('alignment_to_png_report')

    if not name:
        name = '%s.png' % os.path.splitext(os.path.basename(in_fn))[0]
    png_fn = os.path.join(out_dir, name)
    _make_plot(data, png_fn, bounds, dpi, nolegend)
    plot_group = PlotGroup(Constants.PLOT_GROUP_ID,
                           plots=[Plot('alignment_to_png_plot',
                                       os.path.basename(png_fn))])
    report.add_plotgroup(plot_group)
    return report


def make_rainbow_plot(in_fn, png_name, reference=None):
    t1 = time.time()
    data = _read_in_file(in_fn, reference)
    _make_plot(data, png_name)
    t2 = time.time()
    log.info("Plot generated in {s:.2f} sec".format(s=t2 - t1))
