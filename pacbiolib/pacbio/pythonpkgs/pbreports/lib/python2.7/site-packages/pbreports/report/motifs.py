
"""
Create Kinetics report from motifs.gff.gz and motif_summary.csv
"""


import os
import os.path as op
from pprint import pformat
import sys
import gzip
import csv
import logging
import collections
import argparse
import operator

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import Report, PlotGroup, Plot, Table, Column
from pbcommand.models import TaskTypes, FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcommand.utils import setup_log
from pbcore.io.GffIO import GffReader

import pbreports.plot.helper as PH
from pbreports.io.specs import *

__version__ = '2.0'

log = logging.getLogger()

spec = load_spec("motifs")


class Constants(object):
    # XXX TCI stuff
    TOOL_ID = "pbreports.tasks.motifs_report"
    DRIVER_EXE = "python -m pbreports.report.motifs --resolved-tool-contract"
    MAX_MOTIFS_ID = "pbreports.task_options.max_motifs_plot"
    MAX_MOTIFS_DEFAULT = 10

    R_ID = spec.id

    T_ID = "motif_records"
    C_ID = 'motif_id'
    C_POS = 'modified_position'
    C_TYPE = 'modification_type'
    C_PCT_MOTIF = 'percent_motifs_detected'
    C_NMOTIF = 'ndetected_motifs'
    C_NMOTIF_GEN = 'nmotifs_in_genome'
    C_READSCORE = 'mean_readscore'
    C_COV = 'mean_coverage'
    C_PARTNER = 'partner_motif'
    C_IPD = 'mean_ipd_ratio'
    C_GRP = 'group_tag'
    C_OBJ_SCORE = 'objective_score'

    # Plot Groups
    PG_MOD_QV = 'modification_qvs'
    P_MOD_QV = "motifs"

    PG_MOD = 'modifications'
    P_MOD_COV = 'mod_qv_coverage'
    P_MOD_HIST = 'qmod_hist'

    # Images
    I_MOTIFS_QMOD = 'motif_histogram.png'
    I_MOTIFS_QMOD_THUMB = 'motif_histogram_thumb.png'

    I_MOD_HISTOGRAM = 'motif_histogram.png'
    I_KINETICS_SCATTER = 'kinetics_detections.png'
    I_KINETICS_HIST = "kinetics_histogram.png"


class MotifRecord(object):

    def __init__(self, motif_str, center_position, modification_type, fraction,
                 ndetected, genome, group_tag, partner_motif_str, mean_score,
                 mean_ipd_ratio, mean_coverage, objective_score):
        self.motif_str = motif_str
        self.center_position = center_position
        self.modification_type = modification_type
        self.fraction = fraction
        self.ndetected = ndetected
        self.ngenome = genome
        self.group_tag = group_tag
        self.partner_motif_str = partner_motif_str
        self.mean_score = mean_score
        self.mean_ipd_ratio = mean_ipd_ratio
        self.mean_coverage = mean_coverage
        self.objective_score = objective_score

    @staticmethod
    def from_dict(d):
        """Create a new record from the raw CSV header dict"""
        f = lambda x: operator.getitem(d, x)
        r = MotifRecord(f('motifString'), f('centerPos'), f('modificationType'),
                        f('fraction'), f('nDetected'), f(
                            'nGenome'), f('groupTag'),
                        f('partnerMotifString'), f(
                            'meanScore'), f('meanIpdRatio'),
                        f('meanCoverage'), f('objectiveScore'))
        return r

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  m=self.motif_str,
                  t=self.modification_type,
                  s=self.mean_score,
                  o=self.objective_score,
                  c=self.mean_coverage)
        return "<{k} {m} {t} score:{s} coverage:{c} object score:{s} >".format(**_d)


def _motif_summary_type_d():
    d = {}
    d["motifString"] = str
    d["centerPos"] = int
    d["modificationType"] = str
    d["fraction"] = float
    d["nDetected"] = int

    d["nGenome"] = int
    d["groupTag"] = str

    d["partnerMotifString"] = str
    d["meanScore"] = float
    d["meanIpdRatio"] = float
    d["meanCoverage"] = float
    d["objectiveScore"] = float

    return d


def _motif_csv_to_records(csv_file):

    records = []

    with open(csv_file, 'r') as f:
        header = f.readline()
        column_type_d = _motif_summary_type_d()
        # the CSV is written with every item wrapped in "Stuff","More Sfuff",
        column_names = [x.replace('"', '') for x in header.strip().split(',')]
        log.info("Found {n} columns".format(n=len(column_names)))

        for line in f:
            items = line.strip().split(',')
            if len(items) != len(column_names):
                _d = dict(n=len(column_names), x=len(items), l=line)
                raise ValueError(
                    "Expected {n} items, found {x}. Unable to process line {l}".format(**_d))

            values_d = {}
            for column_name, value in zip(column_names, items):
                column_type_func = column_type_d[column_name]
                v = column_type_func(value.replace('"', ''))
                values_d[column_name] = v

            if 'modificationType' in values_d:
                if values_d['modificationType'] == 'modified_base':
                    values_d['modificationType'] = 'unknown'

            record = MotifRecord.from_dict(values_d)
            records.append(record)

    log.info("Found {n} motif records from {f}".format(
        n=len(records), f=csv_file))
    return records


def setTickLabelFontSize(ax, minorTickLabelSize, majorTickLabelSize):
    """Convenience function for changing font size of major and minor ticks"""
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(majorTickLabelSize)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(minorTickLabelSize)


def setAxisLabelFontSize(ax, labelSize):
    """Convenience function for changing font size of x- and y-axis labels."""
    t = ax.get_xaxis().get_label()
    t.set_fontsize(labelSize)
    t = ax.get_yaxis().get_label()
    t.set_fontsize(labelSize)


def _createFigTemplate(dims=(8, 6), facecolor='#ffffff', gridcolor='#e0e0e0'):
    fig = plt.figure(figsize=dims)
    ax = fig.add_subplot(111)
    ax.axesPatch.set_facecolor(facecolor)
    ax.grid(color=gridcolor, linewidth=0.5, linestyle='-')
    ax.set_axisbelow(True)

    fig.patch.set_alpha(0.5)
    ax.patch.set_alpha(0.5)

    setTickLabelFontSize(ax, 12, 12)
    setAxisLabelFontSize(ax, 16)
    return (fig, ax)


def readModificationCsvGz(fn):

    gzFile = gzip.GzipFile(fn)
    reader = csv.reader(gzFile)

    records = []
    header = []
    header = reader.next()

    colIdx = 0
    colMap = {}
    for h in header:
        colMap[h] = colIdx
        colIdx += 1

    # Read csv
    n = 0
    kinHit = collections.namedtuple("kinHit", "base coverage score")
    for row, record in enumerate(reader):
        if int(record[colMap['score']]) > 20:
            tupleRec = kinHit(base=record[colMap['base']], coverage=int(
                record[colMap['coverage']]), score=int(record[colMap['score']]))
            records.append(tupleRec)
            n += 1

    # convert to recarray
    kinRec = [('base', '|S1'), ('coverage', '>i4'),
              ('score', '>i4'), ('color', 'b')]
    kinArr = np.zeros(len(records), dtype=kinRec)
    idx = 0
    for rec in records:
        kinArr['base'][idx] = rec.base
        kinArr['coverage'][idx] = rec.coverage
        kinArr['score'][idx] = rec.score
        idx += 1

    return kinArr


def plotKineticsScatter(kinArr, outputFileName):

    handles = []
    colors = ['red', 'green', 'blue', 'magenta']
    bases = ['A', 'C', 'G', 'T']

    fig, ax = _createFigTemplate(dims=(10, 8))

    for i in xrange(4):
        baseHits = kinArr[kinArr['base'] == bases[i]]

        if baseHits.shape[0] > 0:
            # Add a bit of scatter to avoid ugly aliasing in plot due to
            # integer quantization
            cov = baseHits['coverage'] + 0.25 * \
                np.random.randn(baseHits.shape[0])
            score = baseHits['score'] + 0.25 * \
                np.random.randn(baseHits.shape[0])

            pl = ax.scatter(cov, score, c=colors[i], label=bases[
                            i], lw=0, alpha=0.3, s=12)
            handles.append(pl)

    ax.set_xlabel(get_plot_xlabel(spec, Constants.PG_MOD, Constants.P_MOD_COV))
    ax.set_ylabel(get_plot_ylabel(spec, Constants.PG_MOD, Constants.P_MOD_COV))
    plt.legend(handles, bases, loc='upper left')

    if kinArr.shape[0] > 0:
        ax.set_xlim(0, np.percentile(kinArr['coverage'], 95.0) * 1.4)
        ax.set_ylim(0, np.percentile(kinArr['score'], 99.9) * 1.3)

    fig.savefig(outputFileName, dpi=72)
    plt.close(fig)


def plotKineticsHist(kinArr, outputFileName):

    colors = ['red', 'green', 'blue', 'magenta']
    bases = ['A', 'C', 'G', 'T']

    fig, ax = _createFigTemplate(dims=(10, 8))

    # Check for empty or peculiar modifications report:
    d = kinArr['score']
    if d.size == 0:
        binLim = 1.0
    elif np.isnan(np.sum(d)):
        binLim = np.nanmax(d)
    else:
        binLim = np.percentile(d, 99.9) * 1.2

    ax.set_xlim(0, binLim)
    bins = np.arange(0, binLim, step=binLim / 75)

    for i in xrange(4):
        baseHits = kinArr[kinArr['base'] == bases[i]]
        if baseHits.shape[0] > 0:
            _ = ax.hist(baseHits['score'], color=colors[i], label=bases[
                        i], bins=bins, histtype="step", log=True)

    ax.set_xlabel(get_plot_xlabel(
        spec, Constants.PG_MOD, Constants.P_MOD_HIST))
    ax.set_ylabel(get_plot_ylabel(
        spec, Constants.PG_MOD, Constants.P_MOD_HIST))

    if d.size > 0:
        ax.legend(loc='upper right')

    fig.savefig(outputFileName, dpi=72)
    plt.close(fig)


def addQmodPlot(kinData, outputFolder):

    #chartPng = "kineticDetections.png"
    chartPng = Constants.I_KINETICS_SCATTER

    # Generate modification detection plot
    plotKineticsScatter(kinData, os.path.join(outputFolder, chartPng))

    # Put plot into report
    p = Plot(Constants.P_MOD_COV, image=chartPng)
    #graph = GraphItem()
    #graph.title = 'Modification QV vs. Coverage'
    # graph.addImage(chartPng)
    return p


def addQmodHist(kinData, outputFolder):

    #chartPng = "kineticHistogram.png"
    chartPng = Constants.I_KINETICS_HIST

    image_name = os.path.join(outputFolder, chartPng)
    # Generate modification detection plot
    plotKineticsHist(kinData, image_name)

    # Put plot into report
    #graph = GraphItem()
    p = Plot(Coverage.P_MOD_HIST, image=chartPng)

    #graph.title = 'Modification QVs'
    # graph.addImage(chartPng)
    return p


# The following methods generate a motif histogram

def readMotifFiles(gffFile):

    # Read in the existing motifs.gff
    motReader = GffReader(gffFile)
    y = [{"score": x.score, "attributes": x.attributes} for x in motReader]

    # Convert to a recarray
    idx = 0
    kinRec = [('motif', 'a25'), ('score', '>i4')]
    kinArr = np.zeros(len(y), dtype=kinRec)
    for d in y:
        u = d["attributes"]
        kinArr['score'][idx] = d["score"]
        if "motif" in u:
            kinArr['motif'][idx] = u["motif"]
        else:
            kinArr['motif'][idx] = 'Not Clustered'
        idx += 1

    return kinArr


# used by excludeSparseRegions to locate sparse regions in histogram

def findNextSparseRegion(hist, numBins):

    # find first non-empty bin in histogram:
    for start in range(numBins, -1, -1):
        if hist[start] > 0:
            break

    # find the distance to the next nonzero:
    for i in range(start - 1, -1, -1):
        if hist[i] > 0:
            break

    dist = start - i
    start = start - 1
    return start, dist


# find an upper limit for the x-axis that excludes sparse regions

def excludeSparseRegions(data):

    # Try to catch empty motifs:
    if data.size == 0:
        return 1

    maxBins = int(np.max(data))

    if data.size < 10:
        return maxBins

    # If there are at least five ten points, try to identify possible
    # outlier(s):

    # compute histogram
    hist, binEdges = np.histogram(data, bins=maxBins)

    # create a dictionary of sparse regions in histogram
    d = {}
    ell = len(hist) - 1
    while ell > 1:
        start, dist = findNextSparseRegion(hist, ell)
        d[start] = dist
        ell = ell - dist

    # find the first big gap that we can safely exclude
    foundOKcutoff = False
    F = float(np.sum(hist))
    for key in sorted(d.iterkeys()):

        if d[key] >= 10:
            start = key
            # don't exclude more than top 10% of the data
            s = np.sum(hist[(start + 1):]) / F
            if s < 0.1:
                foundOKcutoff = True
                break

    if not foundOKcutoff:
        start = maxBins

    return start


def plotMotifHist(csvFile, kinArr, max_motifs=10):
    # Use kinArr to determine number of motifs:
    # motifs = np.unique( kinArr['motif'] )

    # Use motif_summary.csv to determine number of motifs
    motifs = []

    with open(csvFile, 'r') as f:
        k = 0
        reader = csv.reader(f, delimiter=',')
        reader.next()
        for row in reader:
            motifs.append(row[0])
            k += 1
            if k >= max_motifs:
                break

    # Check to make sure there exists a 'Not Clustered' site in kinArr:
    if 'Not Clustered' in kinArr['motif']:
        motifs.append('Not Clustered')

    numMotifs = len(motifs)

    # Generate a unique color (RGBA tuple) for each motif
    colors = []
    cm = plt.get_cmap('hsv')
    for i in range(numMotifs - 1):
        colors.append(cm(1. * i / (numMotifs - 1)))
    colors.append('0.75')

    # Try to find an acceptable QV upper bound for each motif and take the
    # maximum of those
    binLim = 1
    for i in xrange(numMotifs):
        baseHits = kinArr[kinArr['motif'] == motifs[i]]
        # Try to locate sparse regions in the histogram for exclusion:
        b = excludeSparseRegions(baseHits['score'])
        binLim = max(binLim, b) + 1

    # Try integer bin boundaries to avoid empty bins:
    bins = range(0, binLim, 1)

    # Plot all histograms on the same plot for now:
    fig, ax = _createFigTemplate(dims=(10, 8))
    ax.set_xlim(0, binLim)

    for i in xrange(numMotifs):
        baseHits = kinArr[kinArr['motif'] == motifs[i]]
        if baseHits.shape[0] > 0:
            pl = ax.hist(baseHits['score'], color=colors[i], label=motifs[
                         i], bins=bins, histtype="step", log=True)

    ax.set_xlabel(get_plot_xlabel(
        spec, Constants.PG_MOD_QV, Constants.P_MOD_QV))
    ax.set_ylabel(get_plot_ylabel(
        spec, Constants.PG_MOD_QV, Constants.P_MOD_QV))

    # Display a legend only if at least one motif was found:
    if numMotifs > 0:
        ax.legend(loc='upper right', fontsize=9)

    return fig, ax


def addQmodMotifHist(csvFile, kinData, outputFolder, dpi=72, max_motifs=10):

    # Apart from passing in motif_summary.csv file name, nearly identical to
    # addQmodHist

    image_name = os.path.join(outputFolder, Constants.I_MOTIFS_QMOD)

    # Generate modification detection plot
    fig, ax = plotMotifHist(csvFile, kinData, max_motifs=max_motifs)

    png, thumbpng = PH.save_figure_with_thumbnail(fig, image_name, dpi=dpi)

    log.info((png, thumbpng))

    plot = Plot(Constants.P_MOD_QV, image=os.path.basename(png),
                thumbnail=os.path.basename(thumbpng))

    pg = PlotGroup(Constants.PG_MOD_QV,
                   title=get_plotgroup_title(spec, Constants.PG_MOD_QV),
                   plots=[plot],
                   thumbnail=os.path.basename(thumbpng))
    return pg


def to_table(motif_records):

    columns = [Column(Constants.C_ID),
               Column(Constants.C_POS),
               Column(Constants.C_TYPE),
               Column(Constants.C_PCT_MOTIF),
               Column(Constants.C_NMOTIF),
               Column(Constants.C_NMOTIF_GEN),
               Column(Constants.C_READSCORE),
               Column(Constants.C_COV),
               Column(Constants.C_PARTNER),
               Column(Constants.C_IPD),
               Column(Constants.C_GRP),
               Column(Constants.C_OBJ_SCORE)]

    # Record attr name ordered by index in columns
    attr_names = ['motif_str', 'center_position',
                  'modification_type',
                  'fraction', 'ndetected',
                  'ngenome', 'mean_score',
                  'mean_coverage', 'partner_motif_str',
                  'mean_ipd_ratio',
                  'group_tag', 'objective_score']

    table = Table(Constants.T_ID, columns=columns)

    for record in motif_records:
        for attr_name, column in zip(attr_names, columns):
            v = getattr(record, attr_name)
            table.add_data_by_column_id(column.id, v)

    return table


def to_motifs_report(gff_file, motif_summary_csv, output_dir, max_motifs=10):

    _d = dict(g=gff_file, c=motif_summary_csv, o=output_dir)
    log.info(
        "starting Motif report generations with: \nGFF:{g}\nCSV:{c}\ndir:{o}".format(**_d))

    # Generate a histogram with lines corresponding to motifs
    kinData = readMotifFiles(gff_file)
    plot_group = addQmodMotifHist(motif_summary_csv, kinData, output_dir,
                                  max_motifs=max_motifs)
    plot_groups = [plot_group]

    motif_records = _motif_csv_to_records(motif_summary_csv)
    table = to_table(motif_records)

    r = Report(Constants.R_ID,
               plotgroups=plot_groups,
               tables=[table])
    log.debug(pformat(r.to_dict(), indent=4))
    return spec.apply_view(r)


def to_mod_report(motif_summary_csv, output_dir):

    # Set up the modifications report
    #report = GraphReportItem()
    #report.title = 'Modifications'
    #graphGroup = GraphGroupItem(title ='Kinetic Detections')

    kinData = readModificationCsvGz(motif_summary_csv)

    p1 = addQmodPlot(kinData, output_dir)
    p2 = addQmodHist(kinData, output_dir)
    plots = [p1, p2]

    pg = PlotGroup(Constants.PG_MOD,
                   title=get_plotgroup_title(spec, Constants.PG_MOD),
                   plots=plots)

    r = Report(Constants.R_ID, plotgroups=[pg])

    return spec.apply_view(r)


def _write_report(r, json_file):
    with open(json_file, 'w') as f:
        f.write(r.to_json())
    log.info("Write report {i} to {f}".format(i=r.id, f=json_file))
    return 0


def _to_motif_report(
        gff_file,
        motif_summary_csv,
        output,
        report_json,
        max_motifs=10):
    report = to_motifs_report(gff_file, motif_summary_csv, output,
                              max_motifs=max_motifs)
    # the json report supplied as the basename,
    return _write_report(report, os.path.join(output, report_json))


def _to_mod_report(args):
    log.debug(args)
    report = to_mod_report(args.motif_summary_csv, args.output)
    return _write_report(report, args.report_json)


def args_runner(args):
    return _to_motif_report(
        gff_file=args.gff_file,
        motif_summary_csv=args.motif_summary_csv,
        output=os.path.dirname(args.report_json),
        report_json=os.path.basename(args.report_json),
        max_motifs=args.maxMotifs)


def resolved_tool_contract_runner(resolved_tool_contract):
    report_json = resolved_tool_contract.task.output_files[0]
    return _to_motif_report(
        gff_file=resolved_tool_contract.task.input_files[0],
        motif_summary_csv=resolved_tool_contract.task.input_files[1],
        output=os.path.dirname(report_json),
        report_json=os.path.basename(report_json),
        max_motifs=resolved_tool_contract.task.options[Constants.MAX_MOTIFS_ID])


def get_parser():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        spec.title,
        __doc__,
        Constants.DRIVER_EXE)

    p.add_input_file_type(FileTypes.GFF, 'gff_file',
                          "GFF file", "Path to motifs.gff.gz")
    p.add_input_file_type(FileTypes.CSV, 'motif_summary_csv',
                          "CSV file", 'Path to Motif summary CSV')
    p.add_output_file_type(FileTypes.REPORT, 'report_json',
                           name="Motifs Report",
                           description="Summary of motif analysis results",
                           default_name="motifs_report")
    p.add_int(Constants.MAX_MOTIFS_ID,
              option_str="maxMotifs",
              default=Constants.MAX_MOTIFS_DEFAULT,
              name="Max. number of motifs in QV plot",
              description="Control number of motifs whose QVs are plotted")
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
