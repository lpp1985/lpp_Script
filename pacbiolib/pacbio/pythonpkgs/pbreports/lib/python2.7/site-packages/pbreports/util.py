
import functools
import argparse
import logging
import glob
import time
import math
import os
import sys

import numpy as np

from copy import deepcopy

from pbcore.util.Process import backticks
from pbcore.io import FastaReader, ReferenceSet
from pbcommand.models import FileTypes
from pbcommand import common_options


log = logging.getLogger(__name__)


class Constants(object):
    DPI_ID = "pbreports.task_options.dpi"
    DUMPDATA_ID = "pbreports.task_options.dumpdata"
    LOG_10 = math.log(10)


def setup_log(alog, file_name=None, level=logging.DEBUG, str_formatter=None):
    """
    Util function for setting up logging.

    Due to how smrtpipe logs, the default behavior is that the stdout
    is where the logging is redirected. If a file name is given the log
    will be written to that file.

    :param log: (log instance) Log instance that handlers and filters will
    be added.
    :param file_name: (str, None), Path to file. If None, stdout will be used.
    :param level: (int) logging level
    """
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)

    if str_formatter is None:
        str_formatter = '[%(levelname)s] %(asctime)-15s [%(name)s %(funcName)s %(lineno)d] %(message)s'

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    alog.addHandler(handler)
    alog.setLevel(level)


def log_timing(func):
    """Simple deco to log the runtime of func"""
    started_at = time.time()

    def wrapper(*args, **kw):
        return func(*args, **kw)

    run_time = time.time() - started_at
    name = func.__name__
    log.info("Func {f} took {s:.2f} sec ({m:.2f} min)".format(
        f=name, s=run_time, m=run_time / 60.0))

    return wrapper


def _nfs_exists_check(ff):
    """
    Central place for all NFS hackery

    Return whether a file or a dir ff exists or not.
    Call ls instead of python os.path.exists to eliminate NFS errors.

    Added try/catch black hole exception cases to help trigger an NFS refresh

    :rtype bool:
    """
    # try to trigger refresh for File case
    try:
        f = open(ff, 'r')
        f.close()
    except Exception:
        pass

    # try to trigger refresh for Directory case
    try:
        _ = os.stat(ff)
        _ = os.listdir(ff)
    except Exception:
        pass

    # Call externally
    # this is taken from Yuan
    cmd = "ls %s" % ff
    _, rcode, _ = backticks(cmd)

    return rcode == 0


def _validate_resource(func, resource):
    """Validate the existence of a file/dir"""
    # Attempt to trigger an NFS metadata refresh
    _ = _nfs_exists_check(resource)

    if func(resource):
        return os.path.abspath(resource)
    else:
        raise IOError("Unable to find {f}".format(f=resource))

validate_file = functools.partial(_validate_resource, os.path.isfile)
validate_dir = functools.partial(_validate_resource, os.path.isdir)
validate_output_dir = functools.partial(_validate_resource, os.path.isdir)


def validate_nonempty_file(resource):
    try:
        resource_path = validate_file(resource)
    except IOError as e:
        raise e
    try:
        with open(resource_path) as handle:
            l = [handle.next() for i in range(2)]
    except StopIteration:
        raise IOError("{f} appears to be empty".format(f=resource))
    return resource_path


def validate_report(report):
    """
    Raise ValueError if report contains path seps
    """
    if not os.path.basename(report) == report:
        raise ValueError(
            "Path separators are not allowed: {r}".format(r=report))
    return report


def validate_fofn(fofn):
    """Validate existence of FOFN and files within the FOFN.

    :param fofn: (str) Path to File of file names.
    :raises: IOError if any file is not found.
    :return: (str) abspath of the input fofn
    """
    _ = _nfs_exists_check(fofn)

    if os.path.isfile(fofn):
        file_names = bas_fofn_to_bas_files(os.path.abspath(fofn))
        log.debug("Found {n} files in FOFN {f}.".format(
            n=len(file_names), f=fofn))
        return os.path.abspath(fofn)
    else:
        raise IOError("Unable to find {f}".format(f=fofn))


def fofn_to_files(fofn):
    """Util func to convert a bas/bax fofn file to a list of bas/bax files."""

    _ = _nfs_exists_check(fofn)

    if os.path.exists(fofn):
        with open(fofn, 'r') as f:
            bas_files = {line.strip() for line in f.readlines()}

        for bas_file in bas_files:
            if not os.path.isfile(bas_file):
                # try one more time to find the file by
                # performing an NFS refresh
                found = _nfs_exists_check(bas_file)
                if not found:
                    raise IOError(
                        "Unable to find bas/bax file '{f}'".format(f=bas_file))

        return list(bas_files)
    else:
        raise IOError("Unable to find FOFN {f}".format(f=fofn))

# For backward compatibility
bas_fofn_to_bas_files = fofn_to_files


def movie_to_cell(movie):
    """
    Parse the cell name from the movie
    """
    # TODO This will need to evolve with movienames
    try:
        return '_'.join(os.path.basename(movie).split('_')[:4])
    except IndexError:
        return movie


def path_to_movie(path):
    return path.split(os.path.sep)[-1].split(".")[0]


def get_fasta_readlengths(fasta_file):
    """
    Get a sorted list of contig lengths
    :return: (tuple) 
    """
    lens = []
    with FastaReader(fasta_file) as f:
        for record in f:
            lens.append(len(record.sequence))
    lens.sort()
    return lens


def accuracy_as_phred_qv(accuracy, max_qv=70):
    """
    Convert fractional accuracy to Phred QV: 0.999 --> 30

    returns: float or numpy array
    """
    if isinstance(accuracy, (float, int)):
        assert 0 <= accuracy <= 1.0
        if accuracy == 1:
            return max_qv
        return -10 * math.log(1 - accuracy) / Constants.LOG_10
    else:
        if isinstance(accuracy, (tuple, list)):
            accuracy = np.array(accuracy)
        error_rate = 1.0 - accuracy
        min_error_rate = 10 ** (-max_qv / 10.0)
        zero_error = error_rate < min_error_rate
        error_rate[zero_error] = min_error_rate
        return -10 * np.log(error_rate) / Constants.LOG_10


def compute_n50_from_file(fasta_file):
    """
    Convenience method to get N50 from a fasta file.
    """
    lens = get_fasta_readlengths(fasta_file)
    return compute_n50(lens)


def get_top_contigs(reference, max_contigs):
    """
    Get a list of contigs sorted in descending order by length
    :param reference: (str) path to reference_dir
    :param max_contigs: (int) max number of contigs to return
    """
    re = openReference(reference)
    return get_top_contigs_from_ref_entry(re, max_contigs)


def get_top_contigs_from_ref_entry(ref_entry, max_contigs):
    """
    Get a list of contigs sorted in descending order by length
    :param reference: ref_entry
    :param max_contigs: (int) max number of contigs to return
    """
    contigs = sorted(list(ref_entry.contigs),
                     key=len,
                     reverse=True)
    return contigs[:max_contigs]


def compute_n50(readlengths):
    """
    :param contig_lengths: list of contig lengths
    see get_fasta_readlengths
    """
    # this might not be the best expected behavior.
    # this could be a np.array
    if len(readlengths) == 0:
        return 0

    sorted_readlengths = np.sort(readlengths)
    n50 = sorted_readlengths[0]
    half_length = np.sum(readlengths) / 2.0

    total = 0
    for readlength in sorted_readlengths:
        if total < half_length:
            n50 = readlength
            total += n50
        else:
            break

    return n50


def add_plot_options(parser):
    parser.add_argument("--dpi", action="store", default=60,
                        help="dot/inch")
    parser.add_argument("--dumpdata", action="store_true",
                        help="Dump csv data for the plots.")
    return parser


def add_base_options(parser):
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug output.")
    parser.add_argument("output", type=validate_output_dir,
                        help="Output directory for associated files")
    parser.add_argument("report", type=validate_report,
                        help="Filename of JSON output report. Should be name only," +
                        "and will be written to output dir")
    return parser


def add_base_options_pbcommand(parser, title="JSON report"):
    """
    Eventual replacement for add_base_options(parser).
    """
    # XXX this is handled elsewhere now
    # common_options.add_debug_option(parser.arg_parser.parser)
    # NOTE 'output' is not part of tool contract!
    parser.arg_parser.parser.add_argument(
        "output",
        type=validate_output_dir,
        help="Output directory for associated files")
    parser.add_output_file_type(FileTypes.REPORT, "report", title,
                                description="Filename of JSON output report",
                                default_name="report")
    return parser


def compose(*funcs):
    """Functional composition
    [f, g, h] will be f(g(h(x)))
    """
    def compose_two(f, g):
        def c(x):
            return f(g(x))
        return c
    return functools.reduce(compose_two, funcs)

add_base_and_plot_options = compose(add_base_options, add_plot_options)


def get_base_parser(version, description):
    p = argparse.ArgumentParser(version=version, description=description)
    p = add_base_options(p)
    return p


def get_base_parser_with_plot_options(version, description):
    p = get_base_parser(version, description)
    p = add_plot_options(p)
    return p


def recfromcsv(fname, **kwargs):
    data = np.recfromcsv(fname, **kwargs)

    # np.recfromcsv will return a singleton if there's but a single row.
    # Beyond being unbelievably stupid, this is unindexable. Turn it into a
    # table with 1 row.
    if len(data.shape) < 1:
        tbl = np.recarray((1,), dtype=data.dtype)
        tbl[0] = data
        return tbl

    return data


def openReference(fname):
    """ Take a ReferenceSet, fasta or reference dir path and return a
    referenceSet.
    """
    if os.path.isdir(fname):
        raise ValueError("{r} is a directory, not a ReferenceSet".format(
                         r=fname))
    ref = ReferenceSet(fname)
    return ref


# FIXME can this be combined with compute_n50?
def compute_n50_from_bins(bins):
    """
    Compute n50 from the numpy array when the index is the length
    and the value is the number of items which have that length (i.e.,
    a histogram with the bin widths set to 1).

    :note: Bin width is assumed to be 1

    """
    total = 0
    for i, j in enumerate(bins):
        for _ in xrange(int(j)):
            total += i
    half_total = total / 2.0
    n50 = 0
    # initialize n50 by finding the first bin that != 0
    for i, bin_value in enumerate(bins):
        if bin_value != 0:
            n50 = i
            break
    rtotal = 0
    for i, bin_value in enumerate(bins):
        if bin_value != 0:
            for _ in xrange(int(bin_value)):
                if rtotal < half_total:
                    n50 = i
                    rtotal += i
                    #log.debug(("N50", n50, rtotal, half_total, total))
                else:
                    return n50
    msg = "Unable to compute n50 from {n} bins with sum {x}".format(
        n=len(bins), x=total)
    # warnings.warn(msg)
    log.warn(msg)
    return 0


def _dist_shaper(bmin, bmax, poolby, dist, trim_to=None):
    """Just change the bins and binlabels! Not the sample means etc.

    Args:
        dist: ([bin, ...], [label, ...])
    """
    try:
        bins, labels = dist
        assert len(bins) > 1, "Need more than 1 bin"
        assert len(bins) == len(labels), "Need same bin, label count"
        # bins and labels
        first = 0
        firstl = labels[0]
        last = 0
        lastl = labels[-1]
        bwidth = labels[1] - labels[0]
        for i, l in enumerate(labels):
            if l <= bmin:
                first = i
                firstl = l
            if l <= bmax:
                last = i
                lastl = l
            else:
                break
        # bound the range:
        newbins = bins[first:last + 1]
        newlabels = labels[first:last + 1]

        # pad the range:
        under = firstl - bmin
        if under > 0:
            lpad = int(under / bwidth)
            newbins = [0] * lpad + newbins
            newlabels = [newlabels[0] - (lpad - i) * bwidth
                         for i in range(lpad)] + newlabels
        under = bmax - lastl
        if under > 0:
            rpad = int(round(under / bwidth))
            newbins = newbins + [0] * rpad
            newlabels = newlabels + [newlabels[-1] + i * bwidth
                                     for i in range(1, rpad + 1)]

        assert (len(newbins) % poolby) == 0, ("pooling factor doesn't "
                                              "divide new "
                                              "nbins evenly")

        # collapse with poolby
        cbins = [sum(newbins[start:start + poolby])
                 for start in xrange(0, len(newbins), poolby)]
        clabels = [newlabels[start]
                   for start in xrange(0, len(newbins), poolby)]

        if not trim_to is None:
            cutoff = len(clabels)
            for i, lab in enumerate(clabels):
                if lab > trim_to:
                    cutoff = i
                    break
            cbins = cbins[:cutoff]
            clabels = clabels[:cutoff]

        bins = cbins
        labels = clabels

    except AssertionError as e:
        log.error("Malformed dist_shaper: " + e.message)
        return dist
    return (bins, labels)


def _cont_dist_shaper(shape_func, dist):
    """Just change the bins and binlabels! Not the sample means etc."""
    dist = deepcopy(dist)
    newbins, newlabels = shape_func((dist.bins, dist.labels))
    dist.bins = newbins
    dist.minBinValue = newlabels[0]
    dist.maxBinValue = newlabels[-1]
    dist.numBins = len(newbins)
    if len(newlabels) > 1:
        dist.binWidth = newlabels[1] - newlabels[0]
    return dist


def dist_shaper(dist_list, nbins=40, trim_excess=False):
    """Produce a function to modify a distribution to have a number of bins
    close to or greater than the number provided, such that all dists in the
    dist list can be transformed to have a consistent appearance.

    .. note:: Distributions must have the same bin width!
    .. note:: Distributions must have the same bin boundary locations (e.g. all
              bins must be integer multiples of binwidth)!

    This does produce quite a few empty bins when you would think they could be
    redistributed, but the rule is that redistribution can only append with
    even pooling! No splitting data allowed, only joining.

    Args:
        dist_list: [([bin, ...], [label, ...]), ...]
    """
    try:
        assert nbins > 0, "nbins must be greater than 0"
        assert len(dist_list) >= 1
        assert len(dist_list[0]) > 1
        bmin = sys.maxint
        bmax = 0
        # prepare for a relaxation of the binwidth requirement (wont happen for
        # a while):
        # we want the largest binwidth, as we have to pad up so that the
        # coarsest distribution has ~ 40 nbins, but they all cover the same
        # range
        bwidth = 0
        for bins, labels in dist_list:
            first = None
            last = None
            if sum(bins) > 0:
                for i, (bv, bl) in enumerate(zip(bins, labels)):
                    if bv != 0:
                        last = bl
                        if first is None:
                            first = bl
                if bmin > first:
                    bmin = first
                if bmax < last:
                    bmax = last
            binWidth = labels[1] - labels[0]
            if binWidth > bwidth:
                bwidth = binWidth
        assert bmin <= bmax
        assert bwidth != 0

        # you have to add on that last bin
        curbins = int(round((bmax - bmin) / bwidth + 1))
        poolby = curbins / float(nbins)

        # poolby should be an int >= 1 before it hits _dist_shaper, so if it is
        # less than that, lets adjust the other values here and now:
        if poolby < 1:
            if trim_excess:
                poolby = 1.0
                curbins = (bmax - bmin) / bwidth + 1
                nbins = curbins
            else:
                # the total to be padded:
                pad = int(round((1 - poolby) * nbins))
                pad_left = min(pad / 2, (bmin / bwidth))
                pad_right = pad - pad_left
                bmin = bmin - (pad_left * bwidth)
                bmax = bmax + (pad_right * bwidth)
                curbins = (bmax - bmin) / bwidth + 1
                poolby = curbins / float(nbins)

        curbins = int(round(curbins))
        poolby = int(math.ceil(poolby))

        # poolby should be an even divisor of curbins before it hits
        # _dist_shaper, so if it isn't, lets adjust the other values here and
        # now:
        actual_bmax = None
        if trim_excess:
            actual_bmax = bmax
        if poolby > 1:
            # we don't want to end up with fewer bins than predicted, so take
            # the max
            curbins = int(max(round(curbins / poolby), nbins)) * int(poolby)
            # round up that multiplier, but remember that it is relative to
            # bmin, not always zero. we only go to the beginning of the last
            # bin, so subtract one before multiplying in the binwidth
            bmax = (curbins + bmin / bwidth - 1) * bwidth

        # you have to add on that last bin
        curbins = int(round((bmax - bmin) / bwidth) + 1)
        poolby = int(curbins / nbins)
        return functools.partial(_dist_shaper, bmin, bmax, poolby,
                                 trim_to=actual_bmax)
    except AssertionError as e:
        log.error(e.message)
        return lambda x: x


def continuous_dist_shaper(dist_list, nbins=40, trim_excess=False):
    generic_dist_list = [(d.bins, d.labels) for d in dist_list]
    shaper = dist_shaper(generic_dist_list, nbins=nbins,
                         trim_excess=trim_excess)
    return functools.partial(_cont_dist_shaper, shaper)
