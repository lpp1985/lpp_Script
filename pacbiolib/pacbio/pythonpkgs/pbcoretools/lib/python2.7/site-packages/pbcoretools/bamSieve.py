
"""
Tool for subsetting a BAM or PacBio DataSet file based on either a whitelist of
hole numbers or a percentage of reads to be randomly selected.
"""

from __future__ import division
from collections import defaultdict
import subprocess
import warnings
import logging
import random
import os.path as op
import re
import sys

from pysam import AlignmentFile

from pbcommand.common_options import (add_log_quiet_option,
    add_log_verbose_option)
from pbcommand.cli import (pacbio_args_runner,
    get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.io import openDataFile, openDataSet, IndexedBamReader, ReadSet

VERSION = "0.1.1"

log = logging.getLogger(__name__)


def _process_zmw_list(zmw_list):
    zmws = set()
    if zmw_list is None:
        return zmws
    elif isinstance(zmw_list, set):
        return zmw_list
    elif isinstance(zmw_list, (list, tuple)):
        return set(zmw_list)
    elif op.isfile(zmw_list):
        base, ext = op.splitext(zmw_list)
        if ext in [".bam", ".xml"]:
            with openDataFile(zmw_list) as ds_zmw:
                for f in ds_zmw.resourceReaders():
                    zmws.update(set(list(f.holeNumber)))
        else:
            with open(zmw_list) as f:
                lines = f.read().splitlines()
                zmws.update(set([int(x) for x in lines]))
    else:
        zmws.update(set([int(x) for x in zmw_list.split(",")]))
    return zmws


def _anonymize_sequence(rec):
    rseq_ = [random.randint(0,3) for i in range(len(rec.query_sequence))]
    rseq = "".join(["ACTG"[i] for i in rseq_])
    rec.query_sequence = rseq
    return rec


def _create_whitelist(ds_in, percentage=None, count=None):
    zmws = set()
    movies = set()
    for i_file, bam in enumerate(ds_in.resourceReaders()):
        movies.update(set([rg["MovieName"] for rg in bam.readGroupTable]))
        zmws.update(set(bam.pbi.holeNumber))
    if len(movies) > 1:
        warnings.warn("The input BAM/dataset contains multiple movies, "+
                      "which may have overlapping ZMWs.")
    if percentage is not None:
        count = int(len(zmws) * percentage / 100.0)
    zmws = list(zmws)
    have_zmws = set()
    whitelist = set()
    k = 0
    while k < count:
        i_zmw = random.randint(0, len(zmws) - 1)
        if not zmws[i_zmw] in have_zmws:
            whitelist.add(zmws[i_zmw])
            have_zmws.add(zmws[i_zmw])
            k += 1
    return whitelist


def _process_bam_whitelist(bam_in, bam_out, whitelist, blacklist,
                           use_barcodes=False, anonymize=False):
    def _is_whitelisted(x):
        if ((len(whitelist) > 0 and x in whitelist) or
            (len(blacklist) > 0 and not x in blacklist)):
            return True
    have_zmws = set()
    have_records = []
    def _add_read(i_rec, zmw):
        rec = bam_in[i_rec]
        if anonymize:
            _anonymize_sequence(rec.peer)
        bam_out.write(rec.peer)
        have_zmws.add(zmw)
        have_records.append(i_rec)
    if use_barcodes:
        for i_rec in range(len(bam_in.holeNumber)):
            bc_fwd = bam_in.bcForward[i_rec]
            bc_rev = bam_in.bcReverse[i_rec]
            if _is_whitelisted(bc_fwd) or _is_whitelisted(bc_rev):
                _add_read(i_rec, bam_in.holeNumber[i_rec])
    else:
        for i_rec, zmw in enumerate(bam_in.holeNumber):
            if _is_whitelisted(zmw):
                _add_read(i_rec, zmw)
    return len(have_records), have_zmws


def filter_reads(input_bam,
                 output_bam,
                 whitelist=None,
                 blacklist=None,
                 percentage=None,
                 count=None,
                 seed=None,
                 ignore_metadata=False,
                 relative=None,
                 anonymize=False,
                 use_barcodes=False):
    if output_bam is None:
        log.error("Must specify output file")
        return 1
    output_bam = op.abspath(output_bam)
    if not op.isdir(op.dirname(output_bam)):
        log.error("Output path '{d}' does not exist.".format(
                  d=op.dirname(output_bam)))
        return 1
    n_specified = 4 - [whitelist, blacklist, percentage, count].count(None)
    if n_specified != 1:
        log.error("You must choose one and only one of the following "+
                  "options: --whitelist, --blacklist, --count, --percentage")
        return 1
    if seed is not None:
        random.seed(seed)
    if whitelist is None and blacklist is None:
        if not 0 < percentage < 100 and not count > 0:
            log.error("No reads selected for output.")
            return 1
    output_ds = None
    if output_bam.endswith(".xml"):
        if not input_bam.endswith(".xml"):
            print "DataSet output only supported for DataSet inputs."
            return 1
        ds_type = output_bam.split(".")[-2]
        ext2 = {
            "subreadset": "subreads",
            "alignmentset": "subreads",
            "consensusreadset": "ccs",
            "consensusalignmentset": "ccs"
        }
        if not ds_type in ext2:
            raise ValueError("Invalid dataset type 't'".format(t=ds_type))
        output_ds = output_bam
        output_bam = ".".join(output_ds.split(".")[:-2] +
                              [ext2[ds_type], "bam"])
    if output_bam == input_bam:
        log.error("Input and output files must not be the same path")
        return 1
    elif not output_bam.endswith(".bam"):
        log.error("Output file name must end in either '.bam' or '.xml'")
        return 1
    n_file_reads = 0
    have_zmws = set()
    scraps_bam = barcode_set = None
    with openDataFile(input_bam) as ds_in:
        if not isinstance(ds_in, ReadSet):
            raise TypeError("{t} is not an allowed dataset type".format(
                            t=type(ds_in).__name__))
        # TODO(nechols)(2016-03-11): refactor this to enable propagation of
        # filtered scraps
        if not ds_in.isIndexed:
            log.error("Input BAM must have accompanying .pbi index")
            return 1
        for ext_res in ds_in.externalResources:
            if ext_res.barcodes is not None:
                assert barcode_set is None or barcode_set == ext_res.barcodes
                barcode_set = barcode_set
        f1 = ds_in.resourceReaders()[0]
        if percentage is not None or count is not None:
            whitelist = _create_whitelist(ds_in, percentage, count)
        # convert these to Python sets
        _whitelist = _process_zmw_list(whitelist)
        _blacklist = _process_zmw_list(blacklist)
        scraps_in = None
        if output_ds is not None and output_ds.endswith(".subreadset.xml"):
            for ext_res in ds_in.externalResources:
                if ext_res.scraps is not None:
                    if use_barcodes:
                        log.warn("Scraps BAM is present but lacks "+
                                 "barcodes - will not be propagated "+
                                 "to output SubreadSet")
                    else:
                        scraps_in = IndexedBamReader(ext_res.scraps)
                    break
        with AlignmentFile(output_bam, 'wb',
                           template=f1.peer) as bam_out:
            for bam_in in ds_in.resourceReaders():
                n_records, have_zmws_ =_process_bam_whitelist(
                    bam_in, bam_out,
                    whitelist=_whitelist,
                    blacklist=_blacklist,
                    use_barcodes=use_barcodes,
                    anonymize=anonymize)
                n_file_reads += n_records
                have_zmws.update(have_zmws_)
        if scraps_in is not None:
            scraps_bam = re.sub("subreads.bam$", "scraps.bam", output_bam)
            with AlignmentFile(scraps_bam, 'wb',
                               template=scraps_in.peer) as scraps_out:
                for ext_res in ds_in.externalResources:
                    if ext_res.scraps is not None:
                        scraps_in_ = IndexedBamReader(ext_res.scraps)
                        n_records, have_zmws_ =_process_bam_whitelist(
                            scraps_in_, scraps_out, _whitelist, _blacklist,
                            use_barcodes=use_barcodes,
                            anonymize=anonymize)
                        have_zmws.update(have_zmws_)
    if n_file_reads == 0:
        log.error("No reads written")
        return 1
    log.info("{n} records from {z} ZMWs written".format(
        n=n_file_reads, z=len(have_zmws)))
    def _run_pbindex(bam_file):
        try:
            rc = subprocess.call(["pbindex", bam_file])
        except OSError as e:
            if e.errno == 2:
                log.warn("pbindex not present, will not create .pbi file")
            else:
                raise
    _run_pbindex(output_bam)
    if output_ds is not None:
        with openDataSet(input_bam) as ds_in:
            ds_out = ds_in.__class__(output_bam)
            if scraps_bam is not None:
                _run_pbindex(scraps_bam)
                ds_out.externalResources[0].scraps = scraps_bam
                # XXX it doesn't pick up the .pbi file - sort of annoying
                # but since the pbcore API doesn't provide a read for the
                # scraps automatically anyway, the impact is minimal
            if barcode_set is not None:
                ds_out.externalResources[0].barcodes = barcode_set
            if not ignore_metadata:
                ds_out.metadata = ds_in.metadata
                ds_out.updateCounts()
            if relative:
                ds_out.makePathsRelative(op.dirname(output_ds))
            ds_out.write(output_ds)
            log.info("wrote {t} XML to {x}".format(
                     t=ds_out.__class__.__name__, x=output_ds))
    return 0


def _iter_bam_files(input_file):
    if input_file.endswith(".xml"):
        with openDataFile(input_file) as ds_in:
            if not ds_in.isIndexed:
                log.warning("Unindexed file(s), this may be very slow")
            for rr in ds_in.resourceReaders():
                yield rr
    else:
        if op.exists(input_file + ".pbi"):
            with IndexedBamReader(input_file) as bam_in:
                yield bam_in
        else:
            with BamReader(input_file) as bam_in:
                yield bam_in


def show_zmws(input_file):
    zmws = []
    for rr in _iter_bam_files(input_file):
        if isinstance(rr, IndexedBamReader):
            zmws.extend(list([int(x) for x in rr.holeNumber]))
        else:
            zmws.extend([int(rec.HoleNumber) for rec in rr])
    print "\n".join([str(x) for x in sorted(list(set(zmws)))])


def run(args):
    if args.show_zmws:
        if [args.whitelist, args.blacklist, args.percentage].count(None) != 3:
            log.warning("Ignoring unused filtering arguments")
        show_zmws(args.input_bam)
        return 0
    return filter_reads(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        whitelist=args.whitelist,
        blacklist=args.blacklist,
        percentage=args.percentage,
        count=args.count,
        seed=args.seed,
        ignore_metadata=args.ignore_metadata,
        relative=args.relative,
        anonymize=args.anonymize,
        use_barcodes=args.barcodes)


def get_parser():
    p = get_default_argparser_with_base_opts(
            version=VERSION,
            description=__doc__,
            default_level="WARN")
    p.add_argument("input_bam",
                   help="Input BAM or DataSet from which reads will be read")
    p.add_argument("output_bam", nargs='?', default=None,
                   help="Output BAM or DataSet to which filtered reads will "
                        "be written")
    p.add_argument("--show-zmws", action="store_true", default=False,
                   help="Print a list of ZMWs and exit")
    p.add_argument("--whitelist", action="store", default=None,
                   help="Comma-separated list of ZMWs, or file containing " +
                        "whitelist of one hole number per line, or " +
                        "BAM/DataSet file from which to extract ZMWs")
    p.add_argument("--blacklist", action="store", default=None,
                   help="Opposite of --whitelist, specifies ZMWs to discard")
    p.add_argument("--percentage", action="store", type=float, default=None,
                   help="If you prefer to recover a percentage of a SMRTcell "
                        "rather than a specific list of reads specify that "
                        "percentage (range 0-100) here")
    p.add_argument("-n", "--count", action="store", type=int, default=None,
                   help="Recover a specific number of ZMWs picked at random")
    p.add_argument("-s", "--seed", action="store", type=int, default=None,
                   help="Random seed for selecting a percentage of reads")
    p.add_argument("--ignore-metadata", action="store_true",
                   help="Discard input DataSet metadata")
    p.add_argument("--relative", action="store_true",
                   help="Make external resource paths relative")
    p.add_argument("--anonymize", action="store_true",
                   help="Randomize sequences for privacy")
    p.add_argument("--barcodes", action="store_true",
                   help="Indicates that the whitelist or blacklist contains "+
                        "barcode indices instead of ZMW numbers")
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=run,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
