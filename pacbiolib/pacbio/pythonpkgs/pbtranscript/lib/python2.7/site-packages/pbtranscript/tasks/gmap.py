
"""
Calls GMAP.
"""

import subprocess  # FIXME use pbcommand wrapper (once this is stable)
import tempfile
import logging
import os.path as op
import os
import sys

import pysam

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes, SymbolTypes, get_pbparser
from pbcommand.utils import setup_log
from pbcore.io import ContigSet, ReferenceSet, AlignmentSet

from pbtranscript.Utils import filter_sam

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbtranscript.tasks.gmap"
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m pbtranscript.tasks.gmap --resolved-tool-contract"


def run_gmap(transcripts_file, reference_file, alignment_file, nproc):
    db_name = "gmap_db"
    ref_fasta_name = transcripts_fasta_name = None
    with ReferenceSet(reference_file) as rs:
        rfs = rs.toExternalFiles()
        assert len(rfs) == 1
        ref_fasta_name = rfs[0]
    # XXX in our setup, i.e. with the old reference layout, the .fasta might
    # be in a 'sequence' directory, with gmap_db one level up, so we try
    # a couple of paths to look for an existing gmap_db
    ref_fasta_dir = op.dirname(ref_fasta_name)
    ref_path = os.getcwd()
    for base_path in [ref_fasta_dir, op.dirname(ref_fasta_dir)]:
        log.info("Looking for gmap_db in %s" % base_path)
        if op.exists(op.join(base_path, db_name)):
            log.info("Found existing gmap_db, gmap_build will not be run")
            ref_path = base_path
            break
    else:
        log.info("gmap_build will need to be run on %s" % ref_fasta_name)
        # workaround for hardcoded paths in gmap_build
        bin_dir = os.environ['_SMRT_GMAP_BIN']  # FIXME whatever key Herb uses
        args1 = [
            "gmap_build",
            "-B", bin_dir,
            "-k", "12",
            "--db=%s" % db_name,
            # XXX what to do about the directory here?  we really need a way
            # to cache or pre-generate these
            #"--dir=%s" % ref_path,  # FIXME spaces?
            "-D", ref_path,
            ref_fasta_name,
        ]
        log.info("ARGUMENTS: {a}".format(a=" ".join(args1)))
        with open("gmap_build.out", "w") as stdout, \
             open("gmap_build.err", "w") as stderr:
            rc = subprocess.call(args1, stdout=stdout, stderr=stderr)
            assert rc == 0, "unexpected exit code {c}".format(c=rc)
    with ContigSet(transcripts_file) as cs:
        cfs = cs.toExternalFiles()
        assert len(cfs) == 1
        transcripts_fasta_name = cfs[0]
    args2 = [
        "gmap",
        "-n", "0",
        "-t", str(nproc),
        "--sam-use-0M",
        "-f", "samse",
        "-D", ref_path,
        "-d", db_name,
        #"--split-output=gmap_tmp",
        transcripts_fasta_name,
    ]
    log.info("ARGUMENTS: %s" % " ".join(args2))
    sam_file = tempfile.NamedTemporaryFile(suffix=".sam", delete=True).name
    with open(sam_file, "w") as stdout, open("gmap.log", "w") as stderr:
        rc = subprocess.call(args2, stdout=stdout, stderr=stderr)
        assert rc == 0, "unexpected exit code {c}".format(c=rc)
    sam_file_2 = tempfile.NamedTemporaryFile(suffix=".sam", delete=True).name
    filter_sam(sam_file, sam_file_2)
    os.remove(sam_file)
    with pysam.AlignmentFile(sam_file_2, "r") as sam_in:
        log.info("Writing alignments to %s" % alignment_file)
        with pysam.AlignmentFile(alignment_file, "wb",
                                 template=sam_in) as bam_out:
            for rec in sam_in:
                bam_out.write(rec)
    os.remove(sam_file_2)
    # FIXME this bam file of course looks nothing like the PacBio standard!
    # (which also makes testing difficult, since we usually run pbvalidate on
    # all outputs)
    #assert subprocess.call(["pbindex", "gmap.aligned.bam"]) == 0
    pysam.index(alignment_file)
    return 0


def get_contract_parser():
    p = get_pbparser(
        tool_id=Constants.TOOL_ID,
        version=Constants.VERSION,
        name=Constants.TOOL_ID,
        description=__doc__,
        driver_exe=Constants.DRIVER_EXE,
        nproc=SymbolTypes.MAX_NPROC)
    p.add_input_file_type(FileTypes.DS_CONTIG, "seq_in",
                          name="ContigSet",
                          description="Input transcripts")
    p.add_input_file_type(FileTypes.DS_REF, "ref_in",
                          name="ReferenceSet",
                          description="Reference genome")
    p.add_output_file_type(FileTypes.BAM, "aln_out",
                           name="Alignments",
                           description="BAM alignments file",
                           default_name="gmap_alignments")
    return p


def args_runner(args):
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    return run_gmap(
        transcripts_file=rtc.task.input_files[0],
        reference_file=rtc.task.input_files[1],
        alignment_file=rtc.task.output_files[0],
        nproc=rtc.task.nproc)


def main(argv=sys.argv[1:]):
    mp = get_contract_parser()
    return pbparser_runner(
        argv=argv,
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
