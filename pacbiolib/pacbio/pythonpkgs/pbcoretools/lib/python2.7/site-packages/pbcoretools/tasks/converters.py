
"""
Tool contract wrappers for miscellaneous quick functions.
"""

from collections import defaultdict
import functools
import tempfile
import logging
import shutil
import gzip
import re
import os.path as op
import os
import sys

from pbcore.io import (SubreadSet, HdfSubreadSet, FastaReader, FastaWriter,
                       FastqReader, FastqWriter, BarcodeSet, ExternalResource,
                       ExternalResources, openDataSet, ContigSet, ReferenceSet,
                       GmapReferenceSet)
from pbcommand.engine import run_cmd
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, SymbolTypes, OutputFileType

log = logging.getLogger(__name__)

TOOL_NAMESPACE = 'pbcoretools'
DRIVER_BASE = "python -m pbcoretools.tasks.converters "

registry = registry_builder(TOOL_NAMESPACE, DRIVER_BASE)

def _run_bax_to_bam(input_file_name, output_file_name):
    base_name = ".".join(output_file_name.split(".")[:-2])
    input_file_name_tmp = input_file_name
    # XXX bax2bam won't write an hdfsubreadset unless the input is XML too
    if input_file_name.endswith(".bax.h5"):
        input_file_name_tmp = tempfile.NamedTemporaryFile(
            suffix=".hdfsubreadset.xml").name
        ds_tmp = HdfSubreadSet(input_file_name)
        ds_tmp.write(input_file_name_tmp)
    args =[
        "bax2bam",
        "--subread",
        "-o", base_name,
        "--output-xml", output_file_name,
        "--xml", input_file_name_tmp
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    with SubreadSet(output_file_name) as ds:
        ds.assertIndexed()
    return 0


def run_bax_to_bam(input_file_name, output_file_name):
    with HdfSubreadSet(input_file_name) as ds_in:
        movies = set()
        for rr in ds_in.resourceReaders():
            movies.add(rr.movieName)
        if len(movies) > 1:
            out_dir = os.path.dirname(output_file_name)
            ds_out_files = []
            for bax_file in ds_in.toExternalFiles():
                output_file_name_tmp = os.path.join(out_dir, ".".join(
                    os.path.basename(bax_file).split(".")[:-2]) +
                    ".hdfsubreadset.xml")
                rc = _run_bax_to_bam(bax_file, output_file_name_tmp)
                if rc != 0:
                    log.error("bax2bam failed")
                    return rc
                ds_out_files.append(output_file_name_tmp)
            ds = SubreadSet(*ds_out_files)
            ds.name = ds_in.name
            if 'Description' in ds_in.objMetadata:
                ds.objMetadata['Description'] = ds_in.objMetadata['Description']
                ds.metadata.merge(ds_in.metadata)
            ds.write(output_file_name)
        else:
            return _run_bax_to_bam(input_file_name, output_file_name)
    return 0


def run_bam_to_bam(subread_set_file, barcode_set_file, output_file_name,
                   nproc=1, score_mode="symmetric"):
    if not score_mode in ["asymmetric", "symmetric"]:
        raise ValueError("Unrecognized score mode '{m}'".format(m=score_mode))
    bc = BarcodeSet(barcode_set_file)
    if len(bc.resourceReaders()) > 1:
        raise NotImplementedError("Multi-FASTA BarcodeSet input is not supported.")
    new_prefix = re.sub(".subreadset.xml$", "", output_file_name)
    args = [
        "bam2bam",
        "-j", str(nproc),
        "-b", str(nproc),
        "-o", new_prefix,
        "--barcodes", barcode_set_file,
        "--scoreMode", score_mode,
        subread_set_file
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    assert op.isfile(output_file_name)
    tmp_out = op.join(op.dirname(output_file_name),
                      "tmp_" + op.basename(output_file_name))
    shutil.move(output_file_name, tmp_out)
    with SubreadSet(tmp_out, strict=True) as ds:
        with SubreadSet(subread_set_file) as ds_in:
            ds.metadata = ds_in.metadata
            ds.name = ds_in.name + " (barcoded)"
        ds.updateCounts()
        ds.newUuid()
        ds.write(output_file_name)
    return 0


def _unzip_fastx(gzip_file_name, fastx_file_name):
    with gzip.open(gzip_file_name, "rb") as gz_in:
        with open(fastx_file_name, "wb") as fastx_out:
            fastx_out.write(gz_in.read())


def archive_files(input_file_names, output_file_name, remove_path=True):
    """
    Create a gzipped tarball from a list of input files.

    :param remove_path: if True, the directory will be removed from the input
                        file names before archiving.  All inputs and the output
                        file must be in the same directory for this to work.
    """
    if remove_path:
        input_file_names = [op.basename(fn) for fn in input_file_names]
    args = ["tar", "-czf", output_file_name] + input_file_names
    log.info("Running '{a}'".format(a=" ".join(args)))
    _cwd = os.getcwd()
    try:
        # we want the files to have no leading path
        os.chdir(op.dirname(output_file_name))
        result = run_cmd(" ".join(args),
                         stdout_fh=sys.stdout,
                         stderr_fh=sys.stderr)
    except Exception:
        raise
    else:
        if result.exit_code != 0:
            return result.exit_code
    finally:
        os.chdir(_cwd)
    assert op.isfile(output_file_name)
    return 0


def _run_bam_to_fastx(program_name, fastx_reader, fastx_writer,
                     input_file_name, output_file_name, tmp_dir=None):
    assert isinstance(program_name, basestring)
    barcode_mode = False
    if output_file_name.endswith(".gz"):
        with openDataSet(input_file_name) as ds_in:
            barcode_mode = ds_in.isBarcoded
    tmp_out_prefix = tempfile.NamedTemporaryFile(dir=tmp_dir).name
    args = [
        program_name,
        "-o", tmp_out_prefix,
        input_file_name,
    ]
    if barcode_mode:
        args.insert(1, "--split-barcodes")
    log.info(" ".join(args))
    result = run_cmd(" ".join(args),
                     stdout_fh=sys.stdout,
                     stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    else:
        base_ext = re.sub("bam2", "", program_name) 
        if not barcode_mode:
            tmp_out = "{p}.{b}.gz".format(p=tmp_out_prefix, b=base_ext)
            assert os.path.isfile(tmp_out), tmp_out
            if output_file_name.endswith(".gz"):
                log.info("cp {t} {f}".format(t=tmp_out, f=output_file_name))
                shutil.copyfile(tmp_out, output_file_name)
            else:
                _unzip_fastx(tmp_out, output_file_name)
            os.remove(tmp_out)
        else:
            suffix = "{f}.gz".format(f=base_ext)
            tmp_out_dir = op.dirname(tmp_out_prefix)
            tc_out_dir = op.dirname(output_file_name)
            barcoded_file_names = []
            # find the barcoded FASTX files and unzip them to the same
            # output directory and file prefix as the ultimate output
            for fn in os.listdir(tmp_out_dir):
                fn = op.join(tmp_out_dir, fn)
                if fn.startswith(tmp_out_prefix) and fn.endswith(suffix):
                    bc_fwd_rev = fn.split(".")[-3].split("_")
                    suffix2 = ".{f}_{r}.{t}".format(
                        f=bc_fwd_rev[0], r=bc_fwd_rev[1], t=base_ext)
                    assert fn == tmp_out_prefix + suffix2 + ".gz"
                    fn_out = re.sub(".gz$", suffix2, output_file_name)
                    fastx_out = op.join(tc_out_dir, fn_out)
                    _unzip_fastx(fn, fastx_out)
                    barcoded_file_names.append(fn_out)
                    os.remove(fn)
            assert len(barcoded_file_names) > 0
            return archive_files(barcoded_file_names, output_file_name)
    return 0


def split_laa_fastq(input_file_name, output_file_base):
    """
    Split an LAA FASTQ file into one file per barcode.
    """
    if op.getsize(input_file_name) == 0:
        return []
    records = defaultdict(list)
    with FastqReader(input_file_name) as fastq_in:
        for rec in fastq_in:
            bc_id = rec.id.split("_")[0]
            records[bc_id].append(rec)
    outputs = []
    for bc_id in sorted(records.keys()):
        ofn = "{b}.{i}.fastq".format(b=output_file_base, i=bc_id)
        with FastqWriter(ofn) as fastq_out:
            for rec in records[bc_id]:
                fastq_out.writeRecord(rec)
        outputs.append(ofn)
    return outputs


def split_laa_fastq_archived(input_file_name, output_file_name):
    """
    Split an LAA FASTQ file into one file per barcode and package as tar.gz.
    """
    base, ext = op.splitext(output_file_name)
    assert (ext == ".gz")
    if base.endswith(".tar"):
        base, ext2 = op.splitext(base)
    fastq_files = [op.basename(fn) for fn in split_laa_fastq(input_file_name, base)]
    if len(fastq_files) == 0: # workaround for empty input
        with open(output_file_name, "wb") as tar_out:
            return 0
    return archive_files(fastq_files, output_file_name)


def run_fasta_to_fofn(input_file_name, output_file_name):
    args = ["echo", input_file_name, ">", output_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                     stderr_fh=sys.stderr)
    return result.exit_code


def run_fasta_to_referenceset(input_file_name, output_file_name):
    args = ["dataset create", "--type ReferenceSet", "--generateIndices",
            output_file_name, input_file_name]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                     stderr_fh=sys.stderr)
    # the '.py' name difference will be resolved in pbdataset/pbcoretools, but
    # for now, work with either
    if result.exit_code == 127:
        args = ["dataset.py create", "--type ReferenceSet",
                "--generateIndices",
                output_file_name, input_file_name]
        log.info(" ".join(args))
        result = run_cmd(" ".join(args), stdout_fh = sys.stdout,
                         stderr_fh=sys.stderr)
    return result.exit_code


def __run_fasta_to_reference(program_name, dataset_class,
                             input_file_name, output_file_name,
                             organism=None, reference_name=None,
                             ploidy="haploid"):
    if reference_name is None or reference_name == "":
        reference_name = op.splitext(op.basename(input_file_name))[0]
    ds_in = ContigSet(input_file_name)
    if len(ds_in.externalResources) > 1:
        raise TypeError("Only a single FASTA file is supported as input.")
    fasta_file_name = ds_in.externalResources[0].resourceId
    output_dir_name = op.dirname(output_file_name)
    args = [
        program_name,
        "--organism", str(organism) if organism != "" else "unknown",
        "--ploidy", str(ploidy) if ploidy != "" else "unknown",
        "--debug",
        fasta_file_name,
        output_dir_name,
        reference_name
    ]
    log.info(" ".join(args))
    result = run_cmd(" ".join(args), stdout_fh=sys.stdout, stderr_fh=sys.stderr)
    if result.exit_code != 0:
        return result.exit_code
    ref_file = op.join(output_dir_name, reference_name,
                       "{t}.xml".format(t=dataset_class.__name__.lower()))
    assert op.isfile(ref_file), ref_file
    with dataset_class(ref_file, strict=True) as ds_ref:
        ds_ref.makePathsAbsolute()
        log.info("saving final {t} to {f}".format(
                 f=output_file_name, t=dataset_class.__name__))
        ds_ref.write(output_file_name)
    return 0


run_fasta_to_reference = functools.partial(__run_fasta_to_reference,
    "fasta-to-reference", ReferenceSet)
run_fasta_to_gmap_reference = functools.partial(__run_fasta_to_reference,
    "fasta-to-gmap-reference", GmapReferenceSet)


run_bam_to_fasta = functools.partial(_run_bam_to_fastx, "bam2fasta",
    FastaReader, FastaWriter)
run_bam_to_fastq = functools.partial(_run_bam_to_fastx, "bam2fastq",
    FastqReader, FastqWriter)


subreads_from_h5_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                            "Subreads", "Subread data in XML dataset",
                                            "Imported SubreadSet", "subreads")
subreads_barcoded_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                             "SubreadSet",
                                             "Barcoded Subreads",
                                             "Barcoded Subreads DataSet XML",
                                             "subreads_barcoded")

@registry("h5_subreads_to_subread", "0.1.0",
          FileTypes.DS_SUBREADS_H5,
          subreads_from_h5_file_type, is_distributed=True, nproc=1)
def run_bax2bam(rtc):
    return run_bax_to_bam(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2bam_barcode", "0.1.0",
          (FileTypes.DS_SUBREADS, FileTypes.DS_BARCODE),
          subreads_barcoded_file_type,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options={"score_mode":"symmetric"})
def run_bam2bam(rtc):
    return run_bam_to_bam(
        subread_set_file=rtc.task.input_files[0],
        barcode_set_file=rtc.task.input_files[1],
        output_file_name=rtc.task.output_files[0],
        nproc=rtc.task.nproc,
        score_mode=rtc.task.options["pbcoretools.task_options.score_mode"])


fasta_file_type = OutputFileType(FileTypes.FASTA.file_type_id, "fasta", "FASTA file",
                                 "Reads in FASTA format", "reads")
fastq_file_type = OutputFileType(FileTypes.FASTQ.file_type_id, "fastq", "FASTQ file",
                                 "Reads in FASTQ format", "reads")
fasta_gzip_file_type = OutputFileType(FileTypes.GZIP.file_type_id, "fasta_gz",
                                      "FASTA file(s)",
                                      "Seqeunce data converted to FASTA Format",
                                      "reads.fasta") # yuck - could be .tar !
fastq_gzip_file_type = OutputFileType(FileTypes.GZIP.file_type_id, "fastq",
                                      "FASTQ file(s)",
                                      "Sequence data converted to FASTQ format",
                                      "reads.fastq")

@registry("bam2fastq", "0.1.0",
          FileTypes.DS_SUBREADS,
          fastq_file_type, is_distributed=True, nproc=1)
def run_bam2fastq(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fasta_archive", "0.1.0",
          FileTypes.DS_SUBREADS,
          fasta_gzip_file_type, is_distributed=True, nproc=1)
def run_bam2fasta_archive(rtc):
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fastq_archive", "0.1.0",
          FileTypes.DS_SUBREADS,
          fastq_gzip_file_type, is_distributed=True, nproc=1)
def run_bam2fastq_archive(rtc):
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


fofn_file_type = OutputFileType(FileTypes.FOFN.file_type_id, "FOFN file",
                                "FOFN file", "List of input files", "files")

@registry("fasta2fofn", "0.1.0",
          FileTypes.FASTA,
          fofn_file_type, is_distributed=False, nproc=1)
def run_fasta2fofn(rtc):
    return run_fasta_to_fofn(rtc.task.input_files[0], rtc.task.output_files[0])


ref_file_type = OutputFileType(FileTypes.DS_REF.file_type_id, "ReferenceSet",
                               "Reference Dataset",
                               "PacBio Reference DataSet XML", "reference")

@registry("fasta2referenceset", "0.1.0",
          FileTypes.FASTA,
          ref_file_type, is_distributed=True, nproc=1)
def run_fasta2referenceset(rtc):
    return run_fasta_to_referenceset(rtc.task.input_files[0],
                                     rtc.task.output_files[0])


@registry("fasta_to_reference", "0.1.0",
          FileTypes.FASTA,
          ref_file_type, is_distributed=True, nproc=1,
          options={
                "organism": "",
                "ploidy": "haploid",
                "reference_name":""
          })
def run_fasta_to_reference_pbscala(rtc):
    return run_fasta_to_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options["pbcoretools.task_options.reference_name"],
        organism=rtc.task.options["pbcoretools.task_options.organism"],
        ploidy=rtc.task.options["pbcoretools.task_options.ploidy"])


gmap_ref_file_type = OutputFileType(FileTypes.DS_GMAP_REF.file_type_id, "GmapReferenceSet",
                               "GmapReferenceSet XML",
                               "PacBio GMAP Reference DataSet XML", "reference")

@registry("fasta_to_gmap_reference", "0.1.0",
          FileTypes.FASTA,
          gmap_ref_file_type, is_distributed=True, nproc=1,
          options={
                "organism": "",
                "ploidy": "haploid",
                "reference_name":""
          })
def _run_fasta_to_gmap_reference(rtc):
    return run_fasta_to_gmap_reference(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        reference_name=rtc.task.options["pbcoretools.task_options.reference_name"],
        organism=rtc.task.options["pbcoretools.task_options.organism"],
        ploidy=rtc.task.options["pbcoretools.task_options.ploidy"])


fasta_ccs_file_type = OutputFileType(FileTypes.GZIP.file_type_id, "fasta_gz",
                                     "Consensus Sequences (FASTA)",
                                     "Consensus sequences generated from CCS2",
                                     "ccs.fasta")
fastq_ccs_file_type = OutputFileType(FileTypes.GZIP.file_type_id, "fastq_gz",
                                     "Consensus Sequences (FASTQ)",
                                     "Consensus sequences generated from CCS2",
                                     "ccs.fastq")

@registry("bam2fastq_ccs", "0.1.0",
          FileTypes.DS_CCS,
          fastq_ccs_file_type, is_distributed=True, nproc=1)
def run_bam2fastq_ccs(rtc):
    """
    Duplicate of run_bam2fastq, but with ConsensusReadSet as input.
    """
    return run_bam_to_fastq(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("bam2fasta_ccs", "0.1.0",
          FileTypes.DS_CCS,
          fasta_ccs_file_type, is_distributed=True, nproc=1)
def run_bam2fasta_ccs(rtc):
    """
    Duplicate of run_bam2fasta, but with ConsensusReadSet as input.
    """
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0])


consensus_gz_ftype = OutputFileType(FileTypes.GZIP.file_type_id,
                                    "fastq_split_gz",
                                    "Consensus Amplicons",
                                    "Consensus amplicons in FASTQ format, split by barcode",
                                    "consensus_fastq")
chimera_gz_ftype = OutputFileType(FileTypes.GZIP.file_type_id,
                                  "fastq_split_gz",
                                  "Chimeric/Noise Sequences by barcode",
                                  "Chimeric and noise sequences in FASTQ format, split by barcode",
                                  "chimera_fastq")

@registry("split_laa_fastq", "0.1.0",
          (FileTypes.FASTQ, FileTypes.FASTQ),
          (consensus_gz_ftype, chimera_gz_ftype),
          is_distributed=True, nproc=1)
def _run_split_laa_fastq(rtc):
    # XXX a bit of a hack to support unique file names for the FASTQ tarballs
    return max(split_laa_fastq_archived(rtc.task.input_files[0],
                                        rtc.task.output_files[0]),
               split_laa_fastq_archived(rtc.task.input_files[1],
                                        rtc.task.output_files[1]))


@registry("contigset2fasta", "0.1.0",
          FileTypes.DS_CONTIG,
          FileTypes.FASTA,
          is_distributed=True, nproc=1)
def contigset_to_fasta(rtc):
    with ContigSet(rtc.task.input_files[0]) as ds_in:
        if len(ds_in.externalResources) != 1:
            raise ValueError("This task assumes that the ContigSet contains "+
                             "only a single FASTA file.")
        file_name = ds_in.externalResources[0].resourceId
        os.symlink(file_name, rtc.task.output_files[0])
    return 0


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
