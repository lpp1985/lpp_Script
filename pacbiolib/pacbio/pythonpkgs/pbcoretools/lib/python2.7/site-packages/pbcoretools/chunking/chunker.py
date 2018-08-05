import argparse
import os
import sys
import logging
import math

import xml.etree.ElementTree as ET

from pbcore.io import SubreadSet, HdfSubreadSet, AlignmentSet
from pbcommand.validators import validate_file, validate_fofn, fofn_to_files
from pbcommand.models.report import Report, Attribute
from pbcommand.models import PipelineChunk, FileTypes
from pbcommand.common_options import add_debug_option
from pbcommand.cli import get_default_argparser
from pbcommand.utils import compose
from pbcommand.cli.utils import main_runner_default, subparser_builder

import pbcoretools.chunking.chunk_utils as CU


log = logging.getLogger(__name__)

__version__ = "0.1.0"


def chunker_by_max_nchunks(alist, chunk_size):
    """Limit the individual size of each chunk"""
    if chunk_size > len(alist):
        n = len(alist)
    else:
        n = int(len(alist) / chunk_size + 1)
    return [alist[i:i + n] for i in range(0, len(alist), n)]


def chunker_by_max_chunksize(alist, max_nchunks):
    """Limit the max total number of chunks"""
    x = int(math.ceil(len(alist) / max_nchunks + 1))
    return chunker_by_max_nchunks(alist, x)


class Constants(object):
    FOFN_REPORT_ID = "fofn_chunk_report"
    FOFN_ATTRIBUTE_ID = "fofn_nchunks"

    CHUNK_KEY_ALNSET = "$chunk.alignmentset_id"
    CHUNK_KEY_SUBSET = "$chunk.subreadset_id"
    CHUNK_KEY_HDFSET = "$chunk.hdf5subreadset_id"
    CHUNK_KEY_REFSET = "$chunk.reference_id"
    CHUNK_KEY_FOFN = "$chunk.fofn_id"
    CHUNK_KEY_MOVIE_FOFN = "$chunk.movie_id"
    CHUNK_KEY_RGN_FOFN = "$chunk.rgn_id"
    CHUNK_KEY_FASTA = "$chunk.fasta_id"
    CHUNK_KEY_FASTQ = "$chunk.fastq_id"
    CHUNK_KEY_CSV = "$chunk.csv_id"

    CHUNK_KEY_CONTIG_ID = "$chunk.contig_id"
    CHUNK_KEY_CONTIG_NAME = "$chunk.contig_name_id"


def to_report(report_id, attribute_id, value):
    a = Attribute(attribute_id, value)
    r = Report(report_id, attributes=[a])
    return r


def write_report(report_id, attribute_id, value, report_file):
    r = to_report(report_id, attribute_id, value)
    r.write_json(report_file)
    return True


def nchunk_fofn(input_file, max_chunks):
    input_files = fofn_to_files(input_file)
    nchunks = min(len(input_files), max_chunks)
    return nchunks


def fofn_to_chunks(fofn):
    files = fofn_to_files(fofn)
    chunks = []
    for i, f in enumerate(files):
        chunk_id = "chunk-{i}".format(i=i)
        _d = {Constants.CHUNK_KEY_FOFN: f}
        p = PipelineChunk(chunk_id, **_d)
        chunks.append(p)
    return chunks


def add_max_nchunks_option(p):
    p.add_argument('--max-total-chunks', type=int, default=16, help="Chunk into X chunks.")
    return p


def _add_chunk_output_dir_option(p):
    p.add_argument('--output-dir', type=str, required=False, default=os.getcwd(),
                   help="Root directory to write chunked files to")
    return p


def _add_input_file_option(file_id, type_, help_):
    def _f(p):
        p.add_argument(file_id, type=type_, help=help_)
        return p
    return _f


def validate_external_resources(ds):
    not_found_paths = [r for r in ds.toExternalFiles() if not os.path.exists(r)]

    if not_found_paths:
        xs = "\n".join(not_found_paths)
        # raising this make it work better with argparse, instead of having
        # a giant useless stacktrace
        raise argparse.ArgumentTypeError("Invalid Dataset(s). {p} \nUnable to find {n} external Resources\n{x}".format(x=xs, p=ds.fileNames, n=len(not_found_paths)))

    return ds


def validate_external_non_empty_resources(ds):
    if ds.toExternalFiles():
        return validate_external_resources(ds)
    else:
        raise argparse.ArgumentTypeError("Invalid DataSet(s) {p}. \nMust have at least one external resource".format(p=ds.fileNames))


def _validate_dataset(dataset_class):
    def wrapper(path):
        try:
            ds = dataset_class(path)
            # this should also validate index files and nested external resources
            validate_external_non_empty_resources(ds)
            return path
        except Exception as e:
            # If this raises, the error is a mangled from argparse and yields
            # a cryptic error message
            raise argparse.ArgumentTypeError("Invalid DataSet(s) {p} {e}".format(p=path, e=e))

    return wrapper


def _validate_xml(path):
    try:
        p = ET.parse(path)
        _ = p.getroot()
        return path
    except ET.ParseError as e:
        raise argparse.ArgumentTypeError("Invalid XML file {p} {e}".format(p=path, e=e))

validate_xml_file = compose(_validate_xml, validate_file)

validate_subreadset = compose(_validate_dataset(SubreadSet), validate_xml_file)
valdiate_hdfsubreadset = compose(_validate_dataset(HdfSubreadSet), validate_xml_file)
validate_alignmentset = compose(_validate_dataset(AlignmentSet), validate_xml_file)


# These are really 'options', but keeping the naming convention consistent
add_input_fofn_option = _add_input_file_option('input_fofn', validate_fofn, "Path to input.fofn (File of File names)")
add_input_fasta_option = _add_input_file_option('fasta', validate_file, "Path to Fasta file.")
add_input_fasta_reference_option = _add_input_file_option('fasta', validate_file, "Path to PacBio Reference Entry Fasta file.")
add_input_fastq_option = _add_input_file_option('fastq', validate_file, "Path to Fastq file")
add_input_alignmentset_option = _add_input_file_option('alignmentset', validate_alignmentset, "Path to AlignmentSet XML file")
add_input_hdfsubreadset_option = _add_input_file_option('hdfsubreadset', valdiate_hdfsubreadset, "Path to HdfSubreadSet XML file")
add_input_subreadset_option = _add_input_file_option('subreadset', validate_subreadset, "Path to SubreadSet XML file")
add_input_csv_option = _add_input_file_option('csv', validate_file, "Path to CSV")
add_output_chunk_json_report_option = _add_input_file_option('chunk_report_json', str, "Path to chunked JSON output")


def _add_common_chunk_options(p):
    # Order matters!
    add_debug_option(p)
    add_max_nchunks_option(p)
    p = _add_chunk_output_dir_option(p)
    p = add_output_chunk_json_report_option(p)
    return p


def _add_chunk_fofn_options(p):
    p = add_input_fofn_option(p)
    p = _add_common_chunk_options(p)
    return p


def _args_chunk_fofn(args):
    fofn_files = fofn_to_files(args.input_fofn)
    log.info("read in fofn with {n} files.".format(n=len(fofn_files)))
    chunks = CU.write_grouped_fofn_chunks(fofn_files, args.max_total_chunks, args.output_dir, args.chunk_report_json)
    log.debug("Converted {x} Fofn into {n} chunks. Write chunks to {f}".format(n=len(chunks), f=args.chunk_report_json, x=len(fofn_files)))
    return 0


def _add_chunk_fasta_options(p):
    p = add_input_fasta_option(p)
    p = _add_common_chunk_options(p)
    return p


def _args_run_chunk_fasta(args):
    return CU.write_fasta_chunks_to_file(args.chunk_report_json, args.fasta, args.max_total_chunks, args.output_dir, "chunk_fa", 'fasta')


def _add_chunk_fastq_options(p):
    add_input_fastq_option(p)
    _add_common_chunk_options(p)
    return p


def _args_run_chunk_fastq(args):
    return CU.write_fastq_chunks_to_file(args.chunk_report_json, args.fasta, args.max_total_chunks, args.output_dir, "chunk_fq", 'fastq')


def _add_chunk_alignmentset_options(p):
    add_input_alignmentset_option(p)
    add_input_fasta_reference_option(p)
    _add_common_chunk_options(p)
    return p


def _args_run_chunk_alignmentset(args):
    return CU.write_alignmentset_chunks_to_file(args.chunk_report_json,
                                                args.alignmentset,
                                                args.fasta,
                                                args.max_total_chunks,
                                                args.output_dir,
                                                "chunk_alignmentset",
                                                FileTypes.DS_ALIGN.ext)


def _args_run_chunk_subreadset(args):
    return CU.write_subreadset_chunks_to_file(args.chunk_report_json,
                                              args.subreadset,
                                              args.fasta,
                                              args.max_total_chunks,
                                              args.output_dir,
                                              "chunk_subreadset",
                                              FileTypes.DS_SUBREADS.ext)


def _add_chunk_hdfsubreadset_options(p):
    add_input_hdfsubreadset_option(p)
    _add_common_chunk_options(p)
    return p


def _add_chunk_subreadset_options(p):
    add_input_subreadset_option(p)
    add_input_fasta_reference_option(p)
    _add_common_chunk_options(p)
    return p


def _args_run_chunk_hdfsubreadset(args):
    return CU.write_hdfsubreadset_chunks_to_file(args.chunk_report_json,
                                                 args.hdfsubreadset,
                                                 args.max_total_chunks,
                                                 args.output_dir,
                                                 "chunk_hdfsubreadset",
                                                 FileTypes.DS_SUBREADS_H5.ext)


def _add_chunk_csv_options(p):
    p = add_input_csv_option(p)
    p = _add_common_chunk_options(p)
    return p


def _args_run_chunk_csv(args):
    return CU.write_csv_chunks_to_file(args.chunk_report_json, args.csv, args.max_total_chunks, args.output_dir, "chunk", "csv")


def get_parser():
    desc = "Tool to create Chunk json files."
    p = get_default_argparser(__version__, desc)

    sp = p.add_subparsers(help="Subparser Commands")

    def builder(sid_, help_, opt_func_, exe_func_):
        return subparser_builder(sp, sid_, help_, opt_func_, exe_func_)

    builder("fofn", "Create a generic chunk.json from a FOFN.", _add_chunk_fofn_options, _args_chunk_fofn)

    builder("fasta", "Create a chunk.json from a Fasta file", _add_chunk_fasta_options, _args_run_chunk_fasta)

    builder("fastq", "Create a chunk.json from a Fastq file", _add_chunk_fastq_options, _args_run_chunk_fastq)

    builder("alignmentset",
            "Create a chunk.json from an AlignmentSet XML file",
            _add_chunk_alignmentset_options, _args_run_chunk_alignmentset)

    builder("hdfsubreadset",
            "Create a chunk.json from an HdfSubreadSet XML file",
            _add_chunk_hdfsubreadset_options, _args_run_chunk_hdfsubreadset)

    builder("subreadset",
            "Create a chunk.json from an SubreadSet XML file",
            _add_chunk_subreadset_options, _args_run_chunk_subreadset)

    builder("csv", "Create a chunk.json CSV from a CSV file", _add_chunk_csv_options, _args_run_chunk_csv)

    return p


def main(argv=None):
    argv_ = sys.argv if argv is None else argv
    parser = get_parser()
    return main_runner_default(argv_[1:], parser, log)
