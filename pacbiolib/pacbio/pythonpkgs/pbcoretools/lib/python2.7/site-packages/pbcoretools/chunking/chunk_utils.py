import os
import logging
import math
import datetime
import functools
import csv

from pbcore.io import (FastaWriter, FastaReader, FastqReader, FastqWriter,
                       AlignmentSet, HdfSubreadSet, SubreadSet, ReferenceSet,
                       ConsensusReadSet, ContigSet, FastaRecord)
from pbcommand.pb_io.common import write_pipeline_chunks
from pbcommand.pb_io.report import fofn_to_report
from pbcommand.models import PipelineChunk


log = logging.getLogger(__name__)


class Constants(object):
    CHUNK_KEY_HDFSET = "$chunk.hdf5subreadset_id"
    CHUNK_KEY_SUBSET = "$chunk.subreadset_id"
    CHUNK_KEY_CCSSET = "$chunk.ccsset_id"
    CHUNK_KEY_ALNSET = "$chunk.alignmentset_id"
    CHUNK_KEY_CONTIGSET = "$chunk.contigset_id"
    CHUNK_KEY_CCS_ALNSET = "$chunk.ccs_alignmentset_id"
    CHUNK_KEY_FASTA = "$chunk.fasta_id"
    CHUNK_KEY_FASTQ = "$chunk.fastq_id"
    CHUNK_KEY_FOFN = "$chunk.fofn_id"
    CHUNK_KEY_FOFN_REPORT = '$chunk.fofn_report_id'
    CHUNK_KEY_CSV = "$chunk.csv_id"


def write_chunks_to_json(chunks, chunk_file):
    log.debug("Wrote {n} chunks to {f}.".format(n=len(chunks), f=chunk_file))
    write_pipeline_chunks(
        chunks, chunk_file, "Chunks written at {d}".format(d=datetime.datetime.now()))
    return 0


def _to_grouped_items_by_max_total_chunks(items, max_total_chunks):
    """Group items by the max total number of chunks to be created"""
    nitems = len(items)

    max_total_chunks = max(1, min(nitems, max_total_chunks))

    grouped_items = []

    n = int(math.ceil(float(nitems)) / max_total_chunks)
    for i in xrange(max_total_chunks):
        if i != max_total_chunks - 1:
            cs = items[i * n: (i + 1) * n]
        else:
            cs = items[i * n:]
        grouped_items.append(cs)

    return grouped_items


def _to_grouped_items_by_max_size_per_item(items, max_chunks_per_item):
    """Break a list into several chunks where the max size of each item"""
    nitems = len(items)

    chunks = []
    if nitems <= max_chunks_per_item:
        # Each item is it's own group
        for i in items:
            chunks.append(i)
    else:
        n = int(math.ceil(float(nitems) / max_chunks_per_item))
        for i in xrange(n):
            if i != max_chunks_per_item - 1:
                cs = items[i * n:n * (i + 1)]
            else:
                cs = items[i * n:]
            chunks.append(cs)

    return chunks


def __get_nrecords_from_reader(reader):
    n = 0
    for _ in reader:
        n += 1
    return n


def write_chunked_csv(chunk_key, csv_path, max_total_nchunks, dir_name, base_name, ext):
    # This needs to have an ignore emtpy file mode

    with open(csv_path, 'r') as csv_fh:
        reader = csv.DictReader(csv_fh)
        field_names = reader.fieldnames
        nrecords = __get_nrecords_from_reader(reader)

    max_total_nchunks = max(1, min(nrecords, max_total_nchunks))

    n = int(math.ceil(float(nrecords)) / max_total_nchunks)

    nchunks = 0
    with open(csv_path, 'r') as csv_fh:
        reader = csv.DictReader(csv_fh)

        it = iter(reader)
        for i in xrange(max_total_nchunks):

            chunk_id = "_".join([base_name, str(nchunks)])
            chunk_name = ".".join([chunk_id, ext])
            nchunks += 1
            nchunk_records = 0
            csv_chunk_path = os.path.join(dir_name, chunk_name)

            with open(csv_chunk_path, 'w+') as csv_chunk_fh:
                writer = csv.DictWriter(csv_chunk_fh, field_names)
                writer.writeheader()
                if i != max_total_nchunks:
                    for _ in xrange(n):
                        nchunk_records += 1
                        writer.writerow(it.next())
                else:
                    for x in it:
                        nchunk_records += 1
                        writer.writerow(x)

            d = dict(nrecords=nchunk_records)
            d[chunk_key] = os.path.abspath(csv_chunk_path)
            c = PipelineChunk(chunk_id, **d)
            yield c


def write_csv_chunks_to_file(chunk_file, csv_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext):
    chunks = list(write_chunked_csv(Constants.CHUNK_KEY_CSV, csv_path,
                                    max_total_chunks, dir_name, chunk_base_name, chunk_ext))
    return write_chunks_to_json(chunks, chunk_file)


def _get_nrecords_from_fastx(fastx_reader):
    n = 0
    for _ in fastx_reader:
        n += 1
    return n


def write_pbcore_records(pbcore_writer_class, records, file_name):
    n = 0
    with pbcore_writer_class(file_name) as w:
        for record in records:
            w.writeRecord(record)
            n += 1

    log.debug("Completed writing {n} records".format(n=n))


def write_contigset_records(pbcore_writer_class, records, file_name):
    """
    Writes the Chunked fasta files and Writes a ContigSet.xml

    Filename has contigset.xml

    """
    fasta_file_name = ".".join(file_name.split(".")[:-2]) + ".fasta"
    write_pbcore_records(pbcore_writer_class, records, fasta_file_name)
    log.debug("Writing ContigSet XML to {f}".format(f=file_name))
    ds = ContigSet(fasta_file_name)
    ds.write(file_name)


def __to_chunked_fastx_files(write_records_func, pbcore_reader_class, pbcore_writer_class,
                             chunk_key, input_file, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=None):
    """Convert a Fasta/Fasta file to a chunked list of files

    :param write_records_func: Func(writer_class, records, file_name)
    :param pbcore_reader_class Pbcore IO Reader
    :param pbcore_writer_class Pbcore IO Writer
    :param chunk_key: Chunk key to assign to PipelineChunk
    :param input_file: Path to input file


    """
    # grab the number of records so we can chunk it
    with pbcore_reader_class(input_file) as f:
        nrecords = __get_nrecords_from_reader(f)

    max_total_nchunks = max(1, min(nrecords, max_total_nchunks))

    n_per_chunk = int(math.ceil(float(nrecords) / max_total_nchunks))

    log.info("Found {n} total records. Max total chunks {m}. Splitting into chunks of approximately {x} records each".format(n=nrecords, x=n_per_chunk, m=max_total_nchunks))
    nchunks = 0
    with pbcore_reader_class(input_file) as r:
        it = iter(r)
        for i in xrange(max_total_nchunks):
            records = []

            chunk_id = "_".join([base_name, str(nchunks)])
            chunk_name = ".".join([chunk_id, ext])
            nchunks += 1
            fasta_chunk_path = os.path.join(dir_name, chunk_name)

            if i != max_total_nchunks:
                n_left = nrecords - (n_per_chunk * i)
                if n_left < 0 or (n_left == 0 and nchunks != 1):
                    break
                for _ in xrange(min(n_per_chunk, n_left)):
                    records.append(it.next())
            else:
                for x in it:
                    records.append(x)

            write_records_func(pbcore_writer_class, records, fasta_chunk_path)
            total_bases = sum(len(r.sequence) for r in records)
            d = dict(total_bases=total_bases, nrecords=len(records))
            d[chunk_key] = os.path.abspath(fasta_chunk_path)
            if extra_chunk_keys is not None:
                d.update(extra_chunk_keys)
            c = PipelineChunk(chunk_id, **d)
            yield c


_to_chunked_fastx_files = functools.partial(__to_chunked_fastx_files, write_pbcore_records)


def to_chunked_fasta_files(fasta_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=None):
    return _to_chunked_fastx_files(FastaReader, FastaWriter, Constants.CHUNK_KEY_FASTA, fasta_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=extra_chunk_keys)


def to_chunked_fastq_files(fastq_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=None):
    return _to_chunked_fastx_files(FastqReader, FastqWriter, Constants.CHUNK_KEY_FASTQ, fastq_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=extra_chunk_keys)


def to_chunked_contigset_files(dataset_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=None):
    return __to_chunked_fastx_files(write_contigset_records, ContigSet, FastaWriter, Constants.CHUNK_KEY_CONTIGSET, dataset_path, max_total_nchunks, dir_name, base_name, ext, extra_chunk_keys=extra_chunk_keys)


def _write_fasta_chunks_to_file(to_chunk_fastx_file_func, chunk_file, fastx_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext, extra_chunk_keys=None):
    chunks = list(to_chunk_fastx_file_func(
        fastx_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext,
        extra_chunk_keys=extra_chunk_keys))
    write_chunks_to_json(chunks, chunk_file)
    return 0


def write_fasta_chunks_to_file(chunk_file, fasta_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext):
    return _write_fasta_chunks_to_file(to_chunked_fasta_files, chunk_file, fasta_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext)


def write_fastq_chunks_to_file(chunk_file, fasta_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext):
    return _write_fasta_chunks_to_file(to_chunked_fastq_files, chunk_file, fasta_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext)


def write_contigset_chunks_to_file(chunk_file, dataset_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext, extra_chunk_keys=None):
    return _write_fasta_chunks_to_file(to_chunked_contigset_files, chunk_file, dataset_path, max_total_chunks, dir_name, chunk_base_name, chunk_ext, extra_chunk_keys=extra_chunk_keys)


def write_alignmentset_chunks_to_file(chunk_file, alignmentset_path,
                                      reference_path, max_total_chunks,
                                      dir_name, chunk_base_name, chunk_ext):
    chunks = list(to_chunked_alignmentset_files(alignmentset_path,
                                                reference_path,
                                                max_total_chunks,
                                                Constants.CHUNK_KEY_ALNSET,
                                                dir_name, chunk_base_name,
                                                chunk_ext))
    write_chunks_to_json(chunks, chunk_file)
    return 0


def to_chunked_alignmentset_files(alignmentset_path, reference_path,
                                  max_total_nchunks, chunk_key, dir_name,
                                  base_name, ext):
    dset = AlignmentSet(alignmentset_path, strict=True)
    dset_chunks = dset.split(contigs=True, maxChunks=max_total_nchunks,
                             breakContigs=True)

    # sanity checking
    reference_set = ReferenceSet(reference_path, strict=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        d['$chunk.reference_id'] = reference_path
        c = PipelineChunk(chunk_id, **d)
        yield c


def write_subreadset_chunks_to_file(chunk_file, subreadset_path,
                                    reference_path,
                                    max_total_chunks, dir_name,
                                    chunk_base_name, chunk_ext):
    chunks = list(to_chunked_subreadset_files(subreadset_path,
                                              reference_path,
                                              max_total_chunks,
                                              Constants.CHUNK_KEY_SUBSET,
                                              dir_name, chunk_base_name,
                                              chunk_ext))
    write_chunks_to_json(chunks, chunk_file)
    return 0


def write_ccsset_chunks_to_file(chunk_file, ccsset_path,
                                reference_path,
                                max_total_chunks, dir_name,
                                chunk_base_name, chunk_ext):
    chunks = list(to_chunked_ccsset_files(ccsset_path,
                                          reference_path,
                                          max_total_chunks,
                                          Constants.CHUNK_KEY_CCSSET,
                                          dir_name, chunk_base_name,
                                          chunk_ext))
    write_chunks_to_json(chunks, chunk_file)
    return 0


def _to_chunked_dataset_files(dataset_type, dataset_path, reference_path,
                              max_total_nchunks, chunk_key, dir_name,
                              base_name, ext):
    dset = dataset_type(dataset_path, strict=True)
    dset_chunks = dset.split(zmws=True, chunks=max_total_nchunks, ignoreSubDatasets=True)
    d = {}

    # sanity checking
    reference_set = ReferenceSet(reference_path)
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        d['$chunk.reference_id'] = reference_path
        c = PipelineChunk(chunk_id, **d)
        yield c

to_chunked_subreadset_files = functools.partial(_to_chunked_dataset_files,
    SubreadSet)
to_chunked_ccsset_files = functools.partial(_to_chunked_dataset_files,
    ConsensusReadSet)


def _to_barcode_chunked_dataset_files(dataset_type, dataset_path,
                                      max_total_nchunks, chunk_key, dir_name,
                                      base_name, ext):
    """
    Similar to to_chunked_subreadset_files, but chunks reads by barcode lists.
    """
    dset = dataset_type(dataset_path, strict=True)
    dset_chunks = dset.split(chunks=max_total_nchunks, barcodes=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        c = PipelineChunk(chunk_id, **d)
        yield c

to_barcode_chunked_subreadset_files = functools.partial(
    _to_barcode_chunked_dataset_files, SubreadSet)
to_barcode_chunked_ccsset_files = functools.partial(
    _to_barcode_chunked_dataset_files, ConsensusReadSet)


def _write_dataset_barcode_chunks_to_file(chunk_func, chunk_key, chunk_file,
                                          dataset_path, max_total_chunks,
                                          dir_name, chunk_base_name,
                                          chunk_ext):
    """
    Similar to write_subreadset_chunks_to_file, but chunks reads (subread or
    CCS) by barcodes for input to pblaa.
    """
    chunks = list(chunk_func(
        dataset_path=dataset_path,
        max_total_nchunks=max_total_chunks,
        chunk_key=chunk_key,
        dir_name=dir_name,
        base_name=chunk_base_name,
        ext=chunk_ext))
    write_chunks_to_json(chunks, chunk_file)
    return 0

write_subreadset_barcode_chunks_to_file = functools.partial(
    _write_dataset_barcode_chunks_to_file, to_barcode_chunked_subreadset_files,
    Constants.CHUNK_KEY_SUBSET)
write_ccsset_barcode_chunks_to_file = functools.partial(
    _write_dataset_barcode_chunks_to_file, to_barcode_chunked_ccsset_files,
    Constants.CHUNK_KEY_CCSSET)


def _to_zmw_chunked_dataset_files(dataset_type, dataset_path,
                                  max_total_nchunks, chunk_key, dir_name,
                                  base_name, ext, extra_chunk_keys=None):
    """
    Similar to to_chunked_subreadset_files, but chunks reads by ZMW ranges
    for input to pbccs or pbtranscript.
    """
    dset = dataset_type(dataset_path, strict=True)
    dset_chunks = dset.split(chunks=max_total_nchunks, zmws=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        if extra_chunk_keys is not None:
            d.update(extra_chunk_keys)
        c = PipelineChunk(chunk_id, **d)
        yield c

to_zmw_chunked_subreadset_files = functools.partial(
    _to_zmw_chunked_dataset_files, SubreadSet)
to_zmw_chunked_ccsset_files = functools.partial(
    _to_zmw_chunked_dataset_files, ConsensusReadSet)


def _to_bam_chunked_dataset_files(dataset_type, dataset_path,
                                  max_total_nchunks, chunk_key, dir_name,
                                  base_name, ext, extra_chunk_keys=None):
    """
    Similar to to_chunked_subreadset_files, but chunks reads by ZMW ranges
    for input to pbccs or pbtranscript.
    """
    dset = dataset_type(dataset_path, strict=True)
    dset_chunks = dset.split(chunks=max_total_nchunks, zmws=False,
                             ignoreSubDatasets=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        if extra_chunk_keys is not None:
            d.update(extra_chunk_keys)
        c = PipelineChunk(chunk_id, **d)
        yield c

to_bam_chunked_subreadset_files = functools.partial(
    _to_bam_chunked_dataset_files, SubreadSet)


def _write_dataset_chunks_to_file(chunk_func, chunk_key, chunk_file,
                                  dataset_path, max_total_chunks,
                                  dir_name, chunk_base_name, chunk_ext,
                                  extra_chunk_keys=None):
    """
    Similar to write_subreadset_chunks_to_file, but chunks reads (subread or
    CCS) by ZMW ranges for input to pbccs.
    """
    chunks = list(chunk_func(
        dataset_path=dataset_path,
        max_total_nchunks=max_total_chunks,
        chunk_key=chunk_key,
        dir_name=dir_name,
        base_name=chunk_base_name,
        ext=chunk_ext,
        extra_chunk_keys=extra_chunk_keys))
    write_chunks_to_json(chunks, chunk_file)
    return 0

write_subreadset_zmw_chunks_to_file = functools.partial(
    _write_dataset_chunks_to_file, to_zmw_chunked_subreadset_files,
    Constants.CHUNK_KEY_SUBSET)
write_ccsset_zmw_chunks_to_file = functools.partial(
    _write_dataset_chunks_to_file, to_zmw_chunked_ccsset_files,
    Constants.CHUNK_KEY_CCSSET)
write_subreadset_bam_chunks_to_file = functools.partial(
    _write_dataset_chunks_to_file, to_bam_chunked_subreadset_files,
    Constants.CHUNK_KEY_SUBSET)


def write_hdfsubreadset_chunks_to_file(chunk_file, hdfsubreadset_path,
                                       max_total_chunks, dir_name,
                                       chunk_base_name, chunk_ext):
    chunks = list(to_chunked_hdfsubreadset_files(hdfsubreadset_path,
                                                 max_total_chunks,
                                                 Constants.CHUNK_KEY_HDFSET,
                                                 dir_name, chunk_base_name,
                                                 chunk_ext))
    write_chunks_to_json(chunks, chunk_file)
    return 0


def to_chunked_hdfsubreadset_files(hdfsubreadset_path, max_total_nchunks,
                                   chunk_key, dir_name, base_name, ext):
    dset = HdfSubreadSet(hdfsubreadset_path, strict=True)
    dset_chunks = dset.split(chunks=max_total_nchunks, ignoreSubDatasets=True)
    d = {}
    for i, dset in enumerate(dset_chunks):
        chunk_id = '_'.join([base_name, str(i)])
        chunk_name = '.'.join([chunk_id, ext])
        chunk_path = os.path.join(dir_name, chunk_name)
        dset.write(chunk_path)
        d[chunk_key] = os.path.abspath(chunk_path)
        c = PipelineChunk(chunk_id, **d)
        yield c


def write_fofn(paths, fofn_path):
    with open(fofn_path, 'w') as w:
        w.write("\n".join(paths))


def to_chunked_grouped_fofn(fofn_groups, chunk_id_prefix, fofn_chunk_key, report_chunk_key, chunk_dir_name):
    """

    :param fofn_groups: A list of FofnGroups
    :param chunk_id_prefix: Prefix used to create the chunk key and grouped
    Fofn files
    :param fofn_chunk_key: Value of the chunk key to write to the chunk file (e.g., $chunk.my_key)
    :param chunk_dir_name: Directory where the Grouped Fofn files will be
    written to
    :return: list of pipeline chunks
    """

    chunks = []
    for i, fofn_group in enumerate(fofn_groups):
        chunk_id = "_".join([chunk_id_prefix, str(i)])
        fofn_group_name = "".join([chunk_id, ".fofn"])
        fofn_group_path = os.path.join(chunk_dir_name, fofn_group_name)

        write_fofn(fofn_group, fofn_group_path)

        # Write the companion fofn metadata report
        fofn_report_name = "".join([chunk_id, "_report", '.json'])
        fofn_report_path = os.path.join(chunk_dir_name, fofn_report_name)
        fofn_report = fofn_to_report(len(fofn_group))
        fofn_report.write_json(fofn_report_path)

        d = dict(nfofns=len(fofn_group))
        d[fofn_chunk_key] = fofn_group_path
        d[report_chunk_key] = fofn_report_path

        c = PipelineChunk(chunk_id, **d)
        chunks.append(c)

    return chunks


def write_grouped_fofn_chunks(fofn_files, max_total_chunks, chunk_dir_name, chunk_json_path):

    fofn_groups = _to_grouped_items_by_max_total_chunks(
        fofn_files, max_total_chunks)

    chunks = to_chunked_grouped_fofn(
        fofn_groups, 'fofn_group', Constants.CHUNK_KEY_FOFN, Constants.CHUNK_KEY_FOFN_REPORT, chunk_dir_name)

    write_pipeline_chunks(
        chunks, chunk_json_path, "Group Fofn created at {d}".format(d=datetime.datetime.now()))

    return chunks
