"""A simple wrapper for
    pbcoretools.tasks.gather_contigs stuff
and
    pbcoretools.chunking.gather.gather_fastq_contigset
"""
from pbcore.io import (ContigSet)
from pbcoretools.chunking.gather import gather_fastq_contigset
import argparse
import logging
import os
import sys
log = logging.getLogger(__name__)

def __gather_contigset(input_files, output_file, new_resource_file):
    """Copied from pbcoretools.chunking.gather:__gather_contigset()
    """
    skip_empty = True
    if skip_empty:
        _input_files = []
        for file_name in input_files:
            cs = ContigSet(file_name)
            if len(cs.toExternalFiles()) > 0:
                _input_files.append(file_name)
        input_files = _input_files
    tbr = ContigSet(*input_files)
    tbr.consolidate(new_resource_file)
    tbr.newUuid()
    tbr.write(output_file, relPaths=True)
    return output_file
def run_gc_gather(dset_fns, ds_out_fn):
    """Gather contigsets/fasta into 1 contigset/fasta.
    """
    # python -m pbcoretools.tasks.gather_contigs --rtc
    log.info('Gathering {!r} from chunks {!r}'.format(ds_out_fn, dset_fns))
    assert ds_out_fn.endswith('xml')
    new_resource_fn = os.path.splitext(ds_out_fn)[0] + '.fasta'
    __gather_contigset(dset_fns, ds_out_fn, new_resource_fn)
def run_gc_gather_fastq(fastq_fns, fastq_fn):
    """Also write fastq.contigset.xml
    Use the contigset implicitly.
    """
    log.info('gc_gather_fastq({!r}, {!r})'.format(fastq_fns, fastq_fn))
    assert fastq_fns, 'Empty list. gather_fastq_contigset() would produce nothing.'
    gather_fastq_contigset(fastq_fns, fastq_fn)
def run(ifastas, ifastqs, fasta, fastq):
    fasta_dset_fn = fasta
    fasta_dset_fns = open(ifastas).read().strip().split()
    run_gc_gather(fasta_dset_fns, fasta_dset_fn)

    fastq_fn = fastq
    fastq_fns = open(ifastqs).read().strip().split()
    run_gc_gather_fastq(fastq_fns, fastq_fn)
def main(argv=sys.argv):
    description = """Gather results of parallel GC. Both fasta and fastq (for now).
"""
    epilog = """
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ifastas',
        help='File of input fasta contigset XML filenames, whitespace delimited.')
    parser.add_argument('ifastqs',
        help='File of nput fastq filenames, whitespace delimited.')
    parser.add_argument('fasta',
        help='Output fasta contigset XML filename. (fasta resource will be created too.)')
    parser.add_argument('fastq',
        help='Output fastq filename. (contigset of fastq resource might be created too.)')
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.INFO)
    main(sys.argv)
