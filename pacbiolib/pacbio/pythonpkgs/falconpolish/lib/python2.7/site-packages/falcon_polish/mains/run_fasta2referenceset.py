"""A simple wrapper for run_bam_to_fasta
in pbcoretools.tasks.converters
"""
from __future__ import absolute_import
from ..sys import system
#from pbcoretools.tasks.converters import run_fasta_to_reference, run_fasta_to_referenceset
from pbcore.io import (ContigSet, ReferenceSet)
import argparse
import logging
import sys
import os.path as op

log = logging.getLogger(__name__)

def run_fasta_to_referenceset(input_file_name, output_file_name, prog):
    """Copied from pbsmrtpipe/pb_tasks/pacbio.py:run_fasta_to_referenceset()
    """
    args = ['dataset', "create", "--type ReferenceSet", "--generateIndices",
            output_file_name, input_file_name]
    system(" ".join(args))

def run_fasta_to_reference(input_file_name, output_file_name,
                           organism, reference_name,
                           ploidy):
    """Copied from pbcoretools/tasks/converters.py:run_fasta_to_reference()
    """
    ds_in = ContigSet(input_file_name)
    if len(ds_in.externalResources) > 1:
        raise TypeError("Only a single FASTA file is supported as input.")
    fasta_file_name = ds_in.externalResources[0].resourceId
    output_dir_name = op.dirname(output_file_name)
    args = [
        "fasta-to-reference",
        "--organism", organism,
        "--ploidy", ploidy,
        "--debug",
        fasta_file_name,
        output_dir_name,
        reference_name
    ]
    log.info(" ".join(args))
    system(" ".join(args))
    ref_file = op.join(output_dir_name, reference_name, "referenceset.xml")
    assert op.isfile(ref_file)
    with ReferenceSet(ref_file, strict=True) as ds_ref:
        ds_ref.makePathsAbsolute()
        log.info("saving final ReferenceSet to {f!r}".format(f=output_file_name))
        ds_ref.write(output_file_name)

def run(fasta, ref):
    try:
        # This uses Python + BAM library.
        run_fasta_to_referenceset(fasta, ref, 'dataset')
        return
    except Exception:
        log.exception('We will try another someting else.')

    try:
        # This uses Python + BAM library.
        # the '.py' name difference will be resolved in pbdataset/pbcoretools, but
        # for now, work with either
        run_fasta_to_referenceset(fasta, ref, 'dataset.py')
        return
    except Exception:
        log.exception('We will try someting else.')
        raise

    try:
        # This uses pbscala and also runs sawriter.
        reference_name = op.splitext(op.basename(fasta))[0]
        organism = "unknown"
        ploidy = "haploid"
        run_fasta_to_reference(fasta, ref, organism=organism, reference_name=reference_name, ploidy=ploidy)
    except Exception:
        log.exception('Out of ideas.')
        raise


def main(argv=sys.argv):
    description = """Create referenceset XML from fasta.
"""
    epilog = """
There might extra files too. We use fasta-to-reference (from pbscala)
if available (which would also run sawriter).
Otheriwse, we use 'dataset create'.
The fasta might be copied, and the dataset should refer to it absolutely (I think).
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--logging',
    #        help='.ini or .json config file for Python logging module')
    parser.add_argument('fasta',
        help='Input fasta filename.')
    parser.add_argument('ref',
        help='Output referenceset XML filename.')
    args = parser.parse_args(argv[1:])
    log.info('RUNNING run_fasta2reference: {}'.format(repr(args)))
    run(**vars(args))

if __name__ == "__main__":
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    main(sys.argv)
