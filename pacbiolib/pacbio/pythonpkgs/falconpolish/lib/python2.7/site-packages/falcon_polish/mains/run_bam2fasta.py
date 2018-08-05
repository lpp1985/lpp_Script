"""A simple wrapper for run_bam_to_fasta
in pbcoretools.tasks.converters
"""
from pbcoretools.tasks.converters import run_bam_to_fasta
import argparse
import logging
import sys

def main(argv=sys.argv):
    description = """Dump .fasta from .bam (identified from .xml dataset).
"""
    epilog = """
If you have only '.bam', you will first need to run 'dataset create --type ReferenceSet --generateIndices outdataset inbam*'
to generate the dataset.
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--logging',
    #        help='.ini or .json config file for Python logging module')
    parser.add_argument('in_file',
        help='Input dataset XML filename. (Must have .bam resource.)')
    parser.add_argument('out_file',
        help='Output fasta filename.')
    args = parser.parse_args(argv[1:])
    #logging.basicConfig()
    #log.warning('RUNNING filterbam: {}'.format(repr(args)))
    run_bam_to_fasta(args.in_file, args.out_file)

if __name__ == "__main__":
    main(sys.argv)
