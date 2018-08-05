"""A simple wrapper for variantCaller.
Part of GenomicConsensus.

Note: GenomicConsensus outputs .fasta instead of .contigset.xml, so
there is weird code in its resolved_tool_contract_runner() which substitutes
a filename and later converts .fasta to a contigset. That does not exist
for the non-pbcommand entry-point, so we copy that code here.
"""
from __future__ import absolute_import
from ..sys import system, say
from pbcore.io import ContigSet
import argparse
import logging
import sys
import traceback

def run(referenceset, fastq, gff, fasta, contigset, alignmentset, options, log_level):
    #'--log-file foo.log',
    #'--verbose',
    #'--debug', # requires 'ipdb'
    #'-j NWORKERS',
    #'--algorithm quiver',
    #'--diploid', # binary
    #'--minConfidence 40',
    #'--minCoverage 5',
    #'--alignmentSetRefWindows',
    cmd = "variantCaller --log-level {log_level} {options} --referenceFilename {referenceset} -o {fastq} -o {gff} -o {fasta} {alignmentset}"
    system(cmd.format(**locals()))
    try:
        say('Converting fasta {!r} to contigset {!r}'.format(fasta, contigset))
        # Convert to contigset.xml
        import pysam
        pysam.faidx(fasta)
        ds = ContigSet(fasta, strict=True)
        ds.write(contigset, relPaths=True)
        say('Successfully wrapped fasta {!r} in contigset {!r}'.format(fasta, contigset))
    except Exception:
        say(traceback.format_exc())
        say('Skipping conversion to contigset.')

def main(argv=sys.argv):
    description = """Run variantCaller. Result is *both* the consensus fasta *and* the consensus contigset XML. The headers will include enough info to stitch the fasta together when variantCaller is run on "scattered" pieces.
"""
    epilog = """If dependencies are missing, we skip the wrapping of fasta in contigset.
"""
    logging.basicConfig()
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--log-level', default='INFO',
            help='DEBUG or INFO')
    parser.add_argument('--options', default='--alignmentSetRefWindows',
            help='Extra options')
    parser.add_argument('alignmentset',
        help='Input alignmentset filename')
    parser.add_argument('referenceset',
        help='Input referenceset filename')
    parser.add_argument('fastq',
        help='Output polished-fastq filename')
    parser.add_argument('gff',
        help='Output variants-gff filename')
    parser.add_argument('fasta',
        help='Output polished-fasta filename')
    parser.add_argument('contigset',
        help='Output consensus contigset XML filename')
    # TODO(CD): When do we say 'polished', and when 'consensus'?
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    main(sys.argv)
