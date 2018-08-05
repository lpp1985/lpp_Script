"""A simple script to scatter a (filtered) subreadset.

The scatterer in pbcoretools requires a referenceset. (Not sure why.)
So we copy/paste the code here.
"""
from pbcore.io import (SubreadSet)
import argparse
import logging
import os
import sys
log = logging.getLogger(__name__)

def run(subreadset, fofn):
    dir_name = os.getcwd()
    maxChunks = 0
    dset = SubreadSet(subreadset, strict=True)
    fns = dset.toFofn()
    import pprint
    log.info('resources in {!r}:\n{}'.format(subreadset, pprint.pformat(fns)))
    nrecs = len(dset)
    # HG with 70x coverage => 200G bases total
    ts = 50000 # @ 20k/read => 1G bases, ~300MB .gz => ~200 chunks for Human
    ts = 500000 # @ 20k/read => 10G bases, ~3GB .gz => ~20 chunks for Human
    # and we expect about 7-10min per chunk.
    chunks = nrecs // ts
    log.info('num_chunks={:g} ({:g} / {:g})'.format(chunks, nrecs, ts))
    log.info('Splitting with dset.split(zmws=False, chunks={}, ignoreSubDatasets=True, maxChunks={},)'.format(
        chunks, maxChunks))
    dset_chunks = dset.split(zmws=False, chunks=chunks, ignoreSubDatasets=True, maxChunks=maxChunks,
            updateCounts=False,
            #targetSize=1, breakContigs=True
    )

    chunk_fns = []
    for i, dset in enumerate(dset_chunks):
        chunk_name = 'chunk_{:03d}.subreadset.xml'.format(i) # TODO: 02
        chunk_fn = os.path.join(dir_name, chunk_name)
        dset.updateCounts()
        dset.write(chunk_fn, validate=False) # , relPaths=True
        chunk_fns.append(chunk_fn)
    with open(fofn, 'w') as ofs:
        for fn in chunk_fns:
            ofs.write('{}\n'.format(fn))
    log.info('Wrote {} chunks into "{}"'.format(len(dset_chunks), fofn))

def main(argv=sys.argv):
    description = """Scatter subreadsets prior to bam2fasta.
"""
    epilog = """
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--logging',
    #        help='.ini or .json config file for Python logging module')
    parser.add_argument('subreadset',
        help='Input subreadset XML filename. Can be filtered.')
    parser.add_argument('fofn',
        help='File of File Names of scattered alignmentsets. Relative paths.')
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    #from logging_tree import printout
    #printout()
    main(sys.argv)
