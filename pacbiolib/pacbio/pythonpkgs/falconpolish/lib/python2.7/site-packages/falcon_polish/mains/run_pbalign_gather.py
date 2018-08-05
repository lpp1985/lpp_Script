"""A simple wrapper for gather_alignmentset
from pbcoretools.chunking.gather
"""
from pbcoretools.chunking.gather import gather_alignmentset
import argparse
import logging
import os
import sys
log = logging.getLogger(__name__)

def run_gather_unmappeds(iunmappeds, ounmapped):
    # Not sure yet whether unmapped names can be repeated in multiple "chunks".
    # For now, we simple concatentate them all.
    with open(ounmapped, 'w') as ofs:
        unmapped_fns = open(iunmappeds).read().strip().split()
        for fn in unmapped_fns:
            with open(fn) as ifs:
                content = ifs.read()
                ofs.write(content)

def run_gather_alignments(idatasets, odataset):
    odir = os.path.dirname(odataset)
    if odir and not os.path.isdir(odir):
        # I think pbcoretools tend to expect the output dir to exist.
        os.makedirs(odir)
    dset_out_fn = odataset
    dset_fns = open(idatasets).read().strip().split()
    gather_alignmentset(dset_fns, dset_out_fn)

def run(iunmappeds, ounmapped, idatasets, odataset):
    run_gather_unmappeds(iunmappeds, ounmapped)
    run_gather_alignments(idatasets, odataset)

def main(argv=sys.argv):
    description = """Gather results of parallel pbalign (blasr).
"""
    epilog = """
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('iunmappeds',
        help='File of unmapped.txt filenames, whitespace delimited.')
    parser.add_argument('ounmapped',
        help='Unmapped read-names. (Someday we may want the actual reads.)')
    parser.add_argument('idatasets',
        help='File of input alignmentset XML filenames, whitespace delimited.')
    parser.add_argument('odataset',
        help='Output alignmentset XML filename.')
    args = vars(parser.parse_args(argv[1:]))
    run(**args)

if __name__ == "__main__":
    logging.basicConfig()
    log.setLevel(logging.INFO)
    main(sys.argv)
