
# TODO modernize this or delete

"""Create a CSV file of Subread from the *.rgn.h5 file

TODO: An API to  *.rgn.h5 file should replace the raw calls to the *.h5 files.

TODO: Add a P*module level test for checking the that the total number of
subreads written should be checked against the fasta and fastq subread file.

# $Id: //depot/software/assembly/pbpy/bin/makeFilterSubreadSummary.py#9 $
# $Date: 2013/02/04 $
"""
import os
from pprint import pformat
import sys
import re
import time
import logging
import warnings
import argparse

import h5py
import numpy as np

from pbcommand.cli import pacbio_args_runner, \
    get_default_argparser_with_base_opts
from pbcommand.utils import setup_log

from pbreports.io.validators import validate_fofn, bas_fofn_to_bas_files


def to_version(major_version, perforce_str):
    # this should be util func.
    rx = re.compile(r'Change: (\d+)')
    v = rx.search(perforce_str).group(1)
    return '%s.%s' % (str(major_version), v)


MAJOR_VERSION = '1.81'
__revision__ = "$Change: 119000 $"
__version__ = to_version(MAJOR_VERSION, __revision__)

log = logging.getLogger(__name__)


def run(region_files, output_csv):
    """
    Return (bool): based on success of analyzing the rgn files
    """
    log.info("{n} version {v} is running".format(
        n=os.path.basename(__file__), v=__version__))

    # core container of list of lists that have 6 items.
    # FileName, ZMW, start, end, length, PassedFilter as 0/1 (False/True)
    datum = []

    t0 = time.time()

    state = True

    log.info("Found {n} rgn.h5 files to analyze.".format(n=len(region_files)))
    log.info(pformat(region_files))

    for fileName in region_files:
        log.info("Looking for {f}".format(f=fileName))

        if os.path.exists(fileName):
            started_at = time.time()
            data = _get_data_from_rgn_file(fileName)
            datum.extend(data)
            run_time = time.time() - started_at

            log.info("Found {n} different Subreads in {s:.2f} sec ({m:.2f} min)".format(
                n=len(data), s=run_time, m=run_time / 60.0))

        else:
            state = False
            log.error("Unable to find {f}".format(f=fileName))

    # ToDo convert this to use the python stdlib
    with open(output_csv, 'w+') as f:
        log.info("{c} writing {f}".format(
            c=os.path.basename(__file__), f=output_csv))
        header = 'MovieName HoleNumber Start End Length PassedFilter'.split()

        f.write(",".join(header) + "\n")
        for data in datum:
            f.write(",".join([str(i) for i in data]) + "\n")

    run_time = time.time() - t0
    log.info("Completed writing {f} with {n} Subreads in {s:.2f} seconds ({m:.2f} minutes)".format(
        f=output_csv, n=len(datum), s=run_time, m=run_time / 60.0))

    return state


def _get_subreads(name, regions, r_index, hq_index):
    """
    Get the data from '/PulseData/Regions'

    :param name: (str) name of movie
    :param regions: (H5py.Group) Region H5 group
    :param r_index: (int) Region type index
    :param hq_index: (int) HQ Region index/flag

    Format #1
    RegionTypes = [GlobalAccuracy, HQRegion, Insert, Adapter]

    Format #2
    RegionTypes = [Adapter Insert HQRegion]

    """

    # dict of zmw_id -> regions (array of regions)
    zmws = {}

    # zmw_id -> True/False (if passed filter)
    zmw_passed = set([])

    for region in regions:
        zmw_id, rtype, start, end, passed_filter = region
        # log.info((name, region))

        if rtype == r_index:
            rg = [name, zmw_id, start, end, end - start]

            if zmw_id in zmws:
                zmws[zmw_id].append(rg)
            else:
                zmws[zmw_id] = [rg]

        if rtype == hq_index:
            if zmw_id in zmw_passed:
                log.debug(
                    "Found duplicate ZMW {i}. Please check the rgn file version!".format(i=zmw_id))
                continue
            else:
                # "Zero'ed" out regions have start == end == 0
                if start == end and start == 0:
                    # zmw_passed[zmw_id] = False
                    #zmw_passed[zmw_id] = 0
                    pass
                else:
                    # zmw_passed[zmw_id] = True
                    #zmw_passed[zmw_id] = 1
                    zmw_passed.add(zmw_id)
        else:
            continue

    # Store a list of (movie_name, ZMW_ID, start, end, length, 0/1)
    datum = []
    for z_id, rgs in zmws.iteritems():
        for r in rgs:
            line = [x for x in r]
            # 1/0 True/False if passed filter
            passed = 1 if z_id in zmw_passed else 0
            line.append(passed)
            # log.debug(line)
            datum.append(line)

    log.info("Found {n} ZMWs. {p} Passed filtering.".format(
        n=len(zmws), p=len(zmw_passed)))
    return datum


def _get_data_from_rgn_file(fileName):
    """
    Extracts the SubRead length from the rgn.cmp.h5 file

    Returns an Array of arrays of the form:
    [movie name, ZMD id, Start id, End id, length]

    The rgn.h5 File has two supported formats (see below).

     """

    name = _get_movie_name_from_pls_h5(fileName)

    if name is None:
        name = os.path.basename(fileName)
        log.info("Unable to extract movie name from {f}. Setting movie name to {n} from the filename".format(
            f=fileName, n=name))
    else:
        log.info("Found movie name {n} from Region File".format(n=name))

    f = h5py.File(fileName, 'r')

    group_name = '/PulseData/Regions'
    region_type_name = 'RegionTypes'

    for k, v in f[group_name].attrs.items():
        log.debug("Found Metadata {k} = {v} in {n}".format(
            k=k, v=v, n=fileName))

    # Get Version by inspecting the metadata
    if region_type_name in f[group_name].attrs:
        regionTypes = f[group_name].attrs[region_type_name]
        log.info("Found {r} metadata {x}".format(r=region_type_name,
                                                 x=regionTypes))
    else:
        msg = "Unable to find {n} to in rgn.h5 metadata".format(
            n=region_type_name)
        log.error(msg)
        raise KeyError(msg)

    # Create a map of RegionTypeName:Index
    version_dict = dict(zip(regionTypes, range(len(regionTypes))))

    log.debug("Created Insert/HQRegion map {x}".format(x=version_dict))

    insert_region_index = version_dict['Insert']
    hq_region_index = version_dict['HQRegion']

    data = [x for x in _get_subreads(
        name, f[group_name], insert_region_index, hq_region_index)]
    log.info("Found {n} Subreads in {f}".format(n=len(data), f=fileName))

    f.close()

    return data


def _get_movie_name_from_pls_h5(plsH5FileName):
    """
    Extract the Movie name from a pls.h5.

    Returns None if the name can't be found and raises a warning.

    :param plsH5FileName: pls.h5 file name
    :return: (str) movie name (or None)

    :note: This needs to be pulled from an API layer, not a random function :(
    """
    movieName = None

    groupName = '/ScanData/RunInfo'

    f = h5py.File(plsH5FileName, 'r')
    try:
        runInfoGroup = f[groupName]

        attrs = {k: v for k, v in runInfoGroup.attrs.iteritems()}

        if 'MovieName' in attrs:
            # FIXME There are some new Astro movie names which are listed as
            # lists of movie names and often only one movie name.
            # This may not be a reliable solution!!
            if isinstance(attrs['MovieName'], (np.ndarray, list)):
                movieName = attrs['MovieName'][0]
            elif isinstance(attrs['MovieName'], basestring):
                movieName = attrs['MovieName']
            else:
                msg = "Unsupported Movie Name from pls/bas H5 file"
                warnings.warn(msg)
                log.warn(msg)
        else:
            log.warn(
                "Unable to find MovieName in {g} attrs!".format(g=groupName))
        f.close()
    except KeyError as e:
        log.warning("unable to find group {n} in rgn file due to {e}".format(
            n=groupName, e=e.message))
        f.close()

    return movieName


def args_runner(args):
    log.info("Starting {f} v{v}".format(
        f=os.path.basename(__file__), v=__version__))
    region_files = bas_fofn_to_bas_files(args.region_fofn)
    state = run(region_files, args.output_csv)
    rcode = 0 if state else -1
    return rcode


def get_parser():
    """Util func to create an argparse instance

    Removing explicit usage due to issues with thirdparty argparse (Python < 2.7)

    usage = "usage: %prog --input=inputRgn.Fofn --outputCsv=mySubreadSummary.csv"
    """
    desc = 'Tool for generating a CSV file of the filtered Subreads from a file name of files (FOFN).'
    parser = get_default_argparser_with_base_opts(
        version=__version__, description=__doc__)
    parser.add_argument('region_fofn', type=validate_fofn,
                        help='Input Region FOFN path')
    parser.add_argument('-o', '--output-csv', default=None, dest='output_csv',
                        help='Output File to write summary to')
    return parser


def main(argv=sys.argv[1:]):
    """Main point of Entry"""
    return pacbio_args_runner(
        argv=argv,
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
