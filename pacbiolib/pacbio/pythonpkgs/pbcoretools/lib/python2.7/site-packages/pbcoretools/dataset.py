import sys
import argparse
import logging
import time

from pbcommand.cli import get_default_argparser_with_base_opts, pacbio_args_runner
from pbcommand.utils import setup_log

from pbcoretools import DataSetEntryPoints as EntryPoints
from pbcoretools.version import __VERSION__

log = logging.getLogger(__name__)

def get_subparsers():
    sps = [('create', EntryPoints.create_options),
           ('filter', EntryPoints.filter_options),
           ('merge', EntryPoints.merge_options),
           ('split', EntryPoints.split_options),
           ('validate', EntryPoints.validate_options),
           ('summarize', EntryPoints.summarize_options),
           ('consolidate', EntryPoints.consolidate_options),
           ('loadstats', EntryPoints.loadStatsXml_options),
           ('newuuid', EntryPoints.newUniqueId_options),
           ('loadmetadata', EntryPoints.loadMetadataXml_options),
           ('copyto', EntryPoints.copyTo_options),
           ('absolutize', EntryPoints.absolutize_options),
           ('relativize', EntryPoints.relativize_options),
          ]
    return sps

def add_subparsers(parser, sps):
    subparsers = parser.add_subparsers(
        title='DataSet sub-commands', dest='subparser_name',
        help="Type {command} -h for a command's options")
    for command_name, func in sps:
        subparser = subparsers.add_parser(
            command_name,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        subparser = func(subparser)
    return parser

def get_parser():
    description = 'Run dataset.py by specifying a command.'
    parser = get_default_argparser_with_base_opts(
        version=__VERSION__,
        description=description,
        default_level="WARNING")
    parser.add_argument("--strict", default=False, action='store_true',
                        help="Turn on strict tests, raise all errors")
    parser.add_argument("--skipCounts", default=False, action='store_true',
                        help="Turn on strict tests, raise all errors")
    subparser_list = get_subparsers()
    parser = add_subparsers(parser, subparser_list)
    return parser

def run(args):
    log.info("Starting {f} version {v} dataset manipulator".format(
        f=__file__, v=__VERSION__))
    return args.func(args)

def main(argv=sys.argv):
    """Main point of Entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=run,
        alog=log,
        setup_log_func=setup_log)

if __name__ == '__main__':
    sys.exit(main())
