
"""
Calls the ICE algorithm, which stands for 'Iteratively Clustering and Error
correction', to identify de novo consensus isoforms.
"""

import logging
import sys

from pbcommand.utils import setup_log
from pbcommand.cli.core import pbparser_runner

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser, get_argument_parser, add_subset_arguments)

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.subset"
    DRIVER_EXE = "python -m pbtranscript.tasks.subset --resolved-tool-contract"
    PARSER_DESC = __doc__


def get_contract_parser():
    p = get_base_contract_parser(Constants)
    add_subset_arguments(p)
    return p


def args_runner(args):
    return PBTranscript(args, subCommand="subset").start()


def resolved_tool_contract_to_args(resolved_tool_contract):
    args = [
        "subset",
        resolved_tool_contract.task.input_files[0],
        resolved_tool_contract.task.output_files[0],
    ]
    return get_argument_parser().parse_args(args)


def resolved_tool_contract_runner(resolved_tool_contract):
    args = resolved_tool_contract_to_args(resolved_tool_contract)
    return args_runner(args)


def main(argv=sys.argv[1:]):
    mp = get_contract_parser()
    return pbparser_runner(
        argv=argv,
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
