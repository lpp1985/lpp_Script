
# wraps 'pbtranscript classify'

"""
Classifies reads from a fasta/q file.  For each read, identify whether it is
full length, whether 5', 3' and poly A tail have been found. The input is a
ConsensusRead dataset.
"""

import logging
import os.path as op
import sys

from pbcore.io import FastaRecord, FastaWriter

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.Utils import mkdir
from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser,
                                              get_argument_parser,
                                              add_classify_arguments)

log = logging.getLogger(__name__)


def parse_primer_sequences(primers_str):
    """
    Return a list of primer FastaRecords if primers_str only contains
    valid primers. Otherwise raise a ValueError.
    """
    if isinstance(primers_str, str) or isinstance(primers_str, unicode):
        primer_fasta_records = []
        primers_str = str(primers_str)
        if '>' not in primers_str:
            raise ValueError("Invalid primer header, could not find leading '>'.")
        for str_index, primer_str in enumerate(primers_str.split('>')[1:]):
            lines = [line.strip().translate(None, '\'\" ') for line in primer_str.split('\n')]
            lines = [line for line in lines if len(line) > 0] # remove empty lines
            if len(lines) < 2:
                raise ValueError("Primer %s must have a sequence." % lines[0])
            primer_name = lines[0]
            primer_sequence = ''.join(lines[1:])
            primer_index = int(str_index / 2)
            primer_strand = 'F' if str_index % 2 == 0 else 'R'
            expected_primer_name = "{s}{i}".format(s=primer_strand, i=primer_index)
            if primer_name != expected_primer_name:
                raise ValueError("Primers should be placed in order F0, R0, F1, R1...")
            for base in primer_sequence:
                if base.upper() not in ('A', 'T', 'G', 'C'):
                    raise ValueError("Primer sequence %s must only contain ATGC" % primer_sequence)
            primer_fasta_records.append(FastaRecord(header=primer_name, sequence=primer_sequence))

        return primer_fasta_records

    raise ValueError("Input primers_str %s must be either str or unicode" % type(primers_str))


class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.classify"
    DRIVER_EXE = "python -m %s --resolved-tool-contract" % TOOL_ID
    PARSER_DESC = __doc__


def get_contract_parser():
    """Get PbParser for classify."""
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    add_classify_arguments(p)
    return p


def args_runner(args):
    """Call Args runner"""
    return PBTranscript(args, subCommand="classify").start()


def resolved_tool_contract_to_args(resolved_tool_contract):
    """Convert resolved tool contract to args."""
    rtc = resolved_tool_contract
    args = [
        "--verbose",
        "classify",
        resolved_tool_contract.task.input_files[0],
        resolved_tool_contract.task.output_files[0],
        "--flnc", resolved_tool_contract.task.output_files[1],
        "--nfl", resolved_tool_contract.task.output_files[2],
        "--summary", resolved_tool_contract.task.output_files[3],  # JSON
        "--report", resolved_tool_contract.task.output_files[4],  # CSV
        "--min_seq_len", str(rtc.task.options[Constants.MIN_SEQ_LEN_ID]),
        "--cpus", str(resolved_tool_contract.task.nproc),
        "--outDir", op.dirname(rtc.task.output_files[0]),
        "--ignore-empty-output",
    ]
    if rtc.task.options[Constants.IGNORE_POLYA_ID]:
        args.append("--ignore_polyA")

    primers_str_obj = rtc.task.options[Constants.PRIMER_SEQUENCES_ID]
    primers_str = str(primers_str_obj).strip().translate(None, '\'\" ')
    if primers_str_obj is not None and primers_str not in ('None', ''):
        logging.info("Detected customer primer: %s", primers_str)
        # Save primer sequences to a fasta file under output dir
        primer_fasta_records = parse_primer_sequences(primers_str=primers_str)
        d = op.dirname(resolved_tool_contract.task.output_files[2])
        mkdir(d)
        primer_fn = op.join(d, "customer_primers.fasta")
        with FastaWriter(primer_fn) as writer:
            for record in primer_fasta_records:
                writer.writeRecord(record)
        logging.info("Customer primer sequences written to file %s", primer_fn)
        args.append("-p")
        args.append("%s" % primer_fn)
    else:
        logging.info("No customer primer detected.")

    return get_argument_parser().parse_args(args)


def resolved_tool_contract_runner(resolved_tool_contract):
    args = resolved_tool_contract_to_args(resolved_tool_contract)
    return args_runner(args)


def main(argv=sys.argv[1:]):
    return pbparser_runner(
        argv=argv,
        parser=get_contract_parser(),
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
