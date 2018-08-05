"""CLI Tools to facilitate development testing."""
import logging
import os
import socket
import sys
from pbcommand.cli import get_default_argparser

from pbcore.io import (FastaWriter, FastaReader, ReferenceSet)

import pbcore.io.dataset as DIO

from pbcommand.common_options import add_log_debug_option
from pbcommand.cli.utils import main_runner_default
from pbcommand.models.report import Report, Attribute
import pbcommand.cli.utils as U
import pbsmrtpipe.mock as M
from pbsmrtpipe.utils import compose

__version__ = '0.1.0'

log = logging.getLogger(__name__)


def _add_run_random_fasta_file(p):
    add_log_debug_option(p)
    p.add_argument('--max-records', type=int, default=25000, help="Max number of Fasta record to write.")
    U.add_fasta_output(p)
    return p


def _args_run_to_random_fasta_file(args):
    M.write_random_fasta_records(args.fasta_out, args.max_records)
    return 0


def _add_run_fasta_filter_options(p):
    add_log_debug_option(p)
    p.add_argument('--min-length', type=int, default=150, help='Min Length of Sequence to filter')
    U.add_fasta_input(p)
    U.add_fasta_output(p)
    return p


def run_fasta_filter(fasta_in, fasta_out, min_seq_length):
    with FastaWriter(fasta_out) as w:
        with FastaReader(fasta_in) as r:
            for record in r:
                if len(record.sequence) > min_seq_length:
                    w.writeRecord(record)

    return 0


def _args_run_fasta_filter(args):
    return run_fasta_filter(args.fasta_in, args.fasta_out, args.min_length)


def _add_run_random_fastq_options(p):
    add_log_debug_option(p)
    p.add_argument('--max-records', type=int, default=1000, help="Max number of Fasta record to write.")
    U.add_fastq_output(p)
    return p


def _args_run_random_fastq_file(args):
    M.write_random_fastq_records(args.fastq_out, nrecords=args.max_records)
    return 0


def _add_run_random_fofn_options(p):
    add_log_debug_option(p)
    U.add_output_dir_option(p)
    p.add_argument('--nfofns', type=int, default=10, help="Number of mock/random Fofns to write.")
    U.add_fofn_output(p)
    return p


def __to_random_fofn(contents, path):
    with open(path, 'w+') as w:
        w.write(contents)


def run_random_fofn(output_fofn, output_dir, nfofns):

    fofns = []
    for i in xrange(nfofns):
        name = "random_{i}".format(i=i)
        file_name = ".".join([name, 'fofn'])
        p = os.path.join(output_dir, file_name)
        with open(p, 'w+') as w:
            w.write(name)
        fofns.append(p)

    M.write_fofn(output_fofn, fofns)
    return 0


def _args_run_random_fofn(args):
    return run_random_fofn(args.fofn_out, args.output_dir, args.nfofns)


def _dataset_to_attribute_reports(ds):
    is_valid = all(os.path.exists(p) for p in ds.toExternalFiles())
    datum = [("uuid", ds.uuid, "Unique Id"),
             ("total_records", ds.numRecords, "num Records"),
             ("valid_files", is_valid, "External files exist")]
    attributes = [Attribute(x, y, name=z) for x, y, z in datum]
    return attributes


def dataset_to_report(ds):
    """
    :type ds: DataSet
    :param ds:
    :return:
    """
    attributes = _dataset_to_attribute_reports(ds)
    return Report("ds_report", attributes=attributes, dataset_uuids=[ds.uuid])


def subread_dataset_report(subread_path, report_path):
    ds = DIO.DataSet(subread_path)
    report = dataset_to_report(ds)
    report.write_json(report_path)
    return 0


def _add_run_dataset_report(p):
    add_log_debug_option(p)
    U.add_subread_input(p)
    U.add_report_output(p)
    return p


def _args_run_dataset_report(args):
    return subread_dataset_report(args.subread_ds, args.json_report)


def fasta_to_plot_group(fasta_file, output_dir):
    lengths = []
    with FastaReader(fasta_file) as f:
        for record in f:
            lengths.append(len(record.sequence))

    from pbreports.plot.helper import get_fig_axes
    from pbreports.model.model import PlotGroup, Plot
    fig, ax = get_fig_axes()

    if len(lengths) == 1:
        v = lengths[0]
        hrange = (v -1, v + 1)
        ax.hist(lengths, range=hrange)
    else:
        ax.hist(lengths)

    ax.set_title("Sequence Length Histogram")
    ax.set_xlabel("Sequence Length")

    name = "sequence_length_hist.png"
    png_path = os.path.join(output_dir, name)
    fig.savefig(png_path)
    plots = [Plot("sequence_lengths", name)]
    pg = PlotGroup("reference_hist", "Sequence Lengths", plots=plots)
    return pg


def try_fasta_to_plot_group(fasta_file, output_json):
    output_dir = os.path.dirname(output_json)
    plot_groups = []
    try:
        plot_group = fasta_to_plot_group(fasta_file, output_dir)
        plot_groups = [plot_group]
    except ImportError as ex:
        log.warn("pbreports is not installed. Skipping plot group generation {e}".format(e=ex))
    except Exception as ex:
        log.error("Unhandled error {e}".format(e=ex))

    return plot_groups


def fasta_to_report(fasta_file, output_json):

    nrecords = 0
    with FastaReader(fasta_file) as r:
        for _ in r:
            nrecords += 1

    attr = Attribute("num_records", nrecords, "Number of Records")
    plot_groups = try_fasta_to_plot_group(fasta_file, output_json)
    return Report("fasta_report", attributes=[attr], plotgroups=plot_groups)


def run_fasta_report(fasta_file, output_json):
    report = fasta_to_report(fasta_file, output_json)
    report.write_json(output_json)
    return 0


def run_reference_dataset_report(reference_ds, output_json):
    """

    :param reference_ds:
    :type reference_ds: ReferenceSet

    :param output_json:
    :return:
    """
    output_dir = os.path.dirname(output_json)
    host = socket.getfqdn()

    attributes = _dataset_to_attribute_reports(reference_ds)
    _add = attributes.append

    _add(Attribute("host", host, name="Host"))
    _add(Attribute("task_dir", output_dir, name="Task Directory"))

    fasta_file = reference_ds.toExternalFiles()[0]

    plot_groups = try_fasta_to_plot_group(fasta_file, output_dir)
    report = Report("dev_diagnostic_report",
                    attributes=attributes,
                    plotgroups=plot_groups,
                    dataset_uuids=[reference_ds.uuid])

    report.write_json(output_json)
    return 0


def _args_run_reference_dataset_report(args):
    reference_ds = ReferenceSet(args.reference_ds)
    return run_reference_dataset_report(reference_ds, args.json_report)


def _add_run_reference_dataset_report(p):
    opts = [add_log_debug_option,
            U.add_report_output,
            U.add_ds_reference_input]

    f = compose(*opts)
    return f(p)


def get_main_parser():
    p = get_default_argparser(__version__, "Tool For generating MOCK data for development testing.")

    sp = p.add_subparsers(help="Subparser Commands")

    def _builder(subparser_id, desc, options_func, exe_func):
        U.subparser_builder(sp, subparser_id, desc, options_func, exe_func)

    _builder('fasta', "Generate a random Fasta file", _add_run_random_fasta_file, _args_run_to_random_fasta_file)

    _builder('fastq', "Generate a random Fastq File", _add_run_random_fastq_options, _args_run_random_fastq_file)

    _builder('fofn', "Generate a random FOFN file", _add_run_random_fofn_options, _args_run_random_fofn)

    _builder('filter-fasta', "Filter a Fasta file by sequence length", _add_run_fasta_filter_options, _args_run_fasta_filter)

    _builder("dataset-report", "DataSet Report Generator", _add_run_dataset_report, _args_run_dataset_report)

    _builder("reference-ds-report", "Reference DataSet Report Generator", _add_run_reference_dataset_report, _args_run_reference_dataset_report)

    return p


def main(argv=None):

    argv_ = sys.argv if argv is None else argv
    parser = get_main_parser()

    return main_runner_default(argv_[1:], parser, log)


if __name__ == '__main__':
    sys.exit(main())
