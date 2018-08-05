import sys
import logging

from pbcommand.common_options import add_log_debug_option
from pbcommand.cli import pacbio_args_runner, get_default_argparser
from pbcommand.utils import setup_log
from pbcommand.cli.utils import main_runner_default
from pbcommand.validators import validate_file

import pbsmrtpipe.report_renderer as R
import pbcommand.cli.utils as U
from pbcommand.pb_io import load_report_from_json

log = logging.getLogger(__name__)
slog = logging.getLogger('status.' + __name__)

__version__ = '1.0'


def _validate_report(path):
    p = validate_file(path)
    _ = load_report_from_json(path)
    return p


def _add_output_file_option(p):
    p.add_argument('--output-file', required=True, type=str, help="Path of output html")
    return p


def _add_ccs_js_extras_option(p):
    _d = "Write styled CSS and JS dirs/files"
    p.add_argument('--with-extras', action='store_true', help=_d)
    return p


def _add_report_option(p):
    p.add_argument('report_path', type=_validate_report, help="Path to pbreport JSON report")
    return p


def _args_to_render_report(args):
    f = R.write_report_with_html_extras if args.with_extras else R.write_report_to_html
    report = load_report_from_json(args.report_path)
    return f(report, args.output_file)


def get_parser():
    desc = "Transform pbreport Report to HTML file."
    p = get_default_argparser(__version__, desc)
    add_log_debug_option(p)
    U.add_output_dir_option(p)
    _add_report_option(p)
    _add_output_file_option(p)
    _add_ccs_js_extras_option(p)
    p.set_defaults(func=_args_to_render_report)
    return p


def main(argv=None):

    argv_ = sys.argv if argv is None else argv
    parser = get_parser()
    return pacbio_args_runner(argv_[1:], parser, _args_to_render_report, log, setup_log)
