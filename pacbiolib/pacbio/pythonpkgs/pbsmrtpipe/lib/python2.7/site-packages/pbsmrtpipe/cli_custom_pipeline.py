import logging
import sys

from pbcommand.cli import get_default_argparser_with_base_opts, pacbio_args_runner
from pbcommand.utils import setup_log

from pbsmrtpipe.core import registry_runner

log = logging.getLogger(__name__)

__version__ = "0.1.0"


def _get_parser():
    desc = "Custom PipelineTemplate Registry to write pipeline templates to output directory"
    p = get_default_argparser_with_base_opts(__version__, desc, default_level=logging.ERROR)
    p.add_argument('output_dir', help="Path to output directory")
    p.add_argument('--with-xml', action="store_true", default=False, help="Also Write Pipeline Templates as XML")
    return p


def _registry_runner(registry):
    def args_runner(args):
        import pbsmrtpipe.loader as L
        rtasks = L.load_all_tool_contracts()
        return registry_runner(registry, rtasks, args.output_dir, emit_xml=args.with_xml)
    return args_runner


def registry_runner_main(registry):
    def _main(argv=sys.argv):
        return pacbio_args_runner(argv[1:], _get_parser(), _registry_runner(registry), log, setup_log)
    return _main
