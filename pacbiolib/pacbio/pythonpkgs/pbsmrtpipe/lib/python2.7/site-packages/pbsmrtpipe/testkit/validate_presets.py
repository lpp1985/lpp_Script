
import logging
import os.path as op
import os
import sys

from pbcommand.cli.core import get_default_argparser_with_base_opts, \
    pacbio_args_runner
from pbcommand.utils import setup_log
from pbcommand.validators import validate_dir


__version__ = "0.1.0"
log = logging.getLogger(__name__)


def validate_preset_xml(dir_name):
    from pbsmrtpipe.pb_io import parse_pipeline_preset_xml, validate_raw_task_options
    import pbsmrtpipe.loader as L
    rtasks_d, _, _, pts = L.load_all()
    for file_name in os.listdir(dir_name):
        if file_name.endswith(".xml"):
            p = parse_pipeline_preset_xml(op.join(dir_name, file_name))
            if p.pipeline_id is None:
                raise ValueError("{f} does not have pipeline-id set".format(
                                 f=file_name))
            elif not p.pipeline_id in pts:
                raise ValueError("pipeline-id {i} not recognized".format(
                                 i=p.pipeline_id))
            log.info("validating {f}...".format(f=file_name))
            validate_raw_task_options(rtasks_d, dict(p.task_options))
        else:
            log.warn("Skipping non-XML file {f}".format(f=file_name))
    return 0


def args_runner(args):
    return validate_preset_xml(args.dir_name)


def get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__)
    p.add_argument("dir_name", type=validate_dir)
    return p


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
