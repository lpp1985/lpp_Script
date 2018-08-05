
# FIXME(nechols)(2016-01-20): way too much code duplication

import logging
import sys

from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

from pbcoretools.tasks import gather_report
from pbcoretools.chunking.gather import run_main_gather_report

log = logging.getLogger(__name__)

class Constants(object):
    TOOL_ID = "pbcoretools.tasks.gather_json"
    CHUNK_KEY = "$chunk.report_id"
    VERSION = "0.1.0"
    DRIVER = "python -m pbcoretools.tasks.gather_json --resolved-tool-contract "
    OPT_CHUNK_KEY = 'pbcoretools.task_options.gather_report_chunk_key'
    REPORT_TYPE = FileTypes.JSON


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           gather_report.get_parser(constants=Constants),
                           gather_report.args_runner,
                           gather_report.rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
