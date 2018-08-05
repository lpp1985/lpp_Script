"""
Tool contract wrappers for miscellaneous quick functions involving filtration
"""

import logging
import sys
import re

from pbcoretools.DataSetEntryPoints import parse_filter_list
from pbcore.io import openDataSet
from pbcommand.cli import registry_builder, registry_runner, QuickOpt
from pbcommand.models import FileTypes, OutputFileType

log = logging.getLogger(__name__)

TOOL_NAMESPACE = 'pbcoretools'
DRIVER_BASE = "python -m pbcoretools.tasks.filters "

registry = registry_builder(TOOL_NAMESPACE, DRIVER_BASE)


rl_opt = QuickOpt(0, "Minimum subread length",
                  "Minimum length of subreads")

filters_opt = QuickOpt(
    "",
    "Filters to add to the DataSet",
    "A comma separated list of other filters to add to the DataSet")

subreads_file_type = OutputFileType(FileTypes.DS_SUBREADS.file_type_id,
                                    "SubreadSet", "Filtered SubreadSet XML",
                                    "Filtered SubreadSet XML", "filtered")

def sanitize_read_length(read_length):
    if read_length:
        if not re.search('^-?\d*(\.\d*)?$', str(read_length).strip()):
            raise ValueError('read_length filter value "{v}" is not a '
                             'number'.format(v=read_length))
        try:
            return int(read_length)
        except ValueError:
            return int(float(read_length))

def run_filter_dataset(in_file, out_file, read_length, other_filters):
    dataSet = openDataSet(in_file)
    dataSet.updateCounts() # just in case
    if other_filters and other_filters != "None":
        filters = parse_filter_list(str(other_filters).split(','))
        dataSet.filters.addFilter(**filters)
        log.info("{i} other filters added".format(i=len(filters)))
    rlen = sanitize_read_length(read_length)
    if rlen:
        dataSet.filters.addRequirement(
            length=[('>=', rlen)])
    if rlen or other_filters:
        dataSet.updateCounts()
    dataSet.write(out_file)
    return 0

@registry("filterdataset", "0.1.0",
          FileTypes.DS_SUBREADS,
          subreads_file_type, is_distributed=True, nproc=1,
          options={"read_length":rl_opt,
                   "other_filters":filters_opt})
def run_filterDataSet(rtc):
    return run_filter_dataset(
        rtc.task.input_files[0], rtc.task.output_files[0],
        rtc.task.options["pbcoretools.task_options.read_length"],
        rtc.task.options["pbcoretools.task_options.other_filters"])

if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
