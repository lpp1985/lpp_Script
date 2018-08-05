
"""
Test that is_distributed is correct in resolved tool contracts.
"""

import os
import json
import logging

from pbcommand.pb_io.tool_contract_io import (load_tool_contract_from,
    load_resolved_tool_contract_from)
from pbsmrtpipe.testkit.loader import dtype_and_uuid_from_dataset_xml
from .base import TestBase

log = logging.getLogger(__name__)


class TestResolvedToolContracts(TestBase):

    def test_resolved_tool_contract_is_distributed(self):
        workflow_options_json = os.path.join(self.job_dir, "workflow",
            "options-workflow.json")
        workflow_options = json.load(open(workflow_options_json))
        distributed_mode = workflow_options["pbsmrtpipe.options.distributed_mode"]
        tasks_dir = os.path.join(self.job_dir, "tasks")
        for task_dir in os.listdir(tasks_dir):
            if task_dir.startswith("."):
                continue
            tc_file = os.path.join(tasks_dir, task_dir, "tool-contract.json")
            tc = load_tool_contract_from(tc_file)
            rtc_file = os.path.join(tasks_dir, task_dir,
                "resolved-tool-contract.json")
            rtc = load_resolved_tool_contract_from(rtc_file)
            if distributed_mode and tc.task.is_distributed:
                self.assertTrue(rtc.task.is_distributed,
                    "Resolved tool contract {f} has unexpected is_distributed=False".format(f=rtc_file))
            else:
                self.assertFalse(rtc.task.is_distributed,
                    "Resolved tool contract {f} has unexpected is_distributed=True".format(f=rtc_file))


    def test_rtc_output_files_in_datastore(self):
        """
        Confirm that all output files listed in resolved tool contracts are
        represented in the datastore.
        """
        datastore = None
        p = os.path.join(self.job_dir, "workflow", "datastore.json")
        with open(p, 'r') as r:
            datastore = json.loads(r.read())
        tasks_dir = os.path.join(self.job_dir, "tasks")
        datastore_output_files = {f["path"] for f in datastore["files"]}
        datastore_uuids = {f["uniqueId"] for f in datastore["files"]}
        for task_dir in os.listdir(tasks_dir):
            if task_dir.startswith("."):
                continue
            rtc_file = os.path.join(tasks_dir, task_dir,
                "resolved-tool-contract.json")
            rtc = load_resolved_tool_contract_from(rtc_file)
            for ofn in rtc.task.output_files:
                self.assertTrue(ofn in datastore_output_files,
                                "{o} not found in datastore".format(o=ofn))
