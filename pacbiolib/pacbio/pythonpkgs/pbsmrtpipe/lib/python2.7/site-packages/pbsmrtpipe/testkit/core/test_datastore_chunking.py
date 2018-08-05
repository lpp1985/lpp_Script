
"""
Specialized tests for datastore behavior when chunking is used, including
propagation of sourceId, name, and description fields from the original task.
Note that this test is *not* suitable for pipelines like dev_diagnostic_stress
which run a specific task more than once.
"""

import unittest
import os
import logging
import json
import re

from pbcommand.pb_io.report import load_report_from_json
from pbcommand.models import FileTypes
from pbcore.io import getDataSetUuid

from .base import TestBase


class TestDataStoreSourceIds(TestBase):

    def test_datastore_unchunked_source_ids_unique(self):
        p = os.path.join(self.job_dir, "workflow", "datastore.json")
        with open(p, 'r') as r:
            d = json.loads(r.read())
            source_ids = set()
            for file_info in d['files']:
                if file_info['isChunked']:
                    continue
                self.assertFalse(file_info['sourceId'] in source_ids,
                                 "{i} occurs more than once".format(
                                 i=file_info['sourceId']))
                source_ids.add(file_info['sourceId'])

    def test_datastore_chunk_gather_consistency(self):
        p = os.path.join(self.job_dir, "workflow", "datastore.json")
        with open(p, 'r') as r:
            d = json.loads(r.read())
            chunked_files = {}
            unchunked_files = {}
            for file_info in d['files']:
                if not file_info['isChunked']:
                    unchunked_files[file_info['sourceId']] = file_info
            for file_info in d['files']:
                if file_info['isChunked']:
                    # CHUNK JSON files from scatter tasks are a special case
                    if file_info['fileTypeId'] == "PacBio.FileTypes.CHUNK":
                        continue
                    gathered = unchunked_files.get(file_info['sourceId'], None)
                    self.assertTrue(gathered is not None,
                                    "Can't find gathered file {i}".format(
                                    i=file_info['sourceId']))
                    self.assertEqual(file_info['name'], gathered['name'])
                    self.assertEqual(file_info['description'],
                                     gathered['description'])
