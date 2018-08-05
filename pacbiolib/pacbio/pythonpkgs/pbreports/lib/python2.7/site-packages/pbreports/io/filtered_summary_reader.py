import sys
import os
import logging

import numpy as np

from pbreports.io.csv_reader import CsvReader

log = logging.getLogger(__name__)


class FilteredSummaryReader(CsvReader):

    # In order to support variation in columns present,
    # we need to specify their types
    # Format is { Name: ( dType, f(str->dType) ) }
    _column_dtypes = {"Movie": ("|S64", str),
                      "ReadId": ("|S64", str),
                      "#Bases": (int, int),
                      "Readlength": (int, int),
                      "ReadScore": (float, float),
                      "Productivity": (int, int),
                      "SequencingZMW": (int, int),
                      "PassedFilter": (int, int),
                      "Sandwiches": (int, int),
                      "Whitelisted": (int, int),
                      "SNR": (float, float),
                      "ArtifactScore": ("|S64", str)}

    def __init__(self, path, column_map=None, strict=False):

        if column_map is None:
            column_map = FilteredSummaryReader._column_dtypes
        CsvReader.__init__(self, path, column_map, strict)

    @property
    def num_reads(self):
        """Get the total number of reads, filtered and unfiltered"""
        bool_indices = self._np_array["SequencingZMW"] > 0
        return len(self._np_array[bool_indices])
