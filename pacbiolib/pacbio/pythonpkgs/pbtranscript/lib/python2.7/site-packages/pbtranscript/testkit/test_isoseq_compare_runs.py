
"""
pysiv2 test for pbsmrtpipe, wrapping compare_isoseq_runs
"""

from unittest import SkipTest
import cPickle
import os.path as op

from pbcommand.models import FileTypes
try:
    from pysiv2.custom.base import TestValuesLoader
except ImportError:
    import unittest
    TestValuesLoader = unittest.TestCase

from pbtranscript.testkit.compare_isoseq_runs import Compare_Isoseq_Runs


class TestCompareIsoseqRuns(TestValuesLoader):
    """
    Class for running pysiv2 tests on pbsmrtpipe jobs.
    """

    @classmethod
    def setUpClass(cls):
        super(TestCompareIsoseqRuns, cls).setUpClass()
        cls.rhs_dir = cls.test_values.get("isoseq", {}).get("ref_dir", None)
        cls.rhs_nfl_pickle = op.join(cls.rhs_dir, "output", "map_noFL",
            "nfl.all.partial_uc.pickle")
        assert op.isfile(cls.rhs_nfl_pickle)
        cls.lhs_dir = None
        cls.lhs_nfl_pickle = None
        cls.lhs_nfl_pickle_tmp = None
        for file_id, file_info in cls.datastore.get_file_dict().iteritems():
            if file_info.file_type_id == FileTypes.PICKLE.file_type_id:
                if "pbtranscript.tasks.cluster" in file_info.file_id:
                    cls.lhs_dir = op.dirname(file_info.path)
                elif "pbtranscript.tasks.gather_nfl_pickle" in file_info.file_id:
                    cls.lhs_nfl_pickle = file_info.path
                elif "pbtranscript.tasks.ice_partial" in file_info.file_id:
                    cls.lhs_nfl_pickle_tmp = file_info.path
        if cls.lhs_nfl_pickle is None:
            cls.lhs_nfl_pickle = cls.lhs_nfl_pickle_tmp
        # XXX we can't subclass this instead because of the way test case
        # setup is performed
        cls._peer = Compare_Isoseq_Runs(
            lhs_dir=cls.lhs_dir,
            rhs_dir=cls.rhs_dir,
            min_cluster_matrix_similarity=0.99,
            min_seq_similarity=0.99,
            print_diff=True)

    def test_isoseq_polished_high_quality_count(self):
        self.assertTrue(self._peer.compare_polished_hq_counts())

    def test_isoseq_polished_high_quality_sequences(self):
        self.assertTrue(self._peer.compare_polished_hq_sequences(
            sort_first=True))

    def test_ref_consensus(self):
        self.assertTrue(self._peer.compare_ref_consensus())

    def test_isoseq_compare_init_pickle(self):
        self.assertTrue(self._peer.compare_init_pickle())

    def test_isoseq_compare_nfl_pickle(self):
        if self.lhs_nfl_pickle is None:
            raise SkipTest("Can't find n.f.l. pickle file from job")
        p1 = p2 = None
        with open(self.rhs_nfl_pickle, 'rb') as f:
            p1 = cPickle.load(f)
        with open(self.lhs_nfl_pickle, 'rb') as f:
            p2 = cPickle.load(f)
        p3 = {k:v for k,v in p1['partial_uc'].iteritems() if k != "nohit"}
        p4 = {k:v for k,v in p2['partial_uc'].iteritems() if k != "nohit"}
        msg = "\n".join(["Mismatch between NFL pickles:", "Reference:"] +
            [ "{k}: {v}".format(k=k, v=v) for k,v in p3.iteritems() ] +
            ["!=", "Job:"] +
            [ "{k}: {v}".format(k=k, v=v) for k,v in p4.iteritems() ])
        self.assertEqual(p3, p4, msg)
