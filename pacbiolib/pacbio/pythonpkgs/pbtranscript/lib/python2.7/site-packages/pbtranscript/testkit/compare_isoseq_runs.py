
"""
Compare the results of a pbtranscript run (or pbsmrtpipe IsoSeq job) to an
existing benchmark.
"""

from cPickle import *
import argparse
import difflib
import os.path as op

import numpy as np
from scipy.sparse import lil_matrix

from pbcore.io import *


def filter_clusters_by_size(clusters, min_cluster_size=2):
    """Filter clusters in dict{cluster_id:set(), cluster_id:set(), ...} by size.
    e.g., clusters={0:set("movie/zmw/s_e",
                          "movie/zmw/s_e"),
                    ...
                    18:set("", "")}
          return clusters which have at least $min_cluster_size members.
    """
    return {k:v for (k, v) in clusters.iteritems() if len(v) >= min_cluster_size}

def filter_cluster_seqs_by_full_length_coverage(fa_fn, min_full_length_coverage=2):
    """Filter cluster sequences in input fasta file, by number of supportive
    full length non-chimeric reads in that cluster.
    e.g., fa_fn = "all_quivered_hq.100_30_0.99.fasta"
          $ cat fa_fn
          >c0/f3p1/3650 isosoform=c0;full_length_coverage=3
          agcc..
    returns [FastaRecord("c0/f3p1/3650", "agcc"), ..., FastaRecord()]
    """
    return [record for record in FastaReader(fa_fn)
            if int(record.name.split('/')[1].split('p')[0][1:]) >= min_full_length_coverage]

def jacard_sim(c1, c2):
    """
    c1, c2, are two clusters, each is a list of seq ids, return jacard similarity.
    """
    return len(set(c1).intersection(c2))*1./len(set(c1).union(c2))


def calc_sim_matrix(uc1, uc2):
    """
    Calculate matrix similarity between two sets of clusters.
    uc1: a set of clusters, dict{cluster_id:set(), cluster_id:set(), ...}
    uc2: a set of clusters, dict{cluster_id:set(), cluster_id:set(), ...}

    return (S, keys1, keys2, m, n), where
            S_ij = jacard_sim(uc1_i, uc2_j)
            keys1, keys2 = cluster ids in uc1 and uc2
            m, n = number of clusters in uc1 and uc2
    """
    keys1 = uc1.keys()
    keys1.sort()
    keys2 = uc2.keys()
    keys2.sort()
    m, n = len(keys1), len(keys2)
    S = lil_matrix((m+1, n+1))
    for i in xrange(m):
        for j in xrange(n):
            c1 = uc1[keys1[i]]
            c2 = uc2[keys2[j]]
            sim = jacard_sim(c1, c2)
            S[i, j] = sim
    return S, keys1, keys2, m, n


def get_seq_similarity(lhs_seq, rhs_seq):
    """Return percentage sequence similarity."""
    lhs_tmp_f = 'lhs_%s.fasta' % (lhs_seq.name.split(' ')[0].replace('/', '_'))
    rhs_tmp_f = 'rhs_%s.fasta' % (rhs_seq.name.split(' ')[0].replace('/', '_'))

    with FastaWriter(lhs_tmp_f) as writer:
        writer.writeRecord(lhs_seq.name.split(' ')[0], lhs_seq.sequence)

    with FastaWriter(rhs_tmp_f) as writer:
        writer.writeRecord(rhs_seq.name.split(' ')[0], rhs_seq.sequence)

    cmd = "blasr %s %s -m 4 --bestn 1 2>/dev/null |cut -f 4 -d ' ' " % (lhs_tmp_f, rhs_tmp_f)
    o, c, e = backticks(cmd)
    assert c==0

    backticks('rm %s %s' % (lhs_tmp_f, rhs_tmp_f))
    return float(o[0])


class Compare_Isoseq_Runs(object):
    """Compare two isoseq runs, including
    init.uc.pickle
    output/final.consensus.fasta
    all_quivered_hq_*.fasta
    """
    def __init__(self, lhs_dir, rhs_dir,
                 min_cluster_matrix_similarity,
                 min_seq_similarity,
                 min_flnc_coverage,
                 print_diff=False):
        self.lhs_dir = lhs_dir
        self.rhs_dir = rhs_dir
        self.min_cluster_matrix_similarity = min_cluster_matrix_similarity if min_cluster_matrix_similarity <= 1 else min_cluster_matrix_similarity / 100.0
        self.min_seq_similarity = min_seq_similarity if min_seq_similarity <= 1 else min_seq_similarity / 100.0
        self.min_flnc_coverage = min_flnc_coverage
        self.num_differences = 0
        self.print_diff = print_diff

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if (self.num_differences != 0):
            raise ValueError ("Isoseq results in %s %s are siginifantly different." %
                              (self.lhs_dir, self.rhs_dir))

    def _submit(self, OK):
        """Keep track of comparison results."""
        if not OK:
            self.num_differences += 1

    def _get_fn(self, which_dir, fn):
        """Get file path."""
        ret = fn
        if which_dir is None:
            pass
        elif which_dir == "lhs":
            ret = op.join(self.lhs_dir, fn)
        elif which_dir == "rhs":
            ret = op.join(self.rhs_dir, fn)
        else:
            raise RuntimeError("unknown which_dir %s" % which_dir)

        if not op.exists(ret):
            raise RuntimeError("Unable to find %s" % ret)
        return ret

    def init_pickle_fn(self, which_dir=None):
        """path to initial pickle."""
        return self._get_fn(which_dir, "init.uc.pickle")

    def compare_init_pickle(self, min_cluster_matrix_similarity=0.99):
        """Compare init pickle files."""
        lhs_pickle = load(open(self.init_pickle_fn("lhs")))
        rhs_pickle = load(open(self.init_pickle_fn("rhs")))

        # clusters which are not orphan
        lhs_clusters = filter_clusters_by_size(lhs_pickle, 2)
        rhs_clusters = filter_clusters_by_size(rhs_pickle, 2)

        S, k1, k2, m, n = calc_sim_matrix(lhs_clusters, rhs_clusters)
        OK = (S.sum()/max(m,n) >= min_cluster_matrix_similarity)
        print "Q: Initial clusters in init_pickle similarity >= %f?" % min_cluster_matrix_similarity
        print "A: %s (%f)" % (OK, S.sum()/max(m,n))
        self._submit(OK)
        return OK

    def print_diff_init_pickle(self):
        """Print differences between init_pickle."""
        lhs_pickle = load(open(self.init_pickle_fn("lhs")))
        rhs_pickle = load(open(self.init_pickle_fn("rhs")))

        lhs_clusters = filter_clusters_by_size(lhs_pickle, 2)
        rhs_clusters = filter_clusters_by_size(rhs_pickle, 2)

        print "PRINT_DIFF_INIT_PICKLE STARTED"
        for lhs_cluster in lhs_clusters:
            if (len(lhs_clusters[lhs_cluster]) > 1):
                print "lhs cluster [%s]" % lhs_cluster
                for read in lhs_clusters[lhs_cluster]:
                    print "--> Looking at %s" % (read)
                    for rhs_cluster in rhs_clusters:
                        if read in rhs_clusters[rhs_cluster]:
                            print "-----> read in rhs cluster [%s]" % rhs_cluster
        print "PRINT_DIFF_INIT_PICKLE ENDED"

    def final_consensus_fn(self, which_dir=None):
        """path to output/final.consensus.fasta."""
        try:
            return self._get_fn(which_dir, "output/final.consensus.fasta")
        except IOError, RuntimeError:
            return self._get_fn(which_dir, "output/final.consensus.fa")

    def compare_ref_consensus(self):
        """Compare ref_consensus."""
        lhs_ref_consensus = self.final_consensus_fn("lhs")
        rhs_ref_consensus = self.final_consensus_fn("rhs")

        len1 = len([r for r in FastaReader(lhs_ref_consensus)])
        len2 = len([r for r in FastaReader(rhs_ref_consensus)])
        print "Q: Number of sequences in ref_consensus equalivalent?"
        print "A: %s (%d vs %d)" % (len1==len2, len1, len2)

        lhs_ref_seqs = [r.sequence for r in FastaReader(lhs_ref_consensus)]
        rhs_ref_seqs = [r.sequence for r in FastaReader(rhs_ref_consensus)]

        OK = (lhs_ref_seqs == rhs_ref_seqs)
        print "Q: Unpolished cluster sequences in output/final.consensus.fasta identical?"
        print "A: %s (%d vs %d)" % (OK, len(lhs_ref_seqs), len(rhs_ref_seqs))
        self._submit(OK)
        return OK

    def polished_hq_fn(self, which_dir=None):
        """path to polished hq clusters."""
        return self._get_fn(which_dir, "all_quivered_hq.100_30_0.99.fasta")

    def compare_polished_hq_counts(self, min_flnc_coverage=2):
        """We want to make sure that all polished hq sequences which
        have at least two supportive full length reads almost identical. """
        lhs_hq_fn = self.polished_hq_fn("lhs")
        rhs_hq_fn = self.polished_hq_fn("rhs")

        lhs_filtered_seqs = filter_cluster_seqs_by_full_length_coverage(lhs_hq_fn, min_flnc_coverage)
        rhs_filtered_seqs = filter_cluster_seqs_by_full_length_coverage(rhs_hq_fn, min_flnc_coverage)

        OK = (len(lhs_filtered_seqs) == len(rhs_filtered_seqs))
        print "Q: Number of hq clusters with at least %d full length coverage identical?" % min_flnc_coverage
        print "A: %s (%d vs %d)" % (OK, len(lhs_filtered_seqs), len(rhs_filtered_seqs))
        self._submit(OK)
        return OK

    def compare_polished_hq_sequences(self, min_seq_similarity=0.99, min_flnc_coverage=2,
            sort_first=False):
        lhs_hq_fn = self.polished_hq_fn("lhs")
        rhs_hq_fn = self.polished_hq_fn("rhs")

        lhs_filtered_seqs = filter_cluster_seqs_by_full_length_coverage(lhs_hq_fn, min_flnc_coverage)
        rhs_filtered_seqs = filter_cluster_seqs_by_full_length_coverage(rhs_hq_fn, min_flnc_coverage)
        if sort_first:
            lhs_filtered_seqs.sort(lambda a,b: cmp(a.name, b.name))
            rhs_filtered_seqs.sort(lambda a,b: cmp(a.name, b.name))
        OK = True
        for i in range(0, min(len(lhs_filtered_seqs), len(rhs_filtered_seqs))):
            lhs_seq, rhs_seq = lhs_filtered_seqs[i], rhs_filtered_seqs[i]
            seq_similarity = get_seq_similarity(lhs_seq, rhs_seq)
            if (seq_similarity < min_seq_similarity):
                OK = False
                print "sequence similarity between lhs (%s) and rhs (%s): %f" % \
                        (rhs_filtered_seqs[i].name, lhs_filtered_seqs[i].name,
                         seq_similarity)
        print "Q: Sequence similarity all greater than %f?" % min_seq_similarity
        print "A: %s" % OK
        self._submit(OK)
        return OK

    def run(self):
        """Run comparisons"""
        if self.print_diff:
            self.print_diff_init_pickle()
        self.compare_init_pickle(min_cluster_matrix_similarity=self.min_cluster_matrix_similarity)
        self.compare_ref_consensus()
        self.compare_polished_hq_counts(min_flnc_coverage=self.min_flnc_coverage)
        self.compare_polished_hq_sequences(min_seq_similarity=self.min_seq_similarity,
                                           min_flnc_coverage=self.min_flnc_coverage)


def run(args):
    """Construct an instance of Compare_IsoSeq_Runs and do the comparison."""
    with Compare_Isoseq_Runs(lhs_dir=args.lhs_dir,
                             rhs_dir=args.rhs_dir,
                             min_cluster_matrix_similarity=args.min_cluster_matrix_similarity,
                             min_seq_similarity=args.min_seq_similarity,
                             min_flnc_coverage=args.min_flnc_coverage,
                             print_diff=args.print_diff) as runner:
        runner.run()

def get_parser():
    """Get argument parser."""
    parser = argparse.ArgumentParser(description="Compare results of two isoseq runs.")
    parser.add_argument("lhs_dir", action="store", type=str, help="output of one isoseq run")
    parser.add_argument("rhs_dir", action="store", type=str, help="output of another isoseq run")
    parser.add_argument("--min_cluster_matrix_similarity", action="store", dest="min_cluster_matrix_similarity", type=float, help="minimum cluster matrix similarity.", default=0.99)
    parser.add_argument("--min_seq_similarity", action="store", dest="min_seq_similarity", type=float, help="minimum sequence similarity.", default=0.99)
    parser.add_argument("--min_flnc_coverage", action="store", dest="min_flnc_coverage", type=int, help="minimum full-length-non-chimeric reads coverage.", default=2)
    parser.add_argument("--print_diff", action="store_true", dest="print_diff", help="Print difference for debugging.", default=False)
    return parser

if __name__ == "__main__":
    import sys
    parser = get_parser()
    run(parser.parse_args(sys.argv[1:]))
