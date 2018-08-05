"""Define functions useful for IceInit and IceIterative."""
import os
import os.path as op
import logging
import shutil
import filecmp
import random
import time
from cPickle import dump, load
from collections import defaultdict
import numpy as np
import pysam
from pbcore.util.Process import backticks
from pbcore.io import FastaReader, FastaWriter, FastqWriter, \
        BasH5Reader, ContigSet
from pbtranscript.Utils import realpath, mkdir, execute, \
        write_files_to_fofn, real_upath, \
        get_files_from_file_or_fofn, \
        FILE_FORMATS, guess_file_format
from pbtranscript.RunnerUtils import write_cmd_to_script
from pbtranscript.findECE import findECE
from pbtranscript.io.BasQV import basQVcacher
from pbtranscript.io import BLASRM5Reader, MetaSubreadFastaReader, \
        BamCollection, BamWriter, LA4IceReader
from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper
from pbtranscript.ice_daligner import DalignerRunner
from pbtranscript.ice.ProbModel import ProbFromQV, \
    ProbFromModel, ProbFromFastq

__author__ = 'etseng@pacificbiosciences.com'

# define gcon script for ice.
gcon_py = "python -m pbtranscript.ice_pbdagcon"

# Define data sets for sge sanity check.
dataDir = op.join(op.dirname(op.dirname(op.realpath(__file__))), "data")
GCON_IN_FA = op.join(dataDir, "gcon_in.fasta")
GCON_OUT_FA = op.join(dataDir, "gcon_out.fasta")

random.seed(0)


def check_blasr(required_min_version=5.1):
    """
    blasr version >= required_min_version
    """
    _o, _c, _e = backticks('blasr --version')
    succeed = True
    if _c == 0:
        try:
            v = float('.'.join(_o[0].split()[1].split('.')[0:2]))
            if v < required_min_version:
                succeed = False
        except Exception:
            succeed = False
    else:
        succeed = False

    if not succeed:
        msg = "blasr not installed or version < %s" % (required_min_version)
        logging.error(msg)
        raise RuntimeError(msg)

    return True


def sanity_check_daligner(scriptDir, testDirName="daligner_test_dir"):
    """
    Run daligner on gcon_in.fa, but don't care about results.
    Just make sure it runs.
    """
    scriptDir = realpath(scriptDir)
    testDir = op.join(scriptDir, testDirName)

    mkdir(scriptDir)
    mkdir(testDir)

    testInFa = op.join(testDir, "daligner.fasta")
    if op.exists(testInFa):
        os.remove(testInFa)
    shutil.copy(GCON_IN_FA, testInFa)
    assert op.exists(testInFa)

    runner = DalignerRunner(query_filename=testInFa,
                            target_filename=testInFa,
                            is_FL=True, same_strand_only=True,
                            query_converted=False, target_converted=False,
                            use_sge=False, cpus=4, sge_opts=None)
    runner.run(output_dir=testDir, min_match_len=300, sensitive_mode=False)
    runner.clean_run()

    shutil.rmtree(testDir)
    logging.info("daligner check passed.")
    return True

def sanity_check_gcon():
    """Sanity check gcon."""
    cmd = gcon_py + " --help"

    errmsg = gcon_py + " is not installed."
    execute(cmd=cmd, errmsg=errmsg)
    return gcon_py


def sanity_check_sge(sge_opts, scriptDir, testDirName="gcon_test_dir"):
    """Sanity check if sge can work."""
    scriptDir = realpath(scriptDir)
    testDir = op.join(scriptDir, testDirName)

    if not op.exists(scriptDir):
        os.makedirs(scriptDir)
    if not op.exists(testDir):
        os.makedirs(testDir)

    testSh = op.join(scriptDir, 'test.sh')
    consensusFa = op.join(testDir, "g_consensus.fasta")
    testInFa = op.join(testDir, "gcon_in.fasta")
    if op.exists(testInFa):
        os.remove(testInFa)
    shutil.copy(GCON_IN_FA, testInFa)
    assert op.exists(testInFa)

    cmd = " ".join([gcon_py, real_upath(testInFa),
                    "{testDir}/g_consensus".format(testDir=real_upath(testDir)),
                    "c1"])
    write_cmd_to_script(cmd=cmd, script=testSh)

    assert op.exists(testSh)
    cmd = sge_opts.qsub_cmd(script=real_upath(testSh),
                            num_threads=1, wait_before_exit=True)

    logging.debug("Submitting cmd: " + cmd)
    backticks(cmd)

    if not filecmp.cmp(consensusFa, GCON_OUT_FA):
        errMsg = "Trouble running qsub or output is not as " + \
                 "expected ({0} and {1} must agree). Abort!".format(
                     consensusFa, GCON_OUT_FA)
        logging.error(errMsg)
        return False
    else:
        shutil.rmtree(testDir)
        logging.info("sge and gcon check passed.")
        return True


def eval_blasr_alignment(record, qver_get_func, qvmean_get_func,
                         sID_starts_with_c, qv_prob_threshold, debug=False):
    """
    Takes a BLASRRecord (blasr -m 5) and goes through the
    alignment string
    ex: |||**||||**|||*|*|
    to determine the sequence of 'M' (matches), 'S' (sub), 'I', 'D'

    qver_get_func --- could be either basQV.basQVcacher.get() or
        basQV.basQVcacher.get_smoothed()

    For any non-match, if either or both query/target's QV indicate
    that the event ('S', 'I', 'D') is expected
    (ex: insertion prob >= qv_prob_threshold),
    then it does not count as a penalty.

    Returns: cigar string, binary ECE array

    NOTE: long insertions/deletions are still a difficult problem
    because alignments can be arbitrary right now the quick solution is:
          use probqv get_smoothed
    however, with homopolymers, penalization can still happen unless
    I write code to check specifically for homopolymers, (otherwise the
    cigar_str[-1]=='D' or 'I' sets in). -- Liz
    """
    if debug:
        import pdb
        pdb.set_trace()

    if record.qStrand == '+':
        query_qver_get_func = lambda _name, _pos: qver_get_func(record.qID, _name, min(_pos+record.qStart, record.qLength-1))
        query_qver_get_func1 = lambda _name, _pos: qver_get_func(record.qID, _name, min(_pos+record.qStart+1, record.qLength-1))
    elif record.qStrand == '-':
        query_qver_get_func = lambda _name, _pos: qver_get_func(record.qID, _name, max(0, record.qEnd-1-_pos))
        query_qver_get_func1 = lambda _name, _pos: qver_get_func(record.qID, _name, max(0, record.qEnd-1-_pos-1))
    else:
        raise Exception, "Unknown strand type {0}".format(record.qStrand)

    if record.sStrand == '+':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart, record.sLength-1))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart+1, record.sLength-1))
    elif record.sStrand == '-':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos-1))
    else:
        raise Exception, "Unknown strand type {0}".format(record.sStrand)

    # if mean_qv_for_q|s is not given, always revert back to qv_prob_threshold
    if qvmean_get_func is None:
        mean_qv_for_q = {'D': qv_prob_threshold, 'S': qv_prob_threshold, 'I': qv_prob_threshold}
        mean_qv_for_s = None if sID_starts_with_c else {'D': qv_prob_threshold, 'S': qv_prob_threshold, 'I': qv_prob_threshold}
    else:
        mean_qv_for_q = {'D': qvmean_get_func(record.qID, 'DeletionQV'),
                         'I': qvmean_get_func(record.qID, 'InsertionQV'),
                         'S': qvmean_get_func(record.qID, 'SubstitutionQV')}
        if sID_starts_with_c:
            mean_qv_for_s = None
        else:
            mean_qv_for_s = {'D': qvmean_get_func(record.sID, 'DeletionQV'),
                            'I': qvmean_get_func(record.sID, 'InsertionQV'),
                            'S': qvmean_get_func(record.sID, 'SubstitutionQV')}

    q_index = 0
    s_index = 0
    last_state, last_tracking_nt, homopolymer_so_far = None, None, False
    cigar_str = ''
    # binary array of 0|1 where 1 is a penalty
    ece = np.zeros(len(record.alnStr), dtype=np.int)
    for offset, nt_aln in enumerate(record.alnStr):
        if nt_aln == '|':  # match
            cigar_str += 'M'
            q_index += 1
            s_index += 1
            last_state = 'M'
        elif record.qAln[offset] == '-':  # deletion
            # for deletion, cases where consider a non-match
            # (1) IF last position was not "D"
            #        case 1a: both query and subject has very good prob
            # (2) IF last position was "D" (q_index did not advance)
            #        case 2a: both query and subject has very good prob
            #        case 2b: subject has good prob; query has bad prob AND last-cur pos is NOT homopolymer
            # case 1a and case 2a have the same condition
            # case 2b arises because the only possible explanation would have been query have bad prob,
            #   but it was used to explain the last deletion and the advanced S nucleotide is diff from the last S
            s_is_good = sID_starts_with_c or subject_qver_get_func('InsertionQV', s_index) < mean_qv_for_s
            q_is_good = query_qver_get_func1('DeletionQV', q_index) < mean_qv_for_q
            if last_state != 'D': # entering D state now, record s
                last_tracking_nt = record.sAln[offset]
                homopolymer_so_far = True
                if (s_is_good and q_is_good):
                    ece[offset] = 1
            else:  # already in D state, which means q_index did not advance, hence q_is_good is the same
                homopolymer_so_far = (record.sAln[offset] == last_tracking_nt)
                if (s_is_good and (q_is_good or not homopolymer_so_far)):
                    ece[offset] = 1

            cigar_str += 'D'
            s_index += 1
            last_state = 'D'
        elif record.sAln[offset] == '-':  # insertion
            # for insertion, cases where consider a non-match
            # (1) IF last position was not "I"
            #        case 1a: both query and subject has very good prob
            # (2) IF last position was "I" (s_index did not advance)
            #        case 2a: both query and subject has very good prob
            #        case 2b: query has good prob; subject has bad prob AND last-cur pos is NOT homopolymer
            q_is_good = query_qver_get_func('InsertionQV', q_index) < mean_qv_for_q
            s_is_good = sID_starts_with_c or subject_qver_get_func1('DeletionQV', s_index) < mean_qv_for_s

            if last_state != 'I':
                last_tracking_nt = record.qAln[offset]
                homopolymer_so_far = True
                if (q_is_good and s_is_good):
                    ece[offset] = 1
            else: # already in "I" state, s_index did not advance, s_is_good is the same
                homopolymer_so_far = (record.qAln[offset] == last_tracking_nt)
                if q_is_good and (s_is_good or not homopolymer_so_far):
                    ece[offset] = 1

            cigar_str += 'I'
            q_index += 1
            last_state = 'I'
        else:  # substitution
            cigar_str += 'S'
            if query_qver_get_func('SubstitutionQV', q_index) < mean_qv_for_q and \
               (sID_starts_with_c or subject_qver_get_func('SubstitutionQV', s_index) < mean_qv_for_s):
                ece[offset] = 1
            q_index += 1
            s_index += 1
            last_state = 'S'

    return cigar_str, ece


class HitItem(object):

    """
    Simply define an object class for saving items produced by
    blasr_against_ref or daligner_against_ref.
    """

    def __init__(self, qID, cID, qStart=None, qEnd=None,
                 missed_q=None, missed_t=None,
                 fakecigar=None, ece_arr=None):
        self.qID = qID
        self.cID = cID
        self.qStart = qStart
        self.qEnd = qEnd
        self.missed_q = missed_q
        self.missed_t = missed_t
        self.fakecigar = fakecigar
        self.ece_arr = ece_arr

    def __str__(self):
        return """{qID}/{qStart}_{qEnd} aligns to {cID}""".format(
                qID=self.qID.split(' ')[0], cID=self.cID.split(' ')[0],
                qStart=self.qStart, qEnd=self.qEnd)

def blasr_against_ref(output_filename, is_FL, sID_starts_with_c,
                      qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                      ece_penalty=1, ece_min_len=20, same_strand_only=True,
                      max_missed_start=200, max_missed_end=50):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    qver_get_func --- should be basQV.basQVcacher.get() or
                      .get_smoothed(), or can just pass in
                      lambda (x, y): 1. to ignore QV
    """
    with BLASRM5Reader(output_filename) as reader:
        for r in reader:
            missed_q = r.qStart + r.qLength - r.qEnd
            missed_t = r.sStart + r.sLength - r.sEnd

            if sID_starts_with_c:
                # because all consensus should start with
                # c<cluster_index>
                assert r.sID.startswith('c')
                if r.sID.find('/') > 0:
                    r.sID = r.sID.split('/')[0]
                if r.sID.endswith('_ref'):
                    # probably c<cid>_ref
                    cID = int(r.sID[1:-4])
                else:
                    cID = int(r.sID[1:])
            else:
                cID = r.sID

            # self hit, useless!
            # low identity not allowed
            # opposite strand not allowed!
            if (cID == r.qID or
                    r.identity < 70. or
                    (r.strand == '-' and same_strand_only)):
                yield HitItem(qID=r.qID, cID=cID)
                continue

            # full-length case: allow up to max_missed_start bp of 5' not aligned
            # and max_missed_end bp of 3' not aligned
            # non-full-length case: not really tested...don't use
            if is_FL and (r.sStart > max_missed_start or r.qStart > max_missed_start or
                          (r.sLength - r.sEnd > max_missed_end) or
                          (r.qLength - r.qEnd > max_missed_end)):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                cigar_str, ece_arr = eval_blasr_alignment(
                    record=r,
                    qver_get_func=qver_get_func,
                    qvmean_get_func=qvmean_get_func,
                    sID_starts_with_c=sID_starts_with_c,
                    qv_prob_threshold=qv_prob_threshold)

                if alignment_has_large_nonmatch(ece_arr,
                                                ece_penalty, ece_min_len):
                    yield HitItem(qID=r.qID, cID=cID)
                else:
                    yield HitItem(qID=r.qID, cID=cID,
                                  qStart=r.qStart, qEnd=r.qEnd,
                                  missed_q=missed_q * 1. / r.qLength,
                                  missed_t=missed_t * 1. / r.sLength,
                                  fakecigar=cigar_str,
                                  ece_arr=ece_arr)


def daligner_against_ref(query_dazz_handler, target_dazz_handler, la4ice_filename,
                         is_FL, sID_starts_with_c,
                         qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                         ece_penalty=1, ece_min_len=20, same_strand_only=True, no_qv_or_aln_checking=False,
                         max_missed_start=200, max_missed_end=50):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    Parameters:
      query_dazz_handler - query dazz handler in DalignRunner
      target_dazz_handler - target dazz handler in DalignRunner
      la4ice_filename - la4ice output of DalignRunner
      qver_get_func - returns a list of qvs of (read, qvname)
                      e.g. basQV.basQVcacher.get() or .get_smoothed()
      qvmean_get_func - which returns mean QV of (read, qvname)
    """
    for r in LA4IceReader(la4ice_filename):
        missed_q = r.qStart + r.qLength - r.qEnd
        missed_t = r.sStart + r.sLength - r.sEnd

        r.qID = query_dazz_handler[r.qID].split(' ')[0]
        r.sID = target_dazz_handler[r.sID].split(' ')[0]

        if sID_starts_with_c:
            # because all consensus should start with
            # c<cluster_index>
            assert r.sID.startswith('c')
            if r.sID.find('/') > 0:
                r.sID = r.sID.split('/')[0]
            if r.sID.endswith('_ref'):
                # probably c<cid>_ref
                cID = int(r.sID[1:-4])
            else:
                cID = int(r.sID[1:])
        else:
            cID = r.sID

        # self hit, useless!
        # (identity is removed here, NOT trustworthy using Jason's code calculations)
        # opposite strand not allowed!
        if (cID == r.qID or (r.strand == '-' and same_strand_only)):
            yield HitItem(qID=r.qID, cID=cID)
            continue

        # this is used for partial_uc/nFL reads only
        # simply accepts hits from daligner for the nFL partial hits
        # testing shows that it does not affect much the Quiver consensus calling
        if no_qv_or_aln_checking:
            yield HitItem(qID=r.qID, cID=cID,
                          qStart=r.qStart, qEnd=r.qEnd,
                          missed_q=missed_q * 1. / r.qLength,
                          missed_t=missed_t * 1. / r.sLength,
                          fakecigar=1,
                          ece_arr=1)
            continue

        # full-length case: allow up to 200bp of 5' not aligned
        # and 50bp of 3' not aligned
        if (is_FL and (r.sStart > max_missed_start or r.qStart > max_missed_start or
                       (r.sLength - r.sEnd > max_missed_end) or
                       (r.qLength - r.qEnd > max_missed_end))):
            yield HitItem(qID=r.qID, cID=cID)
        else:
            cigar_str, ece_arr = eval_blasr_alignment(
                record=r,
                qver_get_func=qver_get_func,
                sID_starts_with_c=sID_starts_with_c,
                qv_prob_threshold=qv_prob_threshold,
                qvmean_get_func=qvmean_get_func)
            #else: # don't use QV, just look at alignment

            if alignment_has_large_nonmatch(ece_arr, ece_penalty, ece_min_len):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                yield HitItem(qID=r.qID, cID=cID,
                              qStart=r.qStart, qEnd=r.qEnd,
                              missed_q=missed_q * 1. / r.qLength,
                              missed_t=missed_t * 1. / r.sLength,
                              fakecigar=cigar_str,
                              ece_arr=ece_arr)

def alignment_has_large_nonmatch(ece_arr, penalty, min_len):
    """
    penalty of (-)1: 50%
    penalty of (-)2: 66%
    penalty of (-)4: 80%
    penalty of (-)9: 90%

    Return True when alignment has large non-matches not explained
    by low base QVs (in other words, "reject" as an isoform hit and
    don't put in the same cluster)
    """
    ece_arr = ece_arr * (penalty + 1)
    s = [0] + list(ece_arr - penalty)
    # fix this later to something faster & better
    return len(findECE(s, len(s), min_len, True)) > 0


def possible_merge(r, ece_penalty, ece_min_len,
                   max_missed_start=200, max_missed_end=50):
    """
    r --- BLASRM5Record
    Criteria:
    (1) identity >= 90% and same strand
    (2) check criteria for how much is allowed to differ on the
        5' / 3' ends
    """
    if r.sID == r.qID or r.identity < 90 or r.strand == '-':
        return False
    # intentional here to prevent disrupting future ICE runs
    # MORE lenient on 5' but NOT on 3'
    if ((r.qLength - r.qEnd) > max_missed_end or (r.sLength - r.sEnd) > max_missed_end or
            r.qStart > max_missed_start or r.sStart > max_missed_start):
        return False

    arr = np.array([(x == '*') * 1 for x in r.alnStr])
    if alignment_has_large_nonmatch(ece_arr=arr,
                                    penalty=ece_penalty,
                                    min_len=ece_min_len):
        return False
    return True


def get_the_only_fasta_record(fa):
    """Input fasta file should contain exactly one FastaRecord,
    return the fastas record."""
    rs = [r for r in FastaReader(fa)]
    if len(rs) != 1:
        errMsg = "Cluster fasta file {fa} must contain only one read.".\
            format(fa=fa)
        raise ValueError(errMsg)
    return rs[0]


"""
The following methods was originally created by jchin:
    ~jchin/depot_mp27/jchin/rset_quvier.py
, and then modified by etseng.

Input: input.fasta.fofn (shared),
       per-cluster in.fasta,
       per-cluster g_consensus.fasta

-- input.fasta.fofn should be raw fasta files
   (pls2fasta -maskRegion) of input.fofn

Within each cluster:
1) create in.raw.fasta based on input.fasta.fofn & in.fasta,
   putting in raw (unrolled) fasta of each ZMW
2) blasr (1) to g_consensus.fasta output as SAM

This is faster than using regions.fofn because it still reads
through the whole .bax.h5 files
"""


def is_blank_sam(samfile):
    """
    return True if the SAM file only has @xx header and NO alignment
    """
    with open(samfile) as f:
        for line in f:
            if not line.startswith('@'):
                return False
    return True


def is_blank_bam(bamfile):
    """
    return True if the BAM file only has @xx header and NO alignment
    """
    try:
        with pysam.Samfile(bamfile, 'rb', check_sq=False) as f:
            f.next()
        f.close()
    except StopIteration:
        return True
    return False


def concat_sam(samfiles, outsam_filename):
    """
    Header looks like:
    @HD     VN:1.3.1
    @SQ     SN:c31  LN:3104 M5:ef7d3f84dea9d9face43e6fd5b6336c4
    @RG     ID:2caa54eef6   PU:in.raw_with_partial.fasta       SM:NO_CHIP_ID
    @PG     ID:BLASR        VN:1.3.1.126469 CL:blasr in.raw_with_partial.fasta g_consensus.fasta -nproc 12 -bestn 5 -nCandidates 10 -sam -out out.sam

    NOTE: check for M5 conflicts; manipulate them if it conflicts
    """
    f_sq = open(outsam_filename + '.sq', 'w')
    f_bd = open(outsam_filename + '.bd', 'w')

    rg_line = None
    pg_line = None

    md5_seen = set()

    if len(samfiles) == 0:
        raise ValueError("No sam input files to concatenate.")

    h = open(samfiles[0])
    line = h.readline()
    assert line.startswith('@HD')
    f_sq.write(line)
    line = h.readline()
    assert line.startswith('@SQ')
    line = h.readline()
    assert line.startswith('@RG')
    rg_line = line  # write at the end
    line = h.readline()
    assert line.startswith('@PG')
    pg_line = line  # write at the end
    h.close()

    for f in samfiles:
        with open(f) as h:
            assert h.readline().startswith('@HD')
            line = h.readline()
            assert line.startswith('@SQ')
            # ------- check for MD5 conflicts ----------- #
            m5 = line.strip().split()[-1]
            assert m5.startswith("M5:")
            if m5 not in md5_seen:
                f_sq.write(line)
                md5_seen.add(m5)
            else:
                s = list(m5[3:])
                while True:
                    # create a random m5 string.
                    random.shuffle(s)
                    s = "".join(s)
                    if s not in md5_seen:
                        break
                line = line[:line.find('M5:')] + 'M5:' + s + '\n'
                logging.debug("MD5 conflict: change to {0}".format(s))
                md5_seen.add(s)
                f_sq.write(line)
            # ----- end MD5 checking and writing --------- #
            assert h.readline().startswith('@RG')
            assert h.readline().startswith('@PG')
            for line in h:
                f_bd.write(line)

    f_bd.close()
    f_sq.write(rg_line)
    f_sq.write(pg_line)
    f_sq.close()

    cmd = "cat {0}.sq {0}.bd > {0}".format(real_upath(outsam_filename))
    execute(cmd=cmd,
            errmsg="Failed to concat sam files! Abort.",
            errcls=IOError)

    os.remove(f_sq.name)
    os.remove(f_bd.name)


def concat_bam_header(in_fns, out_fn=None):
    """Concat bam headers in in_fns and save to out_fn,
    return the merged bam header as an instance of BamHeader"""
    from pbtranscript.io import BamHeader
    h = BamHeader(ignore_pg=True)
    for in_fn in in_fns:
        s = pysam.Samfile(in_fn, 'rb')
        h.add(s.header)
        s.close()
    if out_fn is not None:
        o = pysam.Samfile(out_fn, 'wh', header=h.header)
        o.close()
    return h


def concat_bam(in_fns, out_fn):
    """Concat input bam files to an output bam file.
    Note that each input bam has ONLY one reference sequence.
    """
    # construct sam header
    h = concat_bam_header(in_fns)
    o = BamWriter(out_fn, header=h)
    for index, in_fn in enumerate(in_fns):
        s = pysam.Samfile(in_fn, 'rb')
        for r in s:
            r.tid = index # Overwrite tid !!!
            o.write(r)
        s.close()
    o.close()


def convert_fofn_to_fasta(fofn_filename, out_filename, fasta_out_dir,
                          force_overwrite=False):
    """
    For each .bax.h5 file, create .bax.h5.fasta file and save paths to
    out_filename, which should usually be 'input.fasta.fofn'
    Modified: 09/14/2015, both ends of subreads in fasta files will
    be trimmed in IceQuiver (trim_and_write_raw_file) instead of here.
    """
    logging.info("Converting fofn {fofn} to fasta.".format(fofn=fofn_filename))
    in_fns = get_files_from_file_or_fofn(fofn_filename)
    out_fns = []
    mkdir(fasta_out_dir)
    for in_fn in in_fns:
        logging.debug("converting h5 file: {f}.".format(f=in_fn))
        if not (in_fn.endswith('.bax.h5') or in_fn.endswith('.bas.h5')):
            raise ValueError("fofn file {fofn} ".format(fofn=fofn_filename) +
                             "should only contain bax/bas.h5 files.")

        # e.g. m111xxxx.1.bax.h5 ==>
        #      tmp_out_file = m11xxxx.1.bax.h5.fasta.tmp
        #      out_file = m11xxxx.1.bax.h5.fasta
        in_basename = op.basename(in_fn)
        out_file = op.join(fasta_out_dir, in_basename + '.fasta')
        if op.exists(out_file) and not force_overwrite:
            logging.debug("File {0} already exists. skipping.".format(out_file))
        else:
            cmd = "pls2fasta {in_fn} ".format(in_fn=real_upath(in_fn)) + \
                  " {out} ".format(out=real_upath(out_file)) + \
                  "-minSubreadLength 300 -minReadScore 750 -trimByRegion"
            execute(cmd=cmd)
        out_fns.append(out_file)
    write_files_to_fofn(out_fns, out_filename)



def build_sa(input_fasta, out_sa):
    """Generate suffix array of input_fasta"""
    if op.exists(input_fasta):
        cmd = "sawriter {o} {i} -blt 8 -welter ".\
            format(o=real_upath(out_sa), i=real_upath(input_fasta))
        dummy_out, code, dummy_msg = backticks(cmd)
        if code == 0:
            return True
        else:
            # If failed to generate suffix array, warning.
            logging.warn("Unable to create suffix array for {f}.".format(f=input_fasta))
            return False
    else:
        raise IOError("Unable to find fasta file {f}.".format(f=input_fasta))


def trim_subreads_and_write(reader, in_seqids, out_file, trim_len, min_len,
                            ignore_keyerror=False, bam=False):
    """ Extract (dump) raw subreads of every zmws from in_seqeids from reader
    to out_file.
        reader --- provides random access to raw subreads in input file.
                   type = MetaSubreadFastaReader, when input files are in FASTA,
                   and reads are in format <movie>/<holeNumber>/<subread or CCS>.
                   type = BamCollection, when input files are in BAM.
        trim_len --- trim the first and last n bases when input is BAM
        min_len --- minimum read length to write a subread when input is BAM
        in_seqids --- zmw ids to dump
        out_file --- a FASTA file when input files are in FASTA; a BAM file when
                     input files are in BAM.
        return movies seen
    """
    movies = set()
    zmw_seen = set()
    f = None # output open file handler

    if bam:
        assert isinstance(reader, BamCollection)
        f = BamWriter(out_file, reader.header)
    else:
        assert isinstance(reader, MetaSubreadFastaReader)
        f = FastaWriter(out_file)

    for seqid in in_seqids:
        zmw = seqid
        try:
            zmw = '/'.join(seqid.split('/')[0:2])
        except ValueError:
            raise ValueError("%s does not contain a valid pacbio zmw id." % seqid)

        if zmw not in zmw_seen:
            movies.add(zmw.split('/')[0])
            zmw_seen.add(zmw)
            try:
                if bam:
                    for rec in reader[zmw].subreads:
                        if len(rec) >= 2*trim_len + min_len:
                            f.write(rec.Clip(rec.readStart+trim_len,
                                             rec.readEnd-trim_len))
                else:
                    for rec in reader[zmw]:
                        if len(rec) >= 2*trim_len + min_len:
                            try:
                                m, hn, s_e = rec.name.split('/')
                                s, e = [int(x) for x in s_e.split('_')]
                                new_id = "%s/%s/%d_%d" % (m, hn, s+trim_len, e-trim_len)
                                f.writeRecord(new_id, rec.sequence[trim_len:-trim_len])
                            except ValueError:
                                raise ValueError("%s is not a valid pacbio subread." % rec.name)
            except KeyError:
                if ignore_keyerror:
                    logging.warning("Ignoring {zmw} because the input FASTA/BAM ".
                                    format(zmw=zmw) + " does not contain it.")
                else:
                    raise ValueError("{0} doesn't exist. Abort!".format(zmw))
    f.close()

    return movies


def blasr_for_quiver(query_fn, ref_fasta, out_fn, bam=False,
                     run_cmd=True, blasr_nproc=12):
    """
    query_fn  --- should be in.raw.fasta|bam
    ref_fasta --- reference fasta (ex: g_consensus.fasta) to align to
    out_fn    --- sam|bam output aligning query_fn to ref_fasta

    blasr query_fn ref_fasta -out out_fn -sam -clipping soft
    blasr query_fn ref_fasta -out out_fn -bam
    """
    cmd = "blasr {i} ".format(i=real_upath(query_fn)) + \
          "{r} ".format(r=real_upath(ref_fasta)) + \
          "--nproc {n} ".format(n=blasr_nproc) + \
          "--bestn 5 --nCandidates 10 " + \
          ("--sam --clipping soft " if not bam else "--bam ") + \
          "--out {o} ".format(o=real_upath(out_fn)) + \
          "1>/dev/null 2>/dev/null"
    if run_cmd:
        execute(cmd)
    else:
        logging.debug("CMD: " + cmd)
    return cmd


def num_reads_in_fasta(in_fa):
    """Return the number of reads in the in_fa fasta file."""
    if (not in_fa.endswith(".fa")) and (not in_fa.endswith(".fasta")):
        # if not a fasta file, must be a contigset xml
        if not in_fa.endswith(".xml"):
            raise IOError("%s must be a FASTA or ContigSet file." % in_fa)
        return ContigSet(in_fa).numRecords

    if not op.exists(in_fa):
        raise IOError("fasta file {f} does not exist.".format(f=in_fa))
    cmd = "grep '>' {f} | wc -l ".format(f=in_fa)
    logging.debug("CMD: " + cmd)
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        raise RuntimeError("CMD failed: {cmd}\n{e}".
                           format(cmd=cmd, e=_msg))
    return int(_out[0])


def combine_nfl_pickles(splitted_pickles, out_pickle):
    """Combine splitted nfl pickles to a big pickle."""
    logging.debug("Cominbing {N} nfl pickles: {ps} ".
                  format(N=len(splitted_pickles),
                         ps=",".join(splitted_pickles)) +
                  " into a big pickle {p}.".format(p=out_pickle))

    if len(splitted_pickles) == 1:
        logging.debug("Copying the only given pickle to out_pickle.")
        if realpath(splitted_pickles[0]) != realpath(out_pickle):
            shutil.copyfile(splitted_pickles[0], out_pickle)
    else:
        # Combine all partial outputs
        logging.debug("Merging all pickles.")
        partial_uc = defaultdict(lambda: [])
        nohit = set()
        for pf in splitted_pickles:
            logging.debug("Merging {pf}.".format(pf=pf))
            a = load(open(pf))
            nohit.update(a['nohit'])
            for k, v in a['partial_uc'].iteritems():
                partial_uc[k] += v

        logging.debug("Dumping all to {f}".format(f=out_pickle))
        # Dump to one file
        partial_uc = dict(partial_uc)
        with open(out_pickle, 'w') as f:
            dump({'nohit': nohit, 'partial_uc': partial_uc}, f)
        logging.debug("{f} created.".format(f=out_pickle))


def cid_with_annotation(cid):
    """Given a cluster id, return cluster id with human readable annotation.
    e.g., c0 --> c0 isoform=c0
          c0/89/3888 -> c0/89/3888 isoform=c0;full_length_coverage=89;isoform_length=3888
          c0/f89p190/3888 -> c0/f89p190/3888 isoform=c0;full_length_coverage=89;non_full_length_coverage=190;isoform_length=3888
    """
    fields = cid.split('/')
    short_id, fl_coverage, nfl_coverage, seq_len = None, None, None, None
    if len(fields) != 1 and len(fields) != 3:
        raise ValueError("Not able to process isoform id: {cid}".format(cid=cid))
    short_id = fields[0]
    if len(fields) == 3:
        seq_len = fields[2]
        if "f" in fields[1]:
            if "p" in fields[1]: # f89p190
                fl_coverage = fields[1].split('p')[0][1:]
                nfl_coverage = fields[1].split('p')[1]
            else: # f89
                fl_coverage = fields[1][1:]
        else:
            fl_coverage = fields[1]

    annotations = ["isoform={short_id}".format(short_id=short_id)]
    if fl_coverage is not None:
        annotations.append("full_length_coverage={fl}".format(fl=fl_coverage))
    if nfl_coverage is not None:
        annotations.append("non_full_length_coverage={nfl}".format(nfl=nfl_coverage))
    if seq_len is not None:
        annotations.append("isoform_length={l}".format(l=seq_len))

    return "{cid} {annotation}".format(cid=cid, annotation=";".join(annotations))


class _Parsed_Read_Name(object):
    """ An internal class for parsing read names used in isoseq.
            m...../1234/CCS       --> ccs sequence[:], strand +, get_all True
            m...../1234/0_100     --> subread sequence[0:100], strand +, get_all False
            m...../1234/0_100_CCS --> ccs sequence[0:100], strand +, get_all False
            m...../1234/100_0_CCS --> ccs sequence[100:0], strand -, get_all False
    """
    def __init__(self, seqid):
        try:
            self.movie, self.hn, s_e = seqid.split('/')
            self.hn = int(self.hn)
            self.strand = '+'
            self.s, self.e = 0, 0
            self.is_CCS = True
            if s_e.lower() == "ccs": # m.../1234/CCS
                pass
            elif s_e.endswith('_CCS'): # m.../1234/0_100_CCS
                self.s, self.e = s_e.split('_')[:2]
            else:
                self.s, self.e = s_e.split('_') # m.../1234/0_100
                self.is_CCS = False
            self.s, self.e = int(self.s), int(self.e)
            if self.s > self.e:
                self.s, self.e = self.e, self.s
                self.strand = '-'
        except ValueError:
            raise ValueError("{seqid} is not a valid CCS read".
                             format(seqid=seqid))

    @property
    def get_all(self):
        """Whether this read name represent a full zmw or just part of it."""
        return self.s == self.e == 0

    @property
    def zmwName(self):
        """Return zmw name."""
        return "%s/%d" % (self.movie, self.hn)

    def __repr__(self):
        s_e = "%d_%d" % (self.s, self.e) if self.strand == "+" else \
             "%d_%d" % (self.e, self.s)
        if self.is_CCS:
            return "%s/CCS" % self.zmwName if self.get_all else \
                   "%s/%s_CCS" % (self.zmwName, s_e)
        else:
            return self.zmwName if self.get_all else \
                   "%s/%s" % (self.zmwName, s_e)


def get_qvs_from_bam(reader, parsed_read_name, qv_name):
    """Get qvs from a bam file."""
    if not parsed_read_name.is_CCS:
        raise IOError("get_qvs_from_bam does not support non-ccs bam.")

    zmw = reader[parsed_read_name.zmwName]
    if zmw.ccsRead is None:
        raise IOError("Could not find ccs of zmw %s" % zmw.zmwName)
    qvs = zmw.ccsRead.qv(qv_name)

    if not parsed_read_name.get_all:
        qvs = qvs[parsed_read_name.s:parsed_read_name.e]
    if parsed_read_name.strand == "-":
        qvs = qvs[::-1]
    return qvs


def get_qv_from_bas_handler(bas_handler, parsed_read_name, qv_name):
    """Read QV of type qv_name for movie/hn/s_e from bas_h5 file handler."""
    zmw = bas_handler[parsed_read_name.hn]
    if parsed_read_name.is_CCS:
        if zmw.ccsRead is not None:
            qvs = zmw.ccsRead.qv(qv_name)
        else:  # this is for CCS w/ 0-passed
            qvs = zmw.read().qv(qv_name)
    else:  # subread
        qvs = zmw.read().qv(qv_name)

    if not parsed_read_name.get_all:
        qvs = qvs[parsed_read_name.s:parsed_read_name.e]
    if parsed_read_name.strand == '-':
        qvs = qvs[::-1]
    return qvs


def fafn2fqfn(fafn):
    """return a FASTQ file name which corresponds to a FASTA file."""
    if fafn.find('.') != -1:
        return fafn[:fafn.rfind('.')] + ".fastq"
    else:
        return fafn + ".fastq"


def ice_fa2fq(in_fa, ccs_fofn, out_fq):
    """Convert an input FASTA file to an output FASTQ file,
       reading QVs from the input ccs.h5, ccs.bam or ccs FOFN.
    """
    ccs_fns = get_files_from_file_or_fofn(ccs_fofn)
    fmt = guess_file_format(ccs_fns)

    if fmt == FILE_FORMATS.H5:
        qver = basQVcacher()
        for ccs_fn in ccs_fns:
            qver.add_bash5(ccs_fn)
        bas_handlers = {}
    elif fmt == FILE_FORMATS.BAM:
        qver = BamCollection(*ccs_fns)
    else:
        raise IOError("ice_fa2fq does not support input %s." %
                      ccs_fofn)

    with ContigSetReaderWrapper(in_fa) as reader, \
            FastqWriter(out_fq) as writer:
        for r in reader:
            logging.debug("Getting QVs for {name} ...".format(name=r.name))
            seqid = r.name.split(' ')[0]
            parsed_read_name = _Parsed_Read_Name(seqid)
            if fmt == FILE_FORMATS.H5:
                try:
                    bas_file = qver.bas_files[parsed_read_name.movie][seqid]
                    if bas_file not in bas_handlers:
                        bas_handlers[bas_file] = BasH5Reader(bas_file)
                except KeyError:
                    raise IOError("Could not read {s} from {f}.".
                                  format(s=seqid, f=ccs_fofn))
                qvs = get_qv_from_bas_handler(bas_handler=bas_handlers[bas_file],
                                              parsed_read_name=parsed_read_name,
                                              qv_name="QualityValue")
            elif fmt == FILE_FORMATS.BAM:
                qvs = get_qvs_from_bam(reader=qver,
                                       parsed_read_name=parsed_read_name,
                                       qv_name="QualityValue")
            else:
                assert False

            if len(r.sequence) != len(qvs):
                raise ValueError("Sequence and QVs of {r} should be the same!".
                                 format(r=r.name))
            writer.writeRecord(r.name, r.sequence[:], qvs)

    if fmt == FILE_FORMATS.H5:
        for bas_file, bas_handler in bas_handlers.iteritems():
            logging.debug("Closing {bas_file} ...".format(bas_file=bas_file))
            bas_handler.close()
    elif fmt == FILE_FORMATS.BAM:
        qver.close()


def set_probqv_from_ccs(ccs_fofn, fasta_filename):
    """Set probability and quality values from ccs.h5,
    return probqv, log_info."""
    assert ccs_fofn is not None and fasta_filename is not None
    start_t = time.time()
    probqv = ProbFromQV(input_fofn=ccs_fofn,
                        fasta_filename=fasta_filename)
    msg = "Loading probabilities and QV from " + \
          "{f} + {c} took {t} sec.".format(f=fasta_filename, c=ccs_fofn,
                                           t=(time.time()-start_t))
    return probqv, msg


def set_probqv_from_model():
    """Set probablitiy values from a fixed model,
    return probqv, log_info.
    """
    msg = "Loading predefined probabilities model."
    probqv = ProbFromModel(0.01, 0.07, 0.06)
    return probqv, msg


def set_probqv_from_fq(fastq_filename):
    """Set probability and QVs from FASTQ, return probqv, log_info."""
    assert isinstance(fastq_filename, str)
    start_t = time.time()
    probqv = ProbFromFastq(fastq_filename=fastq_filename)
    msg = "Loading QVs from {f} took {t} sec.".\
          format(f=fastq_filename, t=(time.time()-start_t))
    return probqv, msg


def write_cluster_report(report_fn, uc, partial_uc):
    """
    Write a CSV report to report_fn, each line contains three columns:
        cluster_id, read_id and read_type
    """
    with open(report_fn, 'w') as f:
        f.write("cluster_id,read_id,read_type\n")
        for c in uc.keys():
            for r in uc[c]:
                f.write("c{c},{r},FL\n".format(r=r, c=c))
            if partial_uc is not None and c in partial_uc.keys():
                for r in partial_uc[c]:
                    f.write("c{c},{r},NonFL\n".format(r=r, c=c))

