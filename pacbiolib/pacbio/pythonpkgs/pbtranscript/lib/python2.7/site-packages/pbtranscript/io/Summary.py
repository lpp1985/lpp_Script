"""Define ClassifySummary and ClusterSummary."""

import os.path as op
from collections import OrderedDict

from pbcore.io import ContigSet
from pbcommand.models.report import Report, Attribute


class Summary(object):
    REPORT_ID = None # used for pbcommand report model
    ATTR_LABELS = {}

    def __init__(self):
        pass

    """Supert class for ClassifySummary and ClusterSummary."""
    @property
    def fieldsIDs(self):
        """Specific IDs for Classify and Cluster."""
        raise NotImplementedError("Summary.fieldsIDs() not implemented")

    @property
    def fieldsNames(self):
        """Return all fields names in a list."""
        return [self.ATTR_LABELS[fsid] for fsid in self.fieldsIDs]

    @property
    def fields(self):
        """Return fiels values in a list. Have to match self.fieldsNames"""
        return [getattr(self, name) for name in self.fieldsIDs]

    def __str__(self):
        assert len(self.fieldsNames) == len(self.fields)
        return "\n".join(["{name}={val}".format(name=name, val=val)
                          for name, val in zip(self.fieldsNames, self.fields)])

    def to_report(self, dataset_uuids=()):
        """Convert a summary object to pbcommand.report object."""
        attributes = [Attribute(id_=attribute_id, value=attribute_val, name=attribute_name)
                      for attribute_id, attribute_name, attribute_val
                      in zip(self.fieldsIDs, self.fieldsNames, self.fields)]
        return Report(id_=self.REPORT_ID, attributes=attributes,
                      dataset_uuids=dataset_uuids)

    def write(self, outFile, dataset_uuids=()):
        """Write summary to outFile."""
        with open(outFile, 'w') as writer:
            if outFile.endswith(".json"):
                writer.write(self.to_report(dataset_uuids).to_json())
            else:
                writer.write(self.__str__() + "\n")


class ClassifySummary(Summary):
    REPORT_ID = "pbtranscript_classify"
    ATTR_LABELS = OrderedDict([
        ("num_reads", "Number of consensus reads"),
        ("num_5_seen", "Number of five prime reads"),
        ("num_3_seen", "Number of three prime reads"),
        ("num_polya_seen", "Number of poly-A reads"),
        ("num_filtered_short_reads", "Number of filtered short reads"),
        ("num_nfl", "Number of non-full-length reads"),
        ("num_fl", "Number of full-length reads"),
        ("num_flnc", "Number of full-length non-chimeric reads"),
        ("avg_flnc_len", "Mean full-length non-chimeric read length"),
        ("num_nflnc", "Number of non-full-length non-chimeric reads"),
        ("num_nflc", "Number of non-full-length chimeric reads"),
        ("num_flnc_bases", "Number of full-length non-chimeric bases")
    ])
    ATTR_DESCRIPTIONS = {
        "num_reads": "Total number of CCS reads in the input file; this will be identical to the value in the CCS report",
        "num_5_seen": "The number of CCS reads with a 5-prime signal detected",
        "num_3_seen": "The number of CCS reads with a 3-prime signal detected",
        "num_polya_seen": "The number of CCS reads with poly-A and 3-prime signals detected",
        "num_filtered_short_reads": "The number of CCS reads whose read length is less than the specified minimum sequence length",
        "num_nfl": "The number of non-full-length CCS reads; these are missing the poly-A tail and/or a terminal signal",
        "num_fl": "The number of full-length CCS reads. (Full-length reads are reads which have both prime signals and poly-A detected.)",
        "num_flnc": "The number of full-length CCS reads that are not artificial concatemers",
        "avg_flnc_len": "Mean length of full-length, non-artificial-concatemer CCS reads",
        "num_nflnc": "The number of non-full-length non-artificial-concatemer CCS reads",
        "num_nflc": "The number of non-full-length reads that are artificial concatemers",
        "num_flnc_bases": "Total number of bases in full-length non-artificial-concatemer CCS reads"
    }

    """A ClassifySummary object has all classify summary attributes."""

    def __init__(self):
        Summary.__init__(self)

        self.num_reads = 0  # number of reads from in.fasta
        self.num_5_seen = 0  # number 5' primer seen reads within in.fasta
        self.num_3_seen = 0  # number 3' primer seen reads within in.fasta
        self.num_polya_seen = 0      # number polya seen reads within in.fasta

        # number of filtered short reads whose length < min_seq_len
        self.num_filtered_short_reads = 0

        # number of full-length reads, > min_seq_len,
        # can either be chimeric or non-chimeric
        self.num_fl = 0
        # number of full-length non-chimeric reads, > min_seq_len
        self.num_flnc = 0
        # number of full-length chimeric reads
        self.num_flc = 0
        # total number of bases in full-length non-chimeric reads
        self.num_flnc_bases = 0

        # number of non-full-length reads, > min_seq_len,
        # can either be chimeric or non-chimeric
        self.num_nfl = 0
        # number of non-full-length non-chimeric reads
        self.num_nflnc = None
        # number of non-full-length chimeric reads
        self.num_nflc = None

    @property
    def avg_flnc_len(self):
        """Return average read length of full-length non-chimeric reads."""
        if self.num_flnc is None or self.num_flnc == 0:
            return None
        return int(self.num_flnc_bases / self.num_flnc)

    @property
    def fieldsIDs(self):
        attr = [
            "num_reads",
            "num_5_seen",
            "num_3_seen",
            "num_polya_seen",
            "num_filtered_short_reads",
            "num_nfl",
            "num_fl",
            "num_flnc",
            "num_flnc_bases",
        ]
        if self.num_nflnc is not None and self.num_nflc is not None:
            attr.extend(["num_nflnc", "num_nflc"])

        if self.avg_flnc_len is not None:
            attr.extend(["avg_flnc_len"])
        return attr


class ClusterSummary(Summary):
    REPORT_ID = "pbtranscript_cluster"
    ATTR_LABELS = OrderedDict([
        ("num_consensus_isoforms", "Number of unpolished consensus isoforms"),
        ("num_polished_hq_isoforms", "Number of polished high-quality isoforms"),
        ("num_polished_lq_isoforms", "Number of polished low-quality isoforms"),
        ("avg_consensus_isoform_length", "Mean unpolished consensus isoforms read length"),
        ("num_total_bases", "Total number of bases in unpolished consensus isoforms")
    ])
    ATTR_DESCRIPTIONS = {
        "num_consensus_isoforms": "Total number of consensus isoforms, both high- and low-quality",
        "num_polished_hq_isoforms": "The number of consensus isoforms that have an estimated accuracy above the specified cutoff (0.99 default)",
        "num_polished_lq_isoforms": "The number of consensus isoforms that have an estimated accuracy below the specified cutoff",
        "avg_consensus_isoform_length": "The average length of all consensus isoforms, both high- and low-quality",
        "num_total_bases": "Total number of bases in unpolished consensus isoforms"
    }

    """A ClusterSummary object has all cluster summary attributes."""

    def __init__(self):
        Summary.__init__(self)

        self.num_consensus_isoforms = 0  # number of consensus isoforms.
        # total number of bases in predicted consensus isoforms
        self.num_total_bases = 0
        self.num_polished_hq_isoforms = None
        self.num_polished_lq_isoforms = None

    @property
    def avg_consensus_isoform_length(self):
        """Return average read length of predicted consensus isoforms."""
        if self.num_consensus_isoforms is None or self.num_consensus_isoforms == 0:
            return None
        return int(self.num_total_bases / self.num_consensus_isoforms)

    @property
    def fieldsIDs(self):
        """Return field IDs as a list of strings."""
        fsids = ["num_consensus_isoforms"]
        if self.num_polished_hq_isoforms is not None:
            fsids.extend(["num_polished_hq_isoforms"])
        if self.num_polished_lq_isoforms is not None:
            fsids.extend(["num_polished_lq_isoforms"])
        if self.avg_consensus_isoform_length is not None:
            fsids.extend(["avg_consensus_isoform_length"])
        return fsids


def write_cluster_summary(summary_fn, isoforms_fa, hq_fa=None, lq_fa=None):
    """Extract number of consensus isoforms predicted, and total
    number of bases in all consensuus isoforms from isoforms_fa and write
    the two attributes to summary_fn.

    if hq_fa (polished high-quality isoforms) is not None, report
        the number of polished hq clusters
    if lq_fa (polished high-quality isoforms) is not None, report
        the number of polished hq clusters
    """
    try:
        summary = ClusterSummary()
        dataset_uuids = []
        with ContigSet(isoforms_fa) as reader:
            for r in reader:
                summary.num_consensus_isoforms += 1
                summary.num_total_bases += len(r.sequence[:])
            dataset_uuids.append(reader.uuid)

        if hq_fa is not None and op.getsize(hq_fa) > 0:
            summary.num_polished_hq_isoforms = 0
            with ContigSet(hq_fa) as reader:
                for r in reader:
                    summary.num_polished_hq_isoforms += 1
                dataset_uuids.append(reader.uuid)
        if lq_fa is not None and op.getsize(lq_fa) > 0:
            summary.num_polished_lq_isoforms = 0
            with ContigSet(lq_fa) as reader:
                for r in reader:
                    summary.num_polished_lq_isoforms += 1
                dataset_uuids.append(reader.uuid)
        summary.write(summary_fn, dataset_uuids=dataset_uuids)
    except ZeroDivisionError:
        errMsg = "No consensus isoforms predicted."
        logging.error(errMsg)
        raise RuntimeError(errMsg)
