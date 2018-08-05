"""Summarizes depth of coverage from an AlignmentSet file.

Ported from summarize_coverage.py in pbreports/reports, which was ported from
summarizeCoverage.py in pbpy/bin.
"""

from collections import defaultdict
import functools
import logging
import math
import os
import re
import sys
import time

import numpy

from pbcommand.models import TaskTypes, FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcommand.utils import setup_log
from pbcore.io import GffIO, openDataSet

import pbreports.report.summarize_coverage.interval_tree as interval_tree
from pbreports.util import openReference


log = logging.getLogger(__name__)
__version__ = '0.3.0'


class Constants(object):
    NUM_REGIONS = 1000
    NUM_REGIONS_ID = "pbreports.task_options.num_regions"
    REGION_SIZE = 0
    REGION_SIZE_ID = "pbreports.task_options.region_size"
    FORCE_NUM_REGIONS = False
    FORCE_NUM_REGIONS_ID = "pbreports.task_options.force_num_regions"
    MAX_REGION_SIZE = 100000
    MAX_REGION_SIZE_ID = "pbreports.task_options.max_region_size"
    TOOL_ID = "pbreports.tasks.summarize_coverage"
    MAX_NUM_REGIONS = 40000  # lucky 40000
    BATCH_SIZE = 100000.0


def get_metadata_lines(readers, untruncator):
    """Write the GFF header that contains coverate metadata.
    Replaces SummarizeCmpH5._writeMetaData from pbpy. readers
    is a list of CmpH5Reader and/or pysam.Samfile objects.

    :param readers: list of alignment readers that provide
        referenceInfoTable. In practice these are either CmpH5Readers
        or SamfileAdapters.
    :param untruncator: dict that maps from truncated name to full name.
        If a truncated name does not appear in the dict, then it just
        uses the truncated name.
    """

    metadata_lines = []

    current_time_string = time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime())
    metadata_lines.append("##{k} {v}".format(k="date", v=current_time_string))
    metadata_lines.append("##{k} {v}".format(k="source",
                                             v="PACBIO_AlignmentSummary 1.0"))
    command_line = ' '.join([os.path.basename(__file__)] + sys.argv[1:])
    metadata_lines.append("##{k} {v}".format(k="source-commandline",
                                             v=command_line))

    references = []
    for reader in readers:
        for reference in reader.referenceInfoTable:
            if reference.Length == 0:
                raise Exception("Invalid zero-length contig: {c}"
                                .format(c=reference.FullName))

            ref_key = (reference.FullName, reference.Length)
            if ref_key not in references:
                full_name = untruncator.get(
                    reference.FullName, reference.FullName)
                metadata_lines.append(
                    "##{k} {i} {n}".format(k="sequence-header",
                                           i=full_name,
                                           n=full_name))
                references.append(ref_key)

    for ref_full_name, ref_length in references:
        metadata_lines.append(
            "##{k} {i} 1 {n}".format(
                k="sequence-region",
                i=untruncator.get(ref_full_name, ref_full_name),
                n=ref_length))

    return metadata_lines


def get_region_size(ref_length, num_refs, region_size, num_regions,
                    force_num_regions, max_region_size=0):
    """In the pbpy version of summarizeCoverage, there was some logic spread
    around about the size of regions. This function contains all of that. The
    tricky part is that there are two ways to set region size from the command
    line: regionSize and numRegions. summarizeCoverage also tries to adjust
    these values in certain situations, like if there are many, many
    references.

    :param ref_length: The length of the reference that will be divided into
        regions.
    :param num_refs: The total number of references.
    :param region_size: The desired region size. This can be 0, in which case
        the region size will be determined by num_regions.
    :param num_regions: The desired number of regions per reference. This is
        overridden by a non-zero region_size.
    :param force_num_regions: If False, the number of regions for all
        referenceswill not exceed MAX_NUM_REGIONS. If True, num_regions is
        accepted even if it results in a huge number of regions.
    """

    # Region size overrides everything
    if region_size > 0:
        return region_size

    # And now things are a little more complicated
    if not force_num_regions:
        regions_per_reference = max(
            min(num_regions, Constants.MAX_NUM_REGIONS / num_refs),
            1)
    else:
        regions_per_reference = num_regions
    ugly_region_size = float(ref_length) / regions_per_reference
    pretty_region_size = get_pretty_value(ugly_region_size)

    if max_region_size > 0 and max_region_size < pretty_region_size:
        log.warn("Automatic region size is {r}, above the maximum allowed value - will reduce to {m}".format(
            r=pretty_region_size, m=max_region_size))
        pretty_region_size = max_region_size
    return pretty_region_size


def get_pretty_value(ugly_value):
    """Taken directly from pbpy.

    Returns a number like 1, 2, 5, 20, 50, 100, 200, 500, 1000 which is
    closest to the input value (in log space).

    :param ugly_value: A number for which a nearby, aesthetically pleasing
        number will be found.
    """
    if ugly_value <= 0.0:
        return 0
    v_log10 = math.log(ugly_value) / math.log(10.0)
    iv = math.floor(v_log10)
    mantissa = v_log10 - iv
    targets = [1.0, 2.0, 5.0, 10.0]
    best_dist = 100.0
    best_target = 0
    for t in map(math.log10, targets):
        d = abs(mantissa - t)
        if d < best_dist:
            best_dist = d
            best_target = t
    return max(1, int(round(math.pow(10.0, iv + best_target))))


def project_into_region(intervals, region_start, region_end):
    """This is taken from pbpy, more or less. Names have been changed
    to improve readability.

    JHB's comment:
    Here, I project reads into the range defined by [rangeStart,
    rangeEnd]. Coverage can be most efficiently calculated by first
    obtaining all reads overlapping the range using the
    getOverlappingRanges function then projecting them into the same
    or smaller range.

    :param intervals: A sequence of namedtuples with start and stop
        members.
    :param region_start: Start of the region into which intervals will
        be projected
    :param region_end: End of the region
    """

    region_coverage_arr = numpy.zeros(region_end - region_start,
                                      dtype=numpy.uint32)

    for interval in intervals:
        start = interval.start
        end = interval.stop
        shifted_start = max(region_start, start) - region_start
        shifted_end = min(region_end, end - 1) - region_start
        if shifted_end >= shifted_start:
            region_coverage_arr[shifted_start:(shifted_end + 1)] = \
                region_coverage_arr[shifted_start:(shifted_end + 1)] + 1

    return region_coverage_arr


def get_gaps_from_coverage(coverage_arr):
    """Get the number of contiguous gaps and the number of gap bases
    from the coverage_arr.

    This is taken straight from pbpy but broken out here for ease of
    testing.

    :param coverage_arr: A numpy array produced by project_into_region
    :returns: tuple (n_gaps, tot_gaps) contiguous gaps and total gap bases,
        respectively
    """

    zero_pos_arr = numpy.array(coverage_arr == 0, dtype='i')
    n_gaps = ((numpy.sum(numpy.abs(numpy.diff(zero_pos_arr))) +
               zero_pos_arr[0] + zero_pos_arr[-1]) / 2)
    tot_gaps = numpy.sum(zero_pos_arr)

    return n_gaps, tot_gaps


def get_attributes_from_coverage(coverage_arr):
    """Get the three coverage atttributes for the GFF record; cov, cov2,
    and gaps.
    """

    attributes_list = []
    min_cov = numpy.amin(coverage_arr)
    max_cov = numpy.amax(coverage_arr)
    median_cov = numpy.median(coverage_arr)
    mean_cov = numpy.mean(coverage_arr)
    sd_cov = numpy.std(coverage_arr)

    attributes_list.append(
        ('cov', '%.0f,%.0f,%.0f' % (min_cov, median_cov, max_cov)))
    attributes_list.append(
        ('cov2', '%.3f,%.3f' % (mean_cov, sd_cov)))

    n_gaps, tot_gaps = get_gaps_from_coverage(coverage_arr)
    attributes_list.append(
        ('gaps', '%d,%d' % (n_gaps, tot_gaps)))

    return attributes_list


def build_interval_lists(readers):
    """Create a dictionary with RefGroupId keys and values of
    intervals of alignment starts and ends for that reference.
    """
    interval_lists = defaultdict(list)  # keyed by reference group id
    for reader in readers:
        pbi = reader.pbi
        log.debug("{x}".format(x=reader))
        for ref_id, start, end in zip(pbi.tId, pbi.tStart, pbi.tEnd):
            interval_lists[ref_id].append(interval_tree.Interval(start, end))
    log.debug("Created interval lists for {n} references.".format(
        n=len(interval_lists)))
    return interval_lists


def generate_gff_records(interval_list, readers, ref_id,
                         region_size_func, untruncator):
    """Generator for Gff records for a ref_id.

    :param interval_list: a sequence of interval_tree.Intervals of
        alignments to this reference
    :param reader: CmpH5Reader for SamfileAdapter for file
        containing the alignments
    :param ref_id: ID for this reference
    :param region_size_func: function from reference length to region
        size
    :param untruncator: dict that maps from truncated name to full name.
        If a truncated name does not appear in the dict, then it just
        uses the truncated name.

    :yields: GffIO.Gff3Records
    """
    # Get the appropriate region size for this reference
    for reader in readers:
        try:
            ref_length = reader.referenceInfo(ref_id).Length
            ref_full_name = reader.referenceInfo(ref_id).FullName
            break
        except KeyError:
            pass

    short_name = ref_full_name.split()[0]
    region_size = region_size_func(ref_length)

    if region_size == 0:
        # bug 25079 - /by0 err
        raise ValueError(
            'region_size == 0 for ref_id {r}'.format(r=str(ref_id)))

    log.debug("Chosen region size for reference {i} is {r}"
              .format(i=ref_id, r=region_size))

    log.debug("reference {i} has full name {n} and length {L}"
              .format(i=ref_id, n=ref_full_name, L=ref_length))

    itree = interval_tree.IntervalTree(interval_list)

    # To improve performance, we batch the interval lookups and projections
    # into ranges
    regions_per_batch = int(math.ceil(Constants.BATCH_SIZE / region_size))
    batch_start, batch_end = 0, 0
    batch_coverage_arr = None

    for region_start in xrange(0, ref_length, region_size):
        region_end = region_start + region_size
        # pbpy summarizeCoverage would merge the last region into the
        # penultimate region, so we do that here
        if region_end >= ref_length and region_start > 0:
            continue
        if region_end + region_size >= ref_length:
            region_end = ref_length

        # Check if we need to step to the next batch
        if region_end > batch_end:
            if region_start < batch_end:
                raise ValueError("A region overlaps a batch, which should not "
                                 "happen.")

            batch_start = region_start
            batch_end = region_size * regions_per_batch + batch_end
            if ref_length - region_size <= batch_end:
                batch_end = ref_length
            log.debug("Processing batch ({s}, {e})".format(s=batch_start,
                                                           e=batch_end))

            overlapping_intervals = []
            itree.find_overlapping(batch_start, batch_end,
                                   overlapping_intervals)
            batch_coverage_arr = project_into_region(
                overlapping_intervals, batch_start, batch_end)

        region_start_in_batch = region_start - batch_start
        region_end_in_batch = region_end - batch_start
        region_coverage_arr = batch_coverage_arr[region_start_in_batch:
                                                 region_end_in_batch]

        gff_attributes = get_attributes_from_coverage(region_coverage_arr)

        # Note the region_start + 1. GFF is 1-based and used closed intervals
        # XXX using truncated name (identifier field), see ticket 28667
        gff_record = GffIO.Gff3Record(
            short_name,  # untruncator.get(ref_full_name, ref_full_name),
            region_start + 1, region_end, "region",
            score='0.00', strand='+',
            attributes=gff_attributes)

        yield gff_record


class ReferenceTruncationError(Exception):
    """An error raised when something goes wrong with truncating or expanding
    reference names.
    """
    pass


def get_name_untruncator(repo_path, truncation_regex='\s'):
    """Return a dictionary that maps truncated reference names to full
    reference names.

    :param repo_path: Path to the reference repository that contains the
        full reference names.
    :param truncation_regex: Character at which reference names are truncated.
        For SAM/BAM files, this is whitespace.

    :returns: dict from truncated name to full name

    :raises: ReferenceTruncationError if multiple full names truncate to the
        same name.
    """

    ref_entry = openReference(repo_path)

    truncated_to_full = {}
    for contig in ref_entry.contigs:
        full_ref_name = contig.header
        truncated_ref_name = re.split(truncation_regex, full_ref_name)[0]

        if truncated_ref_name in truncated_to_full:
            msg = ("The full reference '{r}' truncates to '{t}', "
                   "but another reference also truncates to '{t}'."
                   .format(r=full_ref_name, t=truncated_ref_name))
            raise ReferenceTruncationError(msg)

        truncated_to_full[truncated_ref_name] = full_ref_name

    return truncated_to_full


def summarize_coverage(aln_set, aln_summ_gff, ref_set=None,
                       num_regions=Constants.NUM_REGIONS,
                       region_size=Constants.REGION_SIZE,
                       force_num_regions=Constants.FORCE_NUM_REGIONS,
                       max_region_size=Constants.MAX_REGION_SIZE):
    """
    Main point of entry
    """

    if ref_set:
        untruncator = get_name_untruncator(ref_set)
    else:
        # this dict is always used with get(x, x), so when it's empty it will
        # just preserve the original name
        untruncator = {}

    #readers = enumerate_readers(args.alignment_file)
    readers = openDataSet(aln_set).resourceReaders()
    gff_writer = GffIO.GffWriter(aln_summ_gff)

    # First write the metadata. Names of references, command line used, things
    # like that
    metadata_lines = get_metadata_lines(readers, untruncator)
    for metadata_line in metadata_lines:
        gff_writer.writeHeader(metadata_line)
    log.debug("Wrote {n} header lines to {f}"
              .format(n=len(metadata_lines), f=aln_summ_gff))

    # Build lists of intervals for each reference
    interval_lists = build_interval_lists(readers)
    log.debug("Finished creating interval lists for {n} references"
              .format(n=len(interval_lists)))

    # Create a function that gets region size from the reference length by
    # freezing the constant parameters
    get_region_size_frozen = functools.partial(
        get_region_size, num_refs=len(interval_lists),
        region_size=region_size, num_regions=num_regions,
        force_num_regions=force_num_regions,
        max_region_size=max_region_size)

    # Create Gff records and write them
    for ref_group_id in sorted(interval_lists):
        log.debug("Generating coverage GFF records for refGroupID {r}"
                  .format(r=ref_group_id))

        gff_generator = generate_gff_records(
            interval_lists[ref_group_id], readers,
            ref_group_id, get_region_size_frozen,
            untruncator)

        try:

            for gff_record in gff_generator:
                gff_writer.writeRecord(gff_record)

        except ValueError as e:
            log.warn(e)


def args_runner(args):
    summarize_coverage(args.aln_set, args.aln_summ_gff, args.ref_set,
                       args.num_regions, args.region_size,
                       args.force_num_regions)
    return 0


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    summarize_coverage(
        aln_set=rtc.task.input_files[0],
        aln_summ_gff=rtc.task.output_files[0],
        ref_set=rtc.task.input_files[1],
        num_regions=rtc.task.options[Constants.NUM_REGIONS_ID],
        region_size=rtc.task.options[Constants.REGION_SIZE_ID],
        force_num_regions=rtc.task.options[Constants.FORCE_NUM_REGIONS_ID],
        max_region_size=rtc.task.options[Constants.MAX_REGION_SIZE_ID])
    return 0


def add_options_to_parser(p, ds_type=FileTypes.DS_ALIGN):
    p.add_input_file_type(
        ds_type,
        file_id="aln_set",
        name="AlignmentSet",
        description="AlignmentSet")
    p.add_input_file_type(
        FileTypes.DS_REF,
        file_id="ref_set",
        name="Reference dataset",
        description="ReferenceSet or FASTA")
    p.add_output_file_type(
        FileTypes.GFF,
        file_id="aln_summ_gff",
        name="Coverage Summary",
        description="Coverage summary for regions (bins) spanning the reference",
        default_name="alignment_summary")
    p.add_int(
        option_id=Constants.NUM_REGIONS_ID,
        option_str="num_regions",
        default=Constants.NUM_REGIONS,
        name="Number of regions",
        description=("Desired number of genome regions in the summary " +
                     "statistics (used for guidance, not strict). Defaults " +
                     "to 1000"))
    p.add_int(
        option_id=Constants.REGION_SIZE_ID,
        option_str="region_size",
        default=Constants.REGION_SIZE,
        name="Region size",
        description="If supplied, used a fixed genomic region size")
    p.add_int(
        option_id=Constants.MAX_REGION_SIZE_ID,
        option_str="max_region_size",
        default=Constants.MAX_REGION_SIZE,
        name="Maximum region size",
        description="Upper limit for genomic region size (ignored if region_size is set explicitly)")
    p.add_boolean(
        option_id=Constants.FORCE_NUM_REGIONS_ID,
        option_str="force_num_regions",
        default=Constants.FORCE_NUM_REGIONS,
        name="Force the number of regions",
        description=(
            "If supplied, then try to use this number (max value = 40000) of "
            "regions per reference, otherwise the coverage summary report "
            "will optimize the number of regions in the case of many "
            "references.  Not compatible with a fixed region size."))


def _get_parser_core():
    driver_exe = ("python -m "
                  "pbreports.report.summarize_coverage.summarize_coverage "
                  "--resolved-tool-contract ")
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        "Summarize Coverage",
        __doc__,
        driver_exe)
    return p


def get_parser():
    p = _get_parser_core()
    add_options_to_parser(p)
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)

# for 'python -m pbreports.report.sat ...'
if __name__ == "__main__":
    sys.exit(main())
