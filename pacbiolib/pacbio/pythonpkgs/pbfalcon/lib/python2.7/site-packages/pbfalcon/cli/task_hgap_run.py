from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.models import (ResourceTypes, FileTypes, SymbolTypes)
from pbcommand.models.parser import get_pbparser
#from .. import hgap_prepare
from .. import tusks as pbfalcon
import sys
import logging

__version__ = '1.0.0'
log = logging.getLogger(__name__)
TOOL_ID = 'falcon_ns.tasks.task_hgap_run'

# https://github.com/PacificBiosciences/smrtflow/blob/master/smrt-server-analysis/src/main/resources/pipeline-template-view-rules/pipeline_template_view_rules-polished_falcon_fat.json

def add_args_and_options(p):
    # FileType, label, name=title, description
    p.add_input_file_type(FileTypes.JSON, "Label PacBio.FileTypes.json_0", "<FileType id=PacBio.FileTypes.json name=file >", "description for PacBio.FileTypes.json_0")
    p.add_input_file_type(FileTypes.JSON, "Label PacBio.FileTypes.json_1", "<FileType id=PacBio.FileTypes.json name=file >", "description for PacBio.FileTypes.json_1")
    p.add_input_file_type(FileTypes.DS_SUBREADS, "Label PacBio.DataSet.SubreadSet_2", "<DataSetFileType id=PacBio.DataSet.SubreadSet name=file >", "description for PacBio.DataSet.SubreadSet_2")

    # File Type, label, name, description, default file name
    # new outputs
    p.add_output_file_type(FileTypes.FASTA, "preads_id", "preads.fasta", "Tarball of FASTA for preassembly reads", 'preads')
    p.add_output_file_type(FileTypes.FASTA, "polished_fasta_id", "polished.fasta", "FASTA of polished assembly", 'polished')
    p.add_output_file_type(FileTypes.FASTQ, "polished_fastq_id", "polished.fastq", "FASTQ of polished assembly", 'polished')
    p.add_output_file_type(FileTypes.CSV, "polished_csv_id", "polished.csv", "CSV from report of polished assembly", 'polished')
    p.add_output_file_type(FileTypes.DS_ALIGN, "alignmentset_id", "aligned.subreads.alignmentset.xml", "Dataset of BAM files of aligned subreads", 'aligned.subreads')
    p.add_output_file_type(FileTypes.GFF, "gff_id", "alignment.summary.gff", "General Feature Format file for alignment coverage", 'alignment.summary')
    p.add_output_file_type(FileTypes.TXT, "unmapped_id", "unmapped.subreads.txt", "Names of unmapped subreads from pbalign", 'unmapped.subreads')

    # old outputs
    p.add_output_file_type(FileTypes.DS_CONTIG, "contig_id", "polished.contigset.xml", "Contigset of polished FASTA sequences (redundant with polished.fasta)", 'polished.contigset')
    p.add_output_file_type(FileTypes.REPORT, "preassembly_rpt_id", "Preassembly report", "description for <FileType id=PacBio.FileTypes.JsonReport name=report >", 'preassembly_rpt')
    p.add_output_file_type(FileTypes.REPORT, "polished_assembly_rpt_id", "Polished assembly report", "description for <FileType id=PacBio.FileTypes.JsonReport name=report >", 'polished_assembly_rpt')
    p.add_output_file_type(FileTypes.LOG, "out2_id", "Another log output, experimentally", 'description for <FileType id=PacBio.FileTypes.log name=file >', 'out2')
    return p

def get_contract_parser():
    nproc = SymbolTypes.MAX_NPROC
    resource_types = (ResourceTypes.TMP_DIR,)
    # Commandline exe to call "{exe}" /path/to/resolved-tool-contract.json
    driver_exe = "python -m pbfalcon.cli.task_hgap_run --resolved-tool-contract "
    desc = 'pbcommand wrapper for ' + TOOL_ID
    name = 'Tool task_hgap_run'
    p = get_pbparser(TOOL_ID, __version__, name, desc, driver_exe,
            is_distributed=False, nproc=nproc, resource_types=resource_types)
    add_args_and_options(p)
    return p

def run_my_main(input_files, output_files, tmpdir):
    # do stuff. Main should return an int exit code
    rc = pbfalcon.run_hgap(input_files, output_files, tmpdir)
    if rc:
        return rc
    else:
        return 0

def _args_runner(args):
    # this is the args from parser.parse_args()
    # the properties of args are defined as "labels" in the add_args_and_options func.
    # TODO: Convert 'args' to a dict somehow?
    return run_my_main([args.fasta_in], [args.fasta_out], args) # NOT USED and not correct

def _resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    # all options are referenced by globally namespaced id. This allows tools to use other tools options
    # e.g., pbalign to use blasr defined options.
    tempdir = rtc.task.tmpdir_resources[0].path
    # options = rtc.task.options
    return run_my_main(rtc.task.input_files, rtc.task.output_files, tempdir)

def main(argv=sys.argv):
    log.info("Starting {f} version {v} pbcommand example dev app".format(f=__file__, v=__version__))
    p = get_contract_parser()
    return pbparser_runner(argv[1:],
                           p,
                           _args_runner, # argparse runner func
                           _resolved_tool_contract_runner, # tool contract runner func
                           log, # log instance
                           setup_log # setup log func
                           )
if __name__ == '__main__':
    sys.exit(main())
