import logging

from pbsmrtpipe.core import register_pipeline
from pbsmrtpipe.constants import to_pipeline_ns

from .pb_pipelines_sa3 import Constants, Tags, _core_align, _core_gc

log = logging.getLogger(__name__)


def dev_register(relative_id, display_name, tags=(), task_options=None):
    pipeline_id = to_pipeline_ns(relative_id)
    ptags = list(set(tags + (Tags.DENOVO, )))
    return register_pipeline(pipeline_id, display_name, "0.1.0", tags=ptags, task_options=task_options)

def _get_falcon_pipeline(i_cfg, i_fasta_fofn):
    """Basic falcon pipeline components.
    """
    b0 = [
          (i_cfg,        'falcon_ns.tasks.task_falcon_config:0'),
          (i_fasta_fofn, 'falcon_ns.tasks.task_falcon_config:1'),
          ('falcon_ns.tasks.task_falcon_config:0', 'falcon_ns.tasks.task_falcon_make_fofn_abs:0'),
    ]
    br0 = [
          ('falcon_ns.tasks.task_falcon_config:0',        'falcon_ns.tasks.task_falcon0_build_rdb:0'),
          ('falcon_ns.tasks.task_falcon_make_fofn_abs:0', 'falcon_ns.tasks.task_falcon0_build_rdb:1'),
         ]
    br1 = [
          ('falcon_ns.tasks.task_falcon_config:0',     'falcon_ns.tasks.task_falcon0_run_daligner_jobs:0'),
          ('falcon_ns.tasks.task_falcon0_build_rdb:0', 'falcon_ns.tasks.task_falcon0_run_daligner_jobs:1'),
         ]
    rm0 = [('falcon_ns.tasks.task_falcon0_run_daligner_jobs:0', 'falcon_ns.tasks.task_falcon0_rm_las:0')] # rm raw_reads.*.raw_reads.*.las

    # br2: make scripts, LAmerge (e.g., m_00001/merge_00001.sh), LA4Falcon (e.g., preads/c_00001.sh) and db2falcon scripts.
    br2 = [
          ('falcon_ns.tasks.task_falcon_config:0',             'falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:0'),
          ('falcon_ns.tasks.task_falcon0_build_rdb:0',         'falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:1'),
          ('falcon_ns.tasks.task_falcon0_run_daligner_jobs:0', 'falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:2'),
         ]
    # br3: execute LAmerge scripts (e.g., m_00001/merge_00001.sh) to create raw_reads.*.las
    br3 = [('falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:1', 'falcon_ns.tasks.task_falcon0_merge:0')]     # merge.json
    # br4: execute LA4Falcon scripts (e.g., preads/c_00001.sh) to create out.0000*.fasta
    br4 = [('falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:2', 'falcon_ns.tasks.task_falcon0_cons:0'),      # cons.json
           ('falcon_ns.tasks.task_falcon0_merge:0',                    'falcon_ns.tasks.task_falcon0_cons:1')       # merge_done.txt, sentinel
          ]
    rm1 = [('falcon_ns.tasks.task_falcon0_cons:0',   'falcon_ns.tasks.task_falcon1_rm_las:0'),
           ('falcon_ns.tasks.task_falcon0_rm_las:0', 'falcon_ns.tasks.task_falcon1_rm_las:1')] # rm raw_reads.*.las

    bp0 = [
          ('falcon_ns.tasks.task_falcon_config:0',                    'falcon_ns.tasks.task_falcon1_build_pdb:0'),  # config.json
          ('falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:0', 'falcon_ns.tasks.task_falcon1_build_pdb:1'),  # fofn of out.*.fasta
          ('falcon_ns.tasks.task_falcon0_cons:0',                     'falcon_ns.tasks.task_falcon1_build_pdb:2')   # cons_done.txt, sentinel
         ]

    bp1 = [
          ('falcon_ns.tasks.task_falcon_config:0',     'falcon_ns.tasks.task_falcon1_run_daligner_jobs:0'),
          ('falcon_ns.tasks.task_falcon1_build_pdb:0', 'falcon_ns.tasks.task_falcon1_run_daligner_jobs:1'),
         ]
    bp2 = [
          ('falcon_ns.tasks.task_falcon_config:0',             'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:0'),
          ('falcon_ns.tasks.task_falcon1_build_pdb:0',         'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:1'),
          ('falcon_ns.tasks.task_falcon1_run_daligner_jobs:0', 'falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:2'),
         ]
    # bp3: execute LAmerge scripts (e.g., m_00001/merge_00001.sh) to create preads.*.las
    bp3 = [('falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:1', 'falcon_ns.tasks.task_falcon1_merge:0')]     # merge.json
    # bp4: execute db2falcon scripts (e.g., run_db2falcon.sh) to create falcon db.
    bp4 = [('falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:2', 'falcon_ns.tasks.task_falcon1_db2falcon:0'), # db2falcon.json
           ('falcon_ns.tasks.task_falcon1_merge:0',                    'falcon_ns.tasks.task_falcon1_db2falcon:1')  # merge_done.txt, sentinel
          ]
    bf = [
            ('falcon_ns.tasks.task_falcon_config:0',                    'falcon_ns.tasks.task_falcon2_run_asm:0'),  # config.json
            ('falcon_ns.tasks.task_falcon1_run_merge_consensus_jobs:0', 'falcon_ns.tasks.task_falcon2_run_asm:1'),  # fofn of preads.*.las
            ('falcon_ns.tasks.task_falcon1_db2falcon:0',                'falcon_ns.tasks.task_falcon2_run_asm:2')   # db2falcon_done.txt, sentinel
         ]
    rm2 = [('falcon_ns.tasks.task_falcon2_run_asm:0', 'falcon_ns.tasks.task_falcon2_rm_las:0'),
           ('falcon_ns.tasks.task_falcon1_rm_las:0',  'falcon_ns.tasks.task_falcon2_rm_las:1') ] # clean up all *.las regardless

    report_pay = [
          ('falcon_ns.tasks.task_falcon_config:0',
                        'falcon_ns.tasks.task_report_preassembly_yield:0'),
          ('falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs:0',
                        'falcon_ns.tasks.task_report_preassembly_yield:1'),
          ('falcon_ns.tasks.task_falcon0_build_rdb:1',
                        'falcon_ns.tasks.task_report_preassembly_yield:2'),
          ('falcon_ns.tasks.task_falcon0_cons:0',
                        'falcon_ns.tasks.task_report_preassembly_yield:3')  # cons_done.txt, sentinel
    ]
    results = dict()
    results['asm'] = 'falcon_ns.tasks.task_falcon2_run_asm:0'
    return b0 + br0 + br1 + rm0 + br2 + br3 + br4 + rm1 + bp0 + bp1 + bp2 + bp3 + bp4 + bf + rm2 + report_pay, results

def _get_polished_falcon_pipeline():
    subreadset = Constants.ENTRY_DS_SUBREAD

    filt = [(subreadset, 'pbcoretools.tasks.filterdataset:0')]
    btf = [('pbcoretools.tasks.filterdataset:0', 'pbcoretools.tasks.bam2fasta:0')]
    ftfofn = [('pbcoretools.tasks.bam2fasta:0', 'pbcoretools.tasks.fasta2fofn:0')]

    i_fasta_fofn = 'pbcoretools.tasks.fasta2fofn:0'

    gen_cfg = [(i_fasta_fofn, 'falcon_ns.tasks.task_falcon_gen_config:0')]

    i_cfg = 'falcon_ns.tasks.task_falcon_gen_config:0'

    falcon, falcon_results = _get_falcon_pipeline(i_cfg, i_fasta_fofn)

    ref = falcon_results['asm']

    faidx = [(ref, 'pbcoretools.tasks.fasta2referenceset:0')]

    aln = 'pbalign.tasks.pbalign:0'
    ref = 'pbcoretools.tasks.fasta2referenceset:0'

    polish = _core_align(subreadset, ref) + _core_gc(aln,
                                                     ref)
    results = dict()
    results['aln'] = aln
    results['ref'] = ref

    return ((filt + btf + ftfofn + gen_cfg + falcon + faidx + polish), results)

@dev_register("pipe_falcon_with_fofn", "Falcon FOFN Pipeline",
              tags=("local", "chunking", "internal"))
def get_task_falcon_local_pipeline2():
    """Simple falcon local pipeline.
    Use an entry-point for FASTA input.
    """
    return _get_falcon_pipeline('$entry:e_01', '$entry:e_02')[0]

@dev_register("pipe_falcon", "Falcon Pipeline",
              tags=("local", "chunking", "internal"))
def get_task_falcon_local_pipeline1():
    """Simple falcon local pipeline.
    FASTA input comes from config file.
    """
    i_cfg = '$entry:e_01'
    init = [
          (i_cfg, 'falcon_ns.tasks.task_falcon_config_get_fasta:0'),
           ]
    i_fasta_fofn = 'falcon_ns.tasks.task_falcon_config_get_fasta:0' # output from init
    return init + _get_falcon_pipeline(i_cfg, i_fasta_fofn)[0]

@dev_register("polished_falcon", "Polished Falcon Pipeline",
              tags=("chunking", "internal"))
def get_task_polished_falcon_pipeline():
    """Simple polished falcon local pipeline.
    FASTA input comes from the SubreadSet.
    """
    i_cfg = '$entry:e_01'
    subreadset = Constants.ENTRY_DS_SUBREAD

    btf = [(subreadset, 'pbcoretools.tasks.bam2fasta:0')]
    ftfofn = [('pbcoretools.tasks.bam2fasta:0', 'pbcoretools.tasks.fasta2fofn:0')]

    i_fasta_fofn = 'pbcoretools.tasks.fasta2fofn:0'

    falcon, falcon_results = _get_falcon_pipeline(i_cfg, i_fasta_fofn)

    ref = falcon_results['asm']

    faidx = [(ref, 'pbcoretools.tasks.fasta2referenceset:0')]

    ref = 'pbcoretools.tasks.fasta2referenceset:0'

    polish = _core_align(subreadset, ref) + _core_gc('pbalign.tasks.pbalign:0',
                                                     ref)

    return btf + ftfofn + falcon + faidx + polish

# Copied from pb_pipelines_sa3.py
RESEQUENCING_TASK_OPTIONS = {
    "genomic_consensus.task_options.diploid": False,
    "genomic_consensus.task_options.algorithm": "best",
    "pbcoretools.task_options.other_filters": "rq >= 0.7",
    #"pbalign.task_options.algorithm_options": "-minMatch 12 -bestn 10 -minPctSimilarity 70.0 -refineConcordantAlignments",
    #"pbalign.task_options.concordant": True,
}

@dev_register("polished_falcon_lean", "Assembly (HGAP 4) without reports", tags=("internal",),
        task_options=RESEQUENCING_TASK_OPTIONS)
def get_falcon_pipeline_lean():
    """Simple polished falcon local pipeline (sans reports).
    FASTA input comes from the SubreadSet.
    Cfg input is built from preset.xml
    """
    falcon, _ = _get_polished_falcon_pipeline()
    return falcon

@dev_register("polished_falcon_fat", "Assembly (HGAP 4)",
        task_options=RESEQUENCING_TASK_OPTIONS)
def get_falcon_pipeline_fat():
    """Same as polished_falcon_lean, but with reports.
    """
    falcon, results = _get_polished_falcon_pipeline()

    # id's of results from falcon:
    aln = 'pbalign.tasks.pbalign:0'
    ref = 'pbcoretools.tasks.fasta2referenceset:0'

    # summarize the coverage:
    sum_cov = [(aln, "pbreports.tasks.summarize_coverage:0"),
               (ref, "pbreports.tasks.summarize_coverage:1")]

    # gen polished_assembly report:
    # takes alignment summary GFF, polished assembly fastQ
    polished_report = [('pbreports.tasks.summarize_coverage:0', 'pbreports.tasks.polished_assembly:0'),
                       ('genomic_consensus.tasks.variantcaller:2', 'pbreports.tasks.polished_assembly:1')]
    mapping_report = [
        (aln, "pbreports.tasks.mapping_stats_hgap:0"),
        (Constants.ENTRY_DS_SUBREAD, "pbreports.tasks.mapping_stats_hgap:1")
    ]
    coverage_report = [
        (ref, "pbreports.tasks.coverage_report_hgap:0"),
        ("pbreports.tasks.summarize_coverage:0", "pbreports.tasks.coverage_report_hgap:1")
    ]
    fasta_out = [
        ("genomic_consensus.tasks.variantcaller:1", "pbcoretools.tasks.contigset2fasta:0")
    ]

    return falcon + sum_cov + polished_report + mapping_report + coverage_report + fasta_out

def _get_hgap_pypeflow(i_cfg, i_logging_cfg, i_subreadset):
    return [
            (i_cfg,         'falcon_ns.tasks.task_hgap_run:0'),
            (i_logging_cfg, 'falcon_ns.tasks.task_hgap_run:1'),
            (i_subreadset,  'falcon_ns.tasks.task_hgap_run:2'),
           ]

@dev_register("hgap_cmd", "XI- Experimental Assembly (HGAP 5) without reports", tags=("internal",))
def hgap_cmd():
    # from hgap-cfg.json, logging-cfg.json, and subreads-dataset
    """Simple polished HGAP pipeline (sans reports).
    BAM input comes from the SubreadSet.
    hgap-cfg.json comes from $entry:e_01
    logging-cfg.json comes from $entry:e_02
    """
    subreadset = Constants.ENTRY_DS_SUBREAD
    hgap_cfg = '$entry:e_01'
    logging_cfg = '$entry:e_02'
    return _get_hgap_pypeflow(hgap_cfg, logging_cfg, subreadset)

@dev_register("hgap_lean", "X - Experimental Assembly (HGAP 5) without reports", tags=("internal",))
def hgap_lean():
    """GUI polished HGAP pipeline (sans reports).
    BAM input comes from the SubreadSet.
    .cfg inputs are based on pbsmrtpipe options, via task_hgap_prepare
    """
    subreadset = Constants.ENTRY_DS_SUBREAD
    hgap_prepare = [(subreadset,
                   'falcon_ns.tasks.task_hgap_prepare:0')]
    hgap_cfg =     'falcon_ns.tasks.task_hgap_prepare:0'
    logging_cfg =  'falcon_ns.tasks.task_hgap_prepare:1'
    hgap_run = _get_hgap_pypeflow(hgap_cfg, logging_cfg, subreadset)
    return hgap_prepare + hgap_run

@dev_register("hgap_fat", "Assembly (HGAP 5 beta)", tags=("internal",))
def hgap_fat():
    """GUI polished HGAP pipeline.
    BAM input comes from the SubreadSet.
    .cfg inputs are based on pbsmrtpipe options, via task_hgap_prepare
    """
    subreadset = Constants.ENTRY_DS_SUBREAD
    hgap_prepare = [(subreadset,
                   'falcon_ns.tasks.task_hgap_prepare:0')]
    hgap_cfg =     'falcon_ns.tasks.task_hgap_prepare:0'
    logging_cfg =  'falcon_ns.tasks.task_hgap_prepare:1'
    hgap_run = _get_hgap_pypeflow(hgap_cfg, logging_cfg, subreadset)
    return hgap_prepare + hgap_run
