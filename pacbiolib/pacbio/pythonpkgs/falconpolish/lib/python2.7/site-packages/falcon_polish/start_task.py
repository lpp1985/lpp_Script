from falcon_kit.bash import write_script
import contextlib
import cStringIO
import json
import logging
import os
import re
import sys
import tempfile
LOG = logging.getLogger(__name__)

def fn(p): return p

def base_from_done(done_fn):
    """
    >>> base_from_done('/foo/bar/x_done')
    '/foo/bar/x'
    """
    return done_fn[:-5]

@contextlib.contextmanager
def ContentUpdater(fn):
    """Write new content only if differs from old.
    """
    if os.path.exists(fn):
        with open(fn) as f:
            old_content = f.read()
    else:
        old_content = None
    new_writer = cStringIO.StringIO()
    yield new_writer
    new_content = new_writer.getvalue()
    if new_content != old_content:
        with open(fn, 'w') as f:
            f.write(new_content)

def dict2ini(ofs, options_dict):
    ofs.write('[General]\n')
    for key, val in sorted(options_dict.items()):
        ofs.write('{} = {}\n'.format(key, val))
def dict2json(ofs, options_dict):
    content = json.dumps(options_dict, sort_keys=True, indent=4, separators=(',', ': '))
    ofs.write(content + '\n')
def updated_cfg(options_dict):
    opts = dict()
    for key, val in options_dict.iteritems():
        # Drop comments (keys w/ leading ~).
        if key.startswith('~'):
            continue
        # Strip leading and trailing ws, just in case.
        opts[key] = val
    #opts['job_type'] = 'local' # Not needed. We now run falcon PypeTasks within the HGAP refresh loop.
    def add(key, val):
        if not key in opts:
            opts[key] = val
    #add('input_fofn', 'NA') # actually, we do not need this anymore
    add('target', 'assembly')
    add('sge_option_da', 'NA')
    add('sge_option_la', 'NA')
    add('sge_option_pda', 'NA')
    add('sge_option_pla', 'NA')
    add('sge_option_fc', 'NA')
    add('sge_option_cns', 'NA')
    return opts

@contextlib.contextmanager
def run_prepare_falcon(falcon_parameters, i_input_fofn_fn, fc_cfg_fn, fc_json_config_fn):
    wdir = os.path.dirname(fc_json_config_fn)
    #mkdirs(wdir)
    config_falcon = updated_cfg(dict(falcon_parameters))
    config_falcon['input_fofn'] = i_input_fofn_fn
    with ContentUpdater(fc_cfg_fn) as f:
        dict2ini(f, config_falcon)
    with ContentUpdater(fc_json_config_fn) as f:
        dict2json(f, config_falcon)

def task_filterbam(self):
    """
        {
            "other_filters": "rq >= 0.7",
            "read_length": 0
            OR
            "filters": "rq>=.7, length gte 1000, length &lt;= 50000"
        }
    """
    i_dataset_fn = fn(self.dataset)
    o_dataset_fn = fn(self.filtered)
    config = self.parameters.get('pbcoretools.tasks.filterdataset', None) # TODO: Drop this.
    if not config:
        config = self.parameters['pbcoretools']
    other_filters = config.get('other_filters', 'rq >= 0.7')
    read_length = config.get('read_length', 0)
    filters = config.get('filters', None)
    if not filters:
        filters = other_filters + ', length gte {:d}'.format(int(read_length))
    bash = """
set -vex
python -m falcon_polish.mains.run_filterbam {i_dataset_fn} {o_dataset_fn} '{filters}'
""".format(**locals())
    bash_fn = 'run_filterbam.sh'
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_bam2fasta_dexta(self):
    i_dataset_fn = fn(self.dataset)
    o_fasta_done_fn = fn(self.fasta_done)
    o_fasta_fn = base_from_done(o_fasta_done_fn) + '.fasta'
    prefix_basename = os.path.basename(o_fasta_fn)[:-6] #sans .fasta
    prefix = prefix_basename
    actual = '{}.fasta'.format(prefix) # by convention in bam2fasta
    #python -m falcon_polish.mains.run_bam2fasta {i_dataset_fn} {o_fasta_fn}
    bash = """
time bam2fasta -u -o {prefix} {i_dataset_fn}
time dexta -v {actual}
#time mv -f {actual}.dexta {o_fasta_fn}.dexta
touch {o_fasta_done_fn}
""".format(**locals())
    bash_fn = 'run_bam2fasta.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_bam_scatter(self):
    i_dataset_fn = fn(self.dataset)
    o_split_subreadsets_fofn_fn = fn(self.split_subreadsets_fofn)
    bash = """
set -vex
python -m falcon_polish.mains.run_bam_scatter {i_dataset_fn} {o_split_subreadsets_fofn_fn}
""".format(**locals())
    bash_fn = 'run_bam_scatter.sh'
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_prepare_falcon(self):
    """Pre-process FALCON cfg.
    This is super-fast, so it can always run locally.
    """
    i_input_fofn_fn = fn(self.input_fofn)
    fc_cfg_fn = fn(self.fc_cfg)
    fc_json_config_fn = fn(self.fc_json_config)
    falcon_parameters = self.parameters['falcon']
    run_prepare_falcon(falcon_parameters, i_input_fofn_fn, fc_cfg_fn, fc_json_config_fn)
    bash_fn = 'run_prepare_falcon.sh'
    open(bash_fn, 'w').write('') # empty
    self.generated_script_fn = bash_fn
def task_falcon_link(self):
    falcon_asm_done_fn = fn(self.falcon_asm_done)
    done_dir, done_fn = os.path.split(falcon_asm_done_fn)
    falcon_dir = os.path.dirname(os.path.abspath(done_dir))
    falcon_asm_fn = os.path.join(falcon_dir, '2-asm-falcon/p_ctg.fa') # or in abspath(done_dir)
    input_preads_fn = os.path.join(falcon_dir, '0-rawreads/preads/input_preads.fofn')
    length_cutoff_fn = os.path.join(falcon_dir, '0-rawreads/length_cutoff')
    falcon_link_done_fn = fn(self.falcon_link_done)
    o_fasta_fn = 'asm.fasta'
    o_preads_fofn_fn = 'preads.fofn' # for the preassembly report
    o_length_cutoff_fn = 'length_cutoff'
    bash = """
# These from and to symlinks are all by convention.
rm -f {o_fasta_fn}
ln -sf {falcon_asm_fn} {o_fasta_fn}
rm -f {o_preads_fofn_fn}
ln -sf {input_preads_fn} {o_preads_fofn_fn}
rm -f {o_length_cutoff_fn}
ln -sf {length_cutoff_fn} {o_length_cutoff_fn}
touch {falcon_link_done_fn}
ls -ltr
""".format(**locals())
    bash_fn = 'run_falcon_link.sh'
    open(bash_fn, 'w').write(bash)
    self.generated_script_fn = bash_fn
def task_fasta2referenceset(self):
    i_falcon_link_done_fn = fn(self.falcon_link_done)
    idir, done_fn = os.path.split(i_falcon_link_done_fn)
    i_fasta_fn = os.path.join(idir, 'asm.fasta') # by convention
    o_referenceset_fn = fn(self.referenceset)
    bash = """
rm -f {o_referenceset_fn} {i_fasta_fn}.fai
python -m falcon_polish.mains.run_fasta2referenceset {i_fasta_fn} {o_referenceset_fn}
""".format(**locals())
    bash_fn = 'run_fasta2referenceset.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_pbalign_scatter(self):
    """This might have problems if run in /tmp.
    """
    reads_fn = fn(self.dataset)
    referenceset_fn = fn(self.referenceset)
    out_json_fn = fn(self.out_json)
    config = self.parameters.get('pbcoretools.task_options', {})
    max_nchunks = int(config.get('scatter_subread_max_nchunks', '6'))
    bash = r"""
#rm -f {out_json_fn}
python -m pbcoretools.tasks.scatter_subread_reference -v --max_nchunks={max_nchunks} \
        {reads_fn} \
        {referenceset_fn} \
        {out_json_fn}
""".format(**locals())
    bash_fn = 'run_pbalign_scatter.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_pbalign_gather(self):
    o_unmapped_fn = fn(self.o_unmapped)
    o_ds_fn = fn(self.o_ds)
    dos = self.inputs
    unmapped_fns = [fn(v) for k,v in dos.items() if k.startswith('unmapped')]
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('alignmentset')]
    unmapped_fofn_fn = 'unmapped.fofn'
    open(unmapped_fofn_fn, 'w').write('\n'.join(unmapped_fns) + '\n')
    ds_fofn_fn = 'gathered.alignmentsets.fofn'
    open(ds_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')
    bash = r"""
#rm -f {o_ds_fn}
python -m falcon_polish.mains.run_pbalign_gather \
        {unmapped_fofn_fn} \
        {o_unmapped_fn} \
        {ds_fofn_fn} \
        {o_ds_fn}
#pbvalidate {o_ds_fn}
""".format(**locals())
    bash_fn = 'run_pbalign_gather.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_pbalign(self):
    """pbalign will eventually call blasr, like this:
 BlasrService: Align reads to references using blasr.
 BlasrService: Call "blasr /pbi/dept/secondary/siv/testdata/SA3-DS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.all.subreadset.xml /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/falcon_ns.tasks.task_falcon2_run_asm-0/file.fasta -out /scratch/tmpTbV4Ec/wLCUdL.bam  -bam  -bestn 10 -minMatch 12  -nproc 16  -minSubreadLength 50 -minAlnLength 50  -minPctSimilarity 70 -minPctAccuracy 70 -hitPolicy randombest  -randomSeed 1  -minPctSimilarity 70.0 "
 FilterService: Filter alignments using samFilter.
 FilterService: Call "rm -f /scratch/tmpTbV4Ec/aM1Mor.bam && ln -s /scratch/tmpTbV4Ec/wLCUdL.bam /scratch/tmpTbV4Ec/aM1Mor.bam"
 BamPostService: Sort and build index for a bam file.
 BamPostService: Call "samtools sort -m 4G /scratch/tmpTbV4Ec/aM1Mor.bam /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset"
 BamPostService: Call "samtools index /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam.bai"
 BamPostService: Call "pbindex /home/UNIXHOME/cdunn/repo/pb/smrtanalysis-client/smrtanalysis/siv/testkit-jobs/sa3_pipelines/polished-falcon-lambda-007-tiny/job_output/tasks/pbalign.tasks.pbalign-0/aligned.subreads.alignmentset.bam"
] OutputService: Generating the output XML file
    """
    reads_fn = fn(self.dataset)
    referenceset_fn = fn(self.referenceset)
    o_alignmentset_fn = fn(self.alignmentset)
    o_unmapped_fn = fn(self.unmapped)
    tmpdir = tempfile.tempdir
    task_opts = self.parameters['pbalign']
    options = task_opts.get('options', '')
    algorithmOptions = task_opts.get('algorithmOptions', '')

    # Write a file of absolute paths. But that file is relative to CWD:
    abs_alignmentset_fn = 'abs.' + o_alignmentset_fn

    #'--debug', # requires 'ipdb'
    #'--profile', # kinda interesting, but maybe slow?
    #'--algorithmOptions "-minMatch 12 -bestn 10 -minPctSimilarity 70.0"',
    #'--concordant',
    #'--hitPolicy randombest',
    #'--minAccuracy 70.0',
    #'--minLength 50',
    bash = """
o_fn={o_alignmentset_fn}
#rm -f ${{o_fn%.*}}
pbalign --verbose --nproc 16 {options} --tmpDir {tmpdir} --algorithmOptions "{algorithmOptions}" --unaligned {o_unmapped_fn} {reads_fn} {referenceset_fn} {o_alignmentset_fn}
#pbvalidate {o_alignmentset_fn}
dataset relativize {o_alignmentset_fn}
#pbvalidate {o_alignmentset_fn}
""".format(**locals())
    bash_fn = 'run_pbalign.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_gc_scatter(self):
    alignmentset_fn = fn(self.alignmentset)
    referenceset_fn = fn(self.referenceset)
    chunks_fofn_fn = fn(self.out_fofn)
    config = self.parameters.get('pbcoretools.task_options', {})
    max_nchunks = int(config.get('scatter_alignments_reference_max_nchunks', '13'))
    bash = """
python -m falcon_polish.mains.run_gc_scatter \
        --max-nchunks={max_nchunks} \
        {alignmentset_fn} \
        {referenceset_fn} \
        {chunks_fofn_fn}
""".format(**locals())
    bash_fn = 'run_gc_scatter.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_gc_gather(self):
    dos = self.inputs
    ds_out_fn = fn(self.ds_out)
    fastq_out_fn = fn(self.fastq_out)

    fasta_ds_fofn_fn = 'fasta.contigset.fofn'
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('contigset_')]
    open(fasta_ds_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')

    fastq_fofn_fn = 'fastq.fofn'
    dset_fns = [fn(v) for k,v in dos.items() if k.startswith('fastq_')]
    open(fastq_fofn_fn, 'w').write('\n'.join(dset_fns) + '\n')

    bash = r"""
python -m falcon_polish.mains.run_gc_gather \
        {fasta_ds_fofn_fn} \
        {fastq_fofn_fn} \
        {ds_out_fn} \
        {fastq_out_fn}
""".format(**locals())
    bash_fn = 'run_gc_gather.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_genomic_consensus(self):
    alignmentset_fn = fn(self.alignmentset)
    referenceset_fn = fn(self.referenceset)
    polished_fastq_fn = fn(self.polished_fastq)
    variants_gff_fn = fn(self.variants_gff)
    consensus_contigset_fn = fn(self.consensus_contigset)
    task_opts = self.parameters['variantCaller']
    options = task_opts.get('options', '')
    if '--alignmentSetRefWindows' not in options:
        options += ' --alignmentSetRefWindows'
    fasta_fn = re.sub(".contigset.xml", ".fasta", consensus_contigset_fn)
    # Possibly we should escape '{options}'
    bash = """
python -m falcon_polish.mains.run_variantCaller --log-level DEBUG --options '{options}' \
        {alignmentset_fn} \
        {referenceset_fn} \
        {polished_fastq_fn} \
        {variants_gff_fn} \
        {fasta_fn} \
        {consensus_contigset_fn}
""".format(**locals())
    bash_fn = 'run_gc.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_polished_assembly_report(self):
    task_opts = self.parameters['pbreports.tasks.summarize_coverage']
    options = task_opts.get('options', '')
    referenceset_fn = fn(self.referenceset)
    gathered_alignmentset_fn = fn(self.gathered_alignmentset)
    polished_fastq_fn = fn(self.polished_fastq)
    report_fn = fn(self.report_json)
    alignment_summary_gff_fn = 'alignment.summary.gff'
    """
    If necessary, we could call this:
    from pbreports.report.summarize_coverage.summarize_coverage import summarize_coverage
    summarize_coverage(args.aln_set, args.aln_summ_gff, args.ref_set,
                       args.num_regions, args.region_size,
                       args.force_num_regions)
    """
    # https://github.com/PacificBiosciences/pbreports/pull/186
    bash = r"""
python -m pbreports.report.summarize_coverage.summarize_coverage \
        {options} \
        {gathered_alignmentset_fn} \
        {referenceset_fn} \
        {alignment_summary_gff_fn}
python -m pbreports.report.polished_assembly \
        {alignment_summary_gff_fn} \
        {polished_fastq_fn} \
        {report_fn}
""".format(**locals())
    bash_fn = 'run_report.sh'
    #open(bash_fn, 'w').write(bash)
    job_done = bash_fn + '.done'
    hgap_config = self.parameters.get('hgap') # for now, all tasks are tmpdir, or not
    write_script(bash, bash_fn, job_done)
    self.generated_script_fn = bash_fn
def task_foo(self):
    LOG.debug('WARNME1 {!r}'.format(__name__))
    #print repr(self.parameters), repr(self.URL), repr(self.foo1)
    sys.system('touch {}'.format(fn(self.foo2)))
    script_fn = 'noop.sh'
    open(script_fn, 'w').write('echo NOOP raw')
    self.generated_script_fn = script_fn

def task_fastas2fofn(self):
    # Record the fasta filenames in a FOFN, based on a filename convention.
    # We depend on 'done' files, not directly on fastas, so we can
    # delete fastas after we use them, downstream.
    # Note: This is light and quick, so we can do it in-process.
    fofn_fn = fn(self.fofn)
    dos = self.inputs
    fasta_fns = [(base_from_done(fn(v)) + '.dexta') for k,v in dos.items()]
    content = '\n'.join(sorted(fasta_fns)) + '\n'
    # TODO: Do we need ContentUpdater here?
    with ContentUpdater(fofn_fn) as f:
        f.write(content)
    bash_fn = 'run_fastas2fofn.sh'
    open(bash_fn, 'w').write('') # empty
    self.generated_script_fn = bash_fn
