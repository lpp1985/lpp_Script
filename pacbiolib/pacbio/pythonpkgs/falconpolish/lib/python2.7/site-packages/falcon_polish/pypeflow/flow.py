from __future__ import absolute_import
from .. import sys, start_task
from falcon_kit import stats_preassembly

#from pypeflow.pwatcher_bridge import PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase
#from pypeflow.controller import PypeThreadWorkflow # for setNumThreadAllowed()
#from pypeflow.data import makePypeLocalFile, fn
#from pypeflow.task import PypeTask #, PypeThreadTaskBase
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
        makePypeLocalFile, fn, PypeTask)

import falcon_kit.mains.run1
import falcon_kit.run_support

import contextlib
import gzip
import json
import logging
import os
import pprint
import tempfile

log = logging.getLogger(__name__)

#PypeWorkflow = PypeThreadWorkflow
PypeWorkflow = PypeProcWatcherWorkflow

def create_tasks_fasta2DB(split_subreadsets_fofn_pfn, parameters):
    tasks = list()
    next_inputs = dict()
    topdir = os.path.join(os.path.dirname(fn(split_subreadsets_fofn_pfn)), 'run-fastas2fofn')
    # Create the fastas in parallel.
    for i, chunk_fn in enumerate(open(fn(split_subreadsets_fofn_pfn)).read().splitlines()):
        wdir = os.path.join(topdir, 'fasta_job_{:03d}'.format(i)) # TODO: 02
        chunk_pfn = makePypeLocalFile(os.path.join(wdir, chunk_fn))
        fasta_done_fn = os.path.join(wdir, 'chunk_{:03d}_done'.format(i)) # TODO: 02
        # By depending on a sentinel, we are allowed to delete fastas later.
        # Note: i might not match num in chunk_fn, but that is ok
        fasta_done_pfn = makePypeLocalFile(fasta_done_fn)
        make_task = PypeTask(
                inputs = {"dataset": chunk_pfn, },
                outputs =  {"fasta_done": fasta_done_pfn, },
                parameters = parameters,
        )
        task = make_task(start_task.task_bam2fasta_dexta)
        tasks.append(task)
        next_inputs['fasta_{}_done'.format(i)] = fasta_done_pfn
        #fasta_fn = base_from_done(fasta_done_fn) + '.fasta'  # By convention.
    # Create the FOFN of fastas.
    fasta_fofn_fn = os.path.join(topdir, 'fasta.fofn')
    fasta_fofn_pfn = makePypeLocalFile(fasta_fofn_fn)
    make_task = PypeTask(
            inputs = next_inputs,
            outputs =  {"fofn": fasta_fofn_pfn,
            },
            parameters = parameters,
    )
    task = make_task(start_task.task_fastas2fofn)
    tasks.append(task)
    return tasks, fasta_fofn_pfn
def yield_pipeline_chunk_names_from_json(ifs, key):
    d = json.loads(ifs.read())
    for cs in d['chunks']:
        #chunk_id = cs['chunk_id']
        chunk_datum = cs['chunk']
        yield chunk_datum[key]
def create_tasks_pbalign(chunk_json_pfn, referenceset_pfn, parameters):
    """Create a pbalign task for each chunk, plus a gathering task.
    """
    tasks = list()
    gathering = dict()
    chunk_dir = os.path.dirname(fn(chunk_json_pfn))
    for i, subreadset_fn in enumerate(sorted(yield_pipeline_chunk_names_from_json(open(fn(chunk_json_pfn)), '$chunk.subreadset_id'))):
        wdir = 'run-pbalign-{:02d}'.format(i)
        subreadset_fn = os.path.join(chunk_dir, os.path.basename(subreadset_fn))
        subreadset_pfn = makePypeLocalFile(subreadset_fn)
        unmapped_pfn = makePypeLocalFile('{wdir}/unmapped.txt'.format(**locals()))
        alignmentset_pfn = makePypeLocalFile('{wdir}/align.subreads.{i:02d}.alignmentset.xml'.format(**locals()))
        gathering['unmapped_{:02d}'.format(i)] = unmapped_pfn
        gathering['alignmentsets_{:02d}'.format(i)] = alignmentset_pfn
        """Also produces:
        aligned.subreads.i.alignmentset.bam
        aligned.subreads.i.alignmentset.bam.bai
        aligned.subreads.i.alignmentset.bam.pbi
        """
        make_task = PypeTask(
                inputs = {"chunk_json": chunk_json_pfn,
                          "dataset": subreadset_pfn,
                          "referenceset": referenceset_pfn,
                },
                outputs = {"alignmentset": alignmentset_pfn,
                           "unmapped": unmapped_pfn,
                },
                parameters = parameters,
        )
        task = make_task(start_task.task_pbalign)
        tasks.append(task)
    o_alignmentset_pfn = makePypeLocalFile('run-pbalign_gather/aligned.subreads.alignmentset.xml')
    o_unmapped_pfn = makePypeLocalFile('run-pbalign_gather/unmapped.txt')
    make_task = PypeTask(
            inputs = gathering,
            outputs = {"o_ds": o_alignmentset_pfn,
                       "o_unmapped": o_unmapped_pfn,
            },
            parameters = parameters,
    )
    task = make_task(start_task.task_pbalign_gather)
    tasks.append(task)
    return tasks, alignmentset_pfn
def mkdirs(d):
    log.debug('mkdir -p {}'.format(d))
    if not os.path.isdir(d):
        os.makedirs(d)
def create_tasks_gc(fofn_pfn, referenceset_pfn, parameters):
    """Create a gc task for each chunk, plus a gathering task.
    Here is the convoluted workflow:
    1. For each gc instance "chunk":
      A. variantCaller writes .fasta
      B. We create a contigset for the .fasta
    2. We keep the contigset output filenames in a FOFN (from run_gc_scatter)
       and pass that to run_gc_gather().
    3. We read each contigset and add them to a gathered ContigSet.
    4. We "consolidate" their underlying .fasta "resources",
       assuming their filenames match except extenion.
    5. Finally, we write the gathered contigset.
    Whew!
    We also gather fastq here, for convenience.
    """
    tasks = list()
    contigsets = dict()
    fastqs = dict()
    # Assume fofn of gc chunks are all relative to the dir of the fofn.
    for i, alignmentset_bn in enumerate(open(fn(fofn_pfn)).read().split()):
        alignmentset_fn = os.path.join(os.path.dirname(fn(fofn_pfn)), alignmentset_bn)
        wdir = 'run-gc-{:02}'.format(i)
        mkdirs(wdir) # Assume CWD is correct.
        alignmentset_pfn = makePypeLocalFile(alignmentset_fn) # New pfn cuz it was not pfn before.
        polished_fastq_pfn = makePypeLocalFile(os.path.join(wdir, 'consensus.fastq'))
        variants_gff_pfn = makePypeLocalFile(os.path.join(wdir, 'variants.gff'))
        consensus_contigset_pfn = makePypeLocalFile(os.path.join(wdir, 'consensus.contigset.xml'))
        """Also produces:
        consensus.fasta
        consensus.fasta.fai

        And note that these files names are important, as pbcoretools gathering expects
        a particular pattern.
        """
        contigsets['contigset_{:02d}'.format(i)] = consensus_contigset_pfn
        fastqs['fastq_{:02d}'.format(i)] = polished_fastq_pfn
        make_task = PypeTask(
                inputs = {"alignmentset": alignmentset_pfn,
                          "referenceset": referenceset_pfn,},
                outputs = {
                    "polished_fastq": polished_fastq_pfn,
                    "variants_gff": variants_gff_pfn,
                    "consensus_contigset": consensus_contigset_pfn,
                },
                parameters = parameters,
        )
        task = make_task(start_task.task_genomic_consensus)
        tasks.append(task)
    contigset_pfn = makePypeLocalFile('run-gc-gather/contigset.xml')
    gathered_fastq_pfn = makePypeLocalFile('run-gc-gather/gathered.fastq')
    inputs = dict(contigsets)
    inputs.update(fastqs)
    log.debug('inputs to gc_gather:{}'.format(pprint.pformat(contigsets)))
    make_task = PypeTask(
            inputs = inputs,
            outputs = {"ds_out": contigset_pfn,
                       "fastq_out": gathered_fastq_pfn,
            },
            parameters = parameters,
    )
    task = make_task(start_task.task_gc_gather)
    tasks.append(task)
    return tasks, contigset_pfn, gathered_fastq_pfn

def flow(config):
    #import pdb; pdb.set_trace()
    parameters = config
    #exitOnFailure=config['stop_all_jobs_on_failure'] # only matter for parallel jobs
    #wf.refreshTargets(exitOnFailure=exitOnFailure)
    #wf = PypeThreadWorkflow()
    #wf = PypeWorkflow()
    #wf = PypeWorkflow(job_type='local')
    log.debug('config=\n{}'.format(pprint.pformat(config)))
    # Set some defaults on the Workflow.
    concurrent_jobs = 24 # TODO: Configure this.
    wf = PypeWorkflow(
            job_type=config['hgap'].get('job_type'),
            job_queue=config['hgap'].get('job_queue'),
            watcher_type=config['hgap'].get('pwatcher_type', 'blocking'),
            #watcher_directory=config['pwatcher_directory'],
            max_jobs=config['hgap'].get('max_jobs', concurrent_jobs),
    )

    use_tmpdir = config['hgap'].get('use_tmpdir')
    if use_tmpdir:
        log.info('hgap.use_tmpdir={!r}'.format(use_tmpdir))
        if use_tmpdir is not True and '/' in use_tmpdir:
            tempfile.tempdir = use_tmpdir
            log.info('Using tempfile.tempdir={}'.format(tempfile.tempdir))
        else:
            log.info('Keeping tempfile.tempdir={}'.format(tempfile.tempdir))

    dataset_pfn = makePypeLocalFile(config['pbsmrtpipe']['input_files'][0])
    filtered_pfn = makePypeLocalFile('run-filterbam/filtered.subreadset.xml')
    make_task = PypeTask(
            inputs = {"dataset": dataset_pfn, },
            outputs = {"filtered": filtered_pfn, },
            parameters = parameters,
    )
    task = make_task(start_task.task_filterbam)
    wf.addTask(task)

    split_subreadsets_fofn_pfn = makePypeLocalFile('run-bam_scatter/chunked_subreadsets.fofn')
    make_task = PypeTask(
            inputs = {"dataset": filtered_pfn, },
            outputs =  {"split_subreadsets_fofn": split_subreadsets_fofn_pfn, },
            parameters = parameters,
    )
    task = make_task(start_task.task_bam_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    tasks, input_fofn_pfn = create_tasks_fasta2DB(split_subreadsets_fofn_pfn, parameters)
    wf.addTasks(tasks)
    wf.refreshTargets()

    fc_cfg_pfn = makePypeLocalFile('run-falcon/fc.cfg')
    fc_json_config_pfn = makePypeLocalFile('run-falcon/fc.json')
    make_task = PypeTask(
            inputs = {
                      "input_fofn": input_fofn_pfn,
            },
            outputs = {"fc_cfg": fc_cfg_pfn,
                       "fc_json_config": fc_json_config_pfn,
            },
            parameters = parameters,
    )
    task = make_task(start_task.task_prepare_falcon)
    wf.addTask(task)
    wf.refreshTargets()

    input_config_fn = fn(fc_cfg_pfn)
    with sys.cd('run-falcon'):
        falcon_kit.mains.run1.fc_run_logger = falcon_kit.run_support.logger = logging.getLogger('falcon')
        fc_cfg = falcon_kit.run_support.get_dict_from_old_falcon_cfg(
                falcon_kit.run_support.parse_config(input_config_fn))
        # FALCON takes over the workflow for a while.
        # (For debugging, it is still possible to restart just fc_run, if desired.)
        falcon_asm_done_pfn = falcon_kit.mains.run1.run(wf, fc_cfg,
                input_config_fn,
                input_fofn_plf=input_fofn_pfn, # _pfn should be _plf, but oh well
        )
        wf.max_jobs = concurrent_jobs # in case Falcon changed this

    # Here is a hard-linking task to help us attach falcon into the dependency graph.
    falcon_link_done_pfn = makePypeLocalFile('run-falcon_link/falcon_link_done')
    make_task = PypeTask(
            inputs = {"falcon_asm_done": falcon_asm_done_pfn,},
            outputs = {
                       "falcon_link_done": falcon_link_done_pfn,
            },
            parameters = parameters,
    )
    task = make_task(start_task.task_falcon_link)
    wf.addTask(task)

    # The rest of the workflow will operate on datasets, not fasta directly.
    referenceset_pfn = makePypeLocalFile('run-fasta2referenceset/asm.referenceset.xml')
    make_task = PypeTask(
            inputs =  {"falcon_link_done": falcon_link_done_pfn,},
            outputs = {"referenceset": referenceset_pfn,},
            parameters = parameters,
    )
    task = make_task(start_task.task_fasta2referenceset)
    wf.addTask(task)
    wf.refreshTargets()

    # scatter the subreads for pbalign
    """Produces:
    pbalign_chunk.json
    chunk_subreadset_*.subreadset.xml
    """
    pbalign_chunk_json_pfn = makePypeLocalFile('run-pbalign-scatter/pbalign_chunk.json')
    make_task = PypeTask(
            inputs = {"dataset": dataset_pfn,
                      "referenceset": referenceset_pfn,},
            outputs = {"out_json": pbalign_chunk_json_pfn,},
            parameters = parameters,
    )
    task = make_task(start_task.task_pbalign_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    # After scattering, we can specify the pbalign jobs.
    tasks, alignmentset_pfn = create_tasks_pbalign(pbalign_chunk_json_pfn, referenceset_pfn, parameters)
    wf.addTasks(tasks)
    wf.refreshTargets()

    # scatter the alignmentset for genomic_consensus (variantCaller)
    """Produces:
    gc.chunks.fofn
    ???*.congitset.xml ???
    """
    gc_chunks_fofn_pfn = makePypeLocalFile('run-gc_scatter/gc.chunks.fofn')
    make_task = PypeTask(
            inputs = {"alignmentset": alignmentset_pfn,
                      "referenceset": referenceset_pfn,},
            outputs = {"out_fofn": gc_chunks_fofn_pfn,},
            parameters = parameters,
    )
    task = make_task(start_task.task_gc_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    tasks, contigset_pfn, gathered_fastq_pfn = create_tasks_gc(gc_chunks_fofn_pfn, referenceset_pfn, parameters)
    wf.addTasks(tasks)
    wf.refreshTargets()


    # Final report

    polished_assembly_report_json_pfn = makePypeLocalFile('run-polished-assembly-report/polished_assembly_report.json')
    make_task = PypeTask(
            inputs = {"referenceset": referenceset_pfn,
                      "gathered_alignmentset": alignmentset_pfn,
                      "polished_fastq": gathered_fastq_pfn,},
            outputs = {"report_json": polished_assembly_report_json_pfn,},
            parameters = parameters,
    )
    task = make_task(start_task.task_polished_assembly_report)
    wf.addTask(task)

    wf.refreshTargets()

    par_dir = os.path.dirname(fn(polished_assembly_report_json_pfn))
    sys.symlink(os.path.join(par_dir, 'polished_coverage_vs_quality.png'))
    sys.symlink(os.path.join(par_dir, 'polished_coverage_vs_quality_thumb.png'))
    #return
    ##############

    if not os.path.exists('foo.bar1'):
        sys.system('touch foo.bar1')
    foo_fn1 = makePypeLocalFile('foo.bar1')
    foo_fn2 = makePypeLocalFile('foo.bar2')
    make_task = PypeTask(
            inputs = {"foo1": foo_fn1,},
            outputs =  {"foo2": foo_fn2,},
            parameters = parameters,
    )
    task = make_task(start_task.task_foo)
    wf.addTask(task)
    wf.refreshTargets()

"""
also:
genomic_consensus.tasks.gff2vcf-0
genomic_consensus.tasks.gff2bed-0
pbcoretools.tasks.gather_gff-1

/pbi/dept/secondary/siv/smrtlink/smrtlink-alpha/smrtsuite_170220/userdata/jobs_root/000/000114/tasks
"""
