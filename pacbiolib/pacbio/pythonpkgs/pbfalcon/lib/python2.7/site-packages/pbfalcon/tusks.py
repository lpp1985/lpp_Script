from __future__ import absolute_import
from . import report_preassembly
from falcon_polish.sys import symlink, cd, say, system, filesize
from pbcore.io import FastaIO
from pbcommand.engine import run_cmd as pb_run_cmd
from falcon_kit import run_support as support
import falcon_kit.functional
import glob
import hashlib
import itertools
import json
import logging
import os
import re
import StringIO
import sys

log = logging.getLogger(__name__)


def run_cmd(cmd, *args, **kwds):
    say('RUN: %s' %repr(cmd))
    rc = pb_run_cmd(cmd, *args, **kwds)
    say(' RC: %s' %repr(rc))
    if rc.exit_code:
        raise Exception(repr(rc))

def _get_config(fn):
    """Return a dict.
    """
    return support.get_config(support.parse_config(fn))

def _get_config_from_json_fileobj(ifs_json):
    i_json = ifs_json.read()
    say('JSON=\n%s' %i_json[:1024]) # truncated
    return json.loads(i_json)

def run_falcon_config_get_fasta(input_files, output_files):
        i_config_fn, = input_files
        o_fofn_fn, = output_files
        config = _get_config(i_config_fn)
        i_fofn_fn = config['input_fofn']
        if not os.path.isabs(i_fofn_fn):
            i_fofn_fn = os.path.join(os.path.dirname(i_config_fn), i_fofn_fn)
        msg = '%r -> %r' %(i_fofn_fn, o_fofn_fn)
        say(msg)
        with cd(os.path.dirname(i_fofn_fn)):
            return support.make_fofn_abs(i_fofn_fn, o_fofn_fn)
        return 0

def run_falcon_config(input_files, output_files):
        i_config_fn, i_fasta_fofn = input_files
        o_json_fn, = output_files
        log.info('i_config_fn cont: "{}"'.format(open(i_config_fn).read()))
        config = _get_config(i_config_fn)
        config['input_fofn'] = os.path.abspath(i_fasta_fofn)
        config['original_self'] = i_config_fn
        output = json.dumps(config, sort_keys=True, indent=4, separators=(',', ': '))
        out = StringIO.StringIO()
        out.write(output)
        log.info('falcon_config:\n' + output)
        with open(o_json_fn, 'w') as ofs:
            ofs.write(output)
        #return run_cmd('echo hi', open(hello, 'w'), sys.stderr, shell=False)
        return 0

def run_falcon_make_fofn_abs(input_files, output_files):
        i_json_fn, = input_files
        o_fofn_fn, = output_files
        config = _get_config_from_json_fileobj(open(i_json_fn))
        i_fofn_fn = config['input_fofn']
        if not i_fofn_fn.startswith('/'):
            # i_fofn_fn can be relative to the location of the config file.
            original_config_fn = config['original_self']
            i_fofn_fn = os.path.join(os.path.dirname(original_config_fn), i_fofn_fn)
        msg = 'run_falcon_make_fofn_abs(%r -> %r)' %(i_fofn_fn, o_fofn_fn)
        say(msg)
        with cd(os.path.dirname(i_fofn_fn)):
            return support.make_fofn_abs(i_fofn_fn, o_fofn_fn)
        return 0

def read_fns(fofn):
    print('read from fofn:%s' %repr(fofn))
    with open(fofn) as ifs:
        return ifs.read().split()
def write_fns(fofn, fns):
    """Remember the trailing newline, for fasta2DB.
    """
    print('write to fofn:%s' %repr(fofn))
    with open(fofn, 'w') as ofs:
        return ofs.write('\n'.join(list(fns) + ['']))

def run_falcon_build_rdb(input_files, output_files):
    print('output_files: %s' %(repr(output_files)))
    cwd = os.getcwd()
    odir = os.path.realpath(os.path.abspath(os.path.dirname(output_files[0])))
    if True: #debug
        if cwd != odir:
            raise Exception('%r != %r' %(cwd, odir))
    i_json_config_fn, i_fofn_fn = input_files
    print('output_files: %s' %repr(output_files))
    run_daligner_jobs_fn, raw_reads_db_fn, job_done_fn = output_files
    config = _get_config_from_json_fileobj(open(i_json_config_fn))
    script_fn = os.path.join(odir, 'prepare_rdb.sh') # implies run-dir too
    #job_done_fn = os.path.join(odir, 'job.done') # not needed in pbsmrtpipe today tho
    support.build_rdb(i_fofn_fn, config, job_done_fn, script_fn, run_daligner_jobs_fn)
    run_cmd('bash %s' %script_fn, sys.stdout, sys.stderr, shell=False)
    job_descs = falcon_kit.functional.get_daligner_job_descriptions(open(run_daligner_jobs_fn), 'raw_reads')
    # We do not bother to calculate 'single' b/c this is only a sanity-check.
    if not job_descs:
        raise Exception("No daligner jobs generated in '%s' by '%s'." %(run_daligner_jobs_fn, script_fn))
    symlink('raw_reads.db', raw_reads_db_fn)
    return 0

def run_daligner_jobs(input_files, output_files, db_prefix='raw_reads'):
    print('run_daligner_jobs: %s %s' %(repr(input_files), repr(output_files)))
    i_json_config_fn, run_daligner_job_fn, = input_files
    o_fofn_fn, = output_files
    db_dir = os.path.dirname(run_daligner_job_fn)
    cmds = ['pwd', 'ls -al']
    fns = ['.{pre}.bps', '.{pre}.idx', '{pre}.db']
    cmds += [r'\rm -f %s' %fn for fn in fns]
    cmds += ['ln -sf {dir}/%s .' %fn for fn in fns]
    cmd = ';'.join(cmds).format(
            dir=os.path.relpath(db_dir), pre=db_prefix)
    run_cmd(cmd, sys.stdout, sys.stderr, shell=True)
    cwd = os.getcwd()
    config = _get_config_from_json_fileobj(open(i_json_config_fn))
    tasks = create_daligner_tasks(
            run_daligner_job_fn, cwd, db_prefix, db_prefix+'.db', config)
    odirs = []
    for jobd, args in tasks.items():
        with cd(jobd):
            support.run_daligner(**args)
            script_fn = args['script_fn']
            run_cmd('bash %s' %script_fn, sys.stdout, sys.stderr, shell=False)
            odirs.append(os.path.dirname(script_fn))
    write_fns(o_fofn_fn, itertools.chain.from_iterable(glob.glob('%s/*.las' %d) for d in odirs))
    return 0

#    def scripts_daligner(run_jobs_fn, db_prefix, rdb_build_done, pread_aln=False):
def create_daligner_tasks(run_jobs_fn, wd, db_prefix, db_file, config, pread_aln = False):
    tasks = dict() # uid -> parameters-dict

    nblock = support.get_nblock(db_file)

    re_daligner = re.compile(r'\bdaligner\b')

    line_count = 0
    single = (nblock == 1)
    job_descs = falcon_kit.functional.get_daligner_job_descriptions(open(run_jobs_fn), db_prefix, single=single)
    if not job_descs:
        raise Exception("No daligner jobs generated in '%s'." %run_jobs_fn)
    for desc, bash in job_descs.iteritems():
        job_uid = '%08d' %line_count
        line_count += 1
        jobd = os.path.join(wd, "./job_%s" % job_uid)
        support.make_dirs(jobd)
        call = "cd %s; ln -sf ../.%s.bps .; ln -sf ../.%s.idx .; ln -sf ../%s.db ." % (jobd, db_prefix, db_prefix, db_prefix)
        rc = os.system(call)
        if rc:
            raise Exception("Failure in system call: %r -> %d" %(call, rc))
        job_done = os.path.abspath("%s/job_%s_done" %(jobd, job_uid))
        if pread_aln:
            bash = re_daligner.sub("daligner_p", bash)
        script_fn = os.path.join(jobd , "rj_%s.sh"% (job_uid)) # also implies run-dir
        args = {
            'daligner_script': bash,
            'db_prefix': db_prefix,
            'config': config,
            'job_done': job_done,
            'script_fn': script_fn,
        }
        daligner_task = args #make_daligner_task( task_run_daligner )
        tasks[jobd] = daligner_task
    return tasks

def run_merge_consensus_jobs(input_files, output_files, db_prefix='raw_reads', dry_run=False):
    """
    dry_run --- if True, only make scripts, don't actually run them.

    if db_prefix is 'raw_reads', stage 0, make merge scripts and write to o_merge_json_fn,
    then make cons scripts and write them to o_last_json_fn.
    otherwise, db_prefix is 'preads', stage 1, make merge scripts and write to o_merge_json_fn,
    then make db2falcon scripts and write to o_last_json_fn.

    json file contains dict{p_id: dict{'script_fn':script_fn, 'script_dir':script_dir}}
    script_fn and script_dir will be used in scattered tasks later.
    """
    print('run_merge_consensus_jobs: %s %s %s' %(db_prefix, repr(input_files), repr(output_files)))
    i_json_config_fn, run_daligner_job_fn, i_fofn_fn = input_files[0:3]
    o_fofn_fn, o_merge_json_fn, o_last_json_fn = output_files[0:3]

    db_dir = os.path.dirname(run_daligner_job_fn)
    cmds = ['pwd', 'ls -al']
    fns = ['.{pre}.bps', '.{pre}.idx', '{pre}.db']
    cmds += ['rm -f %s' %fn for fn in fns]
    cmds += ['ln -sf {dir}/%s .' %fn for fn in fns]
    cmd = ';'.join(cmds).format(
            dir=os.path.relpath(db_dir), pre=db_prefix)
    run_cmd(cmd, sys.stdout, sys.stderr, shell=True)
    cwd = os.getcwd() # same as dir of o_fofn_fn
    config = _get_config_from_json_fileobj(open(i_json_config_fn))
    # i_fofn_fn has the .las files, so create_merge_tasks does not need to look for them.

    tasks = create_merge_tasks(i_fofn_fn, run_daligner_job_fn, cwd, db_prefix=db_prefix, config=config)
    # For example,
    # tasks      <- {p_id    : (merge_task, cons_task, las_bfn, fasta_bfn)}, where
    # merge_task <- {'config': config, 'merge_subdir':'m_00001', 'script':'LAmerge -v ...'}
    # cons_tasks <- {'db_fn':' abspath/falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0/raw_reads',
    #                'config': config,
    #                'las_fn': 'abspath/falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0/m_00001/raw_reads.1.las',
    #                'out_file_fn':'abspath/falcon_ns.tasks.task_falcon0_run_merge_consensus_jobs-0/preads/out.00001.fasta'}
    # las_bfn    <- 'raw_reads.1.las'
    # fasta_bfn  <- 'out.00001.fasta'

    las_fns = _run_merge_jobs(
            dict((p_id, (argstuple[0], argstuple[2])) for (p_id, argstuple) in tasks.items()),
            dry_run=dry_run, merge_json_fn=o_merge_json_fn)

    if db_prefix == 'raw_reads':
        fasta_fns = _run_consensus_jobs(
            dict((p_id, (argstuple[1], argstuple[3])) for (p_id, argstuple) in tasks.items()),
            dry_run=dry_run, cons_json_fn=o_last_json_fn)
        # Record '*.fasta' from cons jobs in FOFN.
        write_fns(o_fofn_fn, sorted(os.path.abspath(f) for f in fasta_fns))
        assert_nonzero(o_fofn_fn)
        return 0

    # Only reach this line if db_prefix is not 'raw_reads', e.g., is 'preads'

    # Record '*.las' from merge_jobs in FOFN.
    write_fns(o_fofn_fn, sorted(os.path.abspath(f) for f in las_fns))
    assert_nonzero(o_fofn_fn)

    # Generate preads4falcon.fasta from preads.db
    _run_db2falcon_jobs(cwd=cwd, config=config, dry_run=dry_run,
                        db2falcon_json_fn=o_last_json_fn)
    return 0

merged_las_fofn_bfn = 'merged_las.fofn'
#DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def _run_merge_jobs(tasks, dry_run=False, merge_json_fn=None):
    """dry_run --- if True, do not actually run the scripts,
       merge_json_fn --- if not None, write dict{p_id->mege_args} to it
    """
    fns = list()
    json_data = dict()
    for p_id, (merge_args, las_bfn) in tasks.items():
        run_dir = merge_args['merge_subdir']
        job_done = "merge_%05d_done" %p_id
        script_fn = os.path.join(run_dir, "merge_%05d.sh" % (p_id))
        merge_args['job_done'] = job_done
        merge_args['script_fn'] = script_fn
        del merge_args['merge_subdir'] # was just a temporary hack
        # merge_args <- dict{
        # 'job_done'  : 'merge_00001_done',
        # 'script_fn' : 'merge_00001.sh',
        # 'script'    : 'LAmege -v ...',
        # 'config'    : config}
        support.run_las_merge(**merge_args)
        mkdir(run_dir)
        with cd(run_dir):
            if dry_run is False:
                run_cmd('bash %s' %os.path.basename(script_fn), sys.stdout, sys.stderr, shell=False)
        fns.append(os.path.join(run_dir, las_bfn))

        # add script_dir to args for scattered tasks to work in the correct dir
        json_data[p_id] = {'script_dir': os.path.join(os.getcwd(), run_dir), # 'merge_00001.sh'
                           'script_fn': os.path.basename(script_fn)} # '/pbi/.../tasks/falcon_ns.task.task_falcon0_run_merge_consensus_jobs/m_00001',
        json_fn = 'merge_jobs.json' if merge_json_fn is None else merge_json_fn
    # Write dict{p_id: dict{'script_fn':script_fn, 'script_dir':script_dir}} to a json file
    with open(json_fn, 'w') as writer:
        writer.write(json.dumps(json_data) + "\n")

    return fns # *.las, e.g., ['m_00001/raw_reads.1.las', 'm_00002/raw_reads.2.las', 'm_00003/raw_reads.3.las']

def _run_consensus_jobs(tasks, dry_run=False, cons_json_fn=None):
    """dry_run --- if True, do not actually run the scripts
       cons_json_fn --- if not None, write dict{p_id: dict{'script_fn':script_fn, 'script_dir':script_dir}} to it
    """
    fns = list()
    json_data = dict()
    for p_id, (cons_args, fasta_bfn) in tasks.items():
        run_dir = 'preads'
        job_done = "c_%05d_done" %p_id
        script_fn = os.path.join(run_dir, "c_%05d.sh" %(p_id))
        cons_args['job_done'] = job_done
        cons_args['script_fn'] = script_fn
        # cons_args <- dict{
        # 'out_file_fn': abspath to preads/out.00001.fasta
        # 'script_fn'  : c_00001.sh
        # 'job_done'   : c_00001_done
        # 'raw_reads'  : raw_reads
        # 'config'     : config}
        support.run_consensus(**cons_args)
        mkdir(run_dir)
        with cd(run_dir):
            if dry_run is False:
                run_cmd('bash %s' %os.path.basename(script_fn), sys.stdout, sys.stderr, shell=False)
        fns.append(os.path.join(run_dir, fasta_bfn))

        # add script_dir to args for scattered tasks to work in the correct dir
        json_data[p_id] = {'script_fn': os.path.basename(script_fn), # 'c_00001.sh'
                           'script_dir': os.path.join(os.getcwd(), run_dir)} # '/pbi/.../tasks/falcon_ns.tasks.task_falcon0_run_merge_jobs/preads/'

    json_fn = "cons_jobs.json" if cons_json_fn is None else cons_json_fn
    with open(json_fn, 'w') as writer:
        writer.write(json.dumps(json_data) + "\n")
    return fns # *.fasta ['preads/out.0001.fasta', 'preads/out.00002.fasta', 'preads/out.00003.fasta']

def _run_db2falcon_jobs(cwd, config, dry_run, db2falcon_json_fn=None):
    """
    cwd --- current workding directory
    dry_run --- if True, do not actually run the scripts
    db2falcon_json_fn --- if not None, write dict{0: dict('script_fn':script_fn, 'script_dir':script_dir)}
    """
    # Generate preads4falcon.fasta from preads.db
    script_fn = os.path.join(cwd, "run_db2falcon.sh")
    job_done = script_fn + '_done'
    args = {
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
        'preads4falcon_fn': 'preads4falcon.fasta',
        'preads_db': 'preads.db',
    }
    support.run_db2falcon(**args)
    json_fn = "db2falcon.json" if db2falcon_json_fn is None else db2falcon_json_fn

    # add script_dir to args for scattered tasks to work in the correct dir
    json_data = {0: {'script_fn':os.path.basename(script_fn), 'script_dir': cwd}}
    with open(json_fn, 'w') as writer:
        writer.write(json.dumps(json_data) + "\n")

    mkdir(cwd)
    with cd(cwd):
        if dry_run is False:
            run_cmd('bash %s' %os.path.basename(script_fn), sys.stdout, sys.stderr, shell=False)
            assert_nonzero('preads4falcon.fasta')

def run_scripts_in_json(input_files, output_files):
    """
    input_files = ['*.json'] (e.g., merge|cons|db2falcon.json), where
                  *.json <- dict(p_id: dict('script_fn':script_fn, 'script_dir':script_dir))
    output_files = ['*_done.txt'] (e.g., merge_done.txt, cons_done.txt, db2falcon_done.txt)
    execute all script files.
    """
    json_fn = input_files[0]
    txt_fn = output_files[0]

    a = json.load(open(json_fn, 'r'))

    writer = open(txt_fn, 'w')
    for p_id, args in a.iteritems():
        if 'script_fn' not in args:
            raise ValueError("Could not find 'script_fn' in json %r key %r" % (json_fn, p_id))
        script_dir = str(args['script_dir'])
        with cd(script_dir):
            script_fn = str(args['script_fn'])
            run_cmd('bash %s' % script_fn, sys.stdout, sys.stderr, shell=False)
            writer.write(script_fn + "\n")
    writer.close()
    return 0


def create_merge_tasks(i_fofn_fn, run_jobs_fn, wd, db_prefix, config):
    #merge_scripts = bash.scripts_merge(config, db_prefix, run_jobs_fn)
    tasks = {} # pid -> (merge_params, cons_params)
    mjob_data = {}

    with open(run_jobs_fn) as f :
        for l in f:
            l = l.strip().split()
            if l[0] not in ( "LAsort", "LAmerge", "mv" ):
                continue
            if l[0] == "LAsort":
                # We now run this part w/ daligner, but we still need
                # a small script for some book-keeping.
                p_id = int( l[2].split(".")[1] )
                mjob_data.setdefault( p_id, [] )
                #mjob_data[p_id].append(  " ".join(l) ) # Already done w/ daligner!
            if l[0] == "LAmerge":
                l2 = l[2].split(".")
                if l2[1][0] == "L":
                    p_id = int(  l[2].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[2].split(".")[1] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
            if l[0] == "mv":
                l2 = l[1].split(".")
                if l2[1][0] == "L":
                    p_id = int(  l[1].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[1].split(".")[1] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )

    # Could be L1.* or preads.*
    re_las = re.compile(r'\.(\d*)(\.\d*)?\.las$')

    for p_id in mjob_data:
        s_data = mjob_data[p_id]

        support.make_dirs("%s/preads" % (wd) )
        support.make_dirs("%s/las_files" % (wd) )
        merge_subdir = "m_%05d" %p_id
        merge_dir = os.path.join(wd, merge_subdir)
        support.make_dirs(merge_dir)
        #merge_script_file = os.path.abspath( "%s/m_%05d/m_%05d.sh" % (wd, p_id, p_id) )
        merge_script = StringIO.StringIO()
        with cd(merge_dir):
            print("i_fofn_fn=%r" %i_fofn_fn)
            # Since we could be in the gather-task-dir, instead of globbing,
            # we will read the fofn.
            for fn in open(i_fofn_fn).read().splitlines():
                basename = os.path.basename(fn)
                mo = re_las.search(basename)
                if not mo:
                    continue
                left_block = int(mo.group(1))
                if left_block != p_id:
                    # By convention, m_00005 merges L1.5.*.las, etc.
                    continue
                symlink(fn)

        for l in s_data:
            print >> merge_script, l
        las_bfn = '%s.%d.las' %(db_prefix, p_id)
        #print >> merge_script, 'echo %s >| %s' %(las_bfn, merged_las_fofn_bfn)

        #job_done = makePypeLocalFile(os.path.abspath( "%s/m_%05d/m_%05d_done" % (wd, p_id, p_id)  ))
        parameters =  {"script": merge_script.getvalue(),
                       "merge_subdir": merge_subdir,
                       "config": config}
        merge_task = parameters

        fasta_bfn = "out.%05d.fasta" %p_id
        out_file_fn = os.path.abspath("%s/preads/%s" %(wd, fasta_bfn))
        #out_done = makePypeLocalFile(os.path.abspath( "%s/preads/c_%05d_done" % (wd, p_id)  ))
        parameters =  {
                       "db_fn": '{}/{}'.format(os.getcwd(), db_prefix),
                       "las_fn": '{}/{}/{}'.format(os.getcwd(), merge_subdir, las_bfn), # assuming merge ran in merge_dir
                       "out_file_fn": out_file_fn,
                       #"out_done": out_done,
                       "config": config}
        cons_task = parameters
        tasks[p_id] = (merge_task, cons_task, las_bfn, fasta_bfn)
        # tasks <- dict{p_id: (merge_task, cons_task, las_bfn, fasta_bfn)}, where
        # p_id is an integer, e.g., 1
        # merge_task <- dict{'merge_dir': 'm_00001', 'script'='LAmerge -v raw_reads.1 L1.1.1 L1.1.2 L1.1.3', 'script_fn': x}
        # cons_task  <- dict{'db_fn':x, 'las_fn':x, 'out_file_fn':x, 'config':config}
        # las_bfn, e.g. raw_reads.1.las
        # fasta_bfn, e.g., out.00001.fasta

    return tasks

def run_falcon_build_pdb(input_files, output_files):
    print('output_files: %s' %(repr(output_files)))
    cwd = os.getcwd()
    odir = os.path.realpath(os.path.abspath(os.path.dirname(output_files[0])))
    if True: #debug
        if cwd != odir:
            raise Exception('%r != %r' %(cwd, odir))
    i_json_config_fn, i_fofn_fn = input_files[0:2]
    print('output_files: %s' %repr(output_files))
    run_daligner_jobs_fn, = output_files
    config = _get_config_from_json_fileobj(open(i_json_config_fn))
    script_fn = os.path.join(odir, 'prepare_pdb.sh')
    job_done_fn = os.path.join(odir, 'job_done')
    support.build_pdb(i_fofn_fn, config, job_done_fn, script_fn, run_daligner_jobs_fn)
    run_cmd('bash %s' %script_fn, sys.stdout, sys.stderr, shell=False)
    job_descs = falcon_kit.functional.get_daligner_job_descriptions(open(run_daligner_jobs_fn), 'preads')
    # We do not bother to calculate 'single' b/c this is only a sanity-check.
    if not job_descs:
        raise Exception("No daligner jobs generated in '%s' by '%s'." %(run_daligner_jobs_fn, script_fn))
    return 0

def _linewrap_fasta(ifn, ofn):
    """For the pbsmrtpipe validator.
    Not sure whether any tools actually require this.
    """
    n = 0
    with FastaIO.FastaReader(ifn) as fa_in:
        with FastaIO.FastaWriter(ofn) as fa_out:
            for rec in fa_in:
                n += 1
                fa_out.writeRecord(rec)
    return n

def run_falcon_asm(input_files, output_files):
    i_json_config_fn, i_fofn_fn = input_files[0:2]
    o_fasta_fn = output_files[0]
    cwd = os.getcwd()
    pread_dir = os.path.dirname(i_fofn_fn)
    preads4falcon_fasta_fn = os.path.join(pread_dir, 'preads4falcon.fasta')
    db_file = os.path.join(pread_dir, 'preads.db')
    job_done = os.path.join(cwd, 'job_done')
    config = _get_config_from_json_fileobj(open(i_json_config_fn))
    script_fn = os.path.join(cwd ,"run_falcon_asm.sh")
    args = {
        'las_fofn_fn': i_fofn_fn,
        'preads4falcon_fasta_fn': preads4falcon_fasta_fn,
        'db_file_fn': db_file,
        'config': config,
        'job_done': job_done,
        'script_fn': script_fn,
    }
    assert_nonzero(i_fofn_fn)
    assert_nonzero(preads4falcon_fasta_fn)
    support.run_falcon_asm(**args)
    run_cmd('bash %s' %script_fn, sys.stdout, sys.stderr, shell=False)
    p_ctg = 'p_ctg.fa'
    assert_nonzero(p_ctg)
    n_records = _linewrap_fasta(p_ctg, o_fasta_fn)
    if n_records == 0:
        # We already checked 0-length, but maybe this is still possible.
        # Really, we want to detect 0 base-length, but I do not know how yet.
        raise Exception("No records found in primary contigs: '%s'" %os.path.abspath(p_ctg))
    say('Finished run_falcon_asm(%s, %s)' %(repr(input_files), repr(output_files)))
    return 0


def run_rm_las(input_files, output_files, prefix):
    """ Delete all intermediate las files. """
    cmd = "pwd && find .. -type f -name '%s*.las' -delete -print" % prefix
    say(cmd)
    run_cmd(cmd, sys.stdout, sys.stderr)
    with open(output_files[0], 'w') as writer:
        writer.write("#%s" % cmd)
    return 0


def run_hgap(input_files, output_files, tmpdir):
    i_cfg_fn, i_logging_fn, i_subreadset_fn = input_files
    o_preads_fasta_fn, \
    o_polished_fasta_fn, o_polished_fastq_fn, o_polished_csv_fn, \
    o_aligned_subreads_fn, o_alignment_summary_gff_fn, o_unmapped_subreads_txt_fn, \
    o_contigset_fn, o_preass_json_fn, o_polass_json_fn, o_log_fn, = output_files
    # Update the logging-cfg with our log-file.
    logging_cfg = json.loads(open(i_logging_fn).read())
    logging_cfg['handlers']['handler_file_all']['filename'] = o_log_fn
    logging_fn = 'logging.json'
    with open(logging_fn, 'w') as ofs:
        ofs.write(json.dumps(logging_cfg))
    # Update the cfg with our subreadset. (Inside hgap_run?)
    # Run pypeflow.hgap.main.
    cmd = 'TMPDIR={tmpdir} python -m pbfalcon.cli.hgap_run --logging {logging_fn} {i_cfg_fn}'.format(**locals())
    system(cmd)
    # Write Reports
    with open('run-falcon/0-rawreads/report/pre_assembly_stats.json') as stats_ifs: # by convention
        with open(o_preass_json_fn, 'w') as report_ofs:
            report_preassembly.write_report_from_stats(stats_ifs, report_ofs)
    # Symlink expected outputs, by convention.
    symlink('run-falcon/1-preads_ovl/db2falcon/preads4falcon.fasta', o_preads_fasta_fn)
    symlink('run-gc-gather/contigset.fasta', o_polished_fasta_fn)
    symlink('run-gc-gather/gathered.fastq', o_polished_fastq_fn)
    symlink('run-polished-assembly-report/polished_coverage_vs_quality.csv', o_polished_csv_fn)
    symlink('run-polished-assembly-report/alignment.summary.gff', o_alignment_summary_gff_fn)
    symlink('run-pbalign_gather/aligned.subreads.alignmentset.xml', o_aligned_subreads_fn)
    symlink('run-pbalign_gather/unmapped.txt', o_unmapped_subreads_txt_fn)
    symlink('run-gc-gather/contigset.xml', o_contigset_fn)
    symlink('run-polished-assembly-report/polished_assembly_report.json', o_polass_json_fn)
    return 0

def run_report_preassembly_yield(input_files, output_files):
    i_json_config_fn, i_preads_fofn_fn, i_raw_reads_db_fn = input_files[0:3]
    o_json_fn, = output_files
    kwds = {
        'i_json_config_fn': i_json_config_fn,
        'i_raw_reads_db_fn': i_raw_reads_db_fn,
        'i_preads_fofn_fn': i_preads_fofn_fn,
        'o_json_fn': o_json_fn,
    }
    report_preassembly.for_task(**kwds)
    return 0

def assert_nonzero(fn):
    if filesize(fn) == 0:
        raise Exception("0-length filesize for: '%s'" %os.path.abspath(fn))
