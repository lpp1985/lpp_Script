# Much of this was in pbsmrtpipe/tools/chunk_utils.py
from falcon_kit.functional import (get_daligner_job_descriptions, get_script_xformer)
from pbcommand.models import PipelineChunk
from pbcommand.pb_io import write_pipeline_chunks
from falcon_polish.functional import joined_strs
from falcon_polish.sys import symlink
import logging
import os

log = logging.getLogger(__name__)

def symlink_dazzdb(actualdir, db_prefix):
    """Symlink elements of dazzler db.
    For now, 3 files.
    """
    symlink(os.path.join(actualdir, '.%s.bps'%db_prefix))
    symlink(os.path.join(actualdir, '.%s.idx'%db_prefix))
    symlink(os.path.join(actualdir, '%s.db'%db_prefix))

def write_run_daligner_chunks_falcon(
        pread_aln,
        chunk_file,
        config_json_fn,
        run_jobs_fn,
        max_total_nchunks,
        dir_name,
        chunk_base_name,
        chunk_ext,
        chunk_keys):
    db_prefix = 'preads' if pread_aln else 'raw_reads'
    xform_script = get_script_xformer(pread_aln)
    def chunk():
        # cmds is actually a list of small bash scripts, including linefeeds.
        cmds = get_daligner_job_descriptions(open(run_jobs_fn), db_prefix).values()
        if max_total_nchunks < len(cmds):
            log.debug("max_total_nchunks < # daligner cmds: %d < %d" %(
                max_total_nchunks, len(cmds)))
            cmds = joined_strs(cmds, max_total_nchunks)
        symlink_dazzdb(os.path.dirname(run_jobs_fn), db_prefix)
        for i, script in enumerate(cmds):
            chunk_id = '_'.join([chunk_base_name, str(i)])
            chunk_name = '.'.join([chunk_id, chunk_ext])
            chunk_path = os.path.join(dir_name, chunk_name)
            script = xform_script(script)
            open(chunk_path, 'w').write(script)
            d = {}
            d[chunk_keys[1]] = os.path.abspath(chunk_path)
            d[chunk_keys[0]] = config_json_fn
            c = PipelineChunk(chunk_id, **d)
            yield c
    chunks = list(chunk())
    write_pipeline_chunks(chunks, chunk_file, comment=None)
