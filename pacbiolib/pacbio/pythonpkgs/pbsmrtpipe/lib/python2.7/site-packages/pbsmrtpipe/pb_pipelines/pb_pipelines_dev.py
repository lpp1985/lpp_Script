import logging
import itertools

from pbsmrtpipe.core import register_pipeline
from pbsmrtpipe.constants import to_pipeline_ns

from .pb_pipelines_sa3 import Constants, Tags

log = logging.getLogger(__name__)


def dev_register(relative_id, display_name, tags=(), task_options=None):
    pipeline_id = to_pipeline_ns(relative_id)
    ptags = list(set(tags + (Tags.DEV, )))
    return register_pipeline(pipeline_id, display_name, "0.1.0", tags=ptags, task_options=task_options)


def _get_example_options():
    return {"pbsmrtpipe.task_options.dev_message": "Preset Custom Dev Message from register pipeline"}


@dev_register("dev_local", "Dev Local Hello Pipeline", tags=("local", ), task_options=_get_example_options())
def get_dev_local_pipeline():
    """Simple example pipeline"""
    b1 = [('$entry:e_01', 'pbsmrtpipe.tasks.dev_hello_world:0')]

    b2 = [('pbsmrtpipe.tasks.dev_hello_world:0', 'pbsmrtpipe.tasks.dev_hello_worlder:0'),
          ('pbsmrtpipe.tasks.dev_hello_world:0', 'pbsmrtpipe.tasks.dev_hello_garfield:0')]

    b3 = [('pbsmrtpipe.tasks.dev_hello_world:0', 'pbsmrtpipe.tasks.dev_txt_to_fasta:0')]

    b4 = [('pbsmrtpipe.tasks.dev_txt_to_fasta:0', 'pbsmrtpipe.tasks.dev_filter_fasta:0')]

    return b1 + b2 + b3 + b4


@dev_register("dev_01", "Example Dev 01 pipeline", task_options=_get_example_options())
def f():
    """Simplest possible pipeline, a single Task"""
    b = [("$entry:e_01", "pbsmrtpipe.tasks.dev_hello_world:0")]
    return b


@dev_register("dev_01_fail", "Example Dev 01 crashing pipeline")
def f():
    """Simple pipeline that deliberately and immediately crashes"""
    b = [("$entry:e_01", "pbsmrtpipe.tasks.dev_raises_exception:0")]
    return b


@dev_register("dev_01_ds", "Example Dev 01 Subread DataSet pipeline", tags=(Tags.RPT, ))
def f():
    """Simplest possible pipeline, a single Task"""
    b = [("$entry:eid_subread", "pbsmrtpipe.tasks.dev_subread_report:0")]
    return b


@dev_register("dev_02", "Example Dev 02 pipeline")
def f():
    # Reuse existing pipeline and reference a specific task output
    # in the pipeline. The format is {pipeline_id}:{task_id}:{task_out_index}
    # route the output of the existing pipeline to a new task input
    b = [("pbsmrtpipe.pipelines.dev_01:pbsmrtpipe.tasks.dev_hello_world:0", "pbsmrtpipe.tasks.dev_hello_worlder:0")]

    # Redefine new entry points of existing pipelines
    # The format is {pipeline_id}:{entry_label_id}
    b2 = [("$entry:e_txt", "pbsmrtpipe.pipelines.dev_01:$entry:e_01")]

    return b + b2


@dev_register("dev_03", "Example Dev 03 pipelines", tags=('dev', ))
def f():
    """Reuse of pipeline 02"""
    b = [("pbsmrtpipe.pipelines.dev_02:pbsmrtpipe.tasks.dev_hello_worlder:0", "pbsmrtpipe.tasks.dev_hello_garfield:0")]

    b2 = [("$entry:e_txt2", "pbsmrtpipe.pipelines.dev_02:$entry:e_txt")]

    return b + b2


@dev_register("dev_04", "Example Dev 04")
def f():
    b0 = [("pbsmrtpipe.pipelines.dev_03:pbsmrtpipe.tasks.dev_hello_garfield:0", "pbsmrtpipe.tasks.dev_txt_to_fofn:0")]
    b1 = [("pbsmrtpipe.pipelines.dev_03:pbsmrtpipe.tasks.dev_hello_garfield:0", "pbsmrtpipe.tasks.dev_txt_to_fasta:0")]

    # This needs to be ported
    b2 = [("pbsmrtpipe.tasks.dev_txt_to_fasta:0", "pbsmrtpipe.tasks.dev_static_task_tc:0")]

    b3 = [("$entry:e_txt3", "pbsmrtpipe.pipelines.dev_03:$entry:e_txt2")]

    return b0 + b1 + b3


@dev_register("dev_04_w_static_task", "Pipeline that leverages a python static tasks")
def f():

    # pbsmrtpipe.tasks.dev_static_txt_task needs to be ported
    #b = [("pbsmrtpipe.pipelines.dev_03:pbsmrtpipe.tasks.dev_hello_garfield:0", "pbsmrtpipe.tasks.dev_static_txt_task:0")]

    b1 = [("pbsmrtpipe.pipelines.dev_03:pbsmrtpipe.tasks.dev_hello_garfield:0", "pbsmrtpipe.tasks.dev_txt_to_fasta:0")]

    b2 = [("pbsmrtpipe.tasks.dev_txt_to_fasta:0", "pbsmrtpipe.tasks.dev_filter_fasta:0")]

    # this needs to be ported to the new ToolContract interface
    #b2 = [("pbsmrtpipe.tasks.dev_txt_to_fasta:0", "pbsmrtpipe.tasks.dev_static_task_tc:0")]

    b3 = [("$entry:e_txt3", "pbsmrtpipe.pipelines.dev_03:$entry:e_txt2")]

    return b1 + b2 + b3


@dev_register("dev_local_chunk", "Dev Local Hello Chunkable Pipeline", tags=("local", "chunk"))
def get_dev_local_chunk():
    """Simple example pipeline"""
    b1 = [("$entry:e_01", "pbsmrtpipe.tasks.dev_txt_to_fofn:0")]

    b3 = [
        ("pbsmrtpipe.tasks.dev_txt_to_fofn:0", "pbsmrtpipe.tasks.dev_txt_to_fofn_report:0")]

    # Add Chunk-able task using dev_fofn chunk operator
    b4 = [("pbsmrtpipe.tasks.dev_txt_to_fofn_report:0", "pbsmrtpipe.tasks.dev_fofn_example:0"),
          ("pbsmrtpipe.tasks.dev_txt_to_fofn_report:1", "pbsmrtpipe.tasks.dev_fofn_example:1")]

    # Add a task to the chunked output of the txt
    b5 = [("pbsmrtpipe.tasks.dev_txt_to_fofn_report:0", "pbsmrtpipe.tasks.dev_hello_worlder:0")]

    return b1 + b3 + b4 + b5


@dev_register("dev_local_fasta_chunk", "Dev local Task for Chunking pipelines", tags=('chunking',))
def get_dev_local_pipeline_chunk():

    b1 = [("pbsmrtpipe.pipelines.dev_03:pbsmrtpipe.tasks.dev_hello_garfield:0", "pbsmrtpipe.tasks.dev_txt_to_fasta:0")]

    b2 = [("pbsmrtpipe.tasks.dev_txt_to_fasta:0", "pbsmrtpipe.tasks.dev_filter_fasta:0")]

    b3 = [("pbsmrtpipe.tasks.dev_filter_fasta:0", "pbsmrtpipe.tasks.dev_tc_fasta_report:0")]

    return b1 + b2 + b3


@dev_register("dev_dist", "Dev Hello Distributed Workflow Pipeline", tags=('cluster', ))
def get_dist_dev_pipeline():
    """Simple distributed example pipeline"""
    bs = get_dev_local_pipeline()

    b2 = [('pbsmrtpipe.tasks.dev_hello_world:0', 'pbsmrtpipe.tasks.dev_hello_distributed:0'),
          ('pbsmrtpipe.tasks.dev_hello_worlder:0', 'pbsmrtpipe.tasks.dev_hello_distributed:1')]

    return bs + b2


@dev_register("dev_diagnostic", "Reference Set Report", tags=(Tags.RPT, "reference"))
def get_reference_ds_report():
    """Generate a simple report and plot from Reference DataSet"""

    b = [(Constants.ENTRY_DS_REF, "pbsmrtpipe.tasks.dev_reference_ds_report:0")]

    return b


@dev_register("dev_diagnostic_stress", "Reference Set Report + Random tasks", tags=(Tags.RPT, "reference", "stress"))
def get_reference_ds_report():
    """Generate a simple report and plot from Reference DataSet"""

    max_txt_tasks = 50

    def to_b(i):
        return [("pbsmrtpipe.tasks.rset_to_txt:0", "pbsmrtpipe.tasks.dev_hello_worlder:{i}:0".format(i=i))]

    b1 = [(Constants.ENTRY_DS_REF, "pbsmrtpipe.tasks.dev_reference_ds_report:0")]

    b2 = [(Constants.ENTRY_DS_REF, 'pbsmrtpipe.tasks.rset_to_txt:0')]

    bs = list(itertools.chain(*[to_b(x) for x in xrange(max_txt_tasks)]))

    return b1 + b2 + bs
