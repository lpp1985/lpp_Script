import json
import os
import logging
import pprint
import shutil
import uuid

from pbcommand.models import DataStore, DataStoreFile, FileTypes
from pbcommand.models.report import Attribute, Report, Table, Column, Plot, PlotGroup

import pbsmrtpipe

import pbsmrtpipe.report_renderer as R
import pbsmrtpipe.graph.bgraph as B
import pbsmrtpipe.graph.bgraph_utils as BU
import pbsmrtpipe.pb_io as IO
from pbsmrtpipe.graph.models import VALID_ALL_TASK_NODE_CLASSES
from pbsmrtpipe.models import TaskStates, JobResources, RunnableTask
from pbsmrtpipe.utils import setup_log, setup_internal_logs

log = logging.getLogger(__name__)
slog = logging.getLogger('status.' + __name__)


def _to_table(tid, bg, nodes):
    """Create a table from File nodes or Entry nodes"""
    columns = [Column('id', header="Id"),
               Column('is_resolved', header='Is Resolved'),
               Column('path', header="Path")]

    table = Table(tid, columns=columns)
    for node in nodes:
        table.add_data_by_column_id('id', str(node))
        table.add_data_by_column_id('is_resolved', bg.node[node]['is_resolved'])
        try:
            table.add_data_by_column_id('path', bg.node[node]['path'])
        except KeyError as e:
            slog.error("Failed to get path from {n}".format(n=repr(node)))
            slog.error(e)
            table.add_data_by_column_id('path', "NA")

    return table


def _to_report(bg, job_output_dir, job_id, state, was_successful, run_time, error_message=None):
    """ High Level Report of the workflow state

    Write the output of workflow datastore to pbreports report object

    Workflow summary .dot/svg (collapsed workflow)
    Workflow details .dot/svg (chunked workflow)

    To add:
    - Resolved WorkflowSettings (e.g., nproc, max_workers)
    -

    :type bg: BindingsGraph

    """
    emsg = "" if error_message is None else error_message

    attributes = [Attribute('was_successful', was_successful, name="Was Successful"),
                  Attribute('total_run_time_sec', int(run_time), name="Walltime (sec)"),
                  Attribute('error_message', emsg, name="Error Message"),
                  Attribute('job_id', job_id, name="Job Id"),
                  Attribute('job_state', state, name="Job State"),
                  Attribute('job_output_dir', job_output_dir, name="Job Output Directory"),
                  Attribute('pbsmrtpipe_version', pbsmrtpipe.get_version(), name="pbsmrtpipe Version")]

    columns = [Column('task_id', header='Task id'),
               Column('was_successful', header='Was Successful'),
               Column('state', header="Task State"),
               Column('run_time_sec', header="Run Time (sec)"),
               Column('nproc', header="# of procs")]

    tasks_table = Table('tasks', columns=columns)
    for tnode in bg.all_task_type_nodes():
        tasks_table.add_data_by_column_id('task_id', str(tnode))
        tasks_table.add_data_by_column_id('nproc', bg.node[tnode]['nproc'])
        tasks_table.add_data_by_column_id('state', bg.node[tnode]['state'])
        tasks_table.add_data_by_column_id('was_successful', bg.node[tnode]['state'] == TaskStates.SUCCESSFUL)
        # rt_ = bg.node[tnode]['run_time']
        # rtime = None if rt_ is None else int(rt_)
        tasks_table.add_data_by_column_id('run_time_sec', bg.node[tnode]['run_time'])

    ep_table = _to_table("entry_points", bg, bg.entry_binding_nodes())
    fnodes_table = _to_table("file_node", bg, bg.file_nodes())

    report = Report('pbsmrtpipe', tables=[tasks_table, ep_table, fnodes_table],
                    attributes=attributes)
    return report


def _dict_to_report_table(table_id, key_attr, value_attr, d):
    """
    General {k->v} to create a pbreport Table

    :param table_id: Table id
    :param key_attr: Column id
    :param value_attr: Column id
    :param d: dict
    :return:
    """
    columns = [Column(key_attr, header="Attribute"),
               Column(value_attr, header="Value")]

    table = Table(table_id, columns=columns)
    for k, v in d.iteritems():
        table.add_data_by_column_id(key_attr, k)
        table.add_data_by_column_id(value_attr, v)

    return table


def _workflow_opts_to_table(workflow_opts):
    """
    :type workflow_opts: WorkflowLevelOptions
    :param workflow_opts:
    :return:
    """
    tid = "workflow_options"
    wattr = "workflow_attribute"
    wvalue = "workflow_value"
    return _dict_to_report_table(tid, wattr, wvalue, workflow_opts.to_dict())


def _task_opts_to_table(task_opts):
    tattr = "task_attribute_id"
    vattr = "task_value"
    # FIXME this is not being generated correctly
    return _dict_to_report_table("task_options", tattr, vattr, task_opts)


def _to_workflow_settings_report(bg, workflow_opts, task_opts, state, was_successful):

    tables = [_workflow_opts_to_table(workflow_opts), _task_opts_to_table(task_opts)]
    report = Report("workflow_settings_report", tables=tables)
    return report


def to_task_summary_report(bg):

    cs = [Column("workflow_task_id", header="Task Id"),
          Column("workflow_task_status", header="Status"),
          Column("workflow_task_run_time", header="Task Runtime"),
          Column('workflow_task_nproc', header="Number of Procs"),
          Column("workflow_task_emsg", header="Error Message")]

    t = Table("workflow_task_summary", title="Task Summary", columns=cs)
    for tnode in bg.all_task_type_nodes():
        if isinstance(tnode, VALID_ALL_TASK_NODE_CLASSES):
            t.add_data_by_column_id("workflow_task_id", tnode.idx)
            t.add_data_by_column_id("workflow_task_status", bg.node[tnode]['state'])
            t.add_data_by_column_id("workflow_task_run_time", bg.node[tnode]['run_time'])
            t.add_data_by_column_id("workflow_task_nproc", bg.node[tnode]['nproc'])
            t.add_data_by_column_id("workflow_task_emsg", bg.node[tnode]['error_message'])

    return Report("workflow_task_summary", tables=[t])


def _to_workflow_report(job_resources, bg, workflow_opts, task_opts, state, was_successful, plot_images):
    """
    Copy images to image local directory and return a pbreport Report

    """
    plot_groups = []
    if plot_images:
        plots = []
        for i, plot_image in enumerate(plot_images):
            html_image_abs = os.path.join(job_resources.images, os.path.basename(plot_image))
            shutil.copy(plot_image, html_image_abs)
            # Make the file path relative to images/my-image.png
            html_image = os.path.join(os.path.basename(job_resources.images), os.path.basename(plot_image))
            p = Plot("plot_{i}".format(i=i), html_image)
            plots.append(p)

        pg = PlotGroup("workflow_state_plots", plots=plots)
        plot_groups.append(pg)

    return Report("workflow_report", plotgroups=plot_groups)


def datastore_to_report(ds):
    """

    :type ds: DataStore
    :param ds:
    :return:
    """
    attrs = [Attribute("ds_nfiles", len(ds.files), name="Number of files"),
             Attribute("ds_version", ds.version, name="Datastore version"),
             Attribute("ds_created_at", ds.created_at, name="Created At"),
             Attribute("ds_updated_at", ds.updated_at, name="Updated At")]

    columns_names = [("file_id", "File Id"),
                     ("file_type_obj", "File Type"),
                     ("path", "Path"),
                     ("file_size", "Size"),
                     ("created_at", "Created At"),
                     ("modified_at", "Modified At")]

    to_i = lambda s: "ds_" + s
    columns = [Column(to_i(i), header=h) for i, h in columns_names]
    t = Table("datastore", title="DataStore Summary", columns=columns)

    def _to_relative_path(p):
        return "/".join(p.split("/")[-3:])

    for file_id, ds_file in ds.files.iteritems():
        t.add_data_by_column_id(to_i("file_id"), ds_file.file_id)
        t.add_data_by_column_id(to_i("file_type_obj"), ds_file.file_type_id)
        t.add_data_by_column_id(to_i("path"), _to_relative_path(ds_file.path))
        t.add_data_by_column_id(to_i("file_size"), ds_file.file_size)
        t.add_data_by_column_id(to_i("created_at"), ds_file.created_at)
        t.add_data_by_column_id(to_i("modified_at"), ds_file.modified_at)

    return Report("datastore_report", tables=[t], attributes=attrs)


def _get_images_in_dir(dir_name, formats=(".png", ".svg")):
    # report plots only support
    return [os.path.join(dir_name, i_) for i_ in os.listdir(dir_name) if any(i_.endswith(x) for x in formats)]


def write_main_workflow_report(job_id, job_resources, workflow_opts, task_opts, bg_, state_, was_successful_, run_time_sec):
    """
    Write the main workflow level report.

    :type job_resources: JobResources
    :type workflow_opts: WorkflowLevelOptions
    :type bg_: BindingsGraph
    :type was_successful_: bool
    :type run_time_sec: float

    :param job_id:
    :param job_resources:
    :param workflow_opts:
    :param task_opts:
    :param bg_:
    :param state_:
    :param was_successful_:
    :param run_time_sec:
    :return:
    """
    # workflow_json = os.path.join(job_resources.workflow, 'workflow.json')
    # with open(workflow_json, 'w+') as f:
    #     f.write(json.dumps(json_graph.node_link_data(bg_)))

    report_path = os.path.join(job_resources.workflow, 'report-tasks.json')
    report_ = _to_report(bg_, job_resources.root, job_id, state_, was_successful_, run_time_sec)
    report_.write_json(report_path)
    R.write_report_with_html_extras(report_, os.path.join(job_resources.root, 'index.html'), job_resources.html)

    setting_report = _to_workflow_settings_report(bg_, workflow_opts, task_opts, state_, was_successful_)
    R.write_report_to_html(setting_report, os.path.join(job_resources.html, 'settings.html'))

    setting_report = _to_workflow_report(job_resources, bg_, workflow_opts, task_opts, state_, was_successful_, _get_images_in_dir(job_resources.workflow, formats=(".svg",)))
    R.write_report_to_html(setting_report, os.path.join(job_resources.html, 'workflow.html'))


def write_update_main_workflow_report(job_id, job_resources, bg_, state_, was_successful_, run_time_sec):
    """
    This will only update the index.html with the current state of each task
    """

    report_path = os.path.join(job_resources.workflow, 'report-tasks.json')
    report_ = _to_report(bg_, job_resources.root, job_id, state_, was_successful_, run_time_sec)
    report_.write_json(report_path)

    R.write_report_to_html(report_, os.path.join(job_resources.root, 'index.html'))

    return True


def write_task_manifest(manifest_path, tid, task, resource_types, task_version, python_mode_str, cluster_renderer):
    """

    :type manifest_path: str
    :type tid: str
    :type task: Task
    :type task_version: str
    :type python_mode_str: str
    :type cluster_renderer: ClusterTemplateRender | None

    :return:
    """
    # this should be already set in task, or is this for deferred reasons?
    resources_list_d = [dict(resource_type=r, path=p) for r, p in zip(resource_types, task.resources)]

    task_manifest = os.path.join(manifest_path)

    with open(task_manifest, 'w+') as f:
        # this version should be the global pbsmrtpipe version
        # or the manifest file spec version?
        runnable_task = RunnableTask(task, cluster_renderer)
        f.write(json.dumps(runnable_task.to_dict(), sort_keys=True, indent=2))

    log.debug("wrote task id {i} to {p}".format(i=tid, p=manifest_path))
    return True


def _log_pbsmrptipe_header():
    s = '''

       _                        _         _
      | |                      | |       (_)
 _ __ | |__  ___ _ __ ___  _ __| |_ _ __  _ _ __   ___
| '_ \| '_ \/ __| '_ ` _ \| '__| __| '_ \| | '_ \ / _ \\
| |_) | |_) \__ \ | | | | | |  | |_| |_) | | |_) |  __/
| .__/|_.__/|___/_| |_| |_|_|   \__| .__/|_| .__/ \___|
| |                                | |     | |
|_|                                |_|     |_|

'''
    return s


def job_resource_create_and_setup_logs(job_root_dir, bg, task_opts, workflow_level_opts, ep_d):
    """
    Create job resource dirs and setup log handlers

    :type job_root_dir: str
    :type bg: BindingsGraph
    :type task_opts: dict
    :type workflow_level_opts: WorkflowLevelOptions
    :type ep_d: dict
    """

    job_resources = to_job_resources_and_create_dirs(job_root_dir)

    pb_log_path = os.path.join(job_resources.logs, 'pbsmrtpipe.log')
    master_log_path = os.path.join(job_resources.logs, "master.log")
    master_log_level = logging.INFO
    stdout_level = logging.INFO
    if workflow_level_opts.debug_mode:
        master_log_level = logging.DEBUG
        stdout_level = logging.DEBUG

    setup_internal_logs(master_log_path, master_log_level, pb_log_path, stdout_level)

    log.info("Starting pbsmrtpipe v{v}".format(v=pbsmrtpipe.get_version()))
    log.info("\n" + _log_pbsmrptipe_header())

    BU.write_binding_graph_images(bg, job_resources.workflow)

    write_entry_points_json(job_resources.entry_points_json, ep_d)

    # Need to map entry points to a FileType and store in the DataStore? or
    # does DataStore only represent outputs?
    smrtpipe_log_df = DataStoreFile(str(uuid.uuid4()), "pbsmrtpipe::pbsmrtpipe.log", FileTypes.LOG.file_type_id, pb_log_path, name="Analysis Log", description="pbsmrtpipe log")
    master_log_df = DataStoreFile(str(uuid.uuid4()), "pbsmrtpipe::master.log", FileTypes.LOG.file_type_id, master_log_path, name="Master Log", description="Master log")
    ds = write_and_initialize_data_store_json(job_resources.datastore_json, [smrtpipe_log_df, master_log_df])
    slog.info("successfully initialized datastore.")

    write_workflow_settings(workflow_level_opts, os.path.join(job_resources.workflow, 'options-workflow.json'))
    if workflow_level_opts.system_message is not None:
        slog.info("Command: {m}".format(m=workflow_level_opts.system_message))

    slog.info("Entry Points:")
    slog.info("\n" + pprint.pformat(ep_d, indent=4))

    slog.info("Workflow Options:")
    slog.info("\n" + pprint.pformat(workflow_level_opts.to_dict(), indent=4))

    slog.info("Task Options:")
    slog.info("\n" + pprint.pformat(task_opts, indent=4))

    task_opts_path = os.path.join(job_resources.workflow, 'options-task.json')
    with open(task_opts_path, 'w') as f:
        f.write(json.dumps(task_opts, sort_keys=True, indent=4))

    env_path = os.path.join(job_resources.workflow, '.env.json')
    IO.write_env_to_json(env_path)
    log.info("wrote current env to {e}".format(e=env_path))

    try:
        sa_system, sa_components = IO.get_smrtanalysis_system_and_components_from_env()
        log.info(sa_system)
        for c in sa_components:
            log.info(c)
    except Exception:
        # black hole exception
        log.warn("unable to determine SMRT Analysis version.")
        pass

    slog.info("completed setting up job directory resources and logs in {r}".format(r=job_root_dir))
    return job_resources, ds, master_log_df


def to_job_resources_and_create_dirs(root_job_dir):
    """
    Create the necessary job directories.


    :rtype: JobResources
    """
    if not os.path.exists(root_job_dir):
        os.mkdir(root_job_dir)

    def _to_make_dir(*args):
        p = os.path.join(*args)
        if not os.path.exists(p):
            os.mkdir(p)
        return p

    f = _to_make_dir

    attr_values = [root_job_dir,
                   f(root_job_dir, 'workflow'),
                   f(root_job_dir, 'html'),
                   f(root_job_dir, 'logs'),
                   f(root_job_dir, 'tasks'),
                   f(root_job_dir, 'html', 'css'),
                   f(root_job_dir, 'html', 'js'),
                   f(root_job_dir, 'html', 'images'),
                   os.path.join(root_job_dir, 'workflow', 'datastore.json'),
                   os.path.join(root_job_dir, 'workflow', 'entry-points.json')]

    return JobResources(*attr_values)


def write_entry_points_json(file_name, ep_d):
    with open(file_name, 'w') as f:
        f.write(json.dumps(ep_d))


def write_and_initialize_data_store_json(file_name, ds_files):
    ds = DataStore(ds_files)
    ds.write_json(file_name)
    return ds


def write_workflow_settings(workflow_options, file_name):
    """

    :type workflow_options: WorkflowLevelOptions

    :param workflow_options:
    :param file_name:
    :return:
    """
    # temp version of this.
    with open(file_name, 'w') as f:
        f.write(json.dumps(workflow_options.to_dict(), sort_keys=True, indent=4))

    return True
