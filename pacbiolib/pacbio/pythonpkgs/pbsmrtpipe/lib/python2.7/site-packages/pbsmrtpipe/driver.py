from collections import deque
import logging
import os
import socket
import time
import sys
import random
import multiprocessing
import json
import Queue
import pprint
import traceback
import types
import functools
import uuid
import platform

from pbcommand.pb_io import (write_resolved_tool_contract,
                             write_tool_contract,
                             load_report_from_json)
from pbcommand.pb_io.tool_contract_io import write_resolved_tool_contract_avro
from pbcommand.utils import log_traceback
from pbcommand.models import (FileTypes, DataStoreFile)
from pbsmrtpipe.utils import nfs_exists_check

from pbcore.io import getDataSetUuid

import pbsmrtpipe
import pbsmrtpipe.constants as GlobalConstants
from pbsmrtpipe.exceptions import (PipelineRuntimeError,
                                   PipelineRuntimeKeyboardInterrupt,
                                   WorkflowBaseException,
                                   MalformedChunkOperatorError)
import pbsmrtpipe.pb_io as IO
import pbsmrtpipe.graph.bgraph as B
import pbsmrtpipe.graph.bgraph_utils as BU
import pbsmrtpipe.cluster as C
import pbsmrtpipe.tools.runner as T
import pbsmrtpipe.report_renderer as R
import pbsmrtpipe.driver_utils as DU
import pbsmrtpipe.services as WS
from pbsmrtpipe import opts_graph as GX

from pbsmrtpipe.graph.models import (TaskStates,
                                     TaskBindingNode,
                                     TaskChunkedBindingNode,
                                     EntryOutBindingFileNode,
                                     TaskScatterBindingNode)


from pbsmrtpipe.models import (Pipeline, ToolContractMetaTask, MetaTask,
                               GlobalRegistry, TaskResult, validate_operator,
                               AnalysisLink, RunnableTask,
                               ScatterToolContractMetaTask,
                               GatherToolContractMetaTask)
from pbsmrtpipe.engine import TaskManifestWorker
from pbsmrtpipe.pb_io import WorkflowLevelOptions


log = logging.getLogger(__name__)
slog = logging.getLogger('status.' + __name__)

# logging.basicConfig(level=logging.DEBUG)


class Constants(object):
    SHUTDOWN = "SHUTDOWN"

def _init_bg(bg, ep_d):
    """Resolving/Initializing BindingGraph with supplied EntryPoints"""

    # Help initialize graph/epoints
    B.resolve_entry_points(bg, ep_d)
    # Update Task-esque EntryBindingPoint
    B.resolve_entry_binding_points(bg)

    # initialization File node attributes
    for eid, path in ep_d.iteritems():
        B.resolve_entry_point(bg, eid, path)
        B.resolve_successor_binding_file_path(bg)

    return bg


def _get_scatterable_task_id(chunk_operators_d, task_id):
    """Get the companion scatterable task id from the original meta task id"""
    for operator_id, chunk_operator in chunk_operators_d.iteritems():
        if task_id == chunk_operator.scatter.task_id:
            return chunk_operator.scatter.scatter_task_id
    return None


def _is_task_id_scatterable(chunk_operators_d, task_id):
    return _get_scatterable_task_id(chunk_operators_d, task_id) is not None


def _get_chunk_operator_by_scatter_task_id(scatter_task_id, chunk_operators_d):
    for operator_id, chunk_operator in chunk_operators_d.iteritems():
        if scatter_task_id == chunk_operator.scatter.scatter_task_id:
            return chunk_operator

    raise KeyError("Unable to find chunk operator for scatter task id {i}".format(i=scatter_task_id))


def _status(bg):
    """
    Create a terse status messages from the overall status of the graph

    :type bg: `pbsmrtpipe.bgraph.BindingsGraph`
    :param bg:
    :return:
    """
    ntasks = len(bg.all_task_type_nodes())
    ncompleted_tasks = len(B.get_tasks_by_state(bg, TaskStates.SUCCESSFUL))
    return "Workflow status {n}/{t} completed/total tasks".format(t=ntasks, n=ncompleted_tasks)


def _get_report_uuid(path):
    """Get UUID from the file resource or return None"""
    return load_report_from_json(path).uuid


def _get_or_create_uuid_from_file(path, file_type):
    """
    Extract the uuid from the DataSet or Report, or assign a new UUID

    :param path: Path to file

    :rtype: str
    :return: uuid string
    """
    if file_type.file_type_id == FileTypes.REPORT.file_type_id:
        return _get_report_uuid(path)
    elif file_type.file_type_id in FileTypes.ALL_DATASET_TYPES():
        return getDataSetUuid(path)
    else:
        return uuid.uuid4()


def _get_last_lines_of_stderr(n, stderr_path):
    """Read in the last N-lines of the stderr from a task

    This should be smarter to look for stacktraces and common errors.
    """
    lines = []
    try:
        nfs_exists_check(stderr_path)
        with open(stderr_path, 'r') as f:
            lines = deque(f, n)
            lines = [l.rstrip() for l in lines]
    except Exception as e:
        log.exception("Unable to extract stderr from {p} Error {e}".format(p=stderr_path, e=e.message))

    return lines


def _terminate_all_workers(workers, shutdown_event):
    if workers:
        nworkers = len(workers)
        log.info("terminating {n} workers.".format(n=nworkers))
        shutdown_event.set()
        # this is brutal terminate
        for worker in workers:
            log.info("terminating worker {n}".format(n=worker.name))
            worker.terminate()
            # time.sleep(0.25)

        log.info("successfully terminated {n} workers.".format(n=nworkers))


def _terminate_worker(worker):
    """
    :type worker: TaskManifestWorker
    """
    tid = worker.task_id
    name = worker.name
    pid_ = worker.pid
    try:
        worker.terminate()
    except Exception as e:
        log.error("Failed to terminate worker {n} task-id:{i} Pid {p}. {c} {e}".format(n=name, i=tid, p=pid_, e=e.message, c=e.__class__))


def _are_workers_alive(workers):
    return all(w.is_alive() for w in workers.values())


def _is_chunked_task_node_type(tnode):
    # Keep Gather Tasks as non-Chunked.
    return isinstance(tnode, (TaskChunkedBindingNode, TaskScatterBindingNode))


def _write_terminate_script(output_dir):

    def __writer(fx, sx):
        with open(fx, 'w') as f:
            f.write(sx)

    pid = os.getpid()
    kill_script = os.path.join(output_dir, GlobalConstants.PBSMRTPIPE_PID_KILL_FILE_SCRIPT)
    #FIXME(nechols)(2016-11-29) this is a bit hacky
    term_file = os.path.join(output_dir, GlobalConstants.TERM_FILE)

    # kill using the process group id to make sure that the signal is sent
    # to the children in a non-tty
    sx = """#!/bin/bash

# MK. This is taken from https://github.com/cloud66/apps-build/blob/master/gitlab/killtree.sh
# It most certainly doesn't work on BSD

# it is able to kill a process provided and all child sub-processes of that process (doesn't affect parents)
# execute with: killtree <pid> <sig>
# ie: killtree 1234 INT
killtree() {
    local _pid=$1
    local _sig=${2:-INT}
    # needed to stop quickly forking parent from producing child between child killing and parent killing
    #kill -stop ${_pid}
    echo "kill -stop ${_pid}"
    for _child in $(ps -o pid --no-headers --ppid ${_pid}); do
        killtree ${_child} ${_sig}
    done
    kill -${_sig} ${_pid}
    echo "kill -${_sig} ${_pid}"
}

"""

    sx += "echo {i} > {f}\nkilltree {i} INT".format(i=pid, f=term_file)

    __writer(kill_script, sx)

    return pid


def __exe_workflow(global_registry, ep_d, bg, task_opts, workflow_opts, output_dir,
                   workers, shutdown_event, service_uri_or_none):
    """
    Core runner of a workflow.

    :type bg: BindingsGraph
    :type workflow_opts: WorkflowLevelOptions
    :type output_dir: str
    :type service_uri_or_none: str | None

    :param workers: {taskid:Worker}
    :return:

    :rtype: bool

    The function is doing way too much. This is really terrible.
    """
    # this used for the cluster submission.
    job_id = random.randint(100000, 999999)
    started_at = time.time()

    m_ = "Distributed" if workflow_opts.distributed_mode is not None else "Local"

    # Setup logger, job directory and initialize DS
    slog.info("creating job resources in {o}".format(o=output_dir))
    job_resources, ds, master_log_ds_file = DU.job_resource_create_and_setup_logs(output_dir, bg, task_opts, workflow_opts, ep_d)
    slog.info("successfully created job resources.")

    slog.info("starting to execute {m} workflow with assigned job_id {i}".format(i=job_id, m=m_))
    slog.info("system {m} {x} nproc:{n}".format(m=platform.system(), n=multiprocessing.cpu_count(), x=platform.node()))
    slog.info("exe'ing workflow Cluster renderer {c}".format(c=global_registry.cluster_renderer))
    slog.info("Service URI: {i} {t} (fqdn = {h})".format(i=service_uri_or_none, t=type(service_uri_or_none), h=socket.getfqdn()))
    slog.info("Max number of Chunks  {n} ".format(n=workflow_opts.max_nchunks))
    slog.info("Max number of nproc   {n}".format(n=workflow_opts.max_nproc))
    slog.info("Max number of workers {n}".format(n=workflow_opts.max_nworkers))
    slog.info("tmp dir               {n}".format(n=workflow_opts.tmp_dir))

    # Some Pre-flight checks
    # Help initialize graph/epoints
    B.resolve_entry_points(bg, ep_d)
    # Update Task-esque EntryBindingPoint
    B.resolve_entry_binding_points(bg)

    # initialization File node attributes
    for eid, path in ep_d.iteritems():
        B.resolve_entry_point(bg, eid, path)
        B.resolve_successor_binding_file_path(bg)

    # Mark Chunkable tasks
    B.label_chunkable_tasks(bg, global_registry.chunk_operators)

    slog.info("validating binding graph")
    # Check the degree of the nodes
    B.validate_binding_graph_integrity(bg)
    slog.info("successfully validated binding graph.")

    # Add scattered
    # This will add new nodes to the graph if necessary
    B.apply_chunk_operator(bg, global_registry.chunk_operators, global_registry.tasks, workflow_opts.max_nchunks)

    log.debug(BU.to_binding_graph_summary(bg))

    slog.info("pbsmrtpipe main process pid={i} pgroupid={g} ppid={p}".format(i=os.getpid(), g=os.getpgrp(), p=os.getppid()))
    _write_terminate_script(output_dir)

    # "global" file type id counter {str: int} that will be
    # used to generate ids
    file_type_id_to_count = {file_type.file_type_id: 0 for _, file_type in global_registry.file_types.iteritems()}

    # Create a closure to allow file types to assign unique id
    # this returns a func
    to_resolve_files_func = B.to_resolve_files(file_type_id_to_count)

    # Local vars
    max_total_nproc = workflow_opts.total_max_nproc
    max_nworkers = workflow_opts.max_nworkers
    max_nproc = workflow_opts.max_nproc
    max_nchunks = workflow_opts.max_nchunks
    tmp_dir = workflow_opts.tmp_dir

    q_out = multiprocessing.Queue()

    worker_sleep_time = 1
    # To store all the reports that are displayed in the analysis.html
    # {id:task-id, report_path:path/to/report.json}
    analysis_file_links = []

    # Flag for pipeline execution failure
    has_failed = False
    # Time to sleep between each step the execution loop
    # after the first 1 minute of exe, update the sleep time to 2 sec
    sleep_time = 1
    # Running total of current number of slots/cpu's used
    total_nproc = 0

    # Define a bunch of util funcs to try to make the main driver while loop
    # more understandable. Not the greatest model.

    def _to_run_time():
        return time.time() - started_at

    def write_analysis_report(analysis_file_links_):
        analysis_report_html = os.path.join(job_resources.html, 'analysis.html')
        R.write_analysis_link_report(analysis_file_links_, analysis_report_html)

    def update_analysis_file_links(task_id_, report_path_):
        analysis_link = AnalysisLink(task_id_, report_path_)
        log.info("updating report analysis file links {a}".format(a=analysis_link))
        analysis_file_links.append(analysis_link)
        write_analysis_report(analysis_file_links)

    # factories for getting a Worker instance
    # utils for getting the running func and worker type
    def _to_worker(w_is_distributed, wid, task_id, manifest_path_):
        # the IO loading will forceful set this to None
        # if the cluster manager not defined or cluster_mode is False
        if global_registry.cluster_renderer is None or not w_is_distributed:
            r_func = T.run_task_manifest
        else:
            r_func = T.run_task_manifest_on_cluster
        return TaskManifestWorker(q_out, shutdown_event, worker_sleep_time,
                                  r_func, task_id, manifest_path_, name=wid)

    # Define a bunch of util funcs to try to make the main driver while loop
    # more understandable. Not the greatest model.
    def write_report_(bg_, current_state_, was_successful_):
        return DU.write_update_main_workflow_report(job_id, job_resources, bg_, current_state_, was_successful_, _to_run_time())

    def write_task_summary_report(bg_):
        task_summary_report = DU.to_task_summary_report(bg_)
        p = os.path.join(job_resources.html, 'task_summary.html')
        R.write_report_to_html(task_summary_report, p)

    def services_log_update_progress(source_id_, level_, message_):
        if service_uri_or_none is not None:
            total_log_uri = "{u}/log".format(u=service_uri_or_none)
            WS.log_pbsmrtpipe_progress(total_log_uri, message_, level_, source_id_, ignore_errors=True)

    def services_add_datastore_file(datastore_file_):
        if service_uri_or_none is not None:
            total_ds_uri = "{u}/datastore".format(u=service_uri_or_none)
            log.debug("Adding datastore file to services {d}".format(d=datastore_file_))
            WS.add_datastore_file(total_ds_uri, datastore_file_, ignore_errors=True)

    def _update_analysis_reports_and_datastore(tnode_, task_):
        assert (len(tnode_.meta_task.output_file_display_names) ==
                len(tnode_.meta_task.output_file_descriptions) ==
                len(tnode_.meta_task.output_types) == len(task_.output_files))
        for i_file, (file_type_, path_, name, description) in enumerate(zip(
                tnode_.meta_task.output_types, task_.output_files,
                tnode_.meta_task.output_file_display_names,
                tnode_.meta_task.output_file_descriptions)):
            source_id = "{t}-out-{i}".format(t=task_.task_id, i=i_file)
            if tnode_.meta_task.datastore_source_id is not None:
                source_id = tnode_.meta_task.datastore_source_id
            ds_uuid = _get_or_create_uuid_from_file(path_, file_type_)
            is_chunked_ = _is_chunked_task_node_type(tnode_)
            ds_file_ = DataStoreFile(ds_uuid, source_id, file_type_.file_type_id, path_, is_chunked=is_chunked_, name=name, description=description)
            ds.add(ds_file_)
            ds.write_update_json(job_resources.datastore_json)

            # Update Services
            services_add_datastore_file(ds_file_)

            dsr = DU.datastore_to_report(ds)
            R.write_report_to_html(dsr, os.path.join(job_resources.html, 'datastore.html'))
            if file_type_ == FileTypes.REPORT:
                T.write_task_report(job_resources, task_.task_id, path_, DU._get_images_in_dir(task_.output_dir))
                update_analysis_file_links(tnode_.idx, path_)

    def _log_task_failure_and_call_services(task_result, task_id_):
        """
        log the error messages extracted from TaskResult
        :type task_result: TaskResult
        """
        mx = "Task {i} {m}".format(i=task_id_, m=task_result.error_message)
        slog.error(mx)
        log.error(mx)
        services_log_update_progress("pbsmrtpipe::{i}".format(i=task_id_), WS.LogLevels.ERROR, mx)

    def has_available_slots(n):
        if max_total_nproc is None:
            return True
        return total_nproc + n <= max_total_nproc

    # Misc setup
    write_report_(bg, TaskStates.CREATED, False)
    # write empty analysis reports
    write_analysis_report(analysis_file_links)

    # Add Master log to the datastore file
    services_add_datastore_file(master_log_ds_file)

    BU.write_binding_graph_images(bg, job_resources.workflow)

    # write initial report.
    DU.write_main_workflow_report(job_id, job_resources, workflow_opts, task_opts, bg, TaskStates.RUNNING, False, 0.0)

    exit_code = None
    # For book-keeping
    # task id -> tnode
    tid_to_tnode = {}
    # tnode -> Task instance
    tnode_to_task = {}

    is_workflow_distributable = global_registry.cluster_renderer is not None
    # local loop for adjusting sleep time, this will get reset after each new
    # task is created
    niterations = 0
    dt_ramp = 0.25
    # number of iterations before switching to steady state sleep
    stead_state_n = 50
    # sleep for 5 sec
    dt_stead_state = 4
    term_file = os.path.join(output_dir, GlobalConstants.TERM_FILE)
    try:
        log.debug("Starting execution loop... in process {p}".format(p=os.getpid()))

        while True:
            # After the initial startup, bump up the time to reduce resource usage
            # (since multiple instances will be launched from the services)
            niterations += 1
            if niterations < stead_state_n:
                sleep_time = dt_ramp
            else:
                sleep_time = dt_stead_state

            # Convert Task -> ScatterAble task (emits a Chunk.json file)
            B.apply_scatterable(bg, global_registry.chunk_operators, global_registry.tasks)

            # This will add new TaskBinding nodes to the graph if necessary
            B.apply_chunk_operator(bg, global_registry.chunk_operators, global_registry.tasks, max_nchunks)
            # B.write_binding_graph_images(bg, job_resources.workflow)
            # If a TaskScatteredBindingNode is completed successfully and
            # output chunk.json is resolved, read in the file and
            # generate the new chunked tasks. This mutates the graph
            # significantly.
            B.add_gather_to_completed_task_chunks(bg, global_registry.chunk_operators, global_registry.tasks, job_resources.tasks)

            if not _are_workers_alive(workers):
                for tix_, w_ in workers.iteritems():
                    if not w_.is_alive():
                        log.warn("Worker {i} (pid {p}) is not alive for task {x}. Worker exit code {e}.".format(i=w_.name, p=w_.pid, e=w_.exitcode, x=tix_))
                        #w_.terminate()

            # Check if Any tasks are running or that there still runnable tasks
            is_completed = bg.is_workflow_complete()
            has_task_running = B.has_running_task(bg)

            if not has_task_running:
                if not is_completed:
                    if not B.has_task_in_states(bg, TaskStates.RUNNABLE_STATES()):
                        if not B.has_next_runnable_task(bg):
                            msg = "Unable to find runnable task or any tasks running and workflow is NOT completed."
                            log.error(msg)
                            log.error(BU.to_binding_graph_summary(bg))
                            services_log_update_progress("pbsmrtpipe", WS.LogLevels.ERROR, msg)
                            raise PipelineRuntimeError(msg)

            time.sleep(sleep_time)

            # log.debug("Sleeping for {s}".format(s=sleep_time))
            log.debug("\n" + BU.to_binding_graph_summary(bg))
            # BU.to_binding_graph_task_summary(bg)

            # This should only be triggered after events. The main reason
            # to keep updating it was the html report is up to date with the
            # runtime
            write_report_(bg, TaskStates.RUNNING, is_completed)

            if is_completed:
                msg_ = "Workflow is completed. breaking out."
                log.info(msg_)
                services_log_update_progress("pbsmrtpipe", WS.LogLevels.INFO, msg_)
                break

            try:
                result = q_out.get_nowait()
            except Queue.Empty:
                result = None

            # log.info("Results {r}".format(r=result))
            if isinstance(result, TaskResult):
                niterations = 0
                log.debug("Task result {r}".format(r=result))

                tid_, state_, msg_, run_time_ = result
                tnode_ = tid_to_tnode[tid_]
                task_ = tnode_to_task[tnode_]

                # Process Successful Task Result
                if state_ == TaskStates.SUCCESSFUL:
                    msg_ = "Task was successful {r}".format(r=result)
                    slog.info(msg_)

                    # this will raise if a task output is failed to be resolved
                    B.validate_outputs_and_update_task_to_success(bg, tnode_, run_time_, bg.node[tnode_]['task'].output_files)
                    slog.info("Successfully validated outputs of {t}".format(t=repr(tnode_)))

                    services_log_update_progress("pbsmrtpipe::{i}".format(i=tid_), WS.LogLevels.INFO, msg_)
                    B.update_task_output_file_nodes(bg, tnode_, tnode_to_task[tnode_])
                    B.resolve_successor_binding_file_path(bg)

                    total_nproc -= task_.nproc
                    w_ = workers.pop(tid_)
                    _terminate_worker(w_)

                    # Update Analysis Reports and Register output files to Datastore
                    _update_analysis_reports_and_datastore(tnode_, task_)

                    # BU.write_binding_graph_images(bg, job_resources.workflow)
                else:
                    # Process Non-Successful Task Result
                    B.update_task_state(bg, tnode_, state_)
                    _log_task_failure_and_call_services(result, tid_)

                    # let the remaining running jobs continue
                    w_ = workers.pop(tid_)
                    _terminate_worker(w_)

                    total_nproc -= task_.nproc
                    has_failed = True

                    # BU.write_binding_graph_images(bg, job_resources.workflow)

                _update_msg = _status(bg)
                log.info(_update_msg)
                slog.info(_update_msg)

                s_ = TaskStates.FAILED if has_failed else TaskStates.RUNNING

                write_report_(bg, s_, False)
                write_task_summary_report(bg)

            elif isinstance(result, types.NoneType):
                pass
            else:
                log.error("Unexpected queue result type {t} {r}".format(t=type(result), r=result))

            if has_failed:
                log.error("job has failed. breaking out.")
                # Just kill everything
                break

            # Computational resources are tapped
            if len(workers) >= max_nworkers:
                # don't do anything
                continue

            tnode = B.get_next_runnable_task(bg)

            if tnode is None:
                continue
            elif isinstance(tnode, TaskBindingNode):
                niterations = 0
                # Found a Runnable Task

                # base task_id-instance_id
                tid = '-'.join([tnode.meta_task.task_id, str(tnode.instance_id)])

                task_dir = os.path.join(job_resources.tasks, tid)
                if not os.path.exists(task_dir):
                    os.mkdir(task_dir)

                to_resources_func = B.to_resolve_di_resources(task_dir, root_tmp_dir=workflow_opts.tmp_dir)
                input_files = B.get_task_input_files(bg, tnode)

                # convert metatask -> task
                try:
                    task = GX.meta_task_to_task(tnode.meta_task, input_files, task_opts, task_dir, max_nproc, max_nchunks,
                                                to_resources_func, to_resolve_files_func)
                except Exception as e:
                    slog.error("Failed to convert metatask {i} to task. {m}".format(i=tnode.meta_task.task_id, m=e.message))
                    raise

                # log.debug(task)

                bg.node[tnode]['nproc'] = task.nproc

                if not has_available_slots(task.nproc):
                    # not enough slots to run in
                    continue

                bg.node[tnode]['task'] = task
                tnode_to_task[tnode] = task

                if isinstance(tnode.meta_task, (ToolContractMetaTask, ScatterToolContractMetaTask, GatherToolContractMetaTask)):
                    # the task.options have actually already been resolved here, but using this other
                    # code path for clarity
                    if isinstance(tnode.meta_task, ToolContractMetaTask):
                        rtc = IO.static_meta_task_to_rtc(tnode.meta_task, task, task_opts, task_dir, tmp_dir, max_nproc, is_distributed=is_workflow_distributable)
                    elif isinstance(tnode.meta_task, ScatterToolContractMetaTask):
                        rtc = IO.static_scatter_meta_task_to_rtc(tnode.meta_task, task, task_opts, task_dir, tmp_dir, max_nproc, max_nchunks, tnode.meta_task.chunk_keys, is_distributed=is_workflow_distributable)
                    elif isinstance(tnode.meta_task, GatherToolContractMetaTask):
                        # this should always be a TaskGatherBindingNode which will have a .chunk_key
                        rtc = IO.static_gather_meta_task_to_rtc(tnode.meta_task, task, task_opts, task_dir, tmp_dir, max_nproc, tnode.chunk_key, is_distributed=is_workflow_distributable)
                    else:
                        raise TypeError("Unsupported task type {t}".format(t=tnode.meta_task))

                    # write driver manifest, which calls the resolved-tool-contract.json
                    # there's too many layers of indirection here. Partly due to the pre-tool-contract era
                    # python defined tasks.
                    # Always write the RTC json for debugging purposes
                    tc_path = os.path.join(task_dir, GlobalConstants.TOOL_CONTRACT_JSON)
                    write_tool_contract(tnode.meta_task.tool_contract, tc_path)

                    rtc_json_path = os.path.join(task_dir, GlobalConstants.RESOLVED_TOOL_CONTRACT_JSON)
                    rtc_avro_path = os.path.join(task_dir, GlobalConstants.RESOLVED_TOOL_CONTRACT_AVRO)
                    if rtc.driver.serialization == 'avro':
                        # hack to fix command
                        task.cmds[0] = task.cmds[0].replace('.json', '.avro')
                        write_resolved_tool_contract_avro(rtc, rtc_avro_path)
                    # for debugging
                    write_resolved_tool_contract(rtc, rtc_json_path)

                runnable_task_path = os.path.join(task_dir, GlobalConstants.RUNNABLE_TASK_JSON)
                runnable_task = RunnableTask(task, global_registry.cluster_renderer)
                runnable_task.write_json(runnable_task_path)

                # Create an instance of Worker
                w = _to_worker(tnode.meta_task.is_distributed, "worker-task-{i}".format(i=tid), tid, runnable_task_path)

                workers[tid] = w
                w.start()
                total_nproc += task.nproc
                slog.info("Starting worker {i} ({n} workers running, {m} total proc in use)".format(i=tid, n=len(workers), m=total_nproc))

                # Submit job to be run.
                B.update_task_state(bg, tnode, TaskStates.SUBMITTED)
                msg_ = "Updating task {t} to SUBMITTED".format(t=tid)
                log.debug(msg_)
                tid_to_tnode[tid] = tnode
                services_log_update_progress("pbsmrtpipe::{i}".format(i=tnode.idx), WS.LogLevels.INFO, msg_)
                # BU.write_binding_graph_images(bg, job_resources.workflow)

            elif isinstance(tnode, EntryOutBindingFileNode):
                # Handle EntryPoint types. This is not a particularly elegant design :(
                bg.node[tnode]['nproc'] = 1
                log.info("Marking task as completed {t}".format(t=tnode))
                B.update_task_state_to_success(bg, tnode, 0.0)
                # Update output paths
                mock_file_index = 0
                for fnode in bg.successors(tnode):
                    file_path = "/path/to/mock-file-{i}.txt".format(i=mock_file_index)
                    B.update_file_state_to_resolved(bg, fnode, file_path)
                    mock_file_index += 1
            else:
                raise TypeError("Unsupported node type {t} of '{x}'".format(t=type(tnode), x=tnode))

            # Update state of any files
            B.resolve_successor_binding_file_path(bg)
            # Final call to write empty analysis reports
            write_analysis_report(analysis_file_links)

        # end of while loop
        _terminate_all_workers(workers.values(), shutdown_event)

        if has_failed:
            log.debug("\n" + BU.to_binding_graph_summary(bg))

        was_successful = B.was_workflow_successful(bg)
        s_ = TaskStates.SUCCESSFUL if was_successful else TaskStates.FAILED
        write_report_(bg, s_, was_successful)
        #FIXME(nechols)(2016-11-30) very hacky workaround
        if was_successful:
            exit_code = GlobalConstants.EXIT_SUCCESS
        elif os.path.exists(term_file):
            raise PipelineRuntimeKeyboardInterrupt("Job was terminated, ignoring child process errors")
        else:
            exit_code = GlobalConstants.EXIT_FAILURE

    except PipelineRuntimeKeyboardInterrupt:
        write_report_(bg, TaskStates.KILLED, False)
        write_task_summary_report(bg)
        BU.write_binding_graph_images(bg, job_resources.workflow)
        exit_code = GlobalConstants.EXIT_TERMINATED

    except Exception as e:
        slog.error("Unexpected error {e}. Writing reports and shutting down".format(e=str(e)))
        # update workflow reports to failed
        write_report_(bg, TaskStates.FAILED, False)
        write_task_summary_report(bg)
        services_log_update_progress("pbsmrtpipe", WS.LogLevels.ERROR, "Error {e}".format(e=e))
        BU.write_binding_graph_images(bg, job_resources.workflow)
        raise

    finally:
        write_task_summary_report(bg)
        BU.write_binding_graph_images(bg, job_resources.workflow)

    return exit_code


def _write_final_results_message(state):
    if state is True:
        log.info(" Successfully completed workflow.")
    else:
        msg = " Failed to run workflow."
        log.error(msg)


def _validate_entry_points_or_raise(entry_points_d):
    for entry_id, path in entry_points_d.iteritems():
        if not os.path.exists(path):
            raise IOError("Unable to find entry point {e} path {p}".format(e=entry_id, p=path))

    return True


def _load_io_for_workflow(registered_tasks, registered_pipelines, workflow_template_xml_or_pipeline,
                          entry_points_d, preset_jsons, preset_xmls, rc_preset_or_none, force_distribute=None, force_chunk_mode=None, debug_mode=None):
    """
    Load and resolve input IO layer

    # Load Presets and Workflow Options. Resolve and Merge
    # The Order of loading is
    # - rc, workflow.xml, then preset.xml
    # force_distribute will attempt to override ALL settings (if cluster_manager is defined)

    :returns: A tuple of Workflow Bindings, Workflow Level Options, Task Opts, ClusterRenderer)
    :rtype: (List[(str, str)], WorkflowLevelOpts, {TaskId:value}, ClusterRenderer)
    """

    # Load Presets and Workflow Options. Resolve and Merge
    # The Order of loading is
    # - rc, workflow.xml, then preset.xml

    # A little sanity check
    # Validate that entry points exist

    slog.info("validating entry points.")
    _validate_entry_points_or_raise(entry_points_d)
    slog.info("successfully validated {n} entry points".format(n=len(entry_points_d)))

    wopts = {}
    topts = {}

    if rc_preset_or_none is None:
        rc_preset = IO.load_preset_from_env()
    else:
        rc_preset = IO.parse_pipeline_preset_xml(rc_preset_or_none)

    if isinstance(workflow_template_xml_or_pipeline, Pipeline):
        # Use default values defined in the Pipeline
        builder_record = IO.BuilderRecord(workflow_template_xml_or_pipeline.all_bindings, workflow_template_xml_or_pipeline.task_options, {})
        slog.info("Loaded pipeline Id {p}".format(p=workflow_template_xml_or_pipeline.pipeline_id))
    else:
        slog.info("Loading workflow template.")
        builder_record = IO.parse_pipeline_template_xml(workflow_template_xml_or_pipeline, registered_pipelines)
        slog.info("successfully loaded workflow template.")

    preset_xml_record = preset_json_record = None
    if preset_jsons:
        slog.info("Loading preset(s) {p}".format(p=preset_jsons))
        preset_json_record = IO.parse_pipeline_preset_jsons(preset_jsons)
        slog.info("successfully loaded preset.")
    else:
        slog.info("No JSON preset provided. Skipping preset json loading.")
    if preset_xmls:
        slog.info("Loading preset(s) {p}".format(p=preset_xmls))
        preset_xml_record = IO.parse_pipeline_preset_xmls(preset_xmls)
        slog.info("successfully loaded preset.")
    else:
        slog.info("No XML preset provided. Skipping preset XML loading.")

    if rc_preset is not None:
        topts.update(dict(rc_preset.task_options))
        wopts.update(dict(rc_preset.workflow_options))

    wopts.update(dict(builder_record.workflow_options))
    topts.update(builder_record.task_options)

    if preset_xml_record is not None:
        wopts.update(dict(preset_xml_record.workflow_options))
        topts.update(dict(preset_xml_record.task_options))
    if preset_json_record is not None:
        wopts.update(dict(preset_json_record.workflow_options))
        topts.update(dict(preset_json_record.task_options))

    workflow_level_opts = IO.WorkflowLevelOptions.from_id_dict(wopts)
    if len(sys.argv) > 0:
        # XXX evil, but this gets sys.argv into pbsmrtpipe.log
        workflow_level_opts.system_message = " ".join(sys.argv)

    # override distributed mode only if provided.
    if isinstance(force_distribute, bool):
        workflow_level_opts.distributed_mode = force_distribute
    workflow_level_opts = IO.validate_or_modify_workflow_level_options(workflow_level_opts)

    slog.info("Successfully validated workflow options.")

    slog.info("validating supplied task options.")
    topts = IO.validate_raw_task_options(registered_tasks, topts)
    slog.info("successfully loaded and validated task options.")

    workflow_bindings = builder_record.bindings

    if isinstance(workflow_level_opts.cluster_manager_path, str):
        cluster_render = C.load_cluster_templates(workflow_level_opts.cluster_manager_path)
    else:
        cluster_render = None

    if isinstance(force_chunk_mode, bool):
        workflow_level_opts.chunk_mode = force_chunk_mode

    workflow_level_opts.max_nchunks = min(workflow_level_opts.max_nchunks, GlobalConstants.MAX_NCHUNKS)

    if workflow_level_opts.distributed_mode is False:
        total_max_nproc = multiprocessing.cpu_count() if workflow_level_opts.total_max_nproc is None else workflow_level_opts.total_max_nproc
        workflow_level_opts.total_max_nproc = min(total_max_nproc, multiprocessing.cpu_count())
        workflow_level_opts.max_nproc = min(workflow_level_opts.max_nproc, workflow_level_opts.total_max_nproc)
        slog.info("local-only mode updating       MAX NPROC to {x}".format(x=workflow_level_opts.max_nproc))
        slog.info("local-only mode updating TOTAL MAX NPROC to {x}".format(x=workflow_level_opts.total_max_nproc))

    if debug_mode is True:
        slog.info("overriding debug-mode to True")
        workflow_level_opts.debug_mode = debug_mode

    log.debug("Resolved workflow level options to {d}".format(d=workflow_level_opts))
    log.debug("\n" + pprint.pformat(workflow_level_opts.to_dict(), indent=4))
    log.debug("Initial resolving of loaded preset.xml and pipeline.xml task options:")
    log.debug("\n" + pprint.pformat(topts))

    return workflow_bindings, workflow_level_opts, topts, cluster_render


def _load_io_for_task(registered_tasks, entry_points_d, preset_jsons, preset_xmls, rc_preset_or_none, force_distribute=None, force_chunk_mode=None, debug_mode=None):
    """Grungy loading of the IO and resolving values

    Returns a tuple of (WorkflowLevelOptions, TaskOptions, ClusterRender)
    """
    slog.info("validating entry points. {e}".format(e=entry_points_d))
    _validate_entry_points_or_raise(entry_points_d)
    slog.info("successfully validated {n} entry points".format(n=len(entry_points_d)))

    wopts = {}
    topts = {}

    if rc_preset_or_none is None:
        rc_preset = IO.load_preset_from_env()
    else:
        rc_preset = IO.parse_pipeline_preset_xml(rc_preset_or_none)

    if rc_preset:
        topts.update(dict(rc_preset.task_options))
        wopts.update(dict(rc_preset.workflow_options))

    if preset_xmls:
        preset_record = IO.parse_pipeline_preset_xmls(preset_xmls)
        wopts.update(dict(preset_record.workflow_options))
        topts.update(dict(preset_record.task_options))

    if preset_jsons:
        preset_record = IO.parse_pipeline_preset_jsons(preset_jsons)
        wopts.update(dict(preset_record.workflow_options))
        topts.update(dict(preset_record.task_options))

    workflow_level_opts = IO.WorkflowLevelOptions.from_id_dict(wopts)

    workflow_level_opts = IO.validate_or_modify_workflow_level_options(workflow_level_opts)

    if isinstance(force_chunk_mode, bool):
        workflow_level_opts.chunk_mode = force_chunk_mode

    # Validate
    topts = IO.validate_raw_task_options(registered_tasks, topts)

    log.debug("Resolved task options to {d}".format(d=workflow_level_opts))
    log.debug(pprint.pprint(workflow_level_opts.to_dict(), indent=4))

    if isinstance(workflow_level_opts.cluster_manager_path, str):
        cluster_render = C.load_cluster_templates(workflow_level_opts.cluster_manager_path)
        # override distributed mode
        if isinstance(force_distribute, bool):
            workflow_level_opts.distributed_mode = force_distribute
    else:
        cluster_render = None

    workflow_level_opts.max_nchunks = min(workflow_level_opts.max_nchunks, GlobalConstants.MAX_NCHUNKS)

    if workflow_level_opts.distributed_mode is False:
        slog.info("local-only mode detected setting total NPROC to {x}".format(x=multiprocessing.cpu_count()))
        workflow_level_opts.total_max_nproc = multiprocessing.cpu_count()

    if debug_mode is True:
        slog.info("overriding debug-mode to True")
        workflow_level_opts.debug_mode = debug_mode

    return workflow_level_opts, topts, cluster_render


def exe_workflow(global_registry, entry_points_d, bg, task_opts, workflow_level_opts, output_dir, service_uri):
    """This is the fundamental entry point to running a pbsmrtpipe workflow.

    :rtype: int
    """

    slog.info("Initializing Workflow")

    # Create workers container here so we can catch exceptions and shutdown
    # gracefully
    workers = {}
    manager = multiprocessing.Manager()
    shutdown_event = manager.Event()

    state = False
    try:
        exit_code = __exe_workflow(global_registry, entry_points_d, bg, task_opts,
                               workflow_level_opts, output_dir,
                               workers, shutdown_event, service_uri)
    except Exception as e:
        if isinstance(e, KeyboardInterrupt):
            emsg = "received SIGINT. Attempting to abort gracefully."
        else:
            emsg = "Unexpected exception. shutting down."

        log.exception(emsg)
        exit_code = GlobalConstants.EXCEPTION_TO_EXIT_CODE.get(e.__class__, GlobalConstants.DEFAULT_EXIT_CODE)
        state = False

    finally:
        # Each of these calls needs to be wrapped in a try black hole
        if workers:
            _terminate_all_workers(workers.values(), shutdown_event)
        _write_final_results_message(state)

    return exit_code


def workflow_exception_exitcode_handler(func):
    """Decorator to call core workflow/task run funcs

    Funcs should return positive int exit code

    It will log the run time and handle exception handling and logging/stderr.

    :rtype: int
    """

    @functools.wraps(func)
    def _wrapper(*args, **kwargs):
        started_at = time.time()

        exit_code = GlobalConstants.DEFAULT_EXIT_CODE
        try:
            exit_code = func(*args, **kwargs)
        except Exception as e:
            emsg = "Error executing function {f} with {c}".format(f=func.__name__, c=e.__class__)
            log.exception(emsg)
            slog.error(e)

            type_, value_, traceback_ = sys.exc_info()

            log_traceback(slog, e, traceback_)
            log_traceback(log, e, traceback_)

            exit_code = GlobalConstants.EXCEPTION_TO_EXIT_CODE.get(e.__class__, GlobalConstants.DEFAULT_EXIT_CODE)
            print "Exit code from {e} class {c}".format(e=exit_code, c=e.__class__)

        finally:
            print "Shutting down."
            run_time = time.time() - started_at
            run_time_min = run_time / 60.0
            _m = {
                0: "was Successful",
                7: "was Terminated"
            }.get(exit_code, "Failed")
            _d = dict(s=_m, r=run_time, x=pbsmrtpipe.get_version(), m=run_time_min, c=exit_code)
            msg = "Completed execution pbsmrtpipe v{x}. Workflow {s} in {r:.2f} sec ({m:.2f} min) with exit code {c}".format(**_d)

            slog.info(msg)
            log.info(msg)

        return exit_code

    return _wrapper


@workflow_exception_exitcode_handler
def run_pipeline(registered_pipelines_d, registered_file_types_d, registered_tasks_d,
                 chunk_operators, workflow_template_xml_or_pipeline, entry_points_d,
                 output_dir, preset_jsons, preset_xmls, rc_preset_or_none, service_uri,
                 force_distribute=None, force_chunk_mode=None, debug_mode=None):
    """
    Entry point for running a pipeline

    :param workflow_template_xml_or_pipeline: path to workflow xml or Pipeline instance
    :param entry_points_d:
    :param output_dir:
    :param preset_jsons: list of path to preset json
    :param preset_xmls: list of path to preset xml
    :return: exit code

    :type registered_tasks_d: dict[str, pbsmrtpipe.pb_tasks.core.MetaTask]
    :type registered_file_types_d: dict[str, pbsmrtpipe.pb_tasks.core.FileType]
    :type workflow_template_xml_or_pipeline: str
    :type output_dir: str
    :type preset_jsons: list[str]
    :type preset_xmls: list[str]
    :type service_uri: str | None
    :type force_distribute: None | bool

    :rtype: int
    """
    log.debug(pprint.pformat(entry_points_d))

    workflow_bindings, workflow_level_opts, task_opts, cluster_render = _load_io_for_workflow(registered_tasks_d,
                                                                                              registered_pipelines_d,
                                                                                              workflow_template_xml_or_pipeline,
                                                                                              entry_points_d, preset_jsons, preset_xmls,
                                                                                              rc_preset_or_none,
                                                                                              force_distribute=force_distribute,
                                                                                              force_chunk_mode=force_chunk_mode,
                                                                                              debug_mode=debug_mode)

    slog.info("building graph")
    bg = B.binding_strs_to_binding_graph(registered_tasks_d, workflow_bindings)
    slog.info("successfully loaded graph from bindings.")

    valid_chunk_operators = {}
    # Disabled chunk operators if necessary
    if workflow_level_opts.chunk_mode is False:
        slog.info("Chunk mode is False. Disabling {n} chunk operators.".format(n=len(chunk_operators)))
    else:
        # Validate chunk operators, or skip if malformed.
        for chunk_operator_id, chunk_operator in chunk_operators.iteritems():
            try:
                validate_operator(chunk_operator, registered_tasks_d)
                valid_chunk_operators[chunk_operator_id] = chunk_operator
            except MalformedChunkOperatorError as e:
                log.warn("Invalid chunk operator {i}. {m}".format(i=chunk_operator_id, m=e.message))

    filtered_chunk_operators_d = _filter_chunk_operators(bg, valid_chunk_operators)
    # Container to hold all the resources
    global_registry = GlobalRegistry(registered_tasks_d,
                                     registered_file_types_d,
                                     filtered_chunk_operators_d,
                                     cluster_render)

    return exe_workflow(global_registry, entry_points_d, bg, task_opts,
                        workflow_level_opts, output_dir, service_uri)


def _filter_chunk_operators(bg, chunk_operators_d):
    """Filters all the chunk operators that don't have tasks in the BindingsGraph

    Returns a chunk operator dict
    """
    operators_d = {}
    task_ids = {tnode.meta_task.task_id for tnode in bg.task_nodes() if isinstance(tnode, TaskBindingNode)}
    for op_id, op, in chunk_operators_d.iteritems():
        if op.scatter.task_id in task_ids:
            operators_d[op_id] = op

    dn = len(chunk_operators_d) - len(operators_d)
    if dn != 0:
        log.warn("Filtered {n} chunk operators from registry.".format(n=dn))

    return operators_d


def _task_to_entry_point_ids(meta_task):
    """Generate entry points from a meta-task. $entry:e_0, $entry:e_1, ...

    This is used to automatically create pipeline entry points from the
    positional inputs of the task.
    """
    def _to_e(i_):
        return "{e}e_{i}".format(i=i_, e=GlobalConstants.ENTRY_PREFIX)
    return [_to_e(i) for i in xrange(len(meta_task.input_types))]


def _task_to_binding_strings(meta_task):
    """
    Create binding strings from a meta task instance. Entry points are create
    to correspond with the positional input type.

    e_0, e_1, ...

    :type meta_task: MetaTask
    :param meta_task:
    :return: List of binding strings
    """
    def _to_b(i_):
        return "{t}:{i}".format(t=meta_task.task_id, i=i_)

    entry_point_ids = _task_to_entry_point_ids(meta_task)

    return [(entry_point_id, _to_b(i)) for i, entry_point_id in enumerate(entry_point_ids)]


def _validate_task_entry_points_or_raise(meta_task, entry_points_d):
    """Validate the entry points are consistent with the MetaTask

    :raises: KeyError
    :rtype: bool
    """
    entry_point_ids = _task_to_entry_point_ids(meta_task)

    zs = zip(entry_point_ids, meta_task.input_types)
    for ep_id, in_type in zs:
        if ep_id not in entry_points_d:
            _d = dict(n=ep_id, i=in_type, v=entry_points_d.keys(), t=meta_task.task_id)
            raise KeyError("Task {t} Required Entry point '{n}' for input type {i} not supplied. Supplied values {v}".format(**_d))

    return True


@workflow_exception_exitcode_handler
def run_single_task(registered_file_types_d, registered_tasks_d, chunk_operators,
                    entry_points_d, task_id, output_dir, preset_jsons, preset_xmls, rc_preset_or_none,
                    service_config,
                    force_distribute=None,
                    force_chunk_mode=None,
                    debug_mode=None):
    """
    Entry Point for running a single task

    :param task_id:
    :param output_dir:
    :return:
    """

    print entry_points_d
    meta_task = registered_tasks_d.get(task_id, None)

    if meta_task is None:
        raise KeyError("Unable to find task id '{i}' in registered tasks. Use "
                       "'show-tasks' to get a list of registered tasks.".format(i=task_id))

    workflow_level_opts, task_opts, cluster_render = _load_io_for_task(registered_tasks_d, entry_points_d,
                                                                       preset_jsons, preset_xmls, rc_preset_or_none,
                                                                       force_distribute=force_distribute,
                                                                       force_chunk_mode=force_chunk_mode,
                                                                       debug_mode=debug_mode)

    slog.info("building bindings graph")
    binding_str = _task_to_binding_strings(meta_task)

    bg = B.binding_strs_to_binding_graph(registered_tasks_d, binding_str)
    slog.info("successfully bindings graph for task {i}".format(i=task_id))

    # Validate chunk operators
    valid_chunk_operators = {k: v for k, v in chunk_operators.iteritems() if validate_operator(v, registered_tasks_d)}
    filtered_chunk_operators_d = _filter_chunk_operators(bg, valid_chunk_operators)
    # Container to hold all the resources
    global_registry = GlobalRegistry(registered_tasks_d,
                                     registered_file_types_d,
                                     filtered_chunk_operators_d,
                                     cluster_render)

    return exe_workflow(global_registry, entry_points_d, bg, task_opts,
                        workflow_level_opts, output_dir, service_config)
