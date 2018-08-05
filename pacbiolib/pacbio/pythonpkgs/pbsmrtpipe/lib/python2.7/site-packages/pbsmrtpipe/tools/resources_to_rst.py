import os
import os.path as op
import sys
import logging
import json
import re

from jinja2 import PackageLoader, Environment

from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner, get_default_argparser_with_base_opts
from pbcommand.validators import validate_dir
import pbsmrtpipe.loader as L

_TEMPLATE_DIR_NAME = 'templates'

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), _TEMPLATE_DIR_NAME)

_LOADER = PackageLoader('pbsmrtpipe.tools', package_path='templates')

ENV = Environment(loader=_LOADER)

log = logging.getLogger(__name__)
slog = logging.getLogger('status.' + __name__)

__version__ = "0.2.0"


def make_rst_table(rows, headers=None):
    """
    Construct RST syntax for a generic table.
    """
    _rows = list(rows)
    if headers is not None:
        assert len(headers) == len(rows[0])
        _rows.append(headers)
    widths = [max([len(_rows[j][i]) for j in range(len(_rows))])
              for i in range(len(_rows[0]))]
    format_str = "| " + \
        " | ".join(["%-{:d}s".format(x) for x in widths]) + " |"
    sep_str = "+" + "+".join(["-" * (x + 2) for x in widths]) + "+"
    table = [sep_str]
    if headers is not None:
        table.append(format_str % tuple(headers))
        table.append(re.sub("-", "=", sep_str))
    for row in rows:
        table.append(format_str % tuple(row))
        table.append(sep_str)
    return "\n".join(table)


def load_pipelines_from_dir(dir_name):
    """
    :arg dir_name: Path to pipeline template dir
    :type dir_name: str
    :rtype: list[Pipeline]
    """
    pipelines = []
    if os.path.exists(dir_name):
        for file_name in os.listdir(dir_name):
            if file_name.endswith(".json"):
                try:
                    pipelines.append(json.load(open(os.path.join(dir_name, file_name))))
                except Exception as e:
                    log.warn("Unable to load Resolved Pipeline Template from {}. {}".format(
                        dir_name, str(e)))
    return pipelines


def sanitize(s):
    # this will create problems in the RST table
    return s.replace("\n", " ")


def task_option_to_rst_table_str(task_options):
    """
    Convert the task options to an RST table
    :rtype: str
    """
    if task_options:
        header = ["Name", "ID", "Default Value", "OptionType", "Description"]
        table = []
        for to in task_options:
            raw_description = to['description']
            # this will generate invalid rst tables if |n are present
            description = sanitize(raw_description)
            option_type = to['optionTypeId']
            row = [to['name'], to['id'], str(to['default']), option_type, description]
            table.append(row)
        rst_table = make_rst_table(table, headers=header)
        return rst_table
    else:
        return ""


def convert_pipeline_to_rst(pipeline_d):
    pipeline_id = pipeline_d['id']
    name = sanitize(pipeline_d['name'])

    # Nat already wrote a RST table converter, so we'll use that
    task_options = pipeline_d['taskOptions']
    task_option_table_str = task_option_to_rst_table_str(task_options)

    _d = dict(pipeline=pipeline_d,
              entry_points=pipeline_d['entryPoints'],
              task_options=pipeline_d['taskOptions'],
              task_table_summary=task_option_table_str)

    t = ENV.get_template("pipeline_details_rst.tmpl")
    s = t.render(**_d)

    return [s, pipeline_id, name]


def generate_index(pipeline_ids, title, doc_version):
    _d = dict(title=title,
              version=doc_version,
              pipeline_ids=pipeline_ids)
    return ENV.get_template("pipeline_index_rst.tmpl").render(**_d)


def _write_file(output_file, content):
    with open(output_file, 'w+') as f:
        f.write(content)
    return content


def write_converted_pipelines(converted_pipelines, doc_output_dir, index_rst="index.rst", title="PacBio Pipelines", doc_version="0.1.0"):

    if not os.path.exists(doc_output_dir):
        os.makedirs(doc_output_dir)

    pipeline_ids = set([])
    for pipeline_str, pipeline_id, pipeline_name in converted_pipelines:
        pipeline_rst = op.join(doc_output_dir, pipeline_id + ".rst")
        _write_file(pipeline_rst, pipeline_str)
        pipeline_ids.add(pipeline_id)

    full_index_rst = os.path.join(doc_output_dir, index_rst)
    index_str = generate_index(pipeline_ids, title, doc_version)
    _write_file(full_index_rst, index_str)

    return 0


def _render_template(tmpl_name, output_file, **kwargs):
    tmpl = ENV.get_template(tmpl_name)
    content = tmpl.render(**kwargs)
    _write_file(output_file, content)
    return content


def convert_pipeline_json_files(output_dir, pipeline_dir, title="PacBio Pipelines", doc_version="0.1.0"):
    pipelines = load_pipelines_from_dir(pipeline_dir)
    converted_pipelines = [convert_pipeline_to_rst(p) for p in pipelines]
    write_converted_pipelines(converted_pipelines, output_dir, title=title)

    _d = dict(title=title, version=doc_version)
    # Write the sphinx resources; conf.py and Makefile
    conf_py_path = os.path.join(output_dir, 'conf.py')
    _render_template("conf_py.tmpl", conf_py_path, **_d)

    makefile_path = os.path.join(output_dir, "Makefile")
    _render_template("makefile.tmpl", makefile_path, **_d)

    return 0


def convert_pipeline_json_files_args(args):

    output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
    # see comments above about IO layer
    # pipelines = L.load_resolved_pipeline_template_jsons_from_dir(args.pipeline_dir)
    return convert_pipeline_json_files(output_dir, args.pipeline_dir, title=args.title, doc_version=args.doc_version)


def get_parser():
    desc = "Generate Pipeline documentation from a directory of Resolved Pipeline Templates"
    p = get_default_argparser_with_base_opts(__version__, desc)

    f = p.add_argument

    f("pipeline_dir", type=validate_dir, help="Path to Pipeline Template JSON Dir")
    f('-o', "--output-dir", default="pipeline-docs", help="Path to RST Output Dir")
    f('-t', '--title', default="PacBio Pipelines", help="Title of Pipeline documents")
    f('-d', '--doc-version', default="0.1.0", help="Version of Pipeline documents")
    return p


def main(argv=sys.argv):
    parser = get_parser()
    return pacbio_args_runner(argv[1:], parser, convert_pipeline_json_files_args, log, setup_log)

if __name__ == "__main__":
    sys.exit(main(argv=sys.argv))
