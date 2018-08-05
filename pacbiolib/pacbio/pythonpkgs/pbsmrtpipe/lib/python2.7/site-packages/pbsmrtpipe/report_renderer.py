from collections import namedtuple
import os
import random
import logging

from jinja2 import Environment, PackageLoader
import shutil

import pbsmrtpipe


log = logging.getLogger(__name__)


_TEMPLATE_DIR_NAME = 'html_templates'

TEMPLATE_DIR = os.path.join(os.path.dirname(pbsmrtpipe.__file__), _TEMPLATE_DIR_NAME)

_LOADER = PackageLoader('pbsmrtpipe', package_path='html_templates')

ENV = Environment(loader=_LOADER)


def table_to_jinja(table):
    """
    Convert to jinja-template

    Zip/Transpose the columns/data so easier to iterate in jinja template

    :type Table

    :param table:
    :return:
    """
    nvalues = len(table.columns[0].values)

    headers = [c.id for c in table.columns]
    items = []
    for c in table.columns:
        items.append([c.values[i] for i in xrange(nvalues)])

    rows = zip(*items)

    _d = dict(table_id=table.id, headers=headers, rows=rows)
    return _d


def render_report(report):
    """
    General rendering a pbreport Report model to a string


    :type report: Report
    :param report:

    :rtype: str
    :return:
    """
    import pbcommand.models.report as RM

    template = ENV.get_template("report.html")

    rattrs_d = [a.to_dict() for a in report.attributes if report.attributes]
    plot_groups_d = [pg.to_dict() for pg in report.plotGroups if report.plotGroups]
    tables_d = [table_to_jinja(t) for t in report.tables if report.tables]

    _d = dict(title="Report {i}".format(i=report.id),
              report_id=report.id,
              rattrs=rattrs_d,
              plot_groups=plot_groups_d,
              tables=tables_d)

    return template.render(**_d)


def _get_js_css_root_dir():
    import pbsmrtpipe
    d = os.path.dirname(pbsmrtpipe.__file__)
    root_dir = os.path.join(d, "html_templates")
    return root_dir


def copy_js_css(root_resource_dir, root_output_dir):

    resource_dir_names = ('css', 'js')
    for resource_dir_name in resource_dir_names:
        x = os.path.join(root_resource_dir, resource_dir_name)
        resource_output_dir = os.path.join(root_output_dir, resource_dir_name)
        for f in os.listdir(x):
            p = os.path.join(x, f)
            r = os.path.join(resource_output_dir, f)
            if not os.path.exists(r):
                shutil.copy(p, r)

    _d = dict(r=root_output_dir, o=root_output_dir)
    # log.debug("completed copying files from {r} to {o}.".format(**_d))


def _write_str_report(s, path, mode):
    with open(path, mode) as f:
        f.write(s)


def write_report_to_html(report, output_file):
    """Render the report and write the file"""
    s = render_report(report)
    _write_str_report(s, output_file, 'w')
    return 0


def write_report_with_html_extras(report, output_file, extras_dir):
    """Write the css/js/html to output file directory"""
    d = _get_js_css_root_dir()

    _ = write_report_to_html(report, output_file)
    copy_js_css(d, extras_dir)
    return 0


def render_analysis_link_report(analysis_report_links):
    # links [{name: "Display name", path: "relative_path"}, ...]
    t = ENV.get_template("analysis.html")
    return t.render(file_links=analysis_report_links)


def write_analysis_link_report(analysis_report_links, output_file):
    s = render_analysis_link_report(analysis_report_links)
    with open(output_file, 'w+') as f:
        f.write(s)
