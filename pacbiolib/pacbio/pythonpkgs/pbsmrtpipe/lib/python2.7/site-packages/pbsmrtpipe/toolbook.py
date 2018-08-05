import operator
import itertools
import tempfile
import StringIO
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from IPython.display import display_svg

from pbsmrtpipe.external_tools import dot_file_to_svg, dot_file_to_png


def _to_tmp_file(suffix):
    t = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    t.close()
    return t.name


def display_dot(dot_file):
    svg_file = _to_tmp_file('.svg')
    # print "write dot file to {f}".format(f=dot_file)
    dot_file_to_svg(dot_file, svg_file)
    with open(svg_file, 'r') as f:
        s = f.read()
    display_svg(s, raw=True)


def networkx_graph_to_dot_str(g):
    f = StringIO.StringIO()
    nx.write_dot(g, f)
    contents = f.getvalue()
    f.close()
    return contents


def networkx_graph_to_png(g, image_path):
    s = networkx_graph_to_dot_str(g)
    x = _to_tmp_file(".dot")
    with open(x, 'w') as f:
        f.write(s)
    dot_file_to_png(x, image_path)


def display_dot_str(dot_str):
    f = _to_tmp_file('.dot')
    with open(f, 'w') as x:
        x.write(dot_str)
    display_dot(f)


def display_networkx_graph(g):
    dot_file = _to_tmp_file('.dot')
    write_dot(g, dot_file)
    display_dot(dot_file)
