"""
MIT License

Copyright (c) 2016 Michael Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import sys
import networkx as nx
import InputOutput
from optparse import OptionParser


def add_options(parser):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-x", type="string", dest="clms_file", help="Cross-link search result file from Xi")
    parser.add_option("-o", type="int", dest="offset", help="Offset to deal with differences in PDB and UNIPROT sequences")
    parser.add_option("-n", type="string", dest="id", help="id")
    parser.add_option("-m", type="int", dest="max_links", help="max_links")
    options, args = parser.parse_args()
    return options, args

parser = OptionParser()
options, args = add_options(parser)


def build_corroborating_information_graph(xl_data):
    """
    Builds the corroborating information graph for CLMS data. Description can be found in main paper

    Parameters
    ----------
    xl_data : Cross-link data

    Returns
    -------
    g : Corroborating information graph
    pers : Personalization vector (containing cross-link scores)
    """
    g = nx.Graph()
    index = 1
    pers = {}
    # Add cross-links as nodes
    for i, score in xl_data:
        g.add_node(index, xl=i)
        pers[index] = score
        index += 1
    # Add edges between if cross-links are within delta distance
    for n in g.nodes(data=True):
        for o in g.nodes(data=True):
            if o[0] > n[0]:
                if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=6, double=True):
                    g.add_edge(n[0], o[0])
    # Close loops
    for i in xrange(0, 2):
        to_add = []
        for n in g.nodes(data=True):
            for o in g.nodes(data=True):
                if o[0] > n[0]:
                    if has_shared_neighbors(o, n, g, number_of_neighbors=1):
                        to_add.append((o[0], n[0]))
        for i, j in to_add:
            g.add_edge(i, j)

    return g, pers


def do_page_rank(ci_graph, pers, decoy_dict):
    """
    Perform PageRank on CI graph
    Parameters
    ----------
    ci_graph : Corroborating information graph
    pers : Personalization vector
    decoy_dict : Contains which links are decoy. Only used for writing data, not for calculation

    Returns
    -------
    None
    """
    ranked_nodes = nx.pagerank(ci_graph, max_iter=10000, alpha=0.85, tol=1e-08, personalization=pers)
    for_sorting = [ (score, node) for node, score in ranked_nodes.iteritems()]
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        res_lower = ci_graph.node[n]['xl'][0]
        res_upper = ci_graph.node[n]['xl'][1]
        xl_ranked.append((res_lower, res_upper, score))
    InputOutput.InputOutput.write_contact_file(xl_ranked, options.id + "_PR.txt", upper_distance=20, decoy_dict=decoy_dict)


def has_shared_neighbors(node_1, node_2, graph, number_of_neighbors=1):
    loop = False
    shared_nodes = []
    for e in graph.neighbors(node_1[0]):
        for f in graph.neighbors(node_2[0]):
            if e == f:
                shared_nodes.append(e)
    if len(shared_nodes) >= number_of_neighbors:
        loop = True
    return loop


def is_neighbourhood(tuple_1, tuple_2, delta=1, double=True):
    is_nei = False
    t_1 = return_sorted_tuple(tuple_1)
    t_2 = return_sorted_tuple(tuple_2)
    m_1 = 0
    m_2 = 0

    if abs(t_1[0] - t_2[0]) <= delta:
        m_1 = 1
    if abs(t_1[1] - t_2[1]) <= delta:
        m_2 = 1

    if double:
        if m_1 == 1 and m_2 == 1:
            is_nei = True
    else:
        if m_1 == 1 or m_2 == 1:
            is_nei = True
    return is_nei


def return_sorted_tuple(tuple):
    """Return the sorted tuple, such as the lower number is always first"""
    list_tup = list(tuple)
    list_tup.sort()
    tuple = (list_tup[0], list_tup[1])
    return tuple


def main():
    xl_data, gt_data, decoy_dict = InputOutput.InputOutput.load_xl_data(options.clms_file, options.offset)
    InputOutput.InputOutput.write_contact_file(gt_data, options.id + "_PSM.txt", upper_distance=20, decoy_dict=
    decoy_dict)
    xl_graph, pers = build_corroborating_information_graph(xl_data)
    do_page_rank(xl_graph, pers, decoy_dict)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")



