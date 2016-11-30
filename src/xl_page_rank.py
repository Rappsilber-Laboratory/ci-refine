"""Author: Michael Schneider
"""
import sys
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")

import networkx as nx
import InputOutput
from optparse import OptionParser

parser = OptionParser()


def add_options(parser):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-e", type="string", dest="example", help="An example")
    parser.add_option("-o", type="int", dest="offset", help="offset")
    parser.add_option("-n", type="string", dest="id", help="id")
    parser.add_option("-m", type="int", dest="max_links", help="max_links")
    options, args = parser.parse_args()
    return options, args

options, args = add_options(parser)


def build_xl_graph(xl_data):
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
                if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=6, double=True): # good results with delta=6
                    g.add_edge(n[0],o[0])
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


def do_page_rank(xl_graph, pers, decoy_dict):
    ranked_nodes = nx.pagerank(xl_graph, max_iter=10000, alpha=0.85, tol=1e-08, personalization=pers)
    for_sorting = [ (score, node) for node, score in ranked_nodes.iteritems()]
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        print score, res_lower, res_upper, decoy_dict[(res_lower, res_upper)]
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
    xl_data, gt_data, decoy_dict = InputOutput.InputOutput.load_xl_data(options.example, options.offset)
    InputOutput.InputOutput.write_contact_file(gt_data, options.id + "_PSM.txt", upper_distance=20, decoy_dict=
    decoy_dict)
    xl_graph, pers = build_xl_graph(xl_data)
    do_page_rank(xl_graph, pers, decoy_dict)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")



