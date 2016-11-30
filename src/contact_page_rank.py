"""Author: Michael Schneider
"""
import datetime
import os
import sys
import networkx as nx
import InputOutput
import numpy
import argparse
import cPickle

options = {}


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="contact_file", help="predicted contacts file in CASP format", required=True)
    parser.add_argument("-l", type=int, dest="length", help="number of residues in the protein", required=True)
    parser.add_argument("-p", dest="pdb_id", help="pdb id + chain id (4+1 letters)")
    parser.add_argument("-f", dest="pdb_file", help="native pdb file. used for reference, not calculation",
                        required=False)
    parser.add_argument("-s", dest="psipred_file", help="sequence and secondary structure file in psipred format",
                        required=True)
    parser.add_argument("-t", type=float, dest="top", help="fraction of top probable contacts to use. 0 < x < 1",
                        required=True)
    parser.add_argument("-a", type=float, dest="alpha", help="dampening parameter alpha", default=0.4, required=False)
    parser.add_argument("-o", dest="out_folder", help="output folder", default=default_output_folder())
    options = parser.parse_args()


def build_ce_graph(xl_data, length, shift_dict, sec_struct):
    # Initialize graph datastructure. The score of the contact will be used as node weights and also the personalization
    # vector will be set to the contact scores
    g = nx.Graph()
    index = 1
    pers = {}
    for score, i in xl_data[:length]:
        g.add_node(index, xl=i, weight=score)
        pers[index] = float(score)
        index += 1
    # Iterate over the nodes (contacts) and get the co-occurance probability matrix for the secondary structure of the
    # centered contact
    for n in g.nodes(data=True):
        sec_lower = sec_struct[n[1]['xl'][0]]
        sec_upper = sec_struct[n[1]['xl'][1]]
        sec_struct_shift_dict = shift_dict[(sec_lower, sec_upper)]
        # Iterate over pairs of nodes
        if sec_struct_shift_dict:
            for o in g.nodes(data=True):
                if o[0] != n[0]:
                    # Compute the relative shift between the contacts in i,j position
                    shift_tuple = (n[1]['xl'][0] - o[1]['xl'][0], n[1]['xl'][1] - o[1]['xl'][1])
                    # Some exception handling
                    if (sec_struct_shift_dict.has_key(shift_tuple) and not
                    numpy.isnan(sec_struct_shift_dict[shift_tuple]) and sec_struct_shift_dict[shift_tuple] != 0.0):
                        # If there is already this edge, keep the edge with the lower weight
                        if g.has_edge(n[0], o[0]):
                            old_weight = g.edge[n[0]][o[0]]['weight']
                            if old_weight > sec_struct_shift_dict[shift_tuple]:
                                g.add_edge(n[0], o[0], weight=sec_struct_shift_dict[shift_tuple])
                        # If the edge does not exist, draw the edge
                        else:
                            g.add_edge(n[0], o[0], weight=sec_struct_shift_dict[shift_tuple])
                        
    return g, pers


def do_page_rank(xl_graph, node_weights, input_alpha, input_len):
    """
    This runs pagerank.
    """
    ranked_nodes = nx.pagerank(xl_graph, alpha=input_alpha, personalization=node_weights, max_iter=100, tol=1e-08)

    for_sorting = [(score, node) for node, score in ranked_nodes.iteritems() if node <= input_len * 999]
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        xl_ranked.append((res_lower, res_upper, score))
    return xl_ranked

# helper functions
def default_output_folder():
    return "../results/" + datetime.datetime.today().date().isoformat() + "/"


def output_file_name():
    output_directory = os.path.abspath(options.out_folder)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    i = 0
    output_file_name = "%s_RRPAR_%s_%s__%s" % (options.pdb_id, options.alpha, options.top, i)
    while os.path.exists(os.path.join(output_directory, output_file_name)):
        i += 1
        output_file_name = "%s_RRPAR_%s_%s__%s" % (options.pdb_id, options.alpha, options.top, i)
    return os.path.join(output_directory, output_file_name)


def return_sorted_tuple(tuple):
    """Return the sorted tuple, such as the lower number is always first"""
    list_tup = list(tuple)
    list_tup.sort()
    tuple = (list_tup[0], list_tup[1])
    return tuple


def main():
    parse_arguments()
    sec_struct = InputOutput.InputOutput.parse_psipred(options.psipred_file)
    shift_dict = cPickle.load(open("probabilities/shifts_sigma_0.05.txt", "rb"))
    xl_data = InputOutput.InputOutput.load_restraints_pr(options.contact_file, seq_sep_min=12)
    xl_graph, node_weights = build_ce_graph(xl_data, int(options.length * options.top), shift_dict, sec_struct)
    xl_ranked = do_page_rank(xl_graph, node_weights, options.alpha, options.length)
    InputOutput.InputOutput.write_contact_file(xl_ranked, output_file_name(), upper_distance=8)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
