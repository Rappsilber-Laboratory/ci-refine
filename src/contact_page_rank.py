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
    parser.add_argument("-c", dest="contact_file", help="Predicted contacts file in CASP format", required=True)
    parser.add_argument("-l", type=int, dest="length", help="Number of residues in the protein", required=True)
    parser.add_argument("-p", dest="pdb_id", help="pdb id + chain id (4+1 letters)")
    parser.add_argument("-s", dest="psipred_file", help="sequence and secondary structure file in psipred format",
                        required=True)
    parser.add_argument("-t", type=float, dest="beta", help="fraction of top probable contacts to use. 0 < x < 1",
                        required=True)
    parser.add_argument("-a", type=float, dest="alpha", help="dampening parameter alpha", default=0.4, required=False)
    parser.add_argument("-o", dest="out_folder", help="output folder", default=default_output_folder())
    options = parser.parse_args()


def build_corroborating_information_graph(contact_data, number_of_contacts, shift_probabilities, secondary_structures):
    """
    Builds the corroborating evidence (corroborating information) graph for predicted contacts. See the original
    article for more details on how the graph is build

    Parameters
    ----------
    contact_data : Contact data file
    number_of_contacts : Number of contacts that will be used to build the CE/CI graph
    shift_probabilities : Shift dictionary, containing the co-occurence probability matrices for contacts
    secondary_structures : Dictionary containing the secondary structure assignments

    Returns
    -------
    g : Corroborating information graph
    pers : Personalization vector (containing contact prediction scores)
    """
    # Initialize graph datastructure. The score of the contact will be used as node weights and also the personalization
    # vector will be set to the contact scores
    g = nx.Graph()
    index = 1
    pers = {}
    for score, i in contact_data[:number_of_contacts]:
        g.add_node(index, xl=i, weight=score)
        pers[index] = float(score)
        index += 1
    # Iterate over the nodes (contacts) and get the co-occurance probability matrix for the secondary structure of the
    # centered contact
    for n in g.nodes(data=True):
        sec_lower = secondary_structures[n[1]['xl'][0]]
        sec_upper = secondary_structures[n[1]['xl'][1]]
        sec_struct_shift_dict = shift_probabilities[(sec_lower, sec_upper)]
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


def do_page_rank(corroborating_information_graph, pers, input_alpha):
    """
    Run PageRank on the corroborating information graph 
    
    Parameters
    ----------
    corroborating_information_graph : Here, contacts (nodes) connected by edges (corroborating information)
    pers : Personalization vector 
    input_alpha : Damping factor for PageRank

    Returns
    -------
    re_ranked_contacts : Contacts ranked by PageRank
    """
    ranked_nodes = nx.pagerank(corroborating_information_graph, alpha=input_alpha, personalization=pers,
                               max_iter=100, tol=1e-08)

    for_sorting = [(score, node) for node, score in ranked_nodes.iteritems()]
    for_sorting.sort()
    for_sorting.reverse()

    re_ranked_contacts = []
    for score, n in for_sorting:
        res_lower = corroborating_information_graph.node[n]['xl'][0]
        res_upper = corroborating_information_graph.node[n]['xl'][1]
        re_ranked_contacts.append((res_lower, res_upper, score))
    return re_ranked_contacts


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
    contact_data = InputOutput.InputOutput.load_restraints_pr(options.contact_file, seq_sep_min=12)
    ci_graph, node_weights = build_corroborating_information_graph(contact_data, int(options.length * options.beta),
                                                                   shift_dict, sec_struct)
    xl_ranked = do_page_rank(ci_graph, node_weights, options.alpha)
    InputOutput.InputOutput.write_contact_file(xl_ranked, output_file_name(), upper_distance=8)


if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
