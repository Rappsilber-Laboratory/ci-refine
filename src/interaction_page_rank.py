"""Author: Michael Schneider
"""
import datetime
import os
import sys
import networkx as nx
import InputOutput
import numpy as np
import argparse
import cPickle
#import matplotlib.pyplot as plt
import pandas as pd
import itertools

options = {}


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="interaction_file", help="predicted contacts file in CASP format", required=False)
    parser.add_argument("-a", type=float, dest="alpha", help="dampening parameter alpha", default=0.85, required=False)
    options = parser.parse_args()


def get_go_term_column(df, go_dict, col_string):
    """
    Parameters
    ----------
    df
    go_dict
    col_string

    Returns
    -------
    """
    for index, row in df.iterrows():
        # go = row['Gene ontology (cellular component)']
        go = row[col_string]
        if isinstance(go, str):
            go_list = go.split(";")
        else:
            go_list = []

        if row['Entry'] not in go_dict:
            go_dict[row['Entry']] = []

        #if go_list is None and len(go_dict[row['Entry']]) == 0:
        #    go_dict[row['Entry']] = []# go_list
        #elif
        #else:
        go_dict[row['Entry']] = list(set(list(itertools.chain(go_dict[row['Entry']], go_list))))

            #print go_dict[row['Entry']]
            #if row['Entry'] not in go_dict:
        #    go_dict[row['Entry']] = go_list
        #elif row['Entry'] in go_dict and go_list is not None:
        #    print go_listrow['Entry']
        #
        #elif row['Entry'] in go_dict and go_list is None:
        #    go_dict[row['Entry']] = go_list


def load_go_terms(uniprot_kb_file):
    go_dict = {}
    df = pd.read_table(uniprot_kb_file)
    get_go_term_column(df, go_dict, 'Gene ontology (cellular component)')
    #get_go_term_column(df, go_dict, 'Gene ontology (biological process)')
    return go_dict

def go_interaction_score(uniprot_1, uniprot_2, go_dict):
    if uniprot_1 not in go_dict or uniprot_2 not in go_dict:
        return 0.0
    if len(go_dict[uniprot_1]) is 0 or len(go_dict[uniprot_2]) is 0:
        return 0.0

    else:
        go_score = 0
        go_term_list_1 = go_dict[uniprot_1]
        go_term_list_2 = go_dict[uniprot_2]
        for go_terms_1 in go_term_list_1:
            for go_terms_2 in go_term_list_2:
                if go_terms_1 == go_terms_2:
                    go_score += 1
        return float(go_score) / max(len(go_term_list_1), len(go_term_list_2))


def max_go_interaction(protein_tuple_1, protein_tuple_2, go_dict):
    go_interaction_scores = [go_interaction_score(protein_tuple_1[0], protein_tuple_2[0], go_dict),
                             go_interaction_score(protein_tuple_1[1], protein_tuple_2[0], go_dict),
                             go_interaction_score(protein_tuple_1[0], protein_tuple_2[1], go_dict),
                             go_interaction_score(protein_tuple_1[1], protein_tuple_2[1], go_dict)]
    return np.mean(go_interaction_scores)


def share_interaction(protein_tuple_1, protein_tuple_2):

    if protein_tuple_1[0] in protein_tuple_2 or protein_tuple_1[1] in protein_tuple_2:
    #if any(i in protein_tuple_1[0] for i in protein_tuple_2[0]) or \
    #   any(i in protein_tuple_1[0] for i in protein_tuple_2[1]) or \
    #   any(i in protein_tuple_1[1] for i in protein_tuple_2[0]) or \
    #   any(i in protein_tuple_1[1] for i in protein_tuple_2[1]):
        #print protein_tuple_1, protein_tuple_2, "true"
        #sys.exit()
    #if any(i in protein_tuple_1 for i in protein_tuple_2):
        return True
    else:
        #print protein_tuple_1, protein_tuple_2, "false"
        return False


def build_ce_graph(interactions, go_dict):
    # Initialize graph datastructure. The score of the contact will be used as node weights and also the personalization
    # vector will be set to the interaction scores
    g = nx.Graph()
    index = 1
    pers = {}
    edges = []
    for score, i, j in interactions:
        #print len(i), joutput_interactions_pagerank_eucleadian
       # sys.exit()
        g.add_node(index, xl=(i, j), weight=score)
        pers[index] = float(score)
        index += 1
    # Iterate over the nodes (interactions) and draw an edge if two interactions share a protein
    index = 0
    for n in g.nodes(data=True):
        for o in g.nodes(data=True):
            if o[0] > n[0]:
                # Draw edge if an interaction shares a protein
                protein_tuple_n = (n[1]['xl'][0], n[1]['xl'][1])
                protein_tuple_o = (o[1]['xl'][0], o[1]['xl'][1])
                #if share_interaction(protein_tuple_n, protein_tuple_o):
                go_score = max_go_interaction(protein_tuple_n, protein_tuple_o, go_dict)
                if go_score >= 0.5:
                    g.add_edge(n[0], o[0])#, weight=go_score)
                    #edges.append((n[0], o[0]))
        #print index / float(len(interactions))
        index += 1
    #g.add_edges_from(edges)
    print len(g.edges())
    return g, pers


def load_interactions(infile):
    file = open(infile)
    interactions = []
    for line in file:
        strline = str(line).strip().split()
        interactions.append((1-((abs(float(strline[0]) - float(strline[2]))) / 37.950100000000006), strline[1], strline[3]))
    file.close()
    interactions.sort(reverse=True)
    for i in interactions:
        print i
    sys.exit()
    return interactions


def do_page_rank(xl_graph, node_weights, input_alpha, input_len):
    """
    This runs pagerank.
    """
    ranked_nodes = nx.pagerank(xl_graph, alpha=input_alpha, personalization=node_weights, max_iter=1000, tol=1e-08)

    for_sorting = [(score, node) for node, score in ranked_nodes.iteritems() if node <= input_len * 999]
    for_sorting.sort(reverse=True)
    xl_ranked = []
    for score, n in for_sorting:
        #print score, n
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        xl_ranked.append((score, res_lower, res_upper))
    return xl_ranked


def main():
    parse_arguments()
    with open('interaction_data_exp_3_eu_6_final.pkl', 'rb') as infile:
        interactions = cPickle.load(infile)
    #interactions = load_interactions("../data/protein_interaction_data/interactions_and_centers.txt")
    go_dict = load_go_terms("uniprot_data_experiment_3_cell.txt")
    #print go_dict
    print("Building CI Graph")
    print("Number of interactions:", len(interactions))
    interactions.sort(reverse=True)
    xl_graph, node_weights = build_ce_graph(interactions, go_dict)
    print("Applying PageRank")
    xl_ranked = do_page_rank(xl_graph, node_weights, options.alpha, 99999999999)

    with open('interaction_data_exp_3_eu_6_final_pr.pkl', 'wb') as outfile:
        cPickle.dump(xl_ranked, outfile, cPickle.HIGHEST_PROTOCOL)
    #print xl_ranked

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
