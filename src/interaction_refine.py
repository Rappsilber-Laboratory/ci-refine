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
import numpy as np
import argparse
import cPickle
import pandas as pd
import itertools

options = {}


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="interaction_file", help="File containing the interactions", required=True)
    parser.add_argument("-g", dest="go", help="Go-term file", required=True)
    parser.add_argument("-o", dest="out_file", help="Go-term file", required=True)
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
        go = row[col_string]
        if isinstance(go, str):
            go_list = go.split(";")
        else:
            go_list = []

        if row['Entry'] not in go_dict:
            go_dict[row['Entry']] = []

        go_dict[row['Entry']] = list(set(list(itertools.chain(go_dict[row['Entry']], go_list))))


def build_ci_graph(interactions, go_dict):
    # Initialize graph datastructure. The score of the contact will be used as node weights and also the personalization
    # vector will be set to the interaction scores
    """
    Build CI graph for protein-protein interactions

    Parameters
    ----------
    interactions : Protein-Protein Interactions
    go_dict : Dictionary containing GO-Term assignments

    Returns
    -------
    g : Corroborating information graph
    pers : Personalization vector (containing interaction scores)
    """
    g = nx.Graph()
    index = 1
    pers = {}
    for score, i, j in interactions:
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
                go_score = max_go_interaction(protein_tuple_n, protein_tuple_o, go_dict)
                if go_score >= 0.5:
                    g.add_edge(n[0], o[0])
        index += 1
    print len(g.edges())
    return g, pers


def do_page_rank(ci_graph, node_weights, input_alpha=0.85):
    """
    Perform PageRank on CI graph
    Parameters
    ----------
    ci_graph : Corroborating information graph
    pers : Personalization vector

    Returns
    -------
    None
    """
    ranked_nodes = nx.pagerank(ci_graph, alpha=input_alpha, personalization=node_weights, max_iter=1000, tol=1e-08)

    for_sorting = [(score, node) for node, score in ranked_nodes.iteritems()]
    for_sorting.sort(reverse=True)
    xl_ranked = []
    for score, n in for_sorting:
        res_lower = ci_graph.node[n]['xl'][0]
        res_upper = ci_graph.node[n]['xl'][1]
        xl_ranked.append((score, res_lower, res_upper))
    return xl_ranked


def load_go_terms(uniprot_kb_file):
    go_dict = {}
    df = pd.read_table(uniprot_kb_file)
    get_go_term_column(df, go_dict, 'Gene ontology (cellular component)')
    #get_go_term_column(df, go_dict, 'Gene ontology (biological process)')
    return go_dict


def max_go_interaction(protein_tuple_1, protein_tuple_2, go_dict):
    go_interaction_scores = [go_interaction_score(protein_tuple_1[0], protein_tuple_2[0], go_dict),
                             go_interaction_score(protein_tuple_1[1], protein_tuple_2[0], go_dict),
                             go_interaction_score(protein_tuple_1[0], protein_tuple_2[1], go_dict),
                             go_interaction_score(protein_tuple_1[1], protein_tuple_2[1], go_dict)]
    return np.mean(go_interaction_scores)


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
    return interactions

def main():
    parse_arguments()
    with open(options.interaction_file, 'rb') as infile:
        interactions = cPickle.load(infile)

    go_dict = load_go_terms(options.go)
    print("Building CI Graph")
    print("Number of interactions:", len(interactions))

    interactions.sort(reverse=True)
    xl_graph, node_weights = build_ci_graph(interactions, go_dict)
    print("Applying PageRank")
    xl_ranked = do_page_rank(xl_graph, node_weights, options.alpha)

    with open(options.out_file, 'wb') as outfile:
        cPickle.dump(xl_ranked, outfile, cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
