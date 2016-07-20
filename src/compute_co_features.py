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
sys.path.append("experiments/pagerank_additional_test_sets/")
from contact_page_rank import build_ce_graph, do_page_rank
options = {}
from sklearn.metrics import accuracy_score, precision_score
from sklearn.metrics import roc_auc_score
from evaluation_helper import load_pdb_ids
import itertools
import random

def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="contact_file_list", help="PDB ID list", required=True)
    parser.add_argument("-p", dest="prediction_folder", required=True)
    parser.add_argument("-s", dest="secondary_structure_folder", required=True)
    parser.add_argument("-m", dest="method_name", required=True)
    options = parser.parse_args()

def load_native_map(pdb_id):
    native_dict = {}
    file = open("../data/native_maps/%s.map"%pdb_id)
    for line in file:
        strline = str(line).strip().split()
        native_dict[(int(strline[0]), int(strline[1]))] = float(strline[-1])
    file.close()
    return native_dict


def shift_matrix():
    matrix = []
    for i in xrange(-8, 9):
        row = []
        for j in xrange(-8, 9):
            if i == 0 and j == 0:
                pass
            else:
                row.append((i, j))
                matrix.append((i, j))
    return matrix


def contacts_to_dict(contacts):
    contact_dict = {}
    for score, seq_tuple in contacts:
        contact_dict[seq_tuple] = score
    return contact_dict


def neighborhood_vector(contacts, i, j, shift_dict, sec_struct):

    sec_struct_dict = shift_dict[(sec_struct[i], sec_struct[j])]
    contact_dict = contacts_to_dict(contacts)

    shift_mat = shift_matrix()
    nej_vector = []
    for x in shift_mat:
        #print x, sec_struct_dict[x]
        if contact_dict.has_key((i-x[0], j-x[1])):
            nej_vector.append(contact_dict[(i-x[0], j-x[1])])
        else:
            nej_vector.append(0)
    return nej_vector
#       print sec_struct_dict[x]

def sec_struct_encoding(sec_struct):
    if sec_struct == "H":
        return [1, 0 ,0]
    if sec_struct == "E":
        return [0, 1, 0]
    if sec_struct == "C":
        return [0, 0, 1]

def co_occurence_feature_with_label(contacts, sec_struct, shift_dict, feature_matrix, labels, native_dict):
    num_contacts = len(contacts)
    # This does not compute the features for shift dict range?
    for i in xrange(0, num_contacts):
        for j in xrange(0, num_contacts):
            if i != j:
                shift_tuple = (contacts[i][1][0] - contacts[j][1][0], contacts[i][1][1] - contacts[j][1][1])
                sec_lower = sec_struct[contacts[i][1][0]]
                sec_upper = sec_struct[contacts[i][1][1]]
                sec_lower_j = sec_struct[contacts[j][1][0]]
                sec_upper_j = sec_struct[contacts[j][1][0]]
                sec_struct_shift_dict = shift_dict[(sec_lower, sec_upper)]
                if shift_tuple in sec_struct_shift_dict:
                    #,
                    vector_i = neighborhood_vector(contacts, contacts[i][1][0], contacts[i][1][1], shift_dict, sec_struct)
                    vector_j = neighborhood_vector(contacts, contacts[j][1][0], contacts[j][1][1], shift_dict, sec_struct)
                    feature_vector = list(itertools.chain(vector_i, vector_j))
                                      #sec_struct_encoding(sec_lower)[0],
                                      #sec_struct_encoding(sec_lower)[1],
                                      #sec_struct_encoding(sec_lower)[2],
                                      #sec_struct_encoding(sec_upper)[0],
                                      #sec_struct_encoding(sec_upper)[1],
                                      #sec_struct_encoding(sec_upper)[2]]  # ,
                      # ,
                    #,
                                      #sec_struct_encoding(sec_lower)[0],
                                      #sec_struct_encoding(sec_lower)[1],
                                      #sec_struct_encoding(sec_lower)[2],
                                      #sec_struct_encoding(sec_upper)[0],
                                      #sec_struct_encoding(sec_upper)[1],
                                      #sec_struct_encoding(sec_upper)[2]]
                                      #abs(contacts[i][1][0]-contacts[i][1][1])]
                                      #abs(contacts[j][1][0] - contacts[j][1][1]),
                                      #sec_struct_encoding(sec_lower_j)[0],
                                      #sec_struct_encoding(sec_lower_j)[1],
                                      #sec_struct_encoding(sec_lower_j)[2],
                                      #sec_struct_encoding(sec_upper_j)[0],
                                      #sec_struct_encoding(sec_upper_j)[1],
                                      #sec_struct_encoding(sec_upper_j)[2]]
                    #vector_i = neighborhood_vector(contacts, contacts[i][1][0], contacts[i][1][1], shift_dict, sec_struct)
                    #vector_j = neighborhood_vector(contacts, contacts[j][1][0], contacts[j][1][1], shift_dict, sec_struct)
                    if native_dict.has_key(contacts[i][1]) and native_dict.has_key(contacts[j][1]):
                        labels.append(1)
                    else:
                        labels.append(0)
                    #feature_vector = list(itertools.chain(vector_i, vector_j))
                    feature_matrix.append(feature_vector)


def co_occurence_features(contacts, sec_struct, shift_dict, feature_matrix, contact_pairs):
    num_contacts = len(contacts)
    for i in xrange(0, num_contacts):
        for j in xrange(0, num_contacts):
            if i != j:
                shift_tuple = (contacts[i][1][0] - contacts[j][1][0], contacts[i][1][1] - contacts[j][1][1])
                sec_lower = sec_struct[contacts[i][1][0]]
                sec_upper = sec_struct[contacts[i][1][1]]
                sec_lower_j = sec_struct[contacts[j][1][0]]
                sec_upper_j = sec_struct[contacts[j][1][0]]
                sec_struct_shift_dict = shift_dict[(sec_lower, sec_upper)]
                if shift_tuple in sec_struct_shift_dict:

                    vector_i = neighborhood_vector(contacts, contacts[i][1][0], contacts[i][1][1], shift_dict, sec_struct)
                    vector_j = neighborhood_vector(contacts, contacts[j][1][0], contacts[j][1][1], shift_dict, sec_struct)
                    feature_vector = list(itertools.chain(vector_i, vector_j))
                    feature_matrix.append(feature_vector)
                    contact_pairs.append((contacts[i][1],contacts[j][1]))


def co_occurance_features_easy(contacts, sec_struct, shift_dict, feature_matrix, contact_pairs):
    num_contacts = len(contacts)
    # This does not compute the features for shift dict range?
    for i in xrange(0, num_contacts):
        for j in xrange(0, num_contacts):
            if i != j:
                shift_tuple = (contacts[i][1][0] - contacts[j][1][0], contacts[i][1][1] - contacts[j][1][1])
                sec_lower = sec_struct[contacts[i][1][0]]
                sec_upper = sec_struct[contacts[i][1][1]]
                sec_lower_j = sec_struct[contacts[j][1][0]]
                sec_upper_j = sec_struct[contacts[j][1][0]]
                sec_struct_shift_dict = shift_dict[(sec_lower, sec_upper)]
                if shift_tuple in sec_struct_shift_dict:
                    #contacts[i][0], contacts[j][0]
                    feature_vector = [contacts[i][0], contacts[j][0], shift_tuple[0], shift_tuple[1],
                                      sec_struct_shift_dict[shift_tuple],
                                      i,
                                      j,
                                      abs(i-j)] # ,
                                      #sec_struct_encoding(sec_lower)[0],
                                      #sec_struct_encoding(sec_lower)[1],
                                      #sec_struct_encoding(sec_lower)[2],
                                      #sec_struct_encoding(sec_upper)[0],
                                      #sec_struct_encoding(sec_upper)[1],
                                      #sec_struct_encoding(sec_upper)[2]]  # ,
                                      #abs(contacts[i][1][0] - contacts[i][1][1])]
                                      #abs(contacts[j][1][0] - contacts[j][1][1]),
                                      #sec_struct_encoding(sec_lower_j)[0],
                                      #sec_struct_encoding(sec_lower_j)[1],
                                      #sec_struct_encoding(sec_lower_j)[2],
                                      #sec_struct_encoding(sec_upper_j)[0],
                                      #sec_struct_encoding(sec_upper_j)[1],
                                      #sec_struct_encoding(sec_upper_j)[2]]
                                      #abs(contacts[i][1][0] - contacts[i][1][1]),
                                      #abs(contacts[j][1][0] - contacts[j][1][1])]
                                      #sec_struct_encoding(sec_lower)[0],
                                      #sec_struct_encoding(sec_lower)[1],
                                      #sec_struct_encoding(sec_lower)[2],
                                      #sec_struct_encoding(sec_upper)[0],
                                      #sec_struct_encoding(sec_upper)[1],
                                      #sec_struct_encoding(sec_upper)[2],
                                      #sec_struct_encoding(sec_lower_j)[0],
                                      #sec_struct_encoding(sec_lower_j)[1],
                                      #sec_struct_encoding(sec_lower_j)[2],
                                      #sec_struct_encoding(sec_upper_j)[0],
                                      #sec_struct_encoding(sec_upper_j)[1],
                                      #sec_struct_encoding(sec_upper_j)[2]]
                    # vector_i = neighborhood_vector(contacts, contacts[i][1][0], contacts[i][1][1], shift_dict, sec_struct)
                    # vector_j = neighborhood_vector(contacts, contacts[j][1][0], contacts[j][1][1], shift_dict, sec_struct)
                    # feature_vector = list(itertools.chain(vector_i, vector_j))
                    feature_matrix.append(feature_vector)
                    contact_pairs.append((contacts[i][1], contacts[j][1]))



def evaluate_ranking(xl_ranked, native_dict, num_top_contacts, seq_sep = 24):
    prediction_labels = []
    true_labels = []
    long_range_contacts = []

    for i, j, score in xl_ranked:
        if abs(int(j)-int(i)) >= seq_sep:
            long_range_contacts.append((i, j, score))
    for i, j, score in long_range_contacts[:num_top_contacts]:
        if native_dict.has_key((i, j)):
            prediction_labels.append(1)
            true_labels.append(1)
        else:
            prediction_labels.append(1)
            true_labels.append(0)
    try:
        out_score = precision_score(true_labels, prediction_labels)
    except: 
        out_score = 0
    return out_score


def load_files(pdb_id, length):
    psipred_file = "".join([options.secondary_structure_folder,
                            pdb_id,
                            ".horiz"])
    contact_file = "".join([options.prediction_folder,
                            pdb_id,
                            ".prediction"])
    sec_struct = InputOutput.InputOutput.parse_psipred(psipred_file)
    xl_data = InputOutput.InputOutput.load_restraints_pr(contact_file, seq_sep_min=12,
                                                         max_contacts=int(int(length)*0.2))
    return xl_data, sec_struct



def write_results_csv(value_dict, out_file):
    file = open(out_file, "w")
    file.write("alpha,beta,acc\n")
    for parameters, values in value_dict.iteritems():
        file.write("%.3f,%.3f,%.3f\n"%(parameters[0], parameters[1], numpy.mean(values)))
    file.close()


def main():
    folder = "/scratch/schneider/projects/pagerank_refinement/data/co-occurence_features_extensive/"
    parse_arguments()
    feature_matrix = []
    labels = []
    protein_data = load_pdb_ids(options.contact_file_list)
    shift_dict = cPickle.load(open("probabilities/shifts_sigma_0.05.txt", "rb"))
    for pdb_id, length in protein_data:
        native_map = load_native_map(pdb_id)
        contacts, sec_struct = load_files(pdb_id, length)
        co_occurence_feature_with_label(contacts, sec_struct, shift_dict, feature_matrix, labels, native_map)
        numpy.savetxt("%s%s_co_occurance_labels.dat"%(folder, pdb_id), labels)
        import scipy.sparse

        x = scipy.sparse.csr_matrix(feature_matrix)
        with open('%s%s_sparse_array_features.dat'%(folder, pdb_id), 'wb') as outfile:
            cPickle.dump(x, outfile, cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
