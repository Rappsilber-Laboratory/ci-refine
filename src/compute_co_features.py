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
    contact_dict[
    for score, seq_tuple in contacts:
        print score, seq_tuple


def neighborhood_vector(contacts, i, j, shift_dict, sec_struct):
    sec_struct_dict = shift_dict[(sec_struct[i], sec_struct[j])]
    contacts_to_dict(contacts)

#    print sec_struct_dict
    shift_mat = shift_matrix()
#    for x in shift_mat:
#       print sec_struct_dict[x]


def co_occurence_feature(contacts, sec_struct, shift_dict):
    shift_mat = shift_matrix()
    contacts_to_dict(contacts)
    #for u_lower, i_upper, score in contacts:
    #    for x in shift_mat:
            #print x


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


def load_files(pdb_id):
    psipred_file = "".join([options.secondary_structure_folder,
                            pdb_id,
                            ".horiz"])
    contact_file = "".join([options.prediction_folder,
                            pdb_id,
                            ".prediction"])
    sec_struct = InputOutput.InputOutput.parse_psipred(psipred_file)
    xl_data = InputOutput.InputOutput.load_restraints_pr(contact_file, seq_sep_min=12)
    return xl_data, sec_struct



def write_results_csv(value_dict, out_file):
    file = open(out_file, "w")
    file.write("alpha,beta,acc\n")
    for parameters, values in value_dict.iteritems():
        file.write("%.3f,%.3f,%.3f\n"%(parameters[0], parameters[1], numpy.mean(values)))
    file.close()


def main():
    
    parse_arguments()
    
    protein_data = load_pdb_ids(options.contact_file_list)
    shift_dict = cPickle.load(open("probabilities/shifts_sigma_0.05.txt", "rb"))
    #print shift_dict
    for pdb_id, length in protein_data[0:2]:
        native_map = load_native_map(pdb_id)
        contacts, sec_struct = load_files(pdb_id)
        neighborhood_vector(contacts, 56, 79, shift_dict, sec_struct)        
#        print contacts


if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
