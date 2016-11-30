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


def evaluate_ranking(xl_ranked, native_dict, num_top_contacts, seq_sep = 24):
    prediction_labels = []
    true_labels = []
    long_range_contacts = []

    for i, j, score in xl_ranked:
        if abs(int(j)-int(i)) >= seq_sep:
            long_range_contacts.append((i, j, score))
    for i, j, score in long_range_contacts[:num_top_contacts]:
        if native_dict.has_key((i, j)):
            prediction_labels.append(1.0)
            true_labels.append(1)
        else:
            prediction_labels.append(1.0)
            true_labels.append(0)
    try:
        #print true_labels, prediction_labels
        #out_score = roc_auc_score(true_labels, prediction_labels)
        out_score = precision_score(true_labels, prediction_labels)

    except:
        out_score = 0
    return out_score

def run_pagerank(xl_data, sec_struct, alpha, beta, shift_dict, input_length):
    xl_graph, node_weights = build_ce_graph(xl_data, int(int(input_length) * beta), shift_dict, sec_struct)
    xl_ranked = do_page_rank(xl_graph, node_weights, alpha, int(input_length))
    return xl_ranked

def run_pagerank_classifier(xl_data, sec_struct, alpha, beta, shift_dict, input_length, contact_dict):
    xl_graph, node_weights = build_ce_graph(xl_data, int(int(input_length) * beta), shift_dict, sec_struct, contact_dict)
    xl_ranked = do_page_rank(xl_graph, node_weights, alpha, int(input_length))
    return xl_ranked


def load_files(pdb_id, sec_struct_folder=None, pred_folder=None):

    psipred_file = "".join([sec_struct_folder,
                            pdb_id,
                            ".horiz"])
    contact_file = "".join([pred_folder,
                            pdb_id,
                            ".prediction"])
    sec_struct = InputOutput.InputOutput.parse_psipred(psipred_file)

    xl_data = InputOutput.InputOutput.load_restraints_pr(contact_file, seq_sep_min=12)
    return xl_data, sec_struct


def load_files_bur(pdb_id, sec_struct_folder=None, pred_folder=None, solv_folder=None):

    psipred_file = "".join([sec_struct_folder,
                            pdb_id,
                            ".horiz"])
    contact_file = "".join([pred_folder,
                            pdb_id,
                            ".prediction"])
    solv_file = "".join([solv_folder,
                              pdb_id,
                              ".solv"])

    sec_struct = InputOutput.InputOutput.parse_psipred(psipred_file)
    bur_dict = InputOutput.InputOutput.parse_psipred(solv_file)
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
    #shift_dict = cPickle.load(open("shifts_moep.p", "rb"))
    #shift_dict = cPickle.load(open("shifts_test_double.p", "rb"))
    #shift_dict = cPickle.load(open("probabilities/shifts_predictor.p", "rb"))
    #shift_dict = cPickle.load(open("shifts_test_metapsicov.p", "rb"))
    #alpha_range = numpy.linspace(0.2, 0.65, 10)
    #beta_range = numpy.linspace(1.5, 3.0, 4)
    #Test:
    #alpha_range = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65]
    #beta_range = [2.0, 2.5]
    #beta_range = [1.5]
    alpha_range = [0.6]
    beta_range = [2.5]
    print alpha_range
    print beta_range
    #for keys, values in shift_dict.iteritems():
    #    print keys, values
    all_values = {}
    for alpha in alpha_range:
        for beta in beta_range:
            all_values[(alpha, beta)] = []
    counter = 1
    for pdb_id, length in protein_data:
        contacts, sec_struct = load_files(pdb_id, sec_struct_folder=options.secondary_structure_folder,
                                                  pred_folder=options.prediction_folder)
        for alpha in alpha_range:
            for beta in beta_range:
                xl_ranked = run_pagerank(contacts, sec_struct, alpha, beta, shift_dict, length)
                native_map = load_native_map(pdb_id)
                acc = evaluate_ranking(xl_ranked, native_map, int(float(length) * 0.2))
                all_values[(alpha, beta)].append(acc)
    
        print "Best parameters after %s proteins" % counter
        best_acc = 0
        best_alpha = 0
        best_beta = 0
        for parameters, values in all_values.iteritems():
            if numpy.mean(values) > best_acc:
                best_alpha = parameters[0]
                best_beta = parameters[1]
                best_acc = numpy.mean(values)
    
        print "Accuracy", round(best_acc, 3)
        print "alpha", best_alpha
        print "beta", best_beta
        counter += 1
            
    write_results_csv(all_values, "%s_parameter_tuning.csv"%options.method_name)
if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
