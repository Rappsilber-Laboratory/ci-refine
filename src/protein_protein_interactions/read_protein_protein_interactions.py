import pandas as pd
from matplotlib.pylab import plt
import cPickle
import seaborn as sns
from sklearn.metrics import roc_curve, precision_recall_curve
import sys
import numpy as np
import itertools


def set_plot_parameters():
    """
    Set parameters of matplotlib plotting for font style and size

    Returns
    -------
    Nothing
    """
    f_size = 8
    plt.rcParams.update({'font.size': f_size,
                           #'font.name':'Arial',
                        'xtick.major.size':2,
                        'xtick.major.width':0.1,
                        'ytick.major.size':2,
                        'ytick.major.width':0.1,
                        'xtick.minor.size':0,
                        'xtick.minor.width':0.0,
                        'ytick.minor.size':0,
                        'ytick.minor.width':0.0})


def tresh_at_cutoff(fpr, thresholds, cutoff=0.53):
    """
    Return the score threshold at a fixed FPR limit

    Parameters
    ----------
    fpr : List of input FPRs
    thresholds : List of score thresholds
    fpr_limit : Maximum FPR value at which to return the threshold

    Returns
    -------
    fp, tresh (float) : FPR value, score threshold
    """
    for fp, thresh in zip(fpr, thresholds):
        if fp >= cutoff:
            return fp, thresh


            #return fp, thresh


def interactions_at_tresh(data, tresh):
    """
    Return the interactions above a given threshold

    Parameters
    ----------
    data : Interaction Data
    tresh : Threshold

    Returns
    -------
    count, interactions_above_tresh : Number of interactions and the list of interactions.

    """
    count = 0
    interactions_above_tresh = []
    for d in data:
        if d[0] >= tresh:
            interactions_above_tresh.append((d[0], d[1], d[2]))
            count += 1
    return count, interactions_above_tresh


def sorted_uniprot_tuple(protein_1, protein_2):
    """
    Sort Uniprot tuple so that the ordering is consistent.

    Parameters
    ----------
    protein_1 : Uniprot string 1
    protein_2: Uniprot string 2

    Returns
    -------
    sorted_tuple: Sorted tuple containing the two Uniprot IDs in alphanumerical order
    """
    sort_list = [protein_1, protein_2]
    sort_list.sort()
    sorted_tuple = (sort_list[0], sort_list[1])
    return sorted_tuple


def filter_interactions(interactions, protein_list):
    """

    Parameters
    ----------
    interactions
    core_complex_list

    Returns
    -------

    """
    filtered_interactions = []

    for score, protein_1, protein_2 in interactions:
        if protein_1 in protein_list and protein_2 in protein_list:
            filtered_interactions.append((score, protein_1, protein_2))

    """
        for complex in core_complex_list:
            if protein_1 in complex:
                protein_1_in_corum = True
            if protein_2 in complex:
                protein_2_in_corum = True

        if protein_1_in_corum and protein_2_in_corum:

    """
    return filtered_interactions


def get_interactions(corom_complex_list):
    interactions = {}
    for complex in corom_complex_list:
        for i in range(0, len(complex)):
            for j in range(i, len(complex)):
                interactions[sorted_uniprot_tuple(complex[i], complex[j])] = 1
    return interactions


def get_proteins(corum_complex_list):
    proteins = {}
    for complex in corum_complex_list:
        for protein in complex:
            proteins[protein] = 1
    return proteins


def get_non_interactions(interactions, corum_interactions):
    non_interactions = {}
    proteins = []
    for score, protein_1, protein_2 in interactions:
        proteins.append(protein_1)
        proteins.append(protein_2)
    proteins = list(set(proteins))
    for i in range(0, len(proteins)):
        for j in range(i, len(proteins)):
            uniprot_tuple = sorted_uniprot_tuple(proteins[i], proteins[j])
            if uniprot_tuple not in corum_interactions:
                non_interactions[uniprot_tuple] = 1
    return non_interactions


def get_protein_intersection(interactions_1, interactions_2, level='protein'):
    filtered_interactions = []
    protein_list = []

    interaction_list = []

    for score, p1, p2 in interactions_2:
        protein_list.append(p1)
        protein_list.append(p2)
    protein_list = list(set(protein_list))

    for score, p1, p2 in interactions_2:
        interaction_list.append((p1, p2))
    interaction_list = list(set(interaction_list))

    for score, p1, p2 in interactions_1:
        if level is 'protein': # Filtering on protein level
            if p1 in protein_list and p2 in protein_list:
                filtered_interactions.append((score, p1, p2))
        elif level is 'interactions':
            if (p1, p2) in interaction_list:
                filtered_interactions.append((score, p1, p2))
    return filtered_interactions


def get_interaction_scoring(interactions, interaction_list):

    y_true = []
    y_score = []
    for score, protein_1, protein_2 in interactions:
        protein_tuple = sorted_uniprot_tuple(protein_1, protein_2)

        if protein_tuple in interaction_list:
            y_true.append(1)
            y_score.append(score)
        else:
            y_true.append(0)
            y_score.append(score)

    return y_true, y_score


def read_corum_complexes(corum_csv_file):
    """
    Read in corum database from CSV file

    Parameters
    ----------
    corum_csv_file : Corum CSV file

    Returns
    -------
    corum_complex_list : List of corum complexes
    """
    core_complexes = pd.read_csv(corum_csv_file, sep=';')
    core_complexes.fillna(0.0, inplace=True)
    core_complex_list = []

    for index, row in core_complexes.iterrows():
        if row["subunits (UniProt IDs)"] is not 0.0:
            split_line = row["subunits (UniProt IDs)"].split(',')
            cleaned_line = [i.strip('(').strip(')') for i in split_line]
            core_complex_list.append(cleaned_line)

    return core_complex_list


def load_interaction_data(interaction_data):

    with open(interaction_data, 'rb') as infile:
        interactions = cPickle.load(infile)

    interactions.sort(reverse=True)

    return interactions


def interactions_at_cutoff(interaction_data_1, protein_list, core_interactions, cutoff=0.55):

    filtered_interactions = filter_interactions(interaction_data_1, protein_list)
    y_true, y_score = get_interaction_scoring(filtered_interactions, core_interactions)
    fpr, tpr, treshs = precision_recall_curve(y_true, y_score)
    precision, threshold = tresh_at_cutoff(fpr, treshs, cutoff=cutoff)  # Define cutoff
    num_interactions, interactions_above_tresh = interactions_at_tresh(interaction_data_1, threshold)  # Exp 3 # Target precision 0.55

    return num_interactions, interactions_above_tresh


def get_precision(filtered_interactions, interactions_corum):
    y_true, y_score = get_interaction_scoring(filtered_interactions, interactions_corum)
    hits = [y for y in y_true if y == 1]
    return round(float(len(hits)) / float(len(y_true)), 2)


def main():
    # Load interaction data
    core_complex_list = read_corum_complexes("../../data/protein_interaction_data/allComplexes.csv")
    core_interactions = get_interactions(core_complex_list)
    protein_list = get_proteins(core_complex_list)
    experiment_list_orig = ['../interaction_data_exp_1_eu_6_final.pkl',
                             '../interaction_data_exp_2_eu_6_final.pkl',
                             '../interaction_data_exp_3_eu_6_final.pkl']

    experiment_list_pr = ['../interaction_data_exp_1_eu_6_final_pr.pkl',
                           '../interaction_data_exp_2_eu_6_final_pr.pkl',
                           '../interaction_data_exp_3_eu_6_final_pr.pkl']

    cutoffs = [0.55]
    for precision_cutoff in cutoffs:
        all_interactions = []
        all_interactions_pr = []
        for experiment_orig, experiment_pr in zip(experiment_list_orig, experiment_list_pr):

            interactions = load_interaction_data(experiment_orig)
            interactions_pr = load_interaction_data(experiment_pr)

            num_interactions, interactions_at_cut = interactions_at_cutoff(interactions,
                                                                           protein_list, core_interactions,
                                                                           cutoff=precision_cutoff)

            num_interactions_pr, interactions_at_cut_pr = interactions_at_cutoff(interactions_pr,
                                                                                 protein_list,
                                                                                 core_interactions,
                                                                                 cutoff=precision_cutoff)

            for i in interactions_at_cut:
                all_interactions.append(sorted_uniprot_tuple(i[1], i[2]))

            for i in interactions_at_cut_pr:
                all_interactions_pr.append(sorted_uniprot_tuple(i[1], i[2]))

        all_interactions = list(set(all_interactions))
        all_interactions_with_score = [(0, i[0], i[1]) for i in all_interactions]

        all_interactions_pr = list(set(all_interactions_pr))
        all_interactions_with_score_pr = [(0, i[0], i[1]) for i in all_interactions_pr]

        filtered_interactions = filter_interactions(all_interactions_with_score, protein_list)
        filtered_interactions_pr = filter_interactions(all_interactions_with_score_pr, protein_list)

        filtered_interactions = get_protein_intersection(filtered_interactions, filtered_interactions_pr, level='protein')
        filtered_interactions_pr = get_protein_intersection(filtered_interactions_pr, filtered_interactions, level='protein')

        print "Orig:", get_precision(filtered_interactions, core_interactions), len(all_interactions_with_score)
        print "PR:", get_precision(filtered_interactions_pr, core_interactions), len(all_interactions_with_score_pr)

#main()
