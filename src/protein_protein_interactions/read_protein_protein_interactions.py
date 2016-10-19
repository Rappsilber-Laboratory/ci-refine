import pandas as pd
import re
import numpy as np
from sklearn import metrics
from matplotlib.pylab import plt
import cPickle
import seaborn as sns
import sys


def set_plot_parameters():
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

set_plot_parameters()
sns.set(style="white", palette="muted")


def tresh_at_fpr(fpr, thresholds, fpr_limit = 0.7):
    for fp, thresh in zip(fpr, thresholds):
        if fp >= 0.0067:
            return fp, thresh
    #return fpr[-1], thresholds[-1]


def interactions_at_tresh(data, tresh):
    count = 0
    interactions_above_tresh = []
    for d in data:
        if d[0] >= tresh:
            interactions_above_tresh.append((d[0], d[1], d[2]))
            #y_true_above_tresh.append(y)
            count += 1
    return count, interactions_above_tresh


def sorted_uniprot_tuple(protein_1, protein_2):
    sort_list = [protein_1, protein_2]
    sort_list.sort()
    sorted_tuple = (sort_list[0], sort_list[1])
    return sorted_tuple

def get_interaction_scoring_list(interactions, core_complex_list):
    y_true = []
    scores = []
    for score, protein_1, protein_2 in interactions:
        # print interaction_1, interaction_2
        is_in_core_complex = False
        for complex in core_complex_list:
            # test overlaps
            if any(i in protein_1 for i in complex) and any(i in protein_2 for i in complex):
            #if interaction_1 in complex and interaction_2 in complex:
                # print interaction_1, interaction_2, complex
                #print protein_1, protein_2, complex
                is_in_core_complex = True
                break
            #else:
            #    is_in_core_complex = False
        if is_in_core_complex:
            y_true.append(1)
            scores.append(score)
        else:
            y_true.append(0)
            scores.append(score)

    return y_true, scores


def filter_interactions(interactions, core_complex_list):
    filtered_interactions = []
    for score, protein_1, protein_2 in interactions:
        # print interaction_1, interaction_2
        protein_1_in_corum = False
        protein_2_in_corum = False

        for complex in core_complex_list:
            if protein_1 in complex or protein_1 in complex:
            #if any(i in protein_1 for i in complex):
                protein_1_in_corum = True
            #if any(i in protein_2 for i in complex):
                protein_2_in_corum = True

        if protein_1_in_corum and protein_2_in_corum:
            filtered_interactions.append((score, protein_1, protein_2))
    return filtered_interactions


def get_interactions(corom_complex_list):
    interactions = {}
    for complex in corom_complex_list:
        for i in range(0, len(complex)):
            for j in range(i, len(complex)):
                interactions[sorted_uniprot_tuple(complex[i], complex[j])] = 1
    return interactions


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


def get_interaction_scoring(interactions, interaction_list, non_interaction_list):

    y_true = []
    scores = []

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    tpr = []
    fpr = []
    tresholds = []
    for score, protein_1, protein_2 in interactions:
        is_in_core_complex = False
        protein_tuple = sorted_uniprot_tuple(protein_1, protein_2)

        # True positives
        if protein_tuple in interaction_list:
            tp += 1
        # False positives
        elif protein_tuple not in interaction_list:
            fp += 1

        tn = len(non_interaction_list) - (tp+fp)

        fn = len(interaction_list) - tp

        fpr.append(float(fp)/float((fp+tn)))
        tpr.append(float(tp)/float((tp+fn)))
        tresholds.append(score)

        #elif protein_tuple in interaction_list

    #    for complex in core_complex_list:
    #        if interaction_1 in complex and interaction_2 in complex:
    #            is_in_core_complex = True
    #            break
    #        else:
    #            is_in_core_complex = False

    #    if is_in_core_complex:
    #        y_true.append(1)
    #        scores.append(score)
    #    else:
    #        y_true.append(0)
    #        scores.append(score)
    return fpr, tpr, tresholds

file = open("../../data/protein_interaction_data/interactions_and_centers.txt")
interactions = []
for line in file:
    strline = str(line).strip().split()

    score = (abs(float(strline[0])-float(strline[2])) / 37.950100000000006)
    interactions.append((score, strline[1], strline[3]))
file.close()

with open('../interaction_data_exp_1_eu_8_final.pkl', 'rb') as infile:
    interactions = cPickle.load(infile)

with open('../interaction_data_exp_1_eu_8_final.pkl', 'rb') as infile:
    interactions_pr = cPickle.load(infile)

interactions.sort(reverse=True)
interactions_pr.sort(reverse=True)

core_complexes = pd.read_csv("../../data/protein_interaction_data/allComplexes.csv", sep=';')
core_complexes.fillna(0.0, inplace=True)

core_complex_list = []
for index, row in core_complexes.iterrows():

    if row["subunits (UniProt IDs)"] is not 0.0:
        split_line = row["subunits (UniProt IDs)"].split(',')
        cleaned_line = [i.strip('(').strip(')') for i in split_line]
        core_complex_list.append(cleaned_line)

print len(core_complex_list)
core_interactions = get_interactions(core_complex_list)
non_interactions = get_non_interactions(interactions, core_interactions)
#sys.exit()
#sys.exit()
print "Scoring interactions"

interactions = filter_interactions(interactions, core_complex_list)
interactions_pr = filter_interactions(interactions_pr, core_complex_list)
print len(interactions)
fpr, tpr, treshs = get_interaction_scoring(interactions, core_interactions, non_interactions)
#print tpr
#print treshs
fpr_pr, tpr_pr, treshs_pr = get_interaction_scoring(interactions_pr, core_interactions, non_interactions)


f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
sns.despine()

#fpr, tpr, treshs = metrics.roc_curve(y_true, scores)
#print metrics.roc_auc_score(y_true, scores)
#sys.exit()
#fpr_pr, tpr_pr, treshs_pr = metrics.roc_curve(y_true_pr, scores_pr)
#print metrics.roc_auc_score(y_true_pr, scores_pr)

print fpr
fp, tresh = tresh_at_fpr(fpr, treshs)
num_interactions, interactions_above_tresh = interactions_at_tresh(interactions, tresh)
print "SEC Interactions", num_interactions#, metrics.precision_score(y_true_above_tresh, [1]*len(y_true_above_tresh))

print fpr_pr
fp, tresh = tresh_at_fpr(fpr_pr, treshs_pr)
num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, tresh)
print "PR Interactions", num_interactions#, metrics.precision_score(y_true_above_tresh, [1]*len(y_true_above_tresh))


[i.set_linewidth(0.4) for i in ax.spines.itervalues()]
plt.xlim((0, 1))
plt.ylim((0, 1))
plt.ylabel("TPR")
plt.xlabel("FPR")

plt.plot(fpr, tpr)
plt.plot(fpr_pr, tpr_pr)
plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.svg", bbox_inches="tight", pad_inches=0.02, dpi=300)
plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.pdf", bbox_inches="tight", pad_inches=0.02, dpi=300)
plt.show()

