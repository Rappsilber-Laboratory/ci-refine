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

def get_interaction_scoring(interactions, core_complex_list):
    y_true = []
    scores = []
    for score, interaction_1, interaction_2 in interactions:
        # print interaction_1, interaction_2
        is_in_core_complex = False
        for complex in core_complex_list:
            if interaction_1 in complex and interaction_2 in complex:
                # print interaction_1, interaction_2, complex
                #print interaction_1, interaction_2, complex
                is_in_core_complex = True
                break
            else:
                is_in_core_complex = False
        if is_in_core_complex:
            y_true.append(1)
            scores.append(score)
        else:
            y_true.append(0)
            scores.append(score)
    return y_true, scores



# print interaction_1, interaction_2

    

file = open("../../data/protein_interaction_data/interactions_and_centers.txt")
interactions = []
for line in file:
    strline = str(line).strip().split()

    score = (abs(float(strline[0])-float(strline[2])) / 37.950100000000006)
    interactions.append((score, strline[1], strline[3]))
file.close()



with open('../interaction_data_euclidean_test.pkl', 'rb') as infile:
    interactions = cPickle.load(infile)

with open('../output_interactions_pagerank_eucleadian.pkl', 'rb') as infile:
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
#sys.exit()
print "Scoring interactions"


y_true, scores = get_interaction_scoring_list(interactions, core_complex_list)

y_true_pr, scores_pr = get_interaction_scoring_list(interactions_pr, core_complex_list)


f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
sns.despine()

fpr, tpr, treshs = metrics.roc_curve(y_true, scores)
print metrics.roc_auc_score(y_true, scores)
#sys.exit()
fpr_pr, tpr_pr, treshs_pr = metrics.roc_curve(y_true_pr, scores_pr)
print metrics.roc_auc_score(y_true_pr, scores_pr)

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

