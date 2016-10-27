import pandas as pd
from matplotlib.pylab import plt
import cPickle
import seaborn as sns
from sklearn.metrics import roc_curve, precision_recall_curve
import sys
import numpy as np


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


def tresh_at_fpr(fpr, thresholds, fpr_limit=0.53):
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
        #print fp, thresh
        if fp <= fpr_limit:

            return fp, thresh


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


def get_interaction_scoring_list(interactions, core_complex_list):
    """
    DEPRECATED

    Parameters
    ----------
    interactions
    core_complex_list

    Returns
    -------

    """
    y_true = []
    scores = []
    for score, protein_1, protein_2 in interactions:

        is_in_core_complex = False
        for complex in core_complex_list:
            # test overlaps
            if any(i in protein_1 for i in complex) and any(i in protein_2 for i in complex):
                is_in_core_complex = True
                break
        if is_in_core_complex:
            y_true.append(1)
            scores.append(score)
        else:
            y_true.append(0)
            scores.append(score)

    return y_true, scores


def filter_interactions(interactions, core_complex_list):
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
        protein_1_in_corum = False
        protein_2_in_corum = False
        protein_1_complex = 0
        protein_2_complex = 0
        for complex in core_complex_list:
            if protein_1 in complex:
                protein_1_in_corum = True
                protein_1_complex = complex

            if protein_2 in complex:
                protein_2_in_corum = True
                protein_2_complex = complex
        if protein_1_in_corum and protein_2_in_corum:
            #
            # print protein_1, protein_1_complex
            #print protein_2, protein_2_complex
            #print
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

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    tpr = []
    fpr = []
    tresholds = []
    count = 0

    y_true = []
    y_score = []
    for score, protein_1, protein_2 in interactions:
        protein_tuple = sorted_uniprot_tuple(protein_1, protein_2)
        """
        # True positives
        if protein_tuple in interaction_list:
            tp += 1
        # False positives
        elif protein_tuple not in interaction_list:
            fp += 1

        tn = len(non_interaction_list) - (tp+fp)
        fn = max(0, len(interaction_list) - (tp+fp))

        fpr.append(float(fp)/float((fp+tn)))
        tpr.append(float(tp)/float((tp+fn)))

        tresholds.append(score)
        count += 1
        """
        #print protein_tuple
        #print interaction_list
        if protein_tuple in interaction_list:
            y_true.append(1)
            y_score.append(score)
        else:
            y_true.append(0)
            y_score.append(score)

    return y_true, y_score # fpr, tpr, tresholds


def main():
    # Load interaction data
    with open('../interaction_data_exp_1_eu_6_final.pkl', 'rb') as infile:
        interactions = cPickle.load(infile)

    with open('../interaction_data_exp_1_eu_6_final_pr.pkl', 'rb') as infile:
        interactions_pr = cPickle.load(infile)
    interactions.sort(reverse=True)
    interactions_pr.sort(reverse=True)
    interactions_pr = [(interaction[0], interaction[1], interaction[2]) for interaction in interactions_pr]
    #for i in interactions_pr:
    #    print i
    #x = [interaction[0] for interaction in interactions_pr]
    #print x
    #plt.plot(x)
    #sys.exit()

    # Read Corom complexes
    core_complexes = pd.read_csv("../../data/protein_interaction_data/allComplexes.csv", sep=';')
    core_complexes.fillna(0.0, inplace=True)
    core_complex_list = []
    for index, row in core_complexes.iterrows():
        if row["subunits (UniProt IDs)"] is not 0.0:
            split_line = row["subunits (UniProt IDs)"].split(',')
            cleaned_line = [i.strip('(').strip(')') for i in split_line]
            core_complex_list.append(cleaned_line)

    with open('core_complexes.pkl', 'wb') as outfile:
        cPickle.dump(core_complex_list, outfile, cPickle.HIGHEST_PROTOCOL)
    # Filter interactions by proteins that are actually present in Corum
    filtered_interactions = filter_interactions(interactions, core_complex_list)
    filtered_interactions_pr = filter_interactions(interactions_pr, core_complex_list)

    # Get Interaction and non interaction dictionaries
    core_interactions = get_interactions(core_complex_list) # Is this on filtered proteins?
    non_interactions = get_non_interactions(filtered_interactions, core_interactions) # Is this on filtered proteins?

    print len(non_interactions)
    print len(filtered_interactions)
    print len(filtered_interactions_pr)
    #for i in range(0, len(non_interactions) - len(filtered_interactions)):
    #    filtered_interactions.append((0, "fake1", "fake2"))
    #for i in range(0, len(non_interactions) - len(filtered_interactions_pr)):
    #    filtered_interactions_pr.append((0, "fake1", "fake2"))
    # Compute FPR, TPR curve and thresholds
    y_true, y_score = get_interaction_scoring(filtered_interactions, core_interactions, non_interactions)
    y_true_pr, y_score_pr = get_interaction_scoring(filtered_interactions_pr, core_interactions, non_interactions)

    #fpr, tpr, treshs = get_interaction_scoring(filtered_interactions, core_interactions, non_interactions)
    #fpr_pr, tpr_pr, treshs_pr = get_interaction_scoring(filtered_interactions_pr, core_interactions, non_interactions)

    #fpr, tpr, treshs = roc_curve(y_true, y_score)
    #fpr_pr, tpr_pr, treshs_pr = roc_curve(y_true_pr, y_score_pr)
    fpr, tpr, treshs = precision_recall_curve(y_true, y_score)
    fpr_pr, tpr_pr, treshs_pr = precision_recall_curve(y_true_pr, y_score_pr)

    print("Cutoff List Standard")
    for p, r, t in zip(fpr, tpr, treshs):
        print p, r, t
    print("Cutoff List CI")
    for p, r, t in zip(fpr_pr, tpr_pr, treshs_pr):
        print p, r, t

    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.176) # Exp 1
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.241044043869)  # Exp 1 Target precision 0.6
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.394482660393)  # Exp 1 Target precision 0.65
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.227327282579) # Exp 2
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.341766955407)  # Exp 2 # Target precision 0.6
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions,
    #                                                                        0.378015642113)  # Exp 2 # Target precision 0.65
    # precision 0.53 returns pretty much all interactions...
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions, 0.018069923595)  # Exp 3 # Target precision 0.53
    #num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions, 0.148888315685)  # Exp 3 # Target precision 0.6
    num_interactions, interactions_above_tresh_orig = interactions_at_tresh(interactions,
                                                                            0.20356058527)  # Exp 3 # Target precision 0.65
    print("SEC Interactions", num_interactions)

    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 0.00010118782017)  # Exp 1
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 0.000138949631778)  # Exp 1 # Target precision 0.6
    num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr,
                                                                       0.000188719955985)  # Exp 1 # Target precision 0.62
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 0.000211079165048)  # Exp 1 # Target precision 0.65
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 6.73617417516e-05) # Exp 2
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 9.7943746924e-05)  # Exp 2 # Target precision 0.6
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr,
    #                                                                   0.000111704795997)  # Exp 2 # Target precision 0.62
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr,
    #                                                                   0.000148467824096)  # Exp 2 # Target precision 0.65

    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 6.64360536357e-06)  # Exp 3 # Target Precision 0.53
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr, 2.36669556527e-05)  # Exp 3 # Target Precision 0.6
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr,
    #                                                                   3.03039503386e-05)  # Exp 3 # Target Precision 0.62
    #num_interactions, interactions_above_tresh = interactions_at_tresh(interactions_pr,
    #                                                                   3.67293035622e-05)  # Exp 3 # Target Precision 0.65



    print("PR Interactions", num_interactions)

    interaction_list = []
    for score, protein_1, protein_2 in interactions_above_tresh:
        interaction_list.append(sorted_uniprot_tuple(protein_1, protein_2))

    with open('interactions_exp1_0.62.pkl', 'wb') as outfile:
            cPickle.dump(interaction_list, outfile, cPickle.HIGHEST_PROTOCOL)
   # 2.45933254562e-05



    # Plot results
    set_plot_parameters()
    sns.set(style="white", palette="muted")
    f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
    sns.despine()

    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]

    plt.xlim((0, 1.0))
    plt.ylim((0, 1.0))
    plt.ylabel("TPR")
    plt.xlabel("FPR")

    plt.plot(fpr, tpr)
    plt.plot(fpr_pr, tpr_pr)
    plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.svg", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.png", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.savefig("protein_protein_interaction_pagerank_alpha_no_pers_0.85.pdf", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.show()
main()
