import cPickle
import itertools
from read_protein_protein_interactions import filter_interactions, get_interaction_scoring, get_interactions


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


def get_precision(filtered_interactions, interactions_corum):
    y_true, y_score = get_interaction_scoring(filtered_interactions, interactions_corum, [])

    hits = [y for y in y_true if y == 1]

    return round(float(len(hits)) / float(len(y_true)), 2)


# Data Scored by Euclidean Distance
with open('interactions_exp1_0.55_final_orig.pkl', 'rb') as infile:
    interactions_1 = cPickle.load(infile)

with open('interactions_exp2_0.55_final_orig.pkl', 'rb') as infile:
    interactions_2 = cPickle.load(infile)

with open('interactions_exp3_0.55_final_orig.pkl', 'rb') as infile:
    interactions_3 = cPickle.load(infile)

# Interactions PR
with open('interactions_exp1_0.55_final.pkl', 'rb') as infile:
    interactions_pr_1 = cPickle.load(infile)

with open('interactions_exp2_0.55_final.pkl', 'rb') as infile:
    interactions_pr_2 = cPickle.load(infile)

with open('interactions_exp3_0.55_final.pkl', 'rb') as infile:
    interactions_pr_3 = cPickle.load(infile)


with open('core_complexes.pkl', 'rb') as infile:
    core_complexes = cPickle.load(infile)


interactions_corum = get_interactions(core_complexes)


all_interactions = list(set(list(itertools.chain(interactions_1, interactions_2, interactions_3))))
all_interactions_with_score = [(0, i[0], i[1]) for i in all_interactions]

all_interactions_pr = list(set(list(itertools.chain(interactions_pr_1, interactions_pr_2, interactions_pr_3))))
all_interactions_with_score_pr = [(0, i[0], i[1]) for i in all_interactions_pr]


filtered_interactions = filter_interactions(all_interactions_with_score, core_complexes)
filtered_interactions_pr = filter_interactions(all_interactions_with_score_pr, core_complexes)

filtered_interactions = get_protein_intersection(filtered_interactions, filtered_interactions_pr, level='protein')
filtered_interactions_pr = get_protein_intersection(filtered_interactions_pr, filtered_interactions, level='protein')

print('Euclidean:')
print('Number of interactions: %s' % len(all_interactions))
print('Number of filtered_interactions: %s' % len(filtered_interactions))
print('Precision: %s' % get_precision(filtered_interactions, interactions_corum))
print
print('PageRank:')
print('Number of interactions: %s' % len(all_interactions_pr))
print('Number of filtered_interactions: %s' % len(filtered_interactions_pr))
print('Precision: %s' % get_precision(filtered_interactions_pr, interactions_corum))


