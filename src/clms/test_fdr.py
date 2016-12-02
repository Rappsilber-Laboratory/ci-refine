from src.plots import plot_fdr
from src.metrics import compute_fdr_from_labels, compute_number_of_entries_at_fdr

def load_data_file(file_name):
    file = open(file_name)
    labels = [str(line).strip().split()[-1] for line in file]
    file.close()
    return labels

labels = load_data_file("HSA_13Perc_new_PSM.txt")
fdr, hits =  compute_fdr_from_labels(labels, target_fdr=99)

for f, h in zip(fdr, hits):
    print f, h

labels_pr = load_data_file("HSA_13Perc_new_PR.txt")
fdr_pr, hits_pr =  compute_fdr_from_labels(labels_pr, target_fdr=99)

for f, h in zip(fdr_pr, hits_pr):
    print f, h

plot_fdr(fdr, hits, fdr_pr, hits_pr, out_folder="../results/pagerank_cross_links/")



