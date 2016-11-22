from plots import plot_fdr
from metrics import compute_fdr_from_labels, compute_number_of_entries_at_fdr

def load_data_file(file_name):
    file = open(file_name)
    labels = [str(line).strip().split()[-1] for line in file]
    file.close()
    return labels

def compute_increase_at_fdr():
    fdr_cutoffs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for cutoff in fdr_cutoffs:

#        labels = load_data_file("HSA_%sPerc_new_PSM.txt" % cutoff)
#        fdr, hits = compute_fdr_from_labels(labels)

        labels_pr = load_data_file("HSA_%sPerc_new_PR.txt" % cutoff)
        fdr_pr, hits_pr = compute_fdr_from_labels(labels_pr)

        print cutoff, 399, hits_pr[-1], ((hits_pr[-1]-399)/float(399))*100


compute_increase_at_fdr()



