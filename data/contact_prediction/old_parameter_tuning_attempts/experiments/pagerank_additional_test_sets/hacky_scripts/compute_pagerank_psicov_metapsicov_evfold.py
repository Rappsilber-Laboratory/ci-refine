import os
import sys

from evaluation_helper import load_pdb_ids, refine_contacts_with_pagerank

def main():
    script_folder="/scratch/schneider/projects/pagerank_refinement/src/"
    test_set_folder="/scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/evfold/"
    output_folder="/scratch/schneider/projects/pagerank_refinement/results/pagerank_refinement/psicov_set/evfold/"

    protein_data = load_pdb_ids("/scratch/kstahl/test_sets/psicovsupp/" + "psicov_lengths")

    os.chdir(script_folder)

    for pdb_id, length in protein_data:
        contact_file= "".join([test_set_folder,
                             pdb_id,
                             ".prediction"])

        psipred_file= "".join(["/scratch/kstahl/test_sets/psicovsupp/",
                             "psipred/",
                             pdb_id,
                             ".horiz"])

        refine_contacts_with_pagerank(pdb_id, contact_file, length, psipred_file, str(2.0), output_folder)         
main()
