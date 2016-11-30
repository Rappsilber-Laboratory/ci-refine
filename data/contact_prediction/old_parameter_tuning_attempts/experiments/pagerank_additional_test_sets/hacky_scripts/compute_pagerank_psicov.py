import os
import sys

from evaluation_helper import load_pdb_ids, refine_contacts_with_pagerank

def main():
    script_folder="/scratch/schneider/projects/pagerank_refinement/src/"
    test_set_folder="/scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/"
    output_folder="/scratch/schneider/projects/pagerank_refinement/results/pagerank_refinement/psicov_set/epc-map_new_par/"

    protein_data = load_pdb_ids(test_set_folder + "psicov_lengths")

    os.chdir(script_folder)

    for pdb_id, length in protein_data:
        contact_file= "".join([test_set_folder,
                             "epc-map/",
                             pdb_id,
                             ".prediction"])

        psipred_file= "".join([test_set_folder,
                             "psipred/",
                             pdb_id,
                             ".horiz"])

        refine_contacts_with_pagerank(pdb_id, contact_file, length, psipred_file, str(3.0), output_folder, alpha=0.5)         
main()
