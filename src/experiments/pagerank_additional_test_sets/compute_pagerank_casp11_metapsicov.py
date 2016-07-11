import os
import sys

from evaluation_helper import load_pdb_ids, refine_contacts_with_pagerank

def main():

    script_folder="/scratch/schneider/projects/pagerank_refinement/src/"
    test_set_folder="/scratch/kstahl/test_sets/casp11/"
    output_folder="/scratch/schneider/projects/pagerank_refinement/results/pagerank_refinement/casp11/new_set/metapsicov/"

    protein_data = load_pdb_ids("casp11_targets_filtered")

    os.chdir(script_folder)

    for pdb_id, length in protein_data:
        contact_file= "".join(["/scratch/kstahl/test_sets/casp11_official/casp_results/",
                             pdb_id+"/",
                             pdb_id,
                             "RR021_1"])

        psipred_file= "".join([test_set_folder,
                             "psipred/",
                             pdb_id,
                             ".horiz"])

        refine_contacts_with_pagerank(pdb_id, contact_file, length, psipred_file, str(2.0), output_folder)         
main()
