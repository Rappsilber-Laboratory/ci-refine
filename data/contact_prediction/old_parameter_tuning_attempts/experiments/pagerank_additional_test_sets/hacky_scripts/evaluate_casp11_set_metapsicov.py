import os 
import sys

from evaluation_helper import load_pdb_ids_dom_def, evaluate_contact_file

def main():
    contact_script_folder="/scratch/schneider/projects/pagerank_refinement/data/casp11_targets/"
    test_set_folder="/scratch/kstahl/test_sets/casp11_official/"

    protein_data = load_pdb_ids_dom_def("/scratch/kstahl/test_sets/casp11_official/casp_results/domains_fm_targets_split")

    os.chdir(contact_script_folder)    

    for pdb_id, length, dom_def in protein_data:
        pdb_file = "".join([test_set_folder,
                             "pdb/",
                             pdb_id,
                             ".pdb"])

        contact_file= "".join([test_set_folder,
                             "stage2_official/",
                             pdb_id,
                             ".prediction"])

        dom_ass = dom_def
        
        evaluate_contact_file(pdb_file, pdb_id, contact_file, length, dom_ass)
main()
