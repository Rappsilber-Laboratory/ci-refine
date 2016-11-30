import os 
import sys

from evaluation_helper import load_pdb_ids_dom_def, evaluate_contact_file

def main():
    contact_script_folder="/scratch/schneider/Software/epc-map/rbocon/release/analysis/"
    test_set_folder="/scratch/schneider/projects/pagerank_refinement/data/casp11_targets/"

    protein_data = load_pdb_ids_dom_def("casp11_domains_available_in_tarball")

    os.chdir(contact_script_folder)    

    for pdb_id, length, dom_def in protein_data:
        pdb_file = "".join([test_set_folder,
                             pdb_id[0:5],
                             ".pdb"])

        contact_file= "".join(["/scratch/schneider/projects/pagerank_refinement/results/pagerank_refinement/casp11/new_set/metapsicov/",
                             pdb_id[0:5],
                             "_RRPAR_0.4_2.0__0"])

        dom_ass = dom_def
        
        evaluate_contact_file(pdb_file, pdb_id, contact_file, length, dom_ass)
main()
