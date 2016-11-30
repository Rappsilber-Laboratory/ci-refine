import os 
import sys

from evaluation_helper import load_pdb_ids_dom_def, evaluate_contact_file

def main():

    contact_script_folder="/scratch/schneider/Software/epc-map/rbocon/release/analysis/"
    test_set_folder="/scratch/schneider/projects/pagerank_refinement/data/casp11_targets/"
    protein_data = load_pdb_ids_dom_def("casp11_domains_available_in_tarball")

    top_groups = ["021",
                  "124",
                  "420",
                  "398",
                  "410",
                  "479",
                  "008",
                  "041",
                  "086",
                  "262"]

    for group_id in top_groups:
        cmd = " ".join(["python ../../evaluation/quick_check_csv.py", 
                       "--file", "%s_casp11_PR_beta_3.txt" % group_id])
        os.system(cmd)


main()
