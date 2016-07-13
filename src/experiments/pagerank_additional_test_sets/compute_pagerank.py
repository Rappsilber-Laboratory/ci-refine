import os
import sys
import argparse
from evaluation_helper import load_pdb_ids, refine_contacts_with_pagerank


options = {}

def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", dest="test_set_folder", help="Folder where the predictions are living", required=True)
    parser.add_argument("-a", dest="alpha", help="alpha parameter of pagerank", required=True)
    parser.add_argument("-b", dest="beta", help="beta parameter of pagerank", required=True)
    parser.add_argument("-s", dest="psipred_folder", help="PSIPRED folder", required=True)
    parser.add_argument("-l", dest="pdb_list", help="PDB id list", required=True)
    parser.add_argument("-o", dest="out_folder", help="Output folder", required=True)
    options = parser.parse_args()



def main():
    parse_arguments()
    script_folder="/scratch/schneider/projects/pagerank_refinement/src/"
    test_set_folder=options.test_set_folder
    output_folder=options.out_folder

    protein_data = load_pdb_ids(options.pdb_list)

    os.chdir(script_folder)

    for pdb_id, length in protein_data:
        contact_file= "".join([test_set_folder,
                             pdb_id,
                             ".prediction"])

        psipred_file= "".join([options.psipred_folder,
                             pdb_id,
                             ".horiz"])

        refine_contacts_with_pagerank(pdb_id, contact_file, length, psipred_file, str(options.beta), output_folder, alpha=options.alpha)         
main()
