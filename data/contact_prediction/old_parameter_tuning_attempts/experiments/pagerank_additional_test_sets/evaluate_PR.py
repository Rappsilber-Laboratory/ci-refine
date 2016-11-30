import os 
import sys
import argparse

from evaluation_helper import load_pdb_ids, evaluate_contact_file

options = {}


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", dest="test_set_folder", help="Folder where the predictions are living", required=True)
    parser.add_argument("-a", dest="alpha", help="alpha parameter of pagerank", required=True)
    parser.add_argument("-b", dest="beta", help="beta parameter of pagerank", required=True)
    parser.add_argument("-l", dest="pdb_list", help="PDB id list", required=True)
    parser.add_argument("-p", dest="pdb_folder", help="Folder where predictions are living", required=True)
    parser.add_argument("-o", dest="out_file", help="Output file", required=True)
    options = parser.parse_args()



def main():
    parse_arguments()

    contact_script_folder="/scratch/schneider/Software/epc-map/rbocon/release/analysis/"
    test_set_folder=options.test_set_folder
    protein_data = load_pdb_ids(options.pdb_list)

    os.chdir(contact_script_folder)    

    for pdb_id, length in protein_data:
        pdb_file = "".join([options.pdb_folder,
                             "pdb/",
                             pdb_id,
                             ".pdb"])

        contact_file= "".join([test_set_folder,
                             pdb_id,
                             "_RRPAR_%s_%s__0" % (options.alpha, options.beta)])

        dom_ass = "1-%s" % length
        out_file = options.out_file
        evaluate_contact_file(pdb_file, pdb_id, contact_file, length, dom_ass, out_file)
main()
