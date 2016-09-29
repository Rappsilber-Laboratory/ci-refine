"""Author: Michael Schneider
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
#sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")
#sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/structure/")
import InputOutput
from optparse import OptionParser

import structure.StructureContainer

## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options(parser):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-p", type="string", dest="pdb_list", help="restraint_file")

    options, args = parser.parse_args()
    return options, args

options, args = add_options(parser)


def write_native_map(pdb_id):
    pdb_file = "/scratch/schneider/projects/pagerank_refinement/data/predictor_results/d329/pdb/%s.pdb"%pdb_id
    tmp_struct = structure.StructureContainer.StructureContainer()
    try:
        tmp_struct.load_structure('xxxx', pdb_id[-1], pdb_file, seqsep=1)
    except:
        tmp_struct.load_structure('xxxx', ' ', pdb_file, seqsep=1)

    native_contacts = tmp_struct.get_contact_map().contact_list()
    InputOutput.InputOutput.write_contact_file(native_contacts, "../data/native_maps/%s.map"%pdb_id)

def main():
    file = open(options.pdb_list)
    for line in file:
        strline = str(line).strip().split()
        write_native_map(strline[0])
    file.close()

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
