"""Author: Michael Schneider
   Very important notes
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("experiments/")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/structure/")

#import helper
import git_tag
from optparse import OptionParser
import numpy
import InputOutput
## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-r", type="string", dest="restraint_file", help="restraint_file")
    parser.add_option("-f", type="string", dest="fasta_file", help="Fasta_file")
    options, args = parser.parse_args()
    return options, args 

options, args = add_options(parser)

most_distal = {'A': 'CB',
               'G': 'CA',
               'S': 'OG',
               'T': 'CG2',
               'C': 'SG',
               'V': 'CG2',
               'L': 'CD2',
               'I': 'CD1',
               'M': 'CE',
               'P': 'CG',
               'F': 'CZ',
               'Y': 'OH',
               'W': 'CZ3',
               'D': 'OD1',
               'E': 'OE1',
               'N': 'OD1',
               'Q': 'NE2',
               'H': 'CE1',
               'K': 'NZ',
               'R': 'NH2' }

def main():
    """Generic main function. Executes main functionality of program
    """
    """Some inits"""

    res = InputOutput.InputOutput.load_rosetta_restraints(options.restraint_file)
    sequence = InputOutput.InputOutput.read_fasta(options.fasta_file)

    seq_map = {}
    for index, aa in enumerate(sequence):
        seq_map[index+1] = aa
    restraints = []
    for i, j in res:
        restraints.append((i, most_distal[seq_map[i]], j, most_distal[seq_map[j]],1.0))
    InputOutput.InputOutput.write_restraint_file(restraints, 'test_fa_restraints.txt', backbone=False, upper_distance = 5.0)


    #print res
    #print sequence
if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  
