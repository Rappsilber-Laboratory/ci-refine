"""Author: Michael Schneider
""" 

import os 
import sys 
sys.path.append("../src")
import test.test_stub 
from optparse import OptionParser

## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-e", type="string", dest="example", help="An example")
    options, args = parser.parse_args()
    return options, args 

def main():
   """Generic main function. Executes main functionality of program
   """
   print test.test_stub.stub_function()
if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  
