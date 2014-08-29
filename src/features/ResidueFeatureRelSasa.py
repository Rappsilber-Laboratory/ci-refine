import os 
import sys
import GENERAL
import os 
import sys 
import random
from ResidueFeature import *

class ResidueFeatureRelSasa( ResidueFeature ):
    
    def __init__(self, structure):
        self.feature_value = None 
        self.feature_name = "RelSasa"
        #self.pops_bin = "/scratch/schneider/Software/Stride/stride" 
        self.pops_bin=  "/scratch/schneider/Software/pops/bin/pops"
        self.pops_out = None
        
        """Residue Surface area from Rost and Sander Proteins 20,216-226"""
        self.surface_areas = {'A':106, 'C':135, 'D':163, 'E':194, 'F':197,
                              'G': 84, 'H':184, 'I':169, 'K':205, 'L':164,
                              'M':188, 'N':157, 'P':136, 'Q':198, 'R':248,
                              'S':130, 'T':142, 'V':142, 'W':227, 'Y':222 }
        
        self.three_to_one_names = { 'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                                    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V' }



        os.system(self.pops_bin + " --pdb " + structure + " --residueOut" )
        
        self.pops_out = GENERAL.file2list( "pops.out" )
       # os.system("rm -rf pops.out")
       # os.system("rm -rf sigma.out")       
    def calculate_feature( self, residue_id, structure ):
        """Derived classes need to calculate features here!"""
          
        for line in self.pops_out: 
            line_list = line.split()
            
            if len(line_list) >= 7 and len(line_list[0]) == 3:
                
                residue_type = self.three_to_one_names[line_list[0]]
                res_num = int(line_list[2])
                sasa = float(line_list[5])
                
                if res_num == residue_id:
                    self.feature_value = ( sasa / float( self.surface_areas[ residue_type ]  ) )
      
        if self.feature_value == None:
            print "Warning! Secondary structure for residue " + str(residue_id) + " not calculated! Check residue id!"
        
        
       
   
