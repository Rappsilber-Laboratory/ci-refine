import os 
import sys
import GENERAL
import os 
import sys 
import random
import numpy
from ResidueFeature import *

class ResidueFeatureSecStruct( ResidueFeature ):
    
    def __init__(self, structure):
        
        self.feature_value = None 
        self.feature_name = "SecStruct"
        '''your stride bin here!!!'''
        self.stride_bin = "/scratch/schneider/Software/Stride/stride -h " 
        self.stride_out = None
        rand_val = random.randint(1, 99999)
        os.system(self.stride_bin + " " + structure + " > temp" + str(rand_val ) )
        
        self.stride_out = GENERAL.file2list( "temp" + str(rand_val ) )
        os.system("rm -rf temp" + str(rand_val ))
        
        self.ss_dict = {}
        self.acc_dict = {}
        self.dnr_dict = {}
        
        max_resnum = 0
        
        for line in self.stride_out: 
            if line.find("ASG") == 0:
                line_list = line.split()
                max_resnum = int(line_list[3])
        
        for i in range(1, max_resnum + 1):
            self.ss_dict[i] = ""
            self.acc_dict[i] = 0.0
            self.dnr_dict[i] = 0.0
        
        for residue_id in range(1, max_resnum + 1):
            acc_count =  0
            don_count = 0
            res_acc_dict = {}
            res_don_dict = {}
            for line in self.stride_out: 
                if line.find("ASG") == 0:
                    line_list = line.split()
                    res_num = int(line_list[3])
                    ss_struct = line_list[5]
                    if ss_struct == 'B':
                        ss_struct = 'C'
                    elif ss_struct == 'b':
                        ss_struct = 'C'
                    elif ss_struct == 'G':
                        ss_struct = 'H'
                    elif ss_struct == 'I':
                        ss_struct = 'H'
                    elif ss_struct == 'T':
                        ss_struct = 'C'
                    if res_num == residue_id:
                        self.ss_dict[residue_id] = ss_struct
              
                if line.find("ACC") == 0:
                    line_list = line.split()
                    res_num = int(line_list[3])
                    
                    donoring_residue = int(line_list[8])
                    
                    
                    if res_num == residue_id:
                        res_acc_dict[donoring_residue] = 1
                        #acc_count += 1
           
                if line.find("DNR") == 0:
                    line_list = line.split()
                    res_num = int(line_list[3])
                    accepting_residue = int(line_list[8])
                    if res_num == residue_id:
                        #accepting_residue = int(line_list[8])
                        res_don_dict[accepting_residue] = 1
                     #   don_count += 1
            self.acc_dict[residue_id] = res_acc_dict
            self.dnr_dict[residue_id] = res_don_dict
    
    def calculate_SS_limits( self, residue_id, structure ):
        #print structure
        current_ss = self.ss_dict[residue_id]
        current_residue_id = residue_id
        lower_limit = residue_id 
        upper_limit = residue_id
        """Get lower SS limit"""
        while self.ss_dict[ current_residue_id ] == current_ss:
            if self.ss_dict.has_key( current_residue_id ):
                lower_limit = current_residue_id
                current_residue_id -= 1
                if current_residue_id < 1:
                    break
            else:
                break
        """Get upper SS limit"""
        current_residue_id = residue_id
        while self.ss_dict[ current_residue_id ] == current_ss:
            if self.ss_dict.has_key( current_residue_id ):
                upper_limit = current_residue_id
                current_residue_id += 1
                if current_residue_id > len(self.ss_dict):
                    break
            else:
                break

        
        current_residue_id = residue_id
        diff_vector = structure.get_bioPBD_pose()[residue_id]["CA"].coord - structure.get_centroid()
        #print "DIST TO CENTROID:", numpy.sqrt(numpy.sum(diff_vector * diff_vector))
        if lower_limit == upper_limit:
            return 1, 0.0, lower_limit, upper_limit, numpy.sqrt(numpy.sum(diff_vector * diff_vector)) 
        else:
            return upper_limit - lower_limit + 1,  structure.get_contact_map().get_dm()[lower_limit - 1, upper_limit - 1], lower_limit, upper_limit, numpy.sqrt(numpy.sum(diff_vector * diff_vector))

    
    def calculate_feature( self, residue_id, structure, RelSasaFeature = None ):
        """Derived classes need to calculate features here!"""
        
      #  do = "something"
      #  

        
      #  for line in self.stride_out: 
      #      if line.find("ASG") == 0:
      #          
      #          line_list = line.split()
      #          res_num = int(line_list[3])
      #          ss_struct = line_list[5]
                
      #          if res_num == residue_id:
      #              self.feature_value = ss_struct

        
        if self.feature_value == None:
            print "Warning! Secondary structure for residue " + str(residue_id) + " not calculated! Check residue id!"
        
        if RelSasaFeature == None:
            self.feature_value = (self.ss_dict[residue_id], 
                                  self.acc_dict[residue_id], 
                                  self.dnr_dict[residue_id], 
                                  self.calculate_SS_limits(residue_id, structure)[0], 
                                  self.calculate_SS_limits(residue_id, structure)[1] )
        else:
            lim_tupel = self.calculate_SS_limits(residue_id, structure)
            
            lower_lim = lim_tupel[2]
            upper_lim = lim_tupel[3]
            dist_to_cen = lim_tupel[4]
            aa_len = lim_tupel[0]
            l_3d_len = lim_tupel[1]
            
            bur_count = 0
            acc_count = 0 
            for i in range(lower_lim, upper_lim + 1):
                RelSasaFeature.calculate_feature(i, structure )
                if RelSasaFeature.get_feature() <= 0.25:
                    bur_count += 1
                if RelSasaFeature.get_feature() > 0.25:
                    acc_count += 1
            
            self.feature_value = (self.ss_dict[residue_id], 
                                  self.acc_dict[residue_id], 
                                  self.dnr_dict[residue_id], 
                                  aa_len, 
                                  l_3d_len,
                                  bur_count,
                                  acc_count,
                                  dist_to_cen
                                  )
            
        
        
       
   
