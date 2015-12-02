import os 
import sys 

class ResidueFeature:
    
    def __init__(self):
        
        self.feature_value = None 
        self.feature_name = None 
               
    def calculate_feature( self, residue_id, structure ):
        """Derived classes need to calculate features here!"""
        
        do = "something" 
    
    def set_name( self, name ):
        """Derived classes need to calculate features here!"""
        
        self.feature_name = name 
   
    def get_name( self ):
        """Derived classes need to calculate features here!"""
        
        return self.feature_name
   
    def get_feature(self):
        
        if self.feature_value == None:
            print "Warning: Feature value not initialized!" 
            return False
        
        else:
            return self.feature_value