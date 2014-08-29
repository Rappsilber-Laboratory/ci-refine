import os
import sys

class InputOutput:
    
    def __init__(self):
        pass
    
    @staticmethod
    def read_fasta( fasta_file ):
        file = open( fasta_file,'r' )
        sequence = ''
        for line in file:
            strline = str(line).strip()
            if strline.count(">") == 0:
                sequence += strline
        file.close()
        return sequence
    
    @staticmethod
    def write_restraint_file( restraint_object, restraint_file_name ):
        file = open( restraint_file_name, 'w' )
        for c_lower, atom_lower, c_upper, atom_upper, prob in restraint_object:
            
            if atom_lower == 'G':
                atom_1 = "CA"
            else:
                atom_1 = "CB"
            
            if atom_upper == 'G':
                atom_2 = "CA"
            else:
                atom_2 = "CB"
            
            file.write( "AtomPair %s  %s %s  %s BOUNDED  %s %s %s NOE \n"%(atom_1,c_lower,atom_2, c_upper, 1.5 ,8.0, 2.0) )
        file.close()
    
    @staticmethod
    def write_contact_file(  contacts, contact_file_name, upper_distance = 8 ):
        file = open( contact_file_name, 'w' )
        for c_lower, atom_lower, c_upper, atom_upper, prob in contacts:
            file.write(" ".join(["%s"%(c_lower),
                                 "%s"%(c_upper),
                                 "%s"%(0),
                                 "%s"%(upper_distance),
                                 "%s"%(prob),
                                 "\n"]))
        file.close()
    
    @staticmethod
    def load_restraints( restraint_file, seq_sep_min = 24, seq_sep_max=9999 ):
        file = open(restraint_file,"r")
        res = []
        res_dict = {}
        for line in file:
            strline = str(line).strip().split()
            if len(strline) > 2:
                if strline[0] != "REMARK" and strline[0] != "METHOD" and len(strline[0]) <= 35:
                    if abs(int(strline[0]) - int(strline[1])) >= seq_sep_min and abs(int(strline[0]) - int(strline[1])) < seq_sep_max:
                        if res_dict.has_key((int(strline[0]), int(strline[1]))) == False or res_dict.has_key((int(strline[0]), int(strline[1]))) == False:
                             res.append( (float(strline[-1]),int(strline[0]), int(strline[1]), float(strline[2]), float(strline[3]) ) )
                             res_dict[(int(strline[0]), int(strline[1]))] = 1
                             res_dict[(int(strline[1]), int(strline[0]))] = 1
        file.close()
        res.sort()
        res.reverse()
        return res

    @staticmethod
    def load_restraints_pr( restraint_file, seq_sep_min = 12, seq_sep_max=9999 ):
        file = open(restraint_file,"r")
        res = []
        res_dict = {}
        counter = 0
        for line in file:
            strline = str(line).strip().split()
            if len(strline) > 2:
                if strline[0] != "REMARK" and strline[0] != "METHOD" and len(strline[0]) <= 35:
                    if abs(int(strline[0]) - int(strline[1])) >= seq_sep_min and abs(int(strline[0]) - int(strline[1])) < seq_sep_max:
                        if res_dict.has_key((int(strline[0]), int(strline[1]))) == False or res_dict.has_key((int(strline[0]), int(strline[1])))==False:
                             if counter == 0:
                                norm = float(strline[-1])
                             res.append( ((( int(strline[0]), int(strline[1])) , float(strline[-1]) / norm )))
                             res_dict[(int(strline[0]), int(strline[1]))] = 1
                             res_dict[(int(strline[1]), int(strline[0]))] = 1
                             counter += 1
        file.close()
        #res.sort()
        #res.reverse()
        return res

