import os
import sys

class InputOutput:
    
    def __init__(self):
        pass


    @staticmethod
    def load_xl_data( xl_file ):
        col_names = {}
        max_score = 0
        xls = []
        gt_data = []
        for line in open(xl_file):
            if line.startswith('LinkID') or not line.strip():
                for index, col in enumerate(str(line).split(',')):
                    col_names[col] = index
            else:
                line = line.split(',')

                from_site = int(line[col_names['fromSite']])
                to_site = int(line[col_names['ToSite']])
                score = float(line[col_names['Score']])
                is_decoy = line[col_names['isDecoy']]



                if from_site > 0 and to_site > 0 and abs(from_site-to_site) >= 12 and is_decoy == 'false':
                    if max_score == 0:
                        max_score = score
                    xls.append(((from_site, to_site), score/max_score))
                    gt_data.append((from_site, 'CA', to_site, 'CA', score / max_score))

        return xls, gt_data

    """
        file  = open(xl_file)
        from_site = 0
        to_site = 0
        score = 0
        gt_data = []
        xls = []
        is_decoy = 0
        for line in file:
            strline = str(line).strip().split(',')
            from_site = int(strline[7])-28
            to_site = int(strline[8])-28
            score = float(strline[9])
            is_decoy = strline[10]
            if from_site > 0 and to_site > 0 and abs(from_site-to_site) >= 12 and is_decoy == 'FALSE':
                xls.append(((from_site, to_site), score/30.0))
                gt_data.append((from_site, 'CA', to_site, 'CA', score))
        file.close()
        InputOutput.InputOutput.write_contact_file(gt_data, 'gt', upper_distance = 20)
        return xls
    """


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

