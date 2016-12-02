"""
MIT License

Copyright (c) 2016 Michael Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import sys 

import numpy
from optparse import OptionParser
import src.features.ResidueFeatureSecStruct as ResidueFeatureSecStruct
import src.structure.StructureContainer as StructureContainer
import numpy as np
import cPickle
## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-p", type="string", dest="pdb_id_list", help="An example")

    options, args = parser.parse_args()
    return options, args 

options, args  = add_options( parser )


def shift_matrix():
    matrix = []
    for i in xrange(-8,9):
        row = []
        for j in xrange(-8,9):
            if i == 0 and j == 0:
                pass
            else:
                row.append((i,j))
                matrix.append((i,j))
    return matrix


def calculate_probabities( all_shifts ):
    shift_len = len(all_shifts[0])
    for i in xrange(0,shift_len):
        shift_val = 0.0
        contact_val = 0.0
        for vec in all_shifts:
            shift_val += 1.0
            if vec[i] ==1.0:
                contact_val += 1.0
        print contact_val / shift_val


def get_relative_sec_struct_pos(ss_dict, i):
    anchor = ss_dict[i]
    pos = 1

    for j in xrange(1,20):
        if ss_dict.has_key(i-j):
            if anchor == ss_dict[i-j]:
                pos +=1
            else:
                return pos
        else:
            return pos
    return pos


def get_contact_vectors(structure, sec_struct,sec_struct_types, shift_mat):
    all_contacts = []
    for i in xrange(9, structure.get_number_of_residues()+1-9):
        for j in xrange(i+1, structure.get_number_of_residues()+1-9):
            if abs(i-j) >= 12:
                distance = structure.get_contact_map().get_mapped_distance(i,j)
                if distance <= 8.0:
                    sec_lower = sec_struct.ss_dict[i]
                    sec_upper = sec_struct.ss_dict[j]
                    if sec_lower == sec_struct_types[0] and sec_upper == sec_struct_types[1]:
                        contact_vector = []
                        for i_shift, j_shift in shift_mat:
                            if abs((i+i_shift)-(j+j_shift)) >= 12:
                                dist_shift = structure.get_contact_map().get_mapped_distance(i+i_shift, j+j_shift)

                                if dist_shift <= 8.0:
                                    contact_vector.append(1.0)
                                else:
                                    contact_vector.append(0.0)
                            else:
                                contact_vector.append(0.0)

                        all_contacts.append(contact_vector)
    return np.array(all_contacts)


def add_contacts( structure,  sec_struct_pair_types, shift_mat, sec_struct,sol):
    all_contacts = 0
    for i in xrange(9, structure.get_number_of_residues()+1-9):
        for j in xrange(i+1, structure.get_number_of_residues()+1-9):
            if abs(i-j) >= 12:
                distance = structure.get_contact_map().get_mapped_distance(i,j)
                if distance <= 8.0:
                    sec_lower = sec_struct[i]
                    sec_upper = sec_struct[j]
                    all_shifts = sec_struct_pair_types[(sec_lower,sec_upper)][0]
                    for i_shift, j_shift in shift_mat:
                        if abs((i+i_shift)-(j+j_shift)) >= 1:
                            dist_shift = structure.get_contact_map().get_mapped_distance(i+i_shift,j+j_shift)
                            if dist_shift <= 8.5:
                                if all_shifts.has_key((i_shift, j_shift)):
                                    if dist_shift <= 8.0:
                                        all_shifts[(i_shift, j_shift)] = all_shifts[(i_shift, j_shift)] + 1.0
                                    elif dist_shift <= 50.0:
                                        all_shifts[(i_shift, j_shift)] = all_shifts[(i_shift, j_shift)] + np.exp(-1.0* ((dist_shift-8.0)**2/0.05))
                                else:
                                    all_shifts[(i_shift, j_shift)] = 1.0
                    all_contacts += 1
                    sec_struct_pair_types[(sec_lower,sec_upper)] = ( all_shifts, all_contacts )
    return all_contacts


def clean_sec_structs(sec_struct):
    for i in xrange(2,len(sec_struct)-1):
        if sec_struct[i-1] == sec_struct[i+1]:
            if sec_struct[i-1] == "H" or sec_struct[i-1] == "E":
                if sec_struct[i] == 'C':
                    sec_struct[i] = sec_struct[i-1]


def load_contact_maps(sec_struct_type):
    all_c = []
    shift_mat = shift_matrix()
    file = open(options.pdb_id_list)
    all_contacts = 0
    for line in file:
        pdb_id = str(line).strip().split()[0][0:5]

        pdb_file = "/scratch/schneider/pdb_select_dataset/%s/%s.pdb"%(pdb_id[0:4],pdb_id)
        tmp_struct = StructureContainer.StructureContainer()
        sec_struct = ResidueFeatureSecStruct.ResidueFeatureSecStruct(pdb_file)

        print pdb_id
        try:
            tmp_struct.load_structure('xxxx', pdb_id[-1],pdb_file, seqsep =1)
        except:
            tmp_struct.load_structure('xxxx', ' ',pdb_file, seqsep =1)

        """Buried features"""
        c = get_contact_vectors(tmp_struct, sec_struct,sec_struct_type, shift_mat)
        for i in c:
            all_c.append(i)
    c = numpy.array(all_c)
    file.close()
    return c


def main():
    
    """Generic main function. Executes main functionality of program
    """
    shift_mat = shift_matrix()
    sec_struct_types = ["H","C","E"]
    sec_struct_pair_types = {}
    for i in xrange(0,len(sec_struct_types)):
        for j in xrange(0,len(sec_struct_types)):
            all_shifts = {}
            all_contacts = 0
            sec_struct_pair_types[(sec_struct_types[i], sec_struct_types[j])] = (all_shifts, all_contacts)

    file = open(options.pdb_id_list)
    all_contacts = 0
    for line in file:
        pdb_id = str(line).strip().split()[0][0:5]

        pdb_file = "/scratch/schneider/pdb_select_dataset/%s/%s.pdb"%(pdb_id[0:4],pdb_id)

        tmp_struct = StructureContainer.StructureContainer()
        sec_struct = ResidueFeatureSecStruct.ResidueFeatureSecStruct(pdb_file)
        sec_struct = sec_struct.ss_dict
        clean_sec_structs(sec_struct)
        print pdb_id
        try:
            tmp_struct.load_structure('xxxx', pdb_id[-1],pdb_file, seqsep =1)
        except:
            tmp_struct.load_structure('xxxx', ' ',pdb_file, seqsep =1)
        all_contacts += add_contacts(tmp_struct, sec_struct_pair_types, shift_mat, sec_struct, None )

    for keys, values in sec_struct_pair_types.iteritems():
        all_shifts = values[0]

        sum_prob = 0
        for shifts, counts in all_shifts.iteritems():
            all_shifts[shifts] = counts / float(all_contacts)
            sum_prob += counts / float(all_contacts)

        for shifts, counts in all_shifts.iteritems():
            all_shifts[shifts] = counts / float(sum_prob)
        sec_struct_pair_types[keys] = all_shifts

    for keys, values in sec_struct_pair_types.iteritems():
        print keys, values

    cPickle.dump(sec_struct_pair_types, open("shifts.p", "wb"), protocol=2)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  
