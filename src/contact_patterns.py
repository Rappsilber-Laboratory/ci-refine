"""Author: Michael Schneider
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/structure/")
import networkx as nx
import InputOutput
import numpy
from optparse import OptionParser
import pdb
import ResidueFeatureSecStruct
import ResidueFeatureRelSasa
import StructureContainer
import numpy as np
import pickle
import cPickle
from sklearn import cluster, datasets, neighbors
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
    for i in xrange(-9,10):
        row = []
        for j in xrange(-9,10):
            if i == 0 and j == 0:
                pass
            else:
                row.append((i,j))
                matrix.append((i,j))
    #for i in matrix:
        #print i
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
            #shift_vec.append(vec[i])

def get_relative_sec_struct_pos(ss_dict, i):
    anchor = ss_dict[i]
    pos = 1
    #print pos
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


def cluster_shift_maps(shift_map_vector):
    #from sk
    print shift_map_vector[1]
    print shift_map_vector[2]
    dist = neighbors.DistanceMetric.get_metric('rogerstanimoto')
    print dist.pairwise(shift_map_vector[0:10])
    k_means = cluster.Ward(n_clusters = 5)
    labels =  k_means.fit_predict(shift_map_vector)
    n_clusters = 5
    shift_mat = shift_matrix()
    print labels
    from collections import Counter
    res = Counter(labels)
    print res
    #num_list = []
    #for label in xrange(0,50):
    #    num_list.append((list(k_means.labels_).count(label),label))#

    #num_list.sort()
    #num_list.reverse()

    #print num_list
    true = []
    false = []
    count = 0
    shift_dict = {}
    #for dummy,label in num_list[:n_clusters]:
    for label in xrange(0,n_clusters):
        #print list(k_means.labels_).count(label)
        i_0 = numpy.array([0]*shift_map_vector.shape[1])
        for i,pred in zip(shift_map_vector, labels):
            if pred == label:
                i_0 = i_0 + numpy.array(i)

        for i in xrange(0,i_0.shape[0]):
            i_0[i] = i_0[i] / float( list(labels).count(label) )
        sum_prob = numpy.sum([i for i in i_0])

        for i in xrange(0,i_0.shape[0]):
            i_0[i] = i_0[i] / sum_prob
        a = numpy.array(i_0)

        for i,pred in zip(shift_map_vector, labels):
            if pred == label:
                true.append(numpy.dot(i,a))
            else:
                false.append( numpy.dot(i,a))

        all_values = {}
        for shift, val in zip(shift_mat, a):
            all_values[shift] = val
        shift_dict[count] = all_values
        count+=1
    return shift_dict
    #print sec_struct_pair_types
    #print numpy.mean(true)
    #print numpy.mean(false)

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


    #print bur_dict
    shift_mat = shift_matrix()
    #
    sec_struct_types = ["H","C","E"]
    bur_types = ["B","A"]
    seq_sep_types = [(12,24),(24,9999)]
    sec_struct_pair_types = {}
    for i in xrange(0,len(sec_struct_types)):
        for j in xrange(0,len(sec_struct_types)):
            #for k in xrange(0,len(seq_sep_types)):
                #for l in xrange(0,len(seq_sep_types)):
            all_shifts = {}
            all_contacts = 0
            sec_struct_pair_types[(sec_struct_types[i],sec_struct_types[j])] = (all_shifts,all_contacts)

    print sec_struct_pair_types
    #return 0
    new_stuff = {}
    for keys,values in sec_struct_pair_types.iteritems():
        #if keys == ('H','H'):
        c = load_contact_maps(keys)

        shift_dict = cluster_shift_maps(c)
        print shift_dict
        new_stuff[keys] = shift_dict

    print new_stuff
    cPickle.dump(new_stuff, open( "shifts.p", "wb" ),protocol=2 )

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


    for keys, values in sec_struct_pair_types.iteritems():
        all_shifts = values[0]

        sum_prob = 0
        for shifts, counts in all_shifts.iteritems():
            all_shifts[shifts] = counts / float(all_contacts)
            sum_prob += counts / float(all_contacts)

        for shifts, counts in all_shifts.iteritems():
            all_shifts[shifts] = counts / float(sum_prob)
            #sum_prob += counts / float(all_contacts)
        sec_struct_pair_types[keys] = all_shifts


    for keys, values in sec_struct_pair_types.iteritems():
        print keys, values


    """
    for keys, values in all_shifts.iteritems():
        all_shifts[keys] = values / float(all_contacts)
    for keys, values in all_shifts.iteritems():
        print keys, values
    #calculate_probabities(all_shifts)
   # print all_shifts
                        #print tmp_struct.get_contact_map().get_mapped_distance(i+i_shift,j+j_shift)
                    #print
                #distance = tmp_struct.get_contact_map().get_mapped_distance(i,j)
                #print distance
    pickle.dump(all_shifts, open( "shifts.p", "wb" ) )
    """
    cPickle.dump(sec_struct_pair_types, open( "shifts.p", "wb" ),protocol=2 )

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  
