
import cPickle
import numpy as np 
shift_dict = cPickle.load(open( "/scratch/schneider/projects/pagerank_refinement/src/probabilities/shifts_sigma_0.05.txt", "rb" ))

#print shift_dict[("H", "H")]

def to_matrix(shift_dict):
    matrix = np.zeros((17, 17))
    print len(shift_dict)
    for keys, values in shift_dict.iteritems():
        matrix[keys[0]][keys[1]] = values
    return matrix

ele = ["H", "E", "C"]

for e in ele:
    for f in ele:
        mat = to_matrix(shift_dict[(e, f)])
        print mat.shape
        np.savetxt("%s_%s_matrix.txt"%(e,f), mat)
