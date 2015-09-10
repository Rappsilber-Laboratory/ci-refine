
import cPickle
import numpy as np 
shift_dict = cPickle.load(open( "shifts.p", "rb" ))

#print shift_dict[("H", "H")]

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

def to_matrix(shift_dict):
    matrix = np.zeros((17, 17))
    
    print len(shift_dict)
    for keys, values in shift_dict.iteritems():
        #print keys
        matrix[keys[0]+8][keys[1]+8] = values
        print keys,keys[0]+8, keys[1]+8, matrix[keys[0]+8][keys[1]+8]
    return matrix

ele = ["H", "E", "C"]
s = shift_matrix()
print s

for e in ele:
    for f in ele:
        print e, f
        mat = to_matrix(shift_dict[(e, f)])
        print mat.shape
        np.savetxt("%s_%s_matrix.txt"%(e,f), mat)
        

