import os 
import sys
import numpy

weight_pos = []
weight_neg = []
weight_neg_pos = []

file = open(sys.argv[-1])

for line in file:
    strline = str(line).split()
    weight_pos.append(float(strline[1]))
    weight_neg.append(float(strline[2]))
    weight_neg_pos.append(float(strline[3]))


    
file.close()    

print "Pos/Neg:", numpy.mean(weight_pos)/numpy.mean(weight_neg)
print "Neg_Pos/Neg:", numpy.mean(weight_neg_pos)/numpy.mean(weight_neg)
