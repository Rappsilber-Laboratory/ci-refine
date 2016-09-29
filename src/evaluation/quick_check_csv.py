import numpy as np 
import os 
import scipy.stats

from optparse import OptionParser

## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("--file", type="string", dest="score_file", help="File with scores")
   # parser.add_option("--", type="string", dest="example", help="An example")
    options, args = parser.parse_args()
    return options, args 
options, args  = add_options( parser )  
file = open(options.score_file)
 

acc_10 = []
acc_5 = []
acc_2 = []
acc_1 = []
acc_1_5 = []

acc_10_mid = []
acc_5_mid = []
acc_2_mid = []
acc_1_mid = []
acc_1_5_mid = []

cov_10 = []
cov_5 = []
cov_2 = []
cov_1 = []
cov_1_5 = []

cov_10_mid = []
cov_5_mid = []
cov_2_mid = []
cov_1_mid = []
cov_1_5_mid = []


counter = 0
for line in file:
    strline = str(line).strip().split()

    #acc_10.append(float(strline[2]))
    #if (float(strline[-2])/ float(strline[-3])) <= 1.0  and (float(strline[-2])/ float(strline[-3])) > 0 and float(strline[-3]) < 9999999:
    #if float(strline[-1]) > 150 and float(strline[-1]) < 200:

    if float(strline[8]) <= 1.1:
        #print strline[0]
    #acc_10_mid.append(float(strline[]))
        acc_10.append(float(strline[2]))
        cov_10.append(float(strline[3]))
        acc_5.append(float(strline[8]))
        cov_5.append(float(strline[9]))
        acc_2.append(float(strline[14]))
        cov_2.append(float(strline[15]))
        acc_1.append(float(strline[20]))
        cov_1.append(float(strline[21]))
        acc_1_5.append(float(strline[26]))
        acc_1_5.append(float(strline[27]))


        acc_10_mid.append(float(strline[2+30]))
        cov_10_mid.append(float(strline[3+30]))
        acc_5_mid.append(float(strline[8+30]))
        cov_5_mid.append(float(strline[9+30]))
        acc_2_mid.append(float(strline[14+30]))
        cov_2_mid.append(float(strline[15+30]))
        acc_1_mid.append(float(strline[20+30]))
        cov_1_mid.append(float(strline[21+30]))
        acc_1_5_mid.append(float(strline[26+30]))
        cov_1_5_mid.append(float(strline[27+30]))

        counter +=1
#    else:
#        print strline
#print counter
"""
good_count = np.sum([1 for a in acc_5 if a <=0.20])
out_string = " ".join(["& Long  &",
                       "%.3f(%.3f)/%.3f &"%(np.mean(acc_10), scipy.stats.sem(acc_10), np.mean(cov_10)),
                       "%.3f(%.3f)/%.3f &"%(np.mean(acc_5), scipy.stats.sem(acc_5), np.mean(cov_5)),
                       "%.3f(%.3f)/%.3f \\"%(np.mean(acc_2), scipy.stats.sem(acc_2), np.mean(cov_2))])

print out_string
out_string = " ".join(["& Medium  &",
                       "%.3f(%.3f)/%.3f &"%(np.mean(acc_10_mid), scipy.stats.sem(acc_10_mid), np.mean(cov_10_mid)),
                       "%.3f(%.3f)/%.3f &"%(np.mean(acc_5_mid), scipy.stats.sem(acc_5_mid), np.mean(cov_5_mid)),
                       "%.3f(%.3f)/%.3f \\"%(np.mean(acc_2_mid), scipy.stats.sem(acc_2_mid), np.mean(cov_2_mid))])
print out_string
print 
"""
out_string = "".join(["%.3f," % np.mean(acc_10),
                       "%.3f," % np.mean(acc_5),
                       "%.3f," % np.mean(acc_2),
                       "%.3f," % np.mean(acc_1),
                       "%.3f," % np.mean(acc_1_5),
                       "%.3f," % np.mean(acc_10_mid),
                       "%.3f," % np.mean(acc_5_mid),
                       "%.3f," % np.mean(acc_2_mid),
                       "%.3f," % np.mean(acc_1_mid),
                       "%.3f" % np.mean(acc_1_5_mid)])
print out_string
