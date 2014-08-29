import numpy as np 
import os 
import scipy.stats
#file = open("ciss_test_cali_ciss_gremlin.txt")
#file = open("d329_cali_ciss_gremlin.txt")

#file = open("d329_gremlin_new.txt")
#file = open("dec_func.txt")
#file = open("all_gremlin_new.txt")
#file = open("svmcon_test_gremlin_new.txt")
#file = open("svmcon_test_cali_ciss_gremlin.txt")
#file = open("casp9_10_test_gremlin_new.txt")
#file = open("all_gremlin_new.txt")
#file = open("all_cali_ciss_gremlin.txt")
#file = open("all_ciss_gremlin.txt")
#file = open("casp9_10_test_cali_ciss_gremlin.txt")
#file = open("d329_dec_func_rescaled.txt")
#file = open("ciss_test_cali_ciss_gremlin.txt")
#file = open("svmcon_test_dec_func_rescaled.txt")
#file = open("casp9_10_dec_func_rescaled.txt")
#file = open("ciss_test_cache.txt")
#file = open("ciss_control_control_gremlin.txt")
#file = open("train_tuned_ciss_gremlin.txt")
#file = open("d329_train_tuned_ciss_gremlin.txt.txt")
#file = open("svmcon_test_train_tuned_ciss_gremlin.txt")
#file = open("casp9_10_train_tuned_ciss_gremlin.txt")
#file = open("all_train_tuned_ciss_gremlin.txt")
#file = open("ciss_test_gremlin_new.txt")
#file = open("all_dec_func_rescaled.txt")

#file = open("train_tuned_count_gremlin.txt")
#file = open("ciss_test_train_tuned_ciss_gremlin.txt")
#file = open("d329_train_tuned_count_gremlin.txt.txt")
#file = open("svmcon_test_train_tuned_count_gremlin.txt")
#file = open("casp9_10_train_tuned_count_gremlin.txt")
#file = open("all_train_tuned_count_gremlin.txt")
#file = open("train_tuned_count_gremlin_mid_long.txt")
#file = open("train_tuned_ciss_gremlin_mid_long.txt")
#file = open("d329_train_tuned_ciss_gremlin_mid_long.txt.txt")
#file = open("d329_train_tuned_count_gremlin_mid_long.txt.txt")
#file = open("all_count_gremlin_mid_long.txt")
#file = open("all_ciss_gremlin_mid_long.txt")
#file = open("ciss_set_training_tuned_less_1.txt")
#file = open("d329_training_tuned_less_1.txt")
#file = open("casp9_10_gremlin_calibrated.txt")
#file = open("all_new_gremlin_calibrated.txt")
#file = open("d329_gremlin_calibrated.txt")
#file = open("all_high_res_calib.txt")
#file = open("all_no_gremlin.txt")

#file = open("all_new_gremlin_calibrated.txt") ## LATEST, IMPORTANT?

#file = open("ciss_set_over.txt")
##file = open("epc_pre_casp10_casp10.txt")
#file = open("gremlin_pre_casp10_casp10.txt")
#file = open("epc_pre_casp_casp10.txt")
#file = open("group222_pre_casp10_casp10.txt")
#file = open("epc_graph_cutoff_big.txt")
#file = open("epc_pre_casp_logreg.txt")
#file = open("group222_pre_casp10_casp10.txt")
#file = open("group222_pre_casp10_few_seq.txt")
#file = open("epc_pre_casp_casp10_few_seq.txt")
#file = open("group222_pre_casp10_casp10_hard.txt")
#file = open("epc_only_pre_casp_casp10_few_seq.txt")
#file = open("epc_pre_casp_sample.txt")
#file = open("group222_at_least.txt")
#file = open("epc_pre_casp_all_casp10.txt")

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

acc_10_mid = []
acc_5_mid = []
acc_2_mid = []

cov_10 = []
cov_5 = []
cov_2 = []

cov_10_mid = []
cov_5_mid = []
cov_2_mid = []


counter = 0
for line in file:
    strline = str(line).strip().split()

    #acc_10.append(float(strline[2]))
    #if (float(strline[-2])/ float(strline[-3])) <= 1.0  and (float(strline[-2])/ float(strline[-3])) > 0 and float(strline[-3]) < 9999999:
    #if float(strline[-1]) > 150 and float(strline[-1]) < 200:

    if 1==1:
    #acc_10_mid.append(float(strline[]))
        acc_10.append(float(strline[2]))
        cov_10.append(float(strline[3]))
        acc_5.append(float(strline[8]))
        cov_5.append(float(strline[9]))
        acc_2.append(float(strline[14]))
        cov_2.append(float(strline[15]))

        acc_10_mid.append(float(strline[2+24]))
        cov_10_mid.append(float(strline[3+24]))
        acc_5_mid.append(float(strline[8+24]))
        cov_5_mid.append(float(strline[9+24]))
        acc_2_mid.append(float(strline[14+24]))
        cov_2_mid.append(float(strline[15+24]))
        counter +=1
print counter

good_count = np.sum([1 for a in acc_5 if a <=0.20])
#print good_count
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

#print "& Long  & 0.335(0.023)/0.038 & 0.278(0.019)/0.062 & 0.205(0.014)/0.110 \\"

#print "%.3f" % np.mean(acc_10), "%.3f" % scipy.stats.sem(acc_10), "%.3f" % np.mean(cov_10), "%.3f" % np.mean(acc_10_mid), "%.3f" % scipy.stats.sem(acc_10_mid) ,"%.3f" % np.mean(cov_10_mid)
#print "%.3f" %np.mean(acc_5),"%.3f" % scipy.stats.sem(acc_5), "%.3f" %np.mean(cov_5), "%.3f" %np.mean(acc_5_mid),"%.3f" % scipy.stats.sem(acc_5_mid) , "%.3f" %np.mean(cov_5_mid)
#print "%.3f" %np.mean(acc_2),"%.3f" % scipy.stats.sem(acc_2), "%.3f" %np.mean(cov_2), "%.3f" %np.mean(acc_2_mid),"%.3f" % scipy.stats.sem(acc_2_mid) , "%.3f" %np.mean(cov_2_mid)
print 
#print "%.3f" % np.median(acc_10), "%.3f" % scipy.stats.sem(acc_10), "%.3f" % np.mean(cov_10), "%.3f" % np.median(acc_10_mid), "%.3f" % scipy.stats.sem(acc_10_mid) ,"%.3f" % np.median(cov_10_mid)
#print "%.3f" %np.median(acc_5),"%.3f" % scipy.stats.sem(acc_5), "%.3f" %np.mean(cov_5), "%.3f" %np.median(acc_5_mid),"%.3f" % scipy.stats.sem(acc_5_mid) , "%.3f" %np.median(cov_5_mid)
#print "%.3f" %np.median(acc_2),"%.3f" % scipy.stats.sem(acc_2), "%.3f" %np.mean(cov_2), "%.3f" %np.median(acc_2_mid),"%.3f" % scipy.stats.sem(acc_2_mid) , "%.3f" %np.median(cov_2_mid)

