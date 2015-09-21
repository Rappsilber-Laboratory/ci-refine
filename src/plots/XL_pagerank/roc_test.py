import matplotlib
#matplotlib.use('GTKCairo')
#matplotlib.use('wxAgg')
import numpy as np
import matplotlib.pyplot as pylab
from sklearn.metrics import roc_curve, auc
from matplotlib.ticker import NullFormatter
import sys

sys.path.append("/home/michael/scripts/")
sys.path.append("/scratch/schneider/old_scratch/scripts")



f_size = 9

#pylab.rcParams.update({'font.size': f_size, 'font.name':'Arial'})
pylab.rcParams.update({'font.size': f_size,
                       #'font.name':'Arial',
                       'xtick.major.size':2,
                       'xtick.major.width':0.1,
                       'ytick.major.size':2,
                       'ytick.major.width':0.1,
                       'xtick.minor.size':0,
                       'xtick.minor.width':0.0,
                       'ytick.minor.size':0,
                       'ytick.minor.width':0.0})

def load_XL_result_data(result_file, seq_sep=12, filter_missing=True):
    file = open(result_file)
    y_true = []
    y_predict = []
    y_distances = []
    for line in file:
        strline = str(line).strip().split()
        upper_bound = float(strline[-2])
        real_distance = float(strline[-1])
        score = float(strline[0])
        lower_res = int(strline[1])
        upper_res = int(strline[2])
        #print real_distance
        if abs(upper_res-lower_res) >= seq_sep and real_distance <= 998.9:
            if real_distance <= upper_bound:
                y_true.append(1)
            else:
                y_true.append(0)
            y_distances.append(real_distance)
            y_predict.append(score)
    file.close()
    return y_true, y_predict, y_distances

def accuracy(y_true, y_predict, n_top):
    count = [y for y in y_true[:n_top]]
    return float(np.sum(count))/float(n_top)

    #for y in y_true:
    #    if y ==1


#gdt_old = [0.6, 0.78, 0.65, 0.36, 0.45, 0.31] ### old TJ vs Michael (TJ's implementation comparison) 
#gdt_new = [0.63, 0.96, 0.56, 0.62, 0.58, 0.29]

#gdt_old = [0.53, 0.95, 0.72, 0.75, 0.57, 0.40, 0.62, 0.75, 0.68] ### Michael implementation vs TJ thesis  
#gdt_new = [0.54, 0.96, 0.72, 0.73, 0.70, 0.41, 0.53, 0.63, 0.83]

#energy_old = [ -165.7, -115.7, -150.3, -187.2, -131.96, -202.4]
#energy_new = [ -165.5, -119.6, -150.2, -190.2, -136.2, -205.95]  
#data = GENERAL.file2list("analysis_results/HSA_10Perc_PSM_distances.txt")
#data2 = GENERAL.file2list("analysis_results/HSA_10Perc_PR_distances.txt")

protein = "HSA"

#y_true, y_predict, y_distances = load_XL_result_data("analysis_results/%s_20Perc_PSM_distances.txt"%protein)
#y_true_2, y_predict_2, y_distances_2 = load_XL_result_data("analysis_results/%s_20Perc_PR_distances.txt"%protein)
y_true, y_predict, y_distances = load_XL_result_data("analysis_results/%s_0.25L_1_PSM_distances.txt"%protein)
y_true_2, y_predict_2, y_distances_2 = load_XL_result_data("analysis_results/%s_0.25L_1_PR_distances.txt"%protein)
fpr, tpr, _ = roc_curve(y_true, y_predict)
fpr_2, tpr_2, _ = roc_curve(y_true_2, y_predict_2)
n_top = int(len(y_true)*0.5)
print n_top
print "Accuracy, AUC"
print accuracy(y_true, y_predict, int(len(y_true)*0.5)), accuracy(y_true, y_predict, int(len(y_true)*0.8)), auc(fpr, tpr)
print accuracy(y_true_2, y_predict_2, int(len(y_true)*0.5)), accuracy(y_true, y_predict, int(len(y_true)*0.8)), auc(fpr_2, tpr_2)
#print np.mean(y_distances[:358]), np.mean(y_distances_2[:358])
fig = matplotlib.pyplot.gcf()

#y_true = [0 for d in data if (float(d.split()[-1]) <= float(d.split()[-1]))]
#print y_true
#print y_predict

#y = [float(d.split()[1]) for d in data if len(d.split()) > 1 ]

#x2 = [float(d.split()[0]) for d in data2 if len(d.split()) > 1 ]

#y2 = [float(d.split()[1]) for d in data2 if len(d.split()) > 1 ]
pylab.plot(fpr,tpr, '#3B4CC0' , linewidth=1.0)
pylab.plot(fpr_2,tpr_2,'#B40426' , linewidth=1.0)

pylab.ylabel("True positive rate")
pylab.xlabel("False positive rate")
#plot([-1000,1],[-1000,1])
#fig.set_size_inches(4.0,3.0)
#print fig.patches
#pylab.ylim((-0.005,1.01))

#pylab.legend(('Pagerank', 'PSM'),loc='lower right')

ax = pylab.gca()
ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 0.5, 1])
ax.set_xticklabels(["0", "0.5", "1"])
ax.set_yticklabels(["0", "0.5", "1"])
[i.set_linewidth(0.4) for i in ax.spines.itervalues()]
#ax[0].legend(loc='lower right', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)


#pylab.savefig("rmsd_score_test.svg",bbox_inches="tight",pad_inches=0.02,dpi=300)
pylab.plot([0, 100], [0, 100], '-k', linewidth = 0.3)
pylab.xlim((0,1))

fig.set_size_inches(2.0,2.0)

pylab.ylim((0,1.00))
pylab.savefig("%s_20Perc_ROC_test_test.svg"%protein,bbox_inches="tight",pad_inches=0.02,dpi=300)


#pylab.show()
