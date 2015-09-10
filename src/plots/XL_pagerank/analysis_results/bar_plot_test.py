__author__ = 'knurps'
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import matplotlib

f_size = 9

#pylab.rcParams.update({'font.size': f_size, 'font.name':'Arial'})
plt.rcParams.update({'font.size': f_size,
                       #'font.name':'Arial',
                       'xtick.major.size':2,
                       'xtick.major.width':0.1,
                       'ytick.major.size':2,
                       'ytick.major.width':0.1,
                       'xtick.minor.size':0,
                       'xtick.minor.width':0.0,
                       'ytick.minor.size':0,
                       'ytick.minor.width':0.0})


def load_file(data_file, psm=1):
    x = []
    y1 = []
    y2 = []
    counter = 0
    file = open(data_file)
    average_y1 = []
    average_y2 = []
    for line in file:
        strline = str(line).strip().split()
        protein = strline[0]
        tag = strline[0].split('_')[0]
        acc_05 = float(strline[1])
        auc = float(strline[-1])
        #f psm == 1:

        if protein.split('_')[-1] == 'PSM' and protein.split('_')[1][0:2] == '20':
            print protein,
            x.append(tag + ' PSM')
            y1.append(acc_05)
            y2.append(auc)
            average_y1.append(acc_05)
            #average_y2.append(auc)
            print acc_05,
        #else:

        if protein.split('_')[-1] == 'PR' and protein.split('_')[1][0:2] == '20':
            #print protein
            x.append(tag + ' PR')
            y1.append(acc_05)
            print acc_05
            y2.append(auc)
            x.append("free" + str(counter))
            y1.append(0.0)
            y2.append(0.0)
            #average_y1.append(acc_05)
            average_y2.append(acc_05)
            counter +=1

    x.append('Mean PSM')
    x.append('Mean PR')
    print average_y1
    print np.mean(average_y1)
    print average_y2
    print np.mean(average_y2)
    y2.append(np.mean(average_y1))
    y2.append(np.mean(average_y2))
    file.close()
    return x, y1, y2

fig = matplotlib.pyplot.gcf()

#data = np.array([0.809, 0.539,  0.75,  0.907])

#data_random = load_distance_data("../Tx812_random.txt")
#data=  load_distance_data("../Tx812_fit.txt")
sns.set(style="white", context="paper")

#f, axes = plt.subplots(1, 1, figsize=(2.0, 3.0), sharex=True)
#y1 = np.arange(1, 5)
#ax = sns.barplot(y1, data, palette="BuGn_d", ax=axes)
#sns.set_color_codes("pastel")

palette = ['#3B4CC0', '#B40426','black'] * 9
print palette

dict = {'linewidth':0}
x, y1, y2 = load_file("pr_clms_data_auc", psm=1)
print x, y2
print y2
ax = sns.barplot(y1, x, orient='h', palette=palette, **dict)

#x, y1, y2 = load_file("pr_clms_data", psm=0)
#print x, y1
#ax = sns.barplot(y2, x, orient='h', color='#B40426', alpha=0.5)

sns.despine()
#plt.setp(f.axes)
plt.tight_layout(h_pad=3)
fig.set_size_inches(2.5,3.3)
ax.set(xlim=(0.4, 0.9))

#plt.savefig("clms_bar_AUC_LR.svg",bbox_inches='tight',pad_inches=0.02,dpi=300)
plt.show()
#>>> import seaborn as sns
#>>> sns.set_style("whitegrid")
#>>> tips = sns.load_dataset("tips")
#>>> ax = sns.barplot(x="day", y="total_bill", data=tips)