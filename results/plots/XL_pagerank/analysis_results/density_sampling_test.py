__author__ = 'knurps'
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import matplotlib

f_size = 12

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


def load_file(data_file):
    y = []
    file = open(data_file)
    for line in file:
        strline = str(line).strip().split()
        y.append(float(strline[-1])*100)
    file.close()
    return y

fig = matplotlib.pyplot.gcf()
sns.set(style="white", context="paper")
ax = fig.add_subplot(111)
#data = np.array([0.809, 0.539,  0.75,  0.907])

#data_random = load_distance_data("../Tx812_random.txt")
#data=  load_distance_data("../Tx812_fit.txt")


#f, axes = plt.subplots(1, 1, figsize=(2.0, 3.0), sharex=True)
#y1 = np.arange(1, 5)
#ax = sns.barplot(y1, data, palette="BuGn_d", ax=axes)
#sns.set_color_codes("pastel")

palette = ['#3B4CC0', '#B40426','black'] * 9
print palette

dict = {'linewidth':0}
y = load_file("hsa_density_psm.txt")
y2 =  load_file("hsa_density_pr.txt")
z = np.array(y2) - np.array(y)
print z
counter = 0
sample_set = []
all_y = []
all_yerr = []
for i in z:
    sample_set.append(i)
    counter += 1
    print i
    if counter == 9:
        all_y.append(np.mean(sample_set))
        all_yerr.append(np.std(sample_set))
        counter = 0
        sample_set = []

x = np.arange(0, len(all_y), 1)

plt.errorbar(x, all_y, yerr=all_yerr, fmt='o', color='#B40426', linewidth = 0.8)
plt.axhline(y=11.3, xmin=-1, xmax=7,color='black', linestyle='--', linewidth=0.8)
#ax = sns.barplot(y1, x, orient='h', palette=palette, **dict)

#x, y1, y2 = load_file("pr_clms_data", psm=0)
#print x, y1
#ax = sns.barplot(y2, x, orient='h', color='#B40426', alpha=0.5)

sns.despine()
#plt.setp(f.axes)
plt.tight_layout(h_pad=3)
fig.set_size_inches(4.0,2.5)
ax.set_xticks([0, 1, 2,3,4,5,6])
ax.set(xlim=(-.25, 6.25))
ax.set(ylim=(0, 14))
ax.set_xticklabels(["0.25L", "0.5L", "0.75L", "1L", "1.25L", "1.5L", "1.75L"])

plt.savefig("density_improvement.svg",bbox_inches='tight',pad_inches=0.02,dpi=300)
plt.show()

