__author__ = 'knurps'
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import random
sns.set(style="white", context="paper")
gammas = sns.load_dataset("gammas")
#print gammas
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

fig = matplotlib.pyplot.gcf()
palette = ['#3B4CC0', '#B40426']

df = pd.read_csv("epc_all_cuts_1.5L_0.4_2.0.txt")
#print df
#print gammas
ax = sns.tsplot(data=df, time="top", unit="protein",
           condition="algorithm", value="acc", ci=5, color=palette)

sns.despine()
#plt.setp(f.axes)
plt.tight_layout(h_pad=3)
fig.set_size_inches(3.3,2.3)
#ax.set(xlim=(0.4, 0.9))
ax.set_xticks([0.25, 0.5, 0.75, 1,1.25, 1.5 ])
ax.set_yticks([0.3, 0.4, 0.5,0.6])
ax.set_xticklabels(["0.25L", "0.5L","0.75L", "1L","1.25L", "1.5L"])
plt.ylabel("Accuracy")
plt.xlabel("Number of top contacts")
#ax.set_yticklabels(["0", "0.5", "1"])
[i.set_linewidth(0.4) for i in ax.spines.itervalues()]

plt.savefig("pagerank_improvement.svg",bbox_inches="tight",pad_inches=0.02,dpi=300)
#plt.show()