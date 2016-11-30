__author__ = 'knurps'
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import random
sns.set(style="white", context="paper")
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

#x = {(0,1):2, (0, 1):1}
ss = ["H","E","C"]
#for e in ss:
#    for f in ss:
e = "H"
f = "H"
fig = matplotlib.pyplot.gcf()
palette = ['#B40426', '#3B4CC0']
mat = np.genfromtxt("%s_%s_matrix.txt"%(e,f))
plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')

ax = sns.heatmap(mat, square=True, cmap="coolwarm", linewidths=0.2)
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout(h_pad=3)
fig.set_size_inches(2.0,2.0)
plt.savefig("%s_%s_matrix_plot.svg"%(e,f),bbox_inches='tight',pad_inches=0.02,dpi=300)
#ax.set_xticklabels([)
#ax.set_yticklabels(["0", "0.5", "1"])
#plt.show()
