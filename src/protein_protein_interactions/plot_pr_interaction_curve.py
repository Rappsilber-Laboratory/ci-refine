import os
import sys
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib

sns.set(style="white", context="paper")


def set_plot_parameters():
    """
    Set parameters of matplotlib plotting for font style and size

    Returns
    -------
    Nothing
    """
    f_size = 7
    matplotlib.rcParams.update({'font.size': f_size,
                                #'font.name':'Arial',
                                'xtick.major.size':2,
                        	'xtick.major.width':0.1,
                        	'ytick.major.size':2,
                        	'ytick.major.width':0.1,
                        	'xtick.minor.size':0,
                        	'xtick.minor.width':0.0,
                        	'ytick.minor.size':0,
                        	'ytick.minor.width':0.0,
				'axes.labelsize': f_size+1,
         			'axes.titlesize':6,
         			'xtick.labelsize':f_size,
         			'ytick.labelsize':f_size})

#fig, ax = plt.subplots()
set_plot_parameters()

fig, ax = plt.subplots()

with open("orig_data.txt", 'rb') as infile:
    orig_data = infile.readlines()
    orig_data_num_interactions = [int(line.split()[2]) for line in orig_data]
    orig_data_precision = [float(line.split()[1]) for line in orig_data]

with open("pr_data.txt", 'rb') as infile:
    pr_data = infile.readlines()
    pr_data_num_interactions = [int(line.split()[2]) for line in pr_data]
    pr_data_precision = [float(line.split()[1]) for line in pr_data]

plt.scatter([0.53], [7209], 20, 'k')
plt.xlim((0.5, 0.75))
plt.plot(orig_data_precision, orig_data_num_interactions)
plt.plot(pr_data_precision, pr_data_num_interactions)
sns.despine()

plt.tight_layout(h_pad=3)
fig.set_size_inches(1.2,1.2)

[i.set_linewidth(0.4) for i in ax.spines.itervalues()]

plt.ylabel("Number of interactions")
plt.xlabel("Precision")


plt.savefig("number_of_interactions.svg",bbox_inches="tight",pad_inches=0.02,dpi=300)
plt.savefig("number_of_interactions.png",bbox_inches="tight",pad_inches=0.02,dpi=300)
plt.savefig("number_of_interactions.pdf",bbox_inches="tight",pad_inches=0.02,dpi=300)


#print orig_data


