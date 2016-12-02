import matplotlib.pyplot as plt
import os
import seaborn as sns


def set_plot_parameters():
    f_size = 8
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


def plot_fdr(x, y, x2, y2, out_folder="./"):
    sns.set(style="white", palette="muted")
    fig = plt.gcf()

    plt.plot(x, y, '#3B4CC0' , linewidth=1.0)
    plt.plot(x2, y2,'#B40426' , linewidth=1.0)

    plt.ylabel("Number of cross-links")
    plt.xlabel("False discovery rate")
    ax = plt.gca()

    ax.set_xticks([1, 5, 10])
    ax.set_xticklabels(["%s" % i for i in [1, 5, 10]])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]

    max_20_y = max([b for a, b in zip(x, y) if a <= 10])
    max_20_y2 = max([b2 for a2, b2 in zip(x2, y2) if a2 <= 10])

    plt.ylim((0,max(max_20_y, max_20_y2)+100))
    fig.set_size_inches(2.0, 2.0)
    plt.xlim((0, 10))
    if not os.path.exists("%s/raster/" % out_folder):
        os.mkdir("%s/raster/" % out_folder)
        os.mkdir("%s/vector" % out_folder)

    plt.savefig("%s/vector/psm_scoring.svg" % out_folder, dpi=300)
    plt.savefig("%s/raster/psm_scoring.pdf" % out_folder, dpi=300)
    plt.savefig("%s/raster/psm_scoring.png" % out_folder, dpi=300)
    plt.show()


def plot_input_fdr_dependency(x, y, out_folder="./"):

    sns.set(style="white", palette="muted")
    fig = plt.gcf()

    plt.plot(x, y, '#B40426', linewidth=1.0)
    plt.plot((6, 20), (0, 0), '0.3', linewidth=1.0)
    plt.ylabel("Percent increase of links at 5% FDR")
    plt.xlabel("Input FDR")
    ax = plt.gca()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]

    fig.set_size_inches(2.0, 2.0)

    if not os.path.exists("%s/raster/" % out_folder):
        os.mkdir("%s/raster/" % out_folder)
        os.mkdir("%s/vector" % out_folder)

    plt.savefig("%s/vector/input_fdr_dependency.svg" % out_folder, dpi=300)
    plt.savefig("%s/raster/input_fdr_dependency.pdf" % out_folder, dpi=300)
    plt.savefig("%s/raster/input_fdr_dependency.png" % out_folder, dpi=300)
    plt.show()


def plot_score_distribution(y, y2, out_folder, out_name='score_distribution'):

    sns.set(style="white", palette="muted")

    # Set up the matplotlib figure
    f, ax = plt.subplots(1, 1, figsize=(2.0, 2.0), sharex=True)
    sns.despine()

    # Plot a simple histogram with binsize determined automatically
    sns.distplot(y, hist=False, rug=False, color="b", ax=ax, kde_kws={"lw": 1})

    # Plot a kernel density estimate and rug plot
    sns.distplot(y2, hist=False, rug=False, color="r", ax=ax, kde_kws={"lw": 1})
    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]
    plt.xlim((-0.1, 1.1))

    plt.ylabel("Score")
    plt.xlabel("Counts")

    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
        os.mkdir("%s/raster/" % out_folder)
        os.mkdir("%s/vector" % out_folder)

    plt.savefig("%s/vector/%s.svg" % (out_folder, out_name), bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.savefig("%s/raster/%s.pdf" % (out_folder, out_name), bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.savefig("%s/raster/%s.png" % (out_folder, out_name), bbox_inches="tight", pad_inches=0.02, dpi=300)


def plot_spectrum(spectrum):
    sns.set(style="white", palette="muted")

    # Set up the matplotlib figure
    f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
    sns.despine()
    # Get vector data from spectrum object
    peak_list = spectrum.getPeaks()
    x = peak_list[:, 0] # m/z
    y = peak_list[:, 1] # intensities
    # Plot spectrum
    kwargs = {'linewidth': 0.5}
    plt.vlines(x, ymin=0, ymax=y, **kwargs)
    # Thinner lines for aesthetic value
    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]
    # Axis naming
    plt.ylabel("Intensity")
    plt.xlabel("m/z")


def plot_peaks(x, y):
    sns.set(style="white", palette="muted")
    # Set up the matplotlib figure
    f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
    sns.despine()
    # Plot spectrum
    kwargs = {'linewidth': 0.5}
    plt.vlines(x, ymin=0, ymax=y, **kwargs)
    # Thinner lines for aesthetic value
    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]
    # Axis naming
    plt.ylabel("Intensity")
    plt.xlabel("m/z")
    plt.show()


def plot_psm(spectrum, theoretical_spectrum):

    sns.set(style="white", palette="muted")

    # Set up the matplotlib figure
    f, ax = plt.subplots(1, 1, figsize=(3.0, 3.0), sharex=True)
    sns.despine()
    # Get vector data from spectrum object
    peak_list = spectrum.getPeaks()
    x = peak_list[:, 0] # m/z
    y = peak_list[:, 1] # intensities
    # Get vector data from theoretical spectrum
    x_theo = []
    y_theo = []
    for frag in theoretical_spectrum:
        x_theo.append(frag[1])
        y_theo.append(max(y)/2.0)
    # Plot spectrum
    kwargs = {'linewidth': 0.5}
    plt.vlines(x, ymin=0, ymax=y, **kwargs)
    plt.vlines(x_theo, ymin=0, ymax=y_theo, color='r', **kwargs)
    # Thinner lines for aesthetic value
    [i.set_linewidth(0.4) for i in ax.spines.itervalues()]
    # Axis naming
    plt.ylabel("Intensity")
    plt.xlabel("m/z")
    plt.show()
