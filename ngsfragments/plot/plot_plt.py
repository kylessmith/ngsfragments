import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.signal import argrelmax
from projectframe import ProjectFrame

# Tell matplot lib to use 'editable text'
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Local imports
from ..correct.correction import gaussian_smooth


def plot_table(pf: ProjectFrame,
               key: str,
               font_size: int = 14,
               show: bool = True,
               save: str | None = None,
               ax: matplotlib.axes._axes.Axes = None) -> matplotlib.axes._axes.Axes:
    """
    Plot table of summary metrics

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame object
        key : str
            Sample key
        font_size : int
            Font size
        show : bool
            Whether to show plot
        save : str
            Path to save plot
        ax : Axes
            Axes object to plot on

    Returns
    -------
        ax : Axes
            Axes object
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)

    # Grab metrics
    sam_fn = key
    genome = str(pf.params["ref_genome"])
    bin_size = int(pf.params["cnv_binsize"])
    purity = pf.anno.loc[key,"purity"]
    ploidy = pf.anno.loc[key,"ploidy"]
    clonal = pf.anno.loc[key,"clonal"]
    n_frags = int(pf.anno.loc[key,"n_fragments"])
    median_var = np.median(pf.obs_intervals[key]["cnv_segments"].loc[:,"var"].values)


    # Record text
    labels = ["Sample:", "Tumour purity:",
              "Ploidy:", "Clonality:",
              "Median Variance:",
              "Number of fragments:",
              "Genome:", "Bin size:"]
    numbers = [sam_fn, '{:.3f}'.format(purity),
               '{:.3f}'.format(ploidy),
               '{:.3f}'.format(clonal),
               '{:.3f}'.format(median_var),
               str(n_frags), 
               genome,
               str(bin_size)]
    text = list(zip(labels, numbers))

    # Table formatting
    color = "lightgrey"
    blank = "white"
    cellColours = [[color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank],
                   [color, color],
                   [blank, blank]]

    # Plot table
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(text, loc='center', cellColours=cellColours, edges="closed")
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)

    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()

    return ax


def plot_mean_window(pf, obs, center_x=None, title=None, xlabel=None, ylabel=None, show=True, save=None, ax=None):
    """
    """

    # Get scores
    window_scores = pf.obs_values[obs]["tss_wps"].engine.df.mean(axis=0).values
    
    # Creat Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)
    
    # Plot metrics
    ax.plot(window_scores, color="black", linewidth=1.0)
    sns.despine(trim=True)

    # Plot center
    if center_x is not None:
        ax.axvline(x=center_x, linestyle="--", color="grey", linewidth=1.0, dash_capstyle="round")
        ax.set_xticks([0, center_x, len(window_scores)])
        ax.set_xticklabels([-center_x, 0, len(window_scores)-center_x])
    
    # Set labels
    if title is not None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=12)
    
    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()


def plot_genome_signal(df, signal_label, chrom_shift, line=True, plot_median=True, smooth=False, title=None,
                       ylabel=None, show=True, save=None, ax=None):
    """
    Plot a signal across the genome
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(10,5), tight_layout=True)

    # Iterate over chroms
    chrom_index = {chrom:df.loc[:,"seqnames"].values==chrom for chrom in np.unique(df.loc[:,"seqnames"].values)}
    for chrom in chrom_index:
        signal = df.loc[chrom_index[chrom], :]
        # Smooth signal
        if smooth:
            values = gaussian_smooth(signal.loc[:,signal_label].values, scale=10)
        else:
            values = signal.loc[:,signal_label].values
        shift = chrom_shift[chrom]
        if line:
            ax.plot(signal.loc[:,"start"].values+shift, values, linewidth=0.5, color="black")
        else:
            ax.scatter(signal.loc[:,"start"].values+shift, values, s=1, color="black")
        ax.axvline(x=shift, color="grey", linestyle="--")

    # Plot 0 line
    if plot_median:
        ax.axhline(y=np.median(df.loc[:,signal_label]), color="darkgrey", linewidth=0.5)
    
    # Set titles
    if title != None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Chromosomes", fontsize=12)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=12)
    
    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()


def get_chrom_shift(chrom_lengths, chrom_order=None):
    """
    Get the shift for each chromosome

    Parameters
    ----------
        chrom_lengths : dict
            Dictionary of chromosome lengths
        chrom_order : list
            List of chromosome names in order
    
    Returns
    -------
        chrom_shift : dict
            Dictionary of chromosome shifts
    """

    chrom_shift = {}
    shift = 0

    if chrom_order is None:
        chrom_order = list(chrom_lengths.keys())

    for chrom in chrom_order:
        chrom_shift[chrom] = shift
        shift += chrom_lengths[chrom]

    return chrom_shift


def plot_cnv(pf: ProjectFrame,
             obs: str,
             cnv_column: str = "median",
             title: str | None = None,
             show: bool = True,
             save: str | None = None,
             plot_max: int = 5,
             plot_min: int = -3,
             ax: matplotlib.axes._axes.Axes = None):
    """
    Plot the CNVs

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame object
        obs : str
            Observation to plot
        cnv_column : str
            Column to plot
        title : str
            Title of plot
        show : bool
            Whether to show plot
        save : str
            Path to save plot
        plot_max : int
            Maximum value to plot
        plot_min : int
            Minimum value to plot
        ax : Axes
            Axes object to plot on
        
    Returns
    -------
        ax : Axes
            Axes object
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(10,5), tight_layout=True)

    # Assign variables
    cnvs = pf.obs_intervals[obs]["cnv_segments"]
    bins = pf.intervals["cnv_bins"]

    # Determine whether chroms are
    chroms = cnvs.index.unique_labels

    # Get chrom shift
    chrom_shift = get_chrom_shift(pf.uns["chrom_lengths"], chroms)

    # Determine last position
    last_chrom = chroms[-1]
    last_position = chrom_shift[last_chrom]
    last_position += int(bins.loc[last_chrom,:].index.ends[-1]) + 10

    # Iterate over chroms
    for chrom in chroms:
        chrom_bins = bins.loc[chrom, obs]
        shift = chrom_shift[chrom]
        ax.scatter(chrom_bins.index.starts + shift, chrom_bins.values, s=1, color="black")
        ax.axvline(x=shift, color="grey", linestyle="--", linewidth=1.0, solid_capstyle="round")

    # Plot segments
    cnv_types = {"NEUT":"grey", "GAIN":"red", "HETD":"blue", "AMP":"darkred"}
    for i in range(cnvs.shape[0]):
        segment = cnvs.iloc[i,:]
        cnv_type = segment.df.loc[:,"Corrected_Call"].values[0]
        color =  cnv_types[cnv_type] if cnv_type in cnv_types else "darkred"
        start = segment.index.starts[0] + chrom_shift[segment.index.labels[0]]
        end = segment.index.ends[0] + chrom_shift[segment.index.labels[0]]
        ax.plot([start, end],
                  [segment.loc[:,cnv_column].values, segment.loc[:,cnv_column].values],
                  color=color, solid_capstyle="butt", linewidth=2.5)

    # Plot 0 line
    ax.axhline(y=0, color="darkgrey", linewidth=1.5)

    # Plot chromosome names
    text_x = []
    for chrom in chroms:
        text_x.append(chrom_shift[chrom] + (pf.uns["chrom_lengths"][chrom] / 2))
    ax.set_xticks(text_x)
    ax.set_xticklabels(chroms, rotation=90)
    
    # Set titles
    if title != None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Chromosomes", fontsize=12)
    ax.set_ylabel("log2 ratio", fontsize=12)
    ax.set_ylim((plot_min, plot_max))
    ax.set_xlim((0,last_position))

    # Plot arrows if segments of the chart
    off_chart = np.logical_or(cnvs.df.loc[:,"median"].values >= plot_max,
                              cnvs.df.loc[:,"median"].values <= plot_min)
    if np.sum(off_chart) > 0:
        off_segments = cnvs.iloc[off_chart,:]
        for i in range(off_segments.shape[0]):
            segment = off_segments.iloc[i,:]
            cnv_type = segment.df.loc[:,"Corrected_Call"].values[0]
            color =  cnv_types[cnv_type] if cnv_type in cnv_types else "darkred"
            start = segment.starts[0] + chrom_shift[segment.index.labels[0]]
            start = start / last_position # convert to ratio
            end = segment.index.ends[0] + chrom_shift[segment.index.labels[0]]
            end = end / last_position # convert to ratio
            position = (start + end) / 2.0 # take middle
            median = segment.df.loc[:,"median"].values[0]

            if median >= plot_max:
                ax.annotate('', xy=(position, 1.01), xycoords='axes fraction', xytext=(position, 1.075),
                            arrowprops=dict(arrowstyle="simple", color=color))

            else:
                ax.annotate('', xy=(position, -0.01), xycoords='axes fraction', xytext=(position, -0.075),
                            arrowprops=dict(arrowstyle="simple", color=color))
    
    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()

    return ax


def fragment_distribution(pf: ProjectFrame,
                          obs: str,
                          title: str = None,
                          show: bool = True,
                          save: str = None,
                          ax: matplotlib.axes._axes.Axes = None):
    """
    Plot the fragment distribution
    """

    # Assign length dist
    length_dist = pf.values["length_dist"].loc[obs,:].values
    
    # Determine total distribution
    if isinstance(length_dist, dict):
        max_length = max(len(length_dist[chrom]) for chrom in length_dist)
        fragment_lengths = np.zeros(max_length, dtype=int)
        for chrom in length_dist:
            fragment_lengths[:len(length_dist[chrom])] = length_dist[chrom]
    else:
        fragment_lengths = length_dist

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(7,5), tight_layout=True)
    
    # Plot fragment lengths
    plt.plot(fragment_lengths, color="black", linewidth=0.5)
    sns.despine(trim=True, ax=ax)

    # Plot local maxima
    extra_ticks = []
    for maxima in argrelmax(fragment_lengths, order=35)[0]:
        plt.axvline(x=maxima, color="grey", linestyle="--", dash_capstyle="round", linewidth=0.5)
        extra_ticks.append(maxima)

    # Set xticks
    plt.xticks(list(plt.xticks()[0]) + extra_ticks)
    
    # Set titles
    if title != None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Read length", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    
    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=80)
    if show:
        plt.show()


def cnv_summary(pf: ProjectFrame,
                key: str,
                add_wps: bool = False,
                show: bool = True,
                save: str = None):
    """
    Plot summary metrics

    Params
    ------
        metrics
            dict
        show
            bool
        save
            bool
    
    Returns
    -------
        fig
            matplotlib.pyplot.figure
    """

    # Initialize matplotlib grid
    fig = plt.figure(figsize=(10, 7), constrained_layout=True)        
    gs = fig.add_gridspec(3, 3)

    # Plot table
    table_ax = fig.add_subplot(gs[0, 0])
    table_ax = plot_table(pf,
                          key,
                          font_size = 5,
                          show = False,
                          ax = table_ax)

    # Plot fragment length distribution
    with sns.axes_style("ticks"):
        if add_wps:
            len_ax = fig.add_subplot(gs[0, 1])
        else:
            len_ax = fig.add_subplot(gs[0, 1:])
    len_ax = fragment_distribution(pf, key, title=None, show=False, ax=len_ax)

    # Plot WPS
    if add_wps:
        wps_ax = fig.add_subplot(gs[0, 2])
        wps_ax = plot_mean_window(pf, key, center_x = 1000, title=None, show=False, xlabel="Distance to TSS", ylabel="WPS", ax=wps_ax)

    # Plot CNV plot
    with sns.axes_style("ticks"):
        cnv_ax = fig.add_subplot(gs[1:3,:])
    cnv_ax = plot_cnv(pf, key, title=None, show=False, ax=cnv_ax)

    # Save of display plot
    if save != None:
        plt.savefig(save, transparent=True, dpi=300)
    if show:
        plt.show()

    return fig