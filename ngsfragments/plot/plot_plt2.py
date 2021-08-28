import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, sem
from scipy.signal import argrelmax

# Tell matplot lib to use 'editable text'
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_cnv(frags, bins, cnvs, cnv_column="median", title=None, show=True, save=None, ax=None):
    """
    Plot the CNVs
    """

    # Create Axes object if not given
    if ax == None:
        with sns.axes_style("ticks"):
            fig, ax = plt.subplots(figsize=(10,5), tight_layout=True)

    # Determine whether chroms are
    chroms = cnvs.index.unique_labels

    # Determine last position
    last_chrom = chroms[-1]
    last_position = frags.chrom_shift[last_chrom]
    last_position += int(bins.loc[last_chrom,:].index.extract_ends()[-1]) + 10

    # Iterate over chroms
    for chrom in chroms:
        chrom_bins = bins.loc[chrom, "ratios"]
        shift = frags.chrom_shift[chrom]
        ax.scatter(chrom_bins.index.extract_starts() + shift, chrom_bins.values, s=1, color="black")
        ax.axvline(x=shift, color="grey", linestyle="--", linewidth=1.0, solid_capstyle="round")

    # Plot segments
    cnv_types = {"NEUT":"grey", "GAIN":"red", "HETD":"blue", "AMP":"darkred"}
    for i in range(cnvs.shape[0]):
        segment = cnvs.iloc[i,:]
        cnv_type = segment.df.loc[:,"Corrected_Call"].values[0]
        color =  cnv_types[cnv_type] if cnv_type in cnv_types else "darkred"
        start = segment.index.extract_starts()[0] + frags.chrom_shift[segment.index.extract_labels()[0]]
        end = segment.index.extract_ends()[0] + frags.chrom_shift[segment.index.extract_labels()[0]]
        ax.plot([start, end],
                  [segment.loc[:,"median"].values, segment.loc[:,"median"].values],
                  color=color, solid_capstyle="butt", linewidth=2.5)

    # Plot 0 line
    ax.axhline(y=0, color="darkgrey", linewidth=1.5)

    # Plot chromosome names
    text_x = []
    for chrom in chroms:
        text_x.append(frags.chrom_shift[chrom] + (frags.genome[chrom] / 2))
    ax.set_xticks(text_x)
    ax.set_xticklabels(chroms, rotation=90)
    
    # Set titles
    if title != None:
        ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel("Chromosomes", fontsize=12)
    ax.set_ylabel("log2 ratio", fontsize=12)
    ax.set_ylim((-3,5))
    ax.set_xlim((0,last_position))

    # Plot arrows if segments of the chart
    off_chart = np.logical_or(cnvs.df.loc[:,"median"].values >= 5.0,
                              cnvs.df.loc[:,"median"].values <= -3.0)
    if np.sum(off_chart) > 0:
        off_segments = cnvs.loc[off_chart,:]
        for i in range(off_segments.shape[0]):
            segment = off_segments.iloc[i,:]
            cnv_type = segment.df.loc[:,"Corrected_Call"].values[0]
            color =  cnv_types[cnv_type] if cnv_type in cnv_types else "darkred"
            start = segment.starts[0] + frags.chrom_shift[segment.index.extract_labels()[0]]
            start = start / last_position # convert to ratio
            end = segment.index.extract_ends()[0] + frags.chrom_shift[segment.index.extract_labels()[0]]
            end = end / last_position # convert to ratio
            position = (start + end) / 2.0 # take middle
            median = segment.loc[:,"median"]

            if median >= 5.0:
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


def fragment_distribution(length_dist, title=None, show=True, save=None, ax=None):
    """
    Plot the fragment distribution
    """
    
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
    for maxima in argrelmax(fragment_lengths, order=35)[0]:
        plt.axvline(x=maxima, color="grey", linestyle="--", dash_capstyle="round", linewidth=0.5)
    
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


def window_sort(windows):
    """
    Sort row by those most correlated with the means
    """
    
    # Calculate mean values
    window_means = windows.mean(axis=0)
    
    # Calculate correlations with means
    correlations = np.zeros(windows.shape[0])
    for i in range(windows.shape[0]):
        #stat, pvalue = pearsonr(windows[i,:], window_means)
        stat, pvalue = spearmanr(windows[i,:], window_means)
        correlations[i] = stat
    
    # Set nan to 0
    correlations[np.isnan(correlations)] = 0
        
    # Sort by correlations
    sort_index = np.argsort(correlations)[::-1]
    
    return(sort_index, window_means)

        
def window_values(windows):
    """
    """
    
    # Create color map
    cmap = LinearSegmentedColormap.from_list('mycmap', ['#660220', '#b01b2f', '#d46151', '#f2a585',
                                                        '#fcdbc8', '#f7f7f7', '#d2e5ef', '#94c5dd',
                                                        '#4794c1', '#2668aa', '#083160'][::-1])
                                                        
    # Normalize colors
    norm = matplotlib.colors.Normalize(vmin=0, vmax=2)
    
    # Sort rows
    sort_index, means = window_sort(windows)
    
    # Plot
    plt.plot(means)
    
    
def mean_window(window_scores, center_x=None, title=None, xlabel=None, ylabel=None, show=True, save=None, ax=None):
    """
    """
    
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
    