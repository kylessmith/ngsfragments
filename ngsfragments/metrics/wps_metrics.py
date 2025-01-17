from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
import numpy as np
import pandas as pd
from typing import Dict

# Local imports
from ..fragments import Fragments
from ..io.read_sam import from_sam
from ..segment.cnv import call_cnvs
from ..correct.correct_intervals import calculate_interval_bias
from ..correct.correction import gaussian_smooth
from ..peak_calling.CallPeaks import call_peaks, normalize_signal
from .metric_utils import score_window, nfr


def wps_nfr_scores(wps_windows: IntervalFrame,
                   center_x: int) -> np.ndarray:
    """
    """

    # Normalize values
    values = wps_windows.df.values
    values = ((values.T - np.min(values,axis=1)) / (np.max(values, axis=1) - np.min(values, axis=1))).T

    # Calculate NFR
    nf = np.sum(values[:,center_x - 50:center_x + 50], axis=1)
    n1 = np.sum(values[:,center_x - 200:center_x - 50], axis=1)
    n2 = np.sum(values[:,center_x + 50:center_x + 200], axis=1)

    nfr_score = np.log2(nf) - np.log2((n1 + n2) / 2)

    return nfr_score


def wps_window_scores(intervals: Fragments | LabeledIntervalArray,
                      windows: LabeledIntervalArray,
                      window_size: int,
                    protection: int = 120,
                    min_length: int = 120,
                    max_length: int = 220,
                    smooth: bool = True,
                    verbose: bool = False) -> IntervalFrame:
    """
    """

    # Convert intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Sort windows
    windows = windows[windows.sorted_index()]

    # Iterate over chromosomes
    chroms = intervals.unique_labels
    values = np.zeros((len(windows), window_size))
    i = 0
    for chrom in chroms:
        if verbose: print(chrom, flush=True)

        # Get wps scores
        chrom_windows = windows.get(chrom)
        if len(chrom_windows) > 0:
            wps = intervals.wps(protection, chrom, min_length, max_length)
            if smooth:
                wps[chrom][:] = gaussian_smooth(normalize_signal(wps[chrom].values), scale=10)
            else:
                wps[chrom][:] = normalize_signal(wps[chrom].values)
            wps = wps[chrom]

            # Iterate over windows
            for interval in chrom_windows:
                if (interval.end - interval.start) == (window_size):
                    value = wps.loc[interval.start:interval.end-1].values
                    if len(value) == (window_size):
                        values[i,:] = value
                    else:
                        #values[i,:] = np.nan
                        values[i,:] = 0
                else:
                    #values[i,:] = np.nan
                    values[i,:] = 0

                i += 1

    # Convert to DataFrame
    wps_df = pd.DataFrame(values, index=range(len(windows)), columns=range(window_size))
    iframe = IntervalFrame(intervals = windows, df=wps_df)

    return iframe

    

def wps_windows(intervals: Fragments | LabeledIntervalArray,
                protection: int = 120,
                min_length: int = 120,
                max_length: int = 220,
                upstream: int = 1000,
                downstream: int = 1000,
                feature: str = 'tss',
                genome_version: str = "hg38",
                smooth: bool = True,
                verbose: bool = False) -> pd.DataFrame:
    """
    Determine wps scores within windows

    Parameters
    ----------
        wps : Dict[pd.Series]
            Dictionary of wps scores
        chroms : np.ndarray
            Chromosome array
        starts : np.ndarray
            Start position array
        ends : np.ndarray
            End position array

    Returns
    -------
        pd.DataFrame
            DataFrame of wps scores within windows
    """

    # Get windows
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)

    # Get genes
    if feature == "tss":
        windows = genome.get_intervals("tss", upstream=upstream, downstream=downstream, filter_column="gene_type", filter_selection="protein_coding")
    elif feature == "tes":
        windows = genome.get_intervals("tes", upstream=upstream, downstream=downstream, filter_column="gene_type", filter_selection="protein_coding")
    else:
        raise NotImplementedError("Only tss and tes are supported currently")
    
    # Convert intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Iterate over chromosomes
    chroms = windows.index.unique_labels[np.in1d(windows.index.unique_labels, intervals.unique_labels)]
    values = np.zeros((windows.shape[0], upstream + downstream + 1))
    
    i = 0
    for chrom in chroms:
        if verbose: print(chrom, flush=True)
        chrom_windows = windows.loc[chrom,:]
        wps = intervals.wps(protection, chrom, min_length, max_length)
        if smooth:
            wps[chrom][:] = gaussian_smooth(normalize_signal(wps[chrom].values), scale=10)
        else:
            wps[chrom][:] = normalize_signal(wps[chrom].values)

        # Score windows
        values[i:i+chrom_windows.shape[0],:] = score_window(wps[chrom], chrom_windows, chrom_windows.df.loc[:,"Strand"].values, upstream + downstream)
        i += chrom_windows.shape[0]

    # Convert to DataFrame
    wps_df = pd.DataFrame(values, index=windows.df.loc[:,"gene_name"].values, columns=np.arange(-upstream, downstream + 1))
    wps_df = wps_df.loc[np.sum(pd.isnull(wps_df),axis=1).values == 0,:]
    
    return wps_df


def wps_score_tfs(intervals: Fragments | LabeledIntervalArray,
                protection: int = 120,
                min_length: int = 120,
                max_length: int = 220,
                upstream: int = 1000,
                downstream: int = 1000,
                genome_version: str = "hg38",
                smooth: bool = True,
                verbose: bool = False) -> pd.DataFrame:
    """
    """

    # Get windows
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)
    tfs = genome["GTRD_1000sites"]

    # Expand intervals
    for tf in tfs:
        starts = tfs[tf].index.starts - upstream
        ends = tfs[tf].index.ends + downstream
        labels = tfs[tf].index.labels

        # Construct IntervalFrame
        new_intervals = LabeledIntervalArray()
        new_intervals.add(starts, ends, labels)
        iframe = IntervalFrame(new_intervals)

        tfs[tf] = iframe

    # Convert intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Iterate over chromosomes
    chroms = intervals.unique_labels
    values = pd.DataFrame(np.zeros((len(tfs), upstream + downstream + 1)),
                          index = list(tfs.keys()))
    n = {tf:0 for tf in tfs.keys()}

    for chrom in chroms:
        if verbose: print(chrom, flush=True)
        wps = intervals.wps(protection, chrom, min_length, max_length)
        if smooth:
            wps[chrom][:] = gaussian_smooth(normalize_signal(wps[chrom].values), scale=10)
        else:
            wps[chrom][:] = normalize_signal(wps[chrom].values)
        
        # Scale
        #wps[chrom][:] = (wps[chrom].values - np.min(wps[chrom].values)) / (np.max(wps[chrom].values) - np.min(wps[chrom].values))
        
        # Iterate over tfs
        for tf in tfs:
            chrom_windows = tfs[tf].loc[chrom,:]
            # Score windows
            if chrom_windows.shape[0] > 0:
                scores = score_window(wps[chrom], chrom_windows, None, upstream + downstream)
                values.loc[tf,:] += scores.sum(axis=0)
                n[tf] += scores.shape[0]

    # Calculate average
    for tf in tfs:
        if n[tf] > 0:
            values.loc[tf,:] = values.loc[tf,:].values / n[tf]

    return values


def wps_nfr_tfs(intervals: Fragments | LabeledIntervalArray,
                protection: int = 120,
                min_length: int = 120,
                max_length: int = 220,
                upstream: int = 1000,
                downstream: int = 1000,
                genome_version: str = "hg38",
                smooth: bool = True,
                verbose: bool = False) -> pd.DataFrame:
    """
    """

    # Get windows
    scores = wps_score_tfs(intervals,
                protection,
                min_length,
                max_length,
                upstream,
                downstream,
                genome_version,
                smooth,
                verbose)

    # Calculate NFR
    nfr_scores = np.zeros(scores.shape[0])
    for i in range(scores.shape[0]):
        nfr_scores[i] = nfr(scores.iloc[i,:].values, center_x=1000, scale=True)[0]

    # Convert to DataFrame
    nfr_df = pd.DataFrame(nfr_scores, index=scores.index.values, columns=["NFR"])

    return nfr_df