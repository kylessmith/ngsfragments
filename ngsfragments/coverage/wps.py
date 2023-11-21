import pandas as pd
import numpy as np
from ailist import LabeledIntervalArray
from typing import Dict

# Local imports
from ..fragments import Fragments
from ..peak_calling.CallPeaks import call_peaks, normalize_signal


def wps(intervals: Fragments | LabeledIntervalArray,
        chrom: str = None,
        protection: int = 60,
        min_length: int = None,
        max_length: int = None):
    """
    Calculate Window Protection Score for each position in AIList range

    Parameters
    ----------
        protection : int
            Protection window to use
        min_length : int
            Minimum length of intervals to include [default = None]
        max_length : int
            Maximum length of intervals to include [default = None]
        label : str
            Label for hierarchical indexing

    Returns
    -------
        scores : dict of pandas.Series or pandas.Series
            Position on index and WPS as values

    """

    if isinstance(intervals, Fragments):
        intervals = intervals.frags
    
    scores = intervals.wps(protection, chrom, min_length, max_length)
    
    return scores


def wps_peaks(intervals: Fragments | LabeledIntervalArray | None = None,
              wps_dict: Dict[str,pd.Series] | None = None,
              merge_distance: int = 5,
              min_length: int = 120,
              max_length: int = 220,
              protection: int = 120,
              min_interval_filter: int = None,
              max_interval_filter: int = None,
              normalize_scores: bool = True,
              verbose: bool = False):
    """
    Call WPS peaks

    Parameters
    ----------
        frags : Fragments
            Fragments object
        wps_dict : dict of pd.Series
            WPS scores for each chromosome
        shift : int
            Shift peaks by this amount
        merge_distance : int
            Merge peaks within this distance
        min_length : int
            Minimum length of peaks to call
        max_length : int
            Maximum length of peaks to call

    Returns
    -------
        peaks : LabeledIntervalArray
            Peaks called from WPS scores
    """

    # Cehck inputs
    if intervals is None and wps_dict is None:
        raise ValueError("Either intervals or wps_dict must be provided")

	# Initailize peaks
    peaks = LabeledIntervalArray()

    # If intervals are provided
    if intervals is not None:
        if isinstance(intervals, Fragments):
            intervals = intervals.frags
        
        # Iterate over chromosomes
        for chrom in intervals.unique_labels:
            if verbose: print(chrom)
            wps_dict = wps(intervals,
                            chrom = chrom,
                            protection = protection,
                            min_length = min_interval_filter,
                            max_length = max_interval_filter)
            if normalize_scores:
                normalize_wps(wps_dict)
            if len(wps_dict[chrom]) > 0:
                peaks.append(call_peaks(wps_dict[chrom],
                                        str(chrom),
                                        merge_distance=merge_distance,
                                        min_length=min_length,
                                        max_length=max_length))

    elif wps_dict is not None:
        # Iterate over wps
        for chrom in wps_dict:
            if verbose: print(chrom)
            if len(wps_dict[chrom]) > 0:
                peaks.append(call_peaks(wps_dict[chrom],
                                        str(chrom),
                                        merge_distance=merge_distance,
                                        min_length=min_length,
                                        max_length=max_length))

    return peaks


def normalize_wps(wps_dict: Dict[str,pd.Series],
				  window: int = 1000,
				  smooth: bool = True,
				  smooth_window: int = 21,
				  polyorder: int = 2, 
				  n_threads: int = 1,
				  method: str = "mean",
				  verbose: bool = False):
    """
	Normalize WPS

    Parameters
    ----------
        wps_dict : dict of pd.Series
            WPS scores for each chromosome
        window : int
            Window size to use for normalization
        smooth : bool
            Smooth signal before normalization
        smooth_window : int
            Window size to use for smoothing
        polyorder : int
            Polynomial order to use for smoothing
        n_threads : int
            Number of threads to use for smoothing
        method : str
            Method to use for normalization
        verbose : bool
            Print progress

    Returns
    -------
        None
    """

    for chrom in wps_dict:
        if verbose: print(chrom)
        if method == "median":
            wps_dict[chrom].iloc[:] = normalize_signal(wps_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=False)
        elif method == "mean":
            wps_dict[chrom].iloc[:] = normalize_signal(wps_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=True)
    
    return None