import pandas as pd
import numpy as np
from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
from typing import Dict

# Local imports
from ..fragments import Fragments
from ..peak_calling.CallPeaks import call_peaks, normalize_signal


def coverage(intervals: Fragments | LabeledIntervalArray,
				chrom: str = None,
				min_length: int = None,
				max_length: int = None):
	"""
	Calculate coverage for each position in AIList range

	Parameters
	----------
		frags : Fragments
			Fragments object to calculate coverage for
		chrom : str
			Chromosome to calculate coverage for [default = None]
		min_length : int
			Minimum length of intervals to include [default = None]
		max_length : int
			Maximum length of intervals to include [default = None]g

	Returns
	-------
		scores : dict of pandas.Series or pandas.Series
			Position on index and coverage as values

	"""

	if isinstance(intervals, Fragments):
		intervals = intervals.frags
	
	scores = intervals.coverage(chrom, min_length, max_length)
	
	return scores


def coverage_peaks(intervals: Fragments | LabeledIntervalArray | None = None,
					cov_dict: Dict[str,pd.Series] | None = None,
					merge_distance: int = 5,
					min_length: int = 50,
					max_length: int = 150,
					min_interval_filter: int = None,
					max_interval_filter: int = None,
					normalize_scores: bool = True,
					verbose: bool = False):
    """
    Call coverage peaks

    Parameters
    ----------
        frags : Fragments
            Fragments object
        cov_dict : dict of pd.Series
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
    if intervals is None and cov_dict is None:
        raise ValueError("Either intervals or cov_dict must be provided")

	# Initailize peaks
    peaks = LabeledIntervalArray()

    # If intervals are provided
    if intervals is not None:
        if isinstance(intervals, Fragments):
            intervals = intervals.frags
        
        # Iterate over chromosomes
        for chrom in intervals.unique_labels:
            if verbose: print(chrom)
            cov_dict = coverage(intervals,
								chrom = chrom,
								min_length = min_interval_filter,
								max_length = max_interval_filter)
            if normalize_scores:
                normalize_coverage(cov_dict)
            if len(cov_dict[chrom]) > 0:
                peaks.append(call_peaks(cov_dict[chrom],
                                        str(chrom),
                                        merge_distance=merge_distance,
                                        min_length=min_length,
                                        max_length=max_length))

    elif cov_dict is not None:
        # Iterate over wps
        for chrom in cov_dict:
            if verbose: print(chrom)
            if len(cov_dict[chrom]) > 0:
                peaks.append(call_peaks(cov_dict[chrom],
                                        str(chrom),
                                        merge_distance=merge_distance,
                                        min_length=min_length,
                                        max_length=max_length))

    return peaks


def normalize_coverage(cov_dict: Dict[str,pd.Series],
						window: int = 1000,
						smooth: bool = True,
						smooth_window: int = 21,
						polyorder: int = 2, 
						n_threads: int = 1,
						method: str = "mean",
						verbose: bool = False):
    """
	Normalize coverage
    
    Parameters
    ----------
        cov_dict : dict of pd.Series
            Coverage to normalize
        window : int
            Window size to use for normalization
        smooth : bool
            Smooth coverage before normalization
        smooth_window : int
            Window size to use for smoothing
        polyorder : int
            Order of polynomial to use for smoothing
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

    for chrom in cov_dict:
        if verbose: print(chrom)
        if method == "median":
            cov_dict[chrom].iloc[:] = normalize_signal(cov_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=False)
        elif method == "mean":
            cov_dict[chrom].iloc[:] = normalize_signal(cov_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=True)

    return None


def bin_counts(intervals: Fragments | LabeledIntervalArray,
				   bin_size: int = 100000,
				   min_length: int = None,
				   max_length: int = None) -> IntervalFrame:
		"""
		Calculate coverage of fragments

		Parameters
		----------
			bin_size : int
				Size of each bin
			min_length : int
				Minimum length of each fragment to consider
			max_length : int
				Maximum length of each fragment to consider

		Returns
		-------
			bins : IntervalFrame
				Bins containing number of hits
		"""

		if isinstance(intervals, Fragments):
			intervals = intervals.frags

		# Determine nhits per bin
		bins = intervals.bin_nhits(bin_size, min_length, max_length)

		# Create interval frame
		bins_iframe = IntervalFrame(bins[0])
		bins_iframe.df.loc[:,"counts"] = bins[1]

		return bins_iframe


def predict_sex(coverage: Fragments | IntervalFrame,
                fracReadsInChrYForMale: float = 0.002,
                chrXMedianForMale: float = -0.5,
                useChrY: bool = True):
    """
    Predict Gender from coverage

    Parameters
    ----------
        coverage : Fragments or IntervalFrame
            Coverage to use for prediction
        fracReadsInChrYForMale : float
            Fraction of read in Y chromosome threshold
        chrXMedianForMale : float
            Median chrX log2 ratio threshold
        useChrY : bool
            Whether to use the Y chromosome

    Returns
    -------
        gender_record : dict
            Record of predicted gender and metrics
    """

    # Convert to IntervalFrame
    if isinstance(coverage, Fragments):
        bin_coverage = bin_counts(coverage)
    else:
        bin_coverage = coverage

    # Record chromosome index
    chroms = bin_coverage.index.unique_labels

    # Determine sex chromosome names
    chrXStr = [chrom for chrom in chroms if "X" in chrom][0]
    chrYStr = [chrom for chrom in chroms if "Y" in chrom][0]

    # Calculate ratios
    median = np.median([np.median(bin_coverage.loc[chrom,"counts"].values) for chrom in chroms])
    chrXratios = np.log2(bin_coverage.loc[chrXStr,"counts"].values / median)
    chrXratios[np.isinf(chrXratios)] = 0

    # Calculate metrics
    total_counts = np.sum([np.sum(bin_coverage.loc[chrom,"counts"].values) for chrom in chroms])
    chrXMedian = np.median(chrXratios)
    chrYCov = np.sum(bin_coverage.loc[chrYStr,"counts"].values) / total_counts

    # Predict gender
    if chrXMedian < chrXMedianForMale:
        if useChrY and (chrYCov < fracReadsInChrYForMale):
            gender = "female"
        else:
            gender = "male"
    else:
        gender = "female"

    gender_record = {"gender":gender, "chrYCovRatio":chrYCov, "chrXMedian":chrXMedian}

    return gender_record