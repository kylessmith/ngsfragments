import pandas as pd
import numpy as np
from intervalframe import IntervalFrame
import math
from scipy.stats import norm
from scipy.signal import fftconvolve
from scipy.ndimage import convolve
import statsmodels.api as sm

# Local imports
from .cylowess.cylowess import lowess
    
    
def match_bins(bins : IntervalFrame,
               bin_bias : IntervalFrame):
    """
    Restrict to bins in bin_bias

    Parameters
    ----------
        bins
            pd.DataFrame
        bin_bias
            pd.DataFrame

    Returns
    -------
        bins
            pd.DataFrame
        bin_bias
            pd.DataFrame
    """

    # Find common chromosomes and sort
    chroms = bin_bias.index.unique_labels[np.in1d(bin_bias.index.unique_labels, 
                                                  bins.index.unique_labels)]
    bins = bins.loc[chroms,:]
    bin_bias = bin_bias.loc[chroms,:]

    # Match bins
    bins = bins.exact_match(bin_bias)

    # Match bin_bias
    bin_bias = bin_bias.exact_match(bins)
    
    return bins, bin_bias


def filter_bins(bins: IntervalFrame,
                bin_bias: IntervalFrame,
                blacklist_cutoff: float = 0.1,
                mapping_cutoff: float = 0.5):
    """
    Filter bins by genomic bin bias

    Parameters
    ----------
        bins : IntervalFrame
            IntervalFrame of bins
        bin_bias : IntervalFrame
            IntervalFrame of bin bias
        blacklist_cutoff : float
            Maximum blacklist fraction

    Returns
    -------
        bins
        bin_bias
    """
    
    # Find bins passing filters
    chosen = bin_bias.df.loc[:,"blacklist"].values < blacklist_cutoff

    # Find if mappability is present
    if "mappability" in bin_bias.columns:
        chosen = np.logical_and(chosen, bin_bias.df.loc[:,"mappability"].values > mapping_cutoff)
    
    # Filter bins
    bin_bias = bin_bias.iloc[chosen,:]
    bins = bins.iloc[chosen,:]
    
    return bins, bin_bias


def binned_bias_correct_counts(values: np.ndarray,
                                bin_bias: IntervalFrame,
                                n_bins: int = 50) -> np.ndarray:
    """
    Bin bias values and correct for bias

    Parameters
    ----------
        values : np.ndarray
            Values to correct
        bias_values : IntervalFrame
            Bias values
        n_bins : int
            Number of bins

    Returns
    -------
        corrected_values : np.ndarray
            Corrected values
    """

    # Iterate over columns
    new_values = values.copy().astype(float)
    for column in bin_bias.columns:
        if column in ["n", "blacklist", "mappability"]:
            continue

        # Get bias values
        bias_values = bin_bias.df.loc[:,column].values

        # Bin values
        bins = np.histogram(bias_values, bins=n_bins)[1]
        assignments = np.digitize(bias_values, bins)

        # Calculate median bias values
        medians = np.zeros(n_bins)
        for i in range(n_bins):
            if np.sum(assignments == (i+1)) > 0:
                medians[i] = np.median(new_values[assignments == (i+1)])

        # Calculate LOWESS regression
        res = sm.nonparametric.lowess(medians, bins[:-1], frac=1.0, it=10, return_sorted=False)

        # Calculate corrected values
        for i in range(n_bins):
            new_values[assignments == (i+1)] = new_values[assignments == (i+1)] / res[i]

    return new_values


def correct_counts(values: np.ndarray,
                   bin_bias: IntervalFrame) -> np.ndarray:
    """
    Correct for biases in genome metrics file

    Parameters
    ----------
        values : np.ndarray
            Values to correct
        bin_bias : pd.DataFrame
            Bin bias

    Returns
    -------
        values : np.ndarray
            Corrected values
    """
    
    # Calculate LOWESS regression
    median = np.median(values)

    # Calculate delta
    delta = delta = (values.max() - values.min()) * 0.01

    # Iterate over columns
    for column in bin_bias.columns:
        if column in ["n", "blacklist"]:
            continue
        
        # Calculate LOWESS regression
        prediction = sm.nonparametric.lowess(values,
                                             bin_bias.loc[:,column].values,
                                             frac=0.6666,
                                             return_sorted=False,
                                             delta=delta)
        #prediction = lowess(values,
        #                    bin_bias.loc[:,column].values.astype(float),
        #                    frac=0.6666,
        #                    delta=delta)

        # Determine corrected values
        residuals = (values - prediction)
        corrected_values = residuals + median

        # Change minimum to zero
        corrected_values[corrected_values < 0] = 0

        # Update values
        values = corrected_values

    # Change minimum to zero
    values[values < 0] = 0

    return values


def gaussian_smooth(values: np.ndarray,
                    scale: int = 1,
                    p_min: float = 1e-10):
    """
    Smooth signal with gaussian window

    Params
    ------
        values
            np.array
        scale
            int
        p_min
            float

    Returns
    -------
        convolvution
            np.array
    """
    
    ## Gaussian scale space
    delta_tau = 0.5 * math.log(2)
    sigma = math.exp((scale - 1) * delta_tau)
    X = norm.ppf(p_min, 0, sigma)
    window = norm.pdf(np.arange(int(round(X)), int(-round(X) + 1)), 0, sigma)[np.newaxis]
    windowSize = len(window[0])
    window = window / np.sum(window[0])

    if windowSize < 1250:
        convolvution = convolve(values[np.newaxis], window, mode='constant', cval=0.0)[0]
    else:
        convolvution = fftconvolve(values[np.newaxis], window, mode='same')[0]

    return convolvution


def correct(bin_coverage: IntervalFrame,
            genome_version: str = "hg19",
            bin_size: int = 100000,
            wgbs: bool = False,
            bin_correct: bool = True,
            verbose: bool = False):
    """
    Correct bias in bin counts

    Params
    ------
        bin_coverage
            pd.DataFrame
        bin_bias_h5_fn
            str
        genome_version
            str
        bin_size
            int
        wgbs
            bool
        verbose
            bool

    Returns
    -------
        corrected_values
            pd.DataFrame
    """

    if verbose: print("Reading genome biases...")
    import genome_info
    g = genome_info.GenomeInfo(genome_version)
    bin_bias = g.calculate_bin_bias(bin_size=bin_size)

    if verbose: print("Matching bins...")
    bin_coverage, bin_bias = match_bins(bin_coverage, bin_bias)

    if verbose: print("Filtering bins...")
    bin_coverage, bin_bias = filter_bins(bin_coverage, bin_bias,
                                         blacklist_cutoff=0.1)
    
    if verbose: print("Correcting bins...")
    if bin_correct:
        bin_coverage.loc[:,"corrected_counts"] = binned_bias_correct_counts(bin_coverage.loc[:,"counts"].values, bin_bias)
    else:
        bin_coverage.loc[:,"corrected_counts"] = correct_counts(bin_coverage.loc[:,"counts"].values.astype(float), bin_bias)
    
    return bin_coverage
