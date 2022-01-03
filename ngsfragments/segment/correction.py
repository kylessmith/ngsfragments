import pandas as pd
import numpy as np
from intervalframe import IntervalFrame, write_h5, read_h5
import h5py
import math
from scipy.stats import norm
from scipy.signal import fftconvolve
from scipy.ndimage import convolve
import statsmodels.api as sm
import os

# Set data directory
data_dir = os.path.split(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])[0]
data_dir = os.path.join(data_dir, "data")


def read_bin_bias(bin_bias_h5_fn, bin_size=100000, genome="hg19", verbose=False):
    """
    Read genome bin bias from given h5 file
    
    Parameters
    ----------
        bin_bias_h5_fn : str
            Name of the h5 file
    
    Returns
    -------
        bin_bias : IntervalFrame
            Bias bins for given chromosome
    """

    # Find distributed h5 reference bins
    if bin_bias_h5_fn is None:
        h5_file = genome + "_ucsc_bin" + str(bin_size) + ".h5"
        if verbose: print("Reference genome:", genome)
        bin_bias_h5_fn = os.path.join(data_dir, h5_file)
        if verbose: print("   Found bin_bias:", h5_file)
    
    # Open h5 file
    bin_bias_h5 = h5py.File(bin_bias_h5_fn, "r")
    
    # Read biases
    bin_bias = read_h5.read_h5_intervalframe(bin_bias_h5["iframe"])
    
    # Close h5 file
    bin_bias_h5.close()
    
    return bin_bias
    
    
def match_bins(bins, bin_bias):
    """
    Restrict to bins in bin_bias

    Params
    ------
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
    
    # Match bins
    bins = bins.exact_match(bin_bias)

    # Match bin_bias
    bin_bias = bin_bias.exact_match(bins)
    #print(bins.shape, bin_bias.shape)
    
    return bins, bin_bias


def filter_bins(bins, bin_bias, n_cutoff=0.1, map_cutoff=0.66, blacklist_cutoff=0.1):
    """
    Filter bins by genomic bin bias

    Params
    ------
        bins
        bin_bias
        n_cutoff
        map_cutoff
        blacklist_cutoff

    Returns
    -------
        bins
        bin_bias
    """
    
    # Find bins passing filters
    chosen = np.logical_and(bin_bias.df.loc[:,"n"].values < n_cutoff,
                            bin_bias.df.loc[:,"mappability"].values > map_cutoff)
    chosen = np.logical_and(chosen, bin_bias.df.loc[:,"blacklist"].values < blacklist_cutoff)
    
    # Filter bins
    bin_bias = bin_bias.iloc[chosen,:]
    bins = bins.iloc[chosen,:]
    
    return bins, bin_bias


def correct_counts(values, bin_bias):
    """
    Correct for biases in genome metrics file

    Params
    ------
        values
            np.array
        bin_bias
            pd.DataFrame

    Returns
    -------
        corrected_values
            np.array
    """
    
    # Calculate LOWESS regression
    median = np.median(values)
    
    # Correct for GC content
    gc_prediction = sm.nonparametric.lowess(values,
                                            bin_bias.loc[:,"gc"].values,
                                            frac=0.6666,
                                            return_sorted=False)
    residuals = (values - gc_prediction)
    corrected_values = residuals + median

    # Correct for repeat content
    repeat_prediction = sm.nonparametric.lowess(corrected_values,
                                                bin_bias.loc[:,"repeat"].values,
                                                frac=0.6666,
                                                return_sorted=False)
    residuals = (corrected_values - repeat_prediction)
    corrected_values = residuals + median

    # Correct for mappability content
    prediction = sm.nonparametric.lowess(corrected_values,
                                                bin_bias.loc[:,"mappability"].values,
                                                frac=0.6666,
                                                return_sorted=False)

    # Determine corrected values
    residuals = (corrected_values - prediction)
    corrected_values = residuals + median

    return corrected_values


def gaussian_smooth(values, scale=1, p_min=1e-10):
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


def correct(bin_coverage, bin_bias_h5_fn=None, genome="hg19", bin_size=100000, verbose=False):
    """
    Correct bias in bin counts

    Params
    ------
        bin_coverage
            pd.DataFrame
        bin_bias_h5_fn
            str
        verbose
            bool

    Returns
    -------
        corrected_values
            pd.DataFrame
    """

    if verbose: print("Reading genome biases...")
    bin_bias = read_bin_bias(bin_bias_h5_fn, bin_size, genome, verbose)

    if verbose: print("Matching bins...")
    bin_coverage, bin_bias = match_bins(bin_coverage, bin_bias)

    if verbose: print("Filtering bins...")
    bin_coverage, bin_bias = filter_bins(bin_coverage, bin_bias, n_cutoff=0.1, map_cutoff=0.66, blacklist_cutoff=0.1)
    
    if verbose: print("Correcting bins...")
    bin_coverage.loc[:,"corrected_counts"] = correct_counts(bin_coverage.loc[:,"counts"].values, bin_bias)
    
    return bin_coverage