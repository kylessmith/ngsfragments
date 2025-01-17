import pandas as pd
import numpy as np
from intervalframe import IntervalFrame
import math
from scipy.stats import norm
from scipy.signal import fftconvolve
from scipy.ndimage import convolve
import statsmodels.api as sm
from loess import loess_2d, loess_1d
from sklearn.cluster import KMeans

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
                mapping_cutoff: float = 0.90,
                keep_sex_chroms: bool = False):
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
        mapping_cutoff : float
            Minimum mappability fraction
        keep_sex_chroms : bool
            Keep sex chromosomes

    Returns
    -------
        bins
        bin_bias
    """
    
    # Find bins passing filters
    chosen = bin_bias.df.loc[:,"blacklist"].values < blacklist_cutoff

    # Keep
    if keep_sex_chroms:
        chosen = np.logical_or(chosen, bin_bias.index.labels == "chrX")
        chosen = np.logical_or(chosen, bin_bias.index.labels == "chrY")

    # Find if mappability is present
    if "mappability" in bin_bias.columns:
        chosen = np.logical_and(chosen, bin_bias.df.loc[:,"mappability"].values > mapping_cutoff)
    
    # Filter bins
    bin_bias = bin_bias.iloc[chosen,:]
    bins = bins.iloc[chosen,:]
    
    return bins, bin_bias


def binned_bias_correct_counts(values: np.ndarray,
                                bin_bias: IntervalFrame,
                                n_bins: int = 15) -> np.ndarray:
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

    # Remove nans
    not_nan = ~pd.isnull(values)
    values = values[not_nan]
    bin_bias = bin_bias.iloc[not_nan,:]

    # Iterate over columns
    new_values = values.copy().astype(float)
    for column in bin_bias.columns:
        if column in ["n", "blacklist"]:
            continue

        # Get bias values
        bias_values = bin_bias.df.loc[:,column].values

        # Bin values
        b = np.histogram(bias_values, bins=n_bins)[1]
        labels = np.digitize(bias_values, b)

        # Calculate median bias values
        label_x = pd.Series(new_values).groupby(labels).median()

        # Calculate LOWESS regression
        res = loess_1d.loess_1d(b, label_x.values, frac=1.0)
        #res = sm.nonparametric.lowess(medians, bins[:-1], frac=1.0, it=10, return_sorted=False)

        # Calculate corrected values
        for i, l in enumerate(label_x.index.values):
            new_values[labels==l] = new_values[labels==l] / res[1][i]

    return new_values


def gc_mappability_correct(values: np.ndarray,
                           bin_bias: IntervalFrame) -> np.ndarray:
    """
    Correct for gc and mappability together

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

    # Run 2D loess
    res = loess_2d.loess_2d(x=bin_bias.df.loc[:,"gc"].values, y=bin_bias.df.loc[:,"mappability"].values, z=values)

    # Calculate corrected values
    new_values = values / res[0]

    return new_values


def gc_mappability_bin_correct(values: np.ndarray,
                                bin_bias: IntervalFrame,
                                n_bins: int = 25) -> np.ndarray:
    """
    Correct for gc and mappability together

    Parameters
    ----------
        values : np.ndarray
            Values to correct
        bin_bias : pd.DataFrame
            Bin bias
        n_bins : int
            Number of bins to use

    Returns
    -------
        values : np.ndarray
            Corrected values
    """

    # Remove nans
    n = len(values)
    not_nan = ~pd.isnull(values)
    values = values[not_nan]
    bin_bias = bin_bias.iloc[not_nan,:]

    # Initialize new values
    new_values = values.copy().astype(float)

    # Determine bins
    m_labels = np.round(bin_bias.df.loc[:,"mappability"].values*100).astype(int).astype(str)
    gc_labels = np.round(bin_bias.df.loc[:,"gc"].values*100).astype(int).astype(str)
    labels = np.array([str(gc_labels[i])+"-"+str(m_labels[i]) for i in range(len(m_labels))])

    # Calculate medians for clusters
    label_gc = pd.Series(bin_bias.df.loc[:,"gc"].values * 100).groupby(labels).median()
    label_m = pd.Series(bin_bias.df.loc[:,"mappability"].values * 100).groupby(labels).median()
    label_x = pd.Series(values).groupby(labels).median()

    # Run 2D loess
    res = loess_2d.loess_2d(x=label_gc.values, y=label_m.values, z=label_x.values, frac=0.33, degree=2)
    res = pd.Series(res[0], index=label_x.index.values)
        
    # Calculate corrected values
    for i, l in enumerate(label_x.index.values):
        new_values[labels==l] = new_values[labels==l] / res[l]
    new_values = new_values * np.nanmedian(values)

    # Correct others
    #bin_bias = bin_bias.copy()
    #bin_bias.df = bin_bias.df.drop(columns=["mappability","blacklist"])
    #if bin_bias.shape[1] > 0:
    #    new_values = chr_bias_correct_counts(new_values, bin_bias)

    # Add nans
    final_values = np.zeros(n)
    final_values[:] = np.nan
    final_values[not_nan] = new_values

    return final_values



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
            use_normal: bool = True,
            bin_correct: bool = True,
            keep_sex_chroms: bool = False,
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
                                         blacklist_cutoff=0.1,
                                         mapping_cutoff=0.90,
                                         keep_sex_chroms = keep_sex_chroms)
    
    if verbose: print("Correcting bins...")
    if bin_correct:
        if use_normal:
            normal = g["normal_bin"+str(bin_size)]
            normal = normal.exact_match(bin_coverage)
            bin_coverage = bin_coverage.exact_match(normal)
            bin_bias = bin_bias.exact_match(normal)
            normal_corrected_counts = binned_bias_correct_counts(bin_coverage.loc[:,"counts"].values, bin_bias)
            bin_coverage.loc[:,"corrected_counts"] = binned_bias_correct_counts(normal_corrected_counts, bin_bias)
        else:
            bin_coverage.loc[:,"corrected_counts"] = binned_bias_correct_counts(bin_coverage.loc[:,"counts"].values, bin_bias)
    else:
        if use_normal:
            normal = g["normal_bin"+str(bin_size)]
            normal = normal.exact_match(bin_coverage)
            bin_coverage = bin_coverage.exact_match(normal)
            bin_bias = bin_bias.exact_match(normal)
            normal_corrected_counts = correct_counts(bin_coverage.loc[:,"counts"].values.astype(float), bin_bias)
            bin_coverage.loc[:,"corrected_counts"] = correct_counts(normal_corrected_counts, bin_bias)
        else:
            bin_coverage.loc[:,"corrected_counts"] = correct_counts(bin_coverage.loc[:,"counts"].values.astype(float), bin_bias)
    
    return bin_coverage
