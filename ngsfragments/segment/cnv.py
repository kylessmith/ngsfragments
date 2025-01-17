import pandas as pd
import numpy as np
import os
from typing import List
from intervalframe import IntervalFrame
from projectframe import ProjectFrame
from sklearn import linear_model

# Local imports
from ..fragments import Fragments
from ..correct.correction import correct, gaussian_smooth
from .smooth_cnv.smooth_cnv import smoothCNV
from ..coverage.coverage import predict_sex, bin_counts
from .cnv_utilities import hmm_classify, merge_segments, validate_distributions, train_hmm


def sex_prediction(data: Fragments | IntervalFrame, 
                    bin_size: int = 100000):
    """
    Call copy number variations
    """

    # Determine bin coverage
    if isinstance(data, IntervalFrame):
        bin_coverage = data
    else:
        bin_coverage = bin_counts(data, bin_size=bin_size)

    # Predict sex
    predicted_sex = predict_sex(bin_coverage)

    return predicted_sex


def call_cnvs(data: Fragments | IntervalFrame, 
              bin_size: int = 100000,
              genome_version: str = "hg19",
              method: str = "online_both",
              outlier_smooth: bool = True,
              gauss_smooth: bool = False,
              use_normal: bool = False,
              keep_sex_chroms: bool = False,
              verbose: bool = False,
              **kwargs):
    """
    Call copy number variations
    """

    # Get main chromosomes
    import genome_info
    g = genome_info.GenomeInfo(genome_version)
    autosomes = g["autosomes"]

    # Determine bin coverage
    if verbose: print("Binning...")
    if isinstance(data, IntervalFrame):
        bin_coverage = data
    else:
        bin_coverage = bin_counts(data, bin_size=bin_size)

    # Remove chromosomes
    #if len(remove_chroms) > 0:
    #    bin_coverage = bin_coverage.iloc[~bin_coverage.index.labels.isin(remove_chroms),:]

    # Predict sex
    if verbose: print("Predicting sex...")
    predicted_sex = predict_sex(bin_coverage)

    # Remove 0s
    bin_coverage = bin_coverage.iloc[bin_coverage.loc[:,"counts"].values > 0,:]

    # Correct values
    if verbose: print("Correcting genome biases...")
    bin_coverage = correct(bin_coverage,
                           genome_version=genome_version,
                           bin_size=bin_size,
                           use_normal=use_normal,
                           keep_sex_chroms=keep_sex_chroms,
                           verbose=verbose)

    # Determine  median
    median = np.median(bin_coverage.loc[autosomes,:].df.loc[:,"corrected_counts"].values)

    # Determine ratios
    if verbose: print("Normalizing...")
    bin_coverage.df.loc[:,"ratios"] = 0.0
    bin_coverage.df.loc[:,"ratios"] = np.log2(bin_coverage.df.loc[:, "corrected_counts"].values / median)
    bin_coverage.df.loc[:,"ratios"].values[np.isinf(bin_coverage.df.loc[:,"ratios"].values)] = 0

    # Smooth ratios
    if verbose: print("Smoothing...")
    for chrom in bin_coverage.index.unique_labels:
        if outlier_smooth:
            bin_coverage.loc[chrom,"ratios"] = smoothCNV(np.repeat(chrom, bin_coverage.index.label_counts[chrom]),
                                                         bin_coverage.loc[chrom,"ratios"].values)
        if gauss_smooth:
            bin_coverage.loc[chrom,"ratios"] = gaussian_smooth(bin_coverage.loc[chrom,"ratios"].values,
                                                               scale=1, p_min=1e-10)

    # Remove nans
    if verbose: print("Filtering...", flush=True)
    bin_coverage = bin_coverage.iloc[~pd.isnull(bin_coverage.loc[:,"ratios"].values),:]

    # Correct chr19 if necessary
    if genome_version == "hg38" or genome_version == "hg19":
        if "chr19" in bin_coverage.index.unique_labels:
            if np.median(bin_coverage.df.loc[bin_coverage.index.labels=="chr19","ratios"].values) < -0.05:
                bin_coverage.df.loc[bin_coverage.index.labels=="chr19","ratios"] = bin_coverage.df.loc[bin_coverage.index.labels=="chr19","ratios"].values + 0.07
        
    # Calculate segments
    if verbose: print("Segmenting...")
    cnv_segments = bin_coverage.segment("ratios", method=method, **kwargs)
    # Annotate segments
    cnv_segments.annotate(bin_coverage, "ratios", method="median")

    return bin_coverage, cnv_segments


def process_cnvs(segments,
                 bins,
                 merge=False,
                 normal = [0.1, 0.5, 0.9],
                 ploidy = [2],
                 estimatePloidy=False,
                 minSegmentBins=25,
                 maxCN=7,
                 n_mads = 1.4826,
                 scStates = [1, 3],
                 verbose=False,
                 **kwargs):
    """
    Process copy number segment calls
    """

    # Merge segments
    if merge:
        segments = merge_segments(segments, bins)

    # Classify calls be HMM
    hmm_states = train_hmm(bins,
                           normal=normal,
                           ploidy=ploidy,
                           gender=None,
                           estimatePloidy=estimatePloidy,
                           minSegmentBins=minSegmentBins,
                           maxCN=maxCN,
                           scStates=scStates,
                           verbose=verbose,
                           **kwargs)
    
    segments = hmm_classify(segments, hmm_states)
    validate_distributions(segments, bins, n_mads=n_mads)

    return segments
