import pandas as pd
import numpy as np
import os
from intervalframe import IntervalFrame
from projectframe import ProjectFrame
from sklearn import linear_model

# Local imports
from ..fragments import Fragments
from ..correct.correction import correct, gaussian_smooth
from .smooth_cnv.smooth_cnv import smoothCNV
from ..coverage.coverage import predict_sex, bin_counts
from .cnv_utilities import hmm_classify, merge_segments, validate_distributions, train_hmm


def determine_normal(data: IntervalFrame,
                     bin_size: int = 100000,
                     genome: str = "hg19",
                     method: str = "bcp_online_both",
                     outlier_smooth: bool = True,
                     gauss_smooth: bool = False,
                     verbose: bool = False,
                     **kwargs):
    """
    """

    bin_coverage = data.loc[:,[data.columns.values[0]]]
    bin_coverage.df.columns = ["counts"]
    bin_coverage = correct(bin_coverage, genome=genome, bin_size=bin_size, verbose=verbose)
    #mean_values = bin_coverage.df.loc[:,"corrected_counts"].values.copy()

    normal_corrected = IntervalFrame(bin_coverage.index)
    normal_corrected.df.loc[:,data.columns.values[0]] = bin_coverage.df.loc[:,"corrected_counts"].values

    for sample in data.df.columns.values[1:]:
        bin_coverage = data.loc[:,[sample]]
        bin_coverage.df.columns = ["counts"]
        bin_coverage = correct(bin_coverage, genome=genome, bin_size=bin_size, verbose=verbose)
        #mean_values = mean_values + bin_coverage.df.loc[:,"corrected_counts"].values
        normal_corrected.df.loc[:,sample] = bin_coverage.df.loc[:,"corrected_counts"].values

    #mean_values = mean_values / data.shape[1]

    return normal_corrected


def call_cnvs(data: Fragments | IntervalFrame, 
              norm_data: Fragments | IntervalFrame = None,
              bin_size: int = 100000,
              genome: str = "hg19",
              method: str = "bcp_online_both",
              outlier_smooth: bool = True,
              gauss_smooth: bool = False,
              verbose: bool = False,
              **kwargs):
    """
    Call copy number variations
    """

    # Determine bin coverage
    if verbose: print("Binning...")
    if isinstance(data, IntervalFrame):
        bin_coverage = data
    else:
        bin_coverage = bin_counts(data, bin_size=bin_size)

    # Predict sex
    if verbose: print("Predicting sex...")
    predicted_sex = predict_sex(bin_coverage)

    # Correct values
    if verbose: print("Correcting genome biases...")
    bin_coverage = correct(bin_coverage, genome=genome, bin_size=bin_size, verbose=verbose)

    # Determine  normal
    if verbose: print("Determining normal...")
    if norm_data is not None:
        if isinstance(norm_data, IntervalFrame):
            bin_coverage = bin_coverage.exact_match(norm_data)
            norm_data = norm_data.exact_match(bin_coverage)
            median = determine_normal(norm_data,
                                      bin_size=bin_size,
                                      genome=genome,
                                      method=method,
                                      outlier_smooth=outlier_smooth,
                                      gauss_smooth=gauss_smooth,
                                      verbose=verbose,
                                      **kwargs)
        else:
            norm_bin_coverage = bin_counts(norm_data, bin_size=bin_size)
            norm_bin_coverage = correct(norm_bin_coverage, genome=genome,
                                        bin_size=bin_size, verbose=verbose)
            median = norm_bin_coverage.df.loc[:,"corrected_counts"].values
    else:
        median = np.median(bin_coverage.df.loc[:,"corrected_counts"].values)

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
        
    # Calculate segments
    if verbose: print("Segmenting...")
    cnv_segments = bin_coverage.segment("ratios", method=method, **kwargs)
    # Annotate segments
    cnv_segments.annotate(bin_coverage, "ratios", method="median")

    return bin_coverage, cnv_segments


def call_cnvs_multi(data: IntervalFrame, 
                    norm_data: IntervalFrame = None,
                    bin_size: int = 100000,
                    genome: str = "hg19",
                    method: str = "bcp_online_both",
                    outlier_smooth: bool = True,
                    gauss_smooth: bool = False,
                    hmm_call: bool = True,
                    median_shift: bool = True,
                    verbose: bool = False,
                    **kwargs):
    """
    Call copy number variations
    """

    # Filter to same normal and test bins
    bin_coverage = data
    if norm_data is not None:
        bin_coverage = bin_coverage.exact_match(norm_data)
        norm_data = norm_data.exact_match(bin_coverage)
        median = determine_normal(norm_data,
                                        bin_size=bin_size,
                                        genome=genome,
                                        method=method,
                                        outlier_smooth=outlier_smooth,
                                        gauss_smooth=gauss_smooth,
                                        verbose=verbose,
                                        **kwargs)

    # Log chromosome lengths
    chrom_ranges = bin_coverage.index.label_ranges
    chrom_lengths = {chrom:chrom_ranges[chrom][1] for chrom in chrom_ranges}
    
    # Iterate over samples
    pf = ProjectFrame()
    pf.uns["chrom_lengths"] = chrom_lengths
    for sample in bin_coverage.df.columns.values:
        pf.add_obs(sample)
        sample_coverage = bin_coverage.loc[:,[sample]]
        sample_coverage.df.columns = ["counts"]
    
        # Predict sex
        if verbose: print("Predicting sex...")
        predicted_sex = predict_sex(sample_coverage)
        pf.add_anno("predicted_gender", sample, predicted_sex["gender"])

        # Correct values
        if verbose: print("Correcting genome biases...")
        sample_coverage = correct(sample_coverage, genome=genome, bin_size=bin_size, verbose=verbose)

        if norm_data is None:
            median = np.median(sample_coverage.df.loc[:,"corrected_counts"].values)

        # Determine ratios
        if verbose: print("Normalizing...")
        sample_coverage.df.loc[:,"ratios"] = 0.0
        sample_coverage.df.loc[:,"ratios"] = np.log2(sample_coverage.df.loc[:, "corrected_counts"].values / median)
        sample_coverage.df.loc[:,"ratios"].values[np.isinf(sample_coverage.df.loc[:,"ratios"].values)] = 0

        # Smooth ratios
        if verbose: print("Smoothing...")
        for chrom in sample_coverage.index.unique_labels:
            if outlier_smooth:
                sample_coverage.loc[chrom,"ratios"] = smoothCNV(np.repeat(chrom, sample_coverage.index.label_counts[chrom]),
                                                                sample_coverage.loc[chrom,"ratios"].values)
            if gauss_smooth:
                sample_coverage.loc[chrom,"ratios"] = gaussian_smooth(sample_coverage.loc[chrom,"ratios"].values,
                                                                      scale=1, p_min=1e-10)

        # Remove nans
        if verbose: print("Filtering...", flush=True)
        sample_coverage = sample_coverage.iloc[~pd.isnull(sample_coverage.loc[:,"ratios"].values),:]

        # Median shift
        if median_shift and norm_data is not None:
            sample_coverage.df.loc[:,"ratios"] = sample_coverage.df.loc[:,"ratios"].values - np.median(sample_coverage.df.loc[:,"ratios"].values)
            
        # Calculate segments
        if verbose: print("Segmenting...")
        cnv_segments = sample_coverage.segment("ratios", method=method, **kwargs)
        # Annotate segments
        cnv_segments.annotate2(sample_coverage, "ratios", method="median")

        # Call HMM
        if hmm_call:
            cnv_segments = process_cnvs(cnv_segments, sample_coverage, merge=True)

        # Add to Project frame
        sample_bins = sample_coverage.loc[:,["ratios"]]
        sample_bins.df.columns = [sample]
        pf.add_intervals("cnv_bins", sample_bins)
        pf.add_obs_intervals(sample, "cnv_segments", cnv_segments)

    return pf


def process_cnvs(segments,
                 bins,
                 merge=False,
                 normal = [0.1, 0.5, 0.9],
                 ploidy = [2],
                 estimatePloidy=False,
                 minSegmentBins=25,
                 maxCN=7,
                 n_mads = 1.4826,
                 verbose=False,
                 **kwargs):
    """
    Process copy number segment calls
    """

    # Merge segments
    if merge:
        segments = merge_segments(segments, bins)

    # Classify calls be HMM
    hmm_states = train_hmm(bins, normal=normal, ploidy=ploidy,
                                     gender=None, estimatePloidy=estimatePloidy,
                                     minSegmentBins=minSegmentBins, maxCN=maxCN,
                                     verbose=verbose, **kwargs)
    segments = hmm_classify(segments, hmm_states)
    validate_distributions(segments, bins, n_mads=n_mads)

    return segments
