import pandas as pd
import numpy as np
import os
from intervalframe import IntervalFrame
from .correction import correct, gaussian_smooth
from .smooth_cnv.smooth_cnv import smoothCNV
from ..utilities import predict_sex
from .cnv_utilities import hmm_classify, merge_segments, validate_distributions, train_hmm


def call_cnvs(frags, bin_size=100000, bin_bias_h5_fn=None, genome="hg19", method="bcp_online_both",
              outlier_smooth=True, gauss_smooth=False, verbose=False, **kwargs):
    """
    Call copy number variations
    """

    # Determine bin coverage
    if verbose: print("Binning...")
    bin_coverage = frags.bin_counts(bin_size=bin_size)

    # Predict sex
    if verbose: print("Predicting sex...")
    predicted_sex = predict_sex(bin_coverage)

    # Correct values
    if verbose: print("Correcting genome biases...")
    # Find distributed h5 reference bins
    if bin_bias_h5_fn is None:
        h5_file = genome + "_ucsc_bin" + str(bin_size) + ".h5"
        if verbose: print("Reference genome:", genome)
        bin_bias_h5_fn = os.path.join(frags.data_dir, h5_file)
        if verbose: print("   Found bin_bias:", h5_file)
        #genome_fn = os.path.join(self.data_dir, genome+"_chrom_sizes.txt")
        #if self.verbose: print("   Found chrom_sizes:", genome+"_chrom_sizes.txt")
        #centro_fn = os.path.join(self.data_dir, genome+"_centro.txt")
        #if self.verbose: print("   Found centromeres:", genome+"_centro.txt")
    bin_coverage = correct(bin_coverage, bin_bias_h5_fn, verbose=verbose)

    # Determine ratios
    if verbose: print("Normalizing...")
    bin_coverage.loc[:,"ratios"] = 0.0
    median = np.median(bin_coverage.loc[:,"corrected_counts"].values)
    bin_coverage.loc[:,"ratios"] = np.log2(bin_coverage.loc[:, "corrected_counts"].values / median)
    bin_coverage.loc[:,"ratios"].values[np.isinf(bin_coverage.loc[:,"ratios"].values)] = 0

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


def process_cnvs(segments, bins, merge=False, normal = [0.1, 0.5, 0.9], ploidy = [2], estimatePloidy=False,
                 minSegmentBins=25, maxCN=7, verbose=False, **kwargs):
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
    validate_distributions(segments, bins)

    return segments