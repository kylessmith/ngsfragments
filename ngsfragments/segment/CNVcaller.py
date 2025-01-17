from typing import List
import math
from ailist import LabeledIntervalArray
from intervalframe import IntervalFrame
from projectframe import ProjectFrame
import genome_info
import pandas as pd
import numpy as np
from joblib import delayed, Parallel
import statsmodels.api as sm
from typing import Dict

# Local imports
from ..fragments import Fragments
from ..coverage.coverage import predict_sex, bin_counts
from ..correct.correction import correct, gaussian_smooth, binned_bias_correct_counts, match_bins, filter_bins, gc_mappability_bin_correct
from .smooth_cnv.smooth_cnv import smoothCNV
from .merge_regions.merge_regions import merge_adjacent
from .cnv_utilities import train_hmm, hmm_classify, validate_distributions


def calculate_MAD(vector: np.ndarray):
    """
    Calculate Median Absolute Deviation (MAD)

    Parameters
    ----------
        vector : np.ndarray
            Vector of values

    Returns
    -------
        mad : float
            Median Absolute Deviation
    """

    # Calculate MAD
    mad = np.median(np.abs(vector - np.median(vector)))

    return mad


def iterative_merge(segments: IntervalFrame,
                    bins: IntervalFrame,
                    column: str = "ratios",
                    merge_MAD: float = 1.4826):
    """
    Iteratively merge segments

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        column : str
            Column name

    Returns
    -------
        merged : IntervalFrame
            IntervalFrame of merged segments
    """

    # Initialize merged values
    merged = LabeledIntervalArray()
    merged_values = []

    for chrom in segments.index.unique_labels:
        chrom = str(chrom)
        chrom_segs = segments.loc[chrom,:]

        # Extract column Series
        ratios = bins.loc[chrom, column]

        # Record current interval
        previous_start = chrom_segs.index[0].start
        previous_end = chrom_segs.index[0].end
        segment_median1 = chrom_segs.df.loc[:,"median"].values[0]
        
        segment_mad1 = merge_MAD * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
        segment_bounds1 = (segment_median1 - segment_mad1,
                        segment_median1 + segment_mad1)

        # Iterate over intervals
        for i in range(1, chrom_segs.shape[0]):
            # Record next interval
            start = int(chrom_segs.index[i].start)
            end = int(chrom_segs.index[i].end)
            segment_median2 = float(chrom_segs.df.loc[:,"median"].values[i])

            # Check if medians overlap
            if segment_median2 >= segment_bounds1[0] and segment_median2 <= segment_bounds1[1]:
                # Update current interval
                previous_end = end
                segment_median1 = np.median(ratios.intersect(previous_start, previous_end, chrom).values)
                segment_mad1 = merge_MAD * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
                segment_bounds1 = (segment_median1 - segment_mad1,
                                segment_median1 + segment_mad1)

            else:
                merged.add(previous_start, previous_end, chrom)
                merged_values.append(segment_median1)
                previous_start = chrom_segs.index[i].start
                previous_end = chrom_segs.index[i].end
                segment_median1 = float(chrom_segs.df.loc[:,"median"].values[i])
                segment_mad1 = merge_MAD * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
                segment_bounds1 = (segment_median1 - segment_mad1,
                                segment_median1 + segment_mad1)

        # Add last interval
        merged.add(previous_start, previous_end, chrom)
        merged_values.append(segment_median1)

    # Create IntervalFrame
    merged_iframe = IntervalFrame(merged)
    merged_iframe.loc[:,"median"] = merged_values

    return merged_iframe


def merge_segments(segments: IntervalFrame,
                   bins: IntervalFrame,
                   column: str = "ratios",
                   merge_MAD: float = 1.4826):
    """
    Merge segments

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        column : str
            Column name

    Returns
    -------
        merged : IntervalFrame
            IntervalFrame of merged segments
    """

    # Initial values
    previous_length = segments.shape[0]
    current_length = 1

    # Check if there is 1 segment
    if previous_length == current_length:
        return segments

    # Continue merging until no more changes
    merged = segments.copy()
    while previous_length != current_length:
        previous_length = current_length
        merged = iterative_merge(merged, bins, column=column, merge_MAD=merge_MAD)
        current_length = merged.shape[0]

    return merged


def data_binning(intervals: LabeledIntervalArray,
                  n: int = 10) -> IntervalFrame:
    """
    Bin probes into bins of size n.

    Parameters
    ----------
        intervals : LabeledIntervalArray
            IntervalFrame of probes.
        n : int
            Number of probes per bin.
    
    Returns
    -------
        IntervalFrame
            IntervalFrame of bins.
    """

    # Initialize bins
    bins = LabeledIntervalArray()

    # Iterate through chromosomes
    for chrom in intervals.unique_labels:
        starts = intervals.get(chrom).starts
        ends = intervals.get(chrom).ends
        n_bins = math.ceil(len(starts) / n)

        for i in range(n, len(starts), n):
            bins.add(starts[i-n], ends[i], chrom)

    # Construct IntervalFrame
    bins_iframe = IntervalFrame(intervals=bins)

    return bins_iframe


def merge_data_bins(intervals: LabeledIntervalArray) -> IntervalFrame:
    """
    Merge bins that are adjacent bins smaller than 1kb.

    Parameters
    ----------
        intervals : LabeledIntervalArray
            IntervalFrame of bins.
        
    Returns
    -------
        IntervalFrame
            IntervalFrame of merged bins.
    """

    # Initialize bins
    new_bins = LabeledIntervalArray()

    # Investigate lengths
    starts = intervals.starts
    ends = intervals.ends
    labels = intervals.labels
    lengths = ends - starts

    # Iterate through bins
    i = 0
    while i < len(lengths):
        # If bin is less than 1kb and adjacent to another bin of the same label, merge
        if lengths[i] < 1000 and i < len(lengths) - 1:
            if labels[i] == labels[i+1]:
                new_bins.add(starts[i], ends[i+1], labels[i])
                i += 2
            else:
                new_bins.add(starts[i], ends[i], labels[i])
                i += 1
        
        elif lengths[i] < 1000000:
            new_bins.add(starts[i], ends[i], labels[i])
            i += 1
        else:
            i += 1
    
    # Construct to IntervalFrame
    bins_iframe = IntervalFrame(intervals=new_bins)

    return bins_iframe


def annotate_bins(iframe,
                    bins = None,
                    bin_size=100000,
                    method="mean"):
    """
    """

    # Chosen function
    if method == "mean":
        func = np.mean
    elif method == "sum":
        func = np.sum
    elif method == "nansum":
        func = np.nansum
    elif method == "median":
        func = np.median
    elif method == "std":
        func = np.std

    if bins is None:
        # Get ranges
        ranges = iframe.index.label_ranges
        range_dict = {chrom:ranges[chrom][1] for chrom in ranges}

        # Create bins
        bin_iframe = IntervalFrame.from_dict_range(range_dict, bin_size)
    else:
        bin_iframe = bins

    values = pd.DataFrame(np.zeros((bin_iframe.df.shape[0],iframe.shape[1])),
                          columns = iframe.columns.values)
    values[:] = np.nan
    for i, index in enumerate(iframe.index.iter_intersect(bin_iframe.index, return_intervals=False, return_index=True)):
        value = iframe.df.iloc[index,:].values
        values.values[i,:] = func(value, axis=0)

    bin_iframe.df = values
    
    return bin_iframe


def regress_normal(values: np.ndarray,
                   normal_values: np.ndarray) -> np.ndarray:
    """
    Regress out normal values from values.

    Parameters
    ----------
        values : np.ndarray
            Values to regress out normal values from.
        normal_values : np.ndarray
            Normal values to regress out.
    
    Returns
    -------
        np.ndarray
            Residuals from regression.
    """

    # Remove nans
    n = len(values)
    not_nan = ~pd.isnull(values)
    values = values[not_nan]
    normal_values = normal_values[not_nan]

    # Regress out normal
    result = sm.GLM(values, normal_values, family=sm.families.Gaussian()).fit()
    residuals = result.predict(normal_values)
    ratios = values / residuals

    # Add nans
    new_ratios = np.zeros(n)
    new_ratios[:] = np.nan
    new_ratios[not_nan] = ratios


    return new_ratios


class CNVcaller(object):
    """
    Copy Number Variation caller
    """

    def __init__(self,
                 cnv_binsize: int = 100000,
                hmm_binsize: int = 1000000,
                method: str = "online_both",
                outlier_smooth: str = True,
                gauss_smooth: bool = False,
                bcp_cutoff: float = 0.3,
                normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                ploidy: List[int] = [2, 3, 4],
                estimatePloidy: bool = True,
                minSegmentBins: int = 25,
                maxCN: int = 5,
                hazard: int = 100,
                shuffles: int = 5000,
                p: float = 0.00005,
                use_normal: bool = False,
                keep_sex_chroms: bool = False,
                scStates: List[int] = [1, 3],
                n_mads: float = 1.4826,
                genome_version: str = "hg38",
                n_per_bin: int = 10,
                n_per_bin_hmm: int = 15,
                verbose: bool = False):
        """
        Initialize CNVcaller
        """

        self.pf = ProjectFrame()
        self.pf.params["ref_genome"] = genome_version
        self.genome = genome_info.GenomeInfo(genome_version)

        self.cnv_binsize = cnv_binsize
        self.hmm_binsize = hmm_binsize
        self.method = method
        self.outlier_smooth = outlier_smooth
        self.gauss_smooth = gauss_smooth
        self.bcp_cutoff = bcp_cutoff
        self.normal = normal
        self.ploidy = ploidy
        self.estimatePloidy = estimatePloidy
        self.minSegmentBins = minSegmentBins
        self.maxCN = maxCN
        self.hazard = hazard
        self.shuffles = shuffles
        self.p = p
        self.use_normal = use_normal
        self.keep_sex_chroms = keep_sex_chroms
        self.scStates = scStates
        self.n_mads = n_mads
        self.n_per_bin = n_per_bin
        self.n_per_bin_hmm = n_per_bin_hmm
        self.verbose = verbose


    def bin_data(self,
                 data: Fragments | IntervalFrame | Dict[str, Fragments],
                 use_hmm_bins: bool = False,
                 remove_nan: bool = False,
                 merge: bool = True) -> IntervalFrame:
        """
        Bin data into bins of size n_per_bin
        """

        # Check if data is Fragments or IntervalFrame
        if isinstance(data, Fragments):
            if use_hmm_bins:
                bins = bin_counts(data,
                                  bin_size = self.hmm_binsize)
            else:
                bins = bin_counts(data,
                                bin_size = self.cnv_binsize)
                
            bins.df.columns = [data.sam_file.split(".")[0]]

        elif isinstance(data, IntervalFrame):
            if use_hmm_bins:
                #bins = data_binning(data.index,
                #                    n = self.n_per_bin_hmm)
                bins = merge_adjacent(data.index, n = self.n_per_bin_hmm)
                bins = IntervalFrame(intervals=bins)
            else:
                #bins = data_binning(data.index,
                #                    n = self.n_per_bin)
                bins = merge_adjacent(data.index, n = self.n_per_bin)
                bins = IntervalFrame(intervals=bins)
            if merge:
                bins = merge_data_bins(intervals = bins.index)

            # Sum values in bins
            bins = annotate_bins(data,
                                bins = bins,
                                method = "nansum")
            
        elif isinstance(data, dict):
            # Get bins
            if use_hmm_bins:
                bins = LabeledIntervalArray.create_bin(self.genome["chrom_sizes"], bin_size=self.hmm_binsize)
            else:
                bins = LabeledIntervalArray.create_bin(self.genome["chrom_sizes"], bin_size=self.cnv_binsize)
            bins = IntervalFrame(intervals=bins)

            # Determine nhits per Fragments
            for key in data:
                bins.df.loc[:,key] = data[key].frags.nhits_from_LabeledIntervalArray(bins.index)
            
        # Remove nans
        if remove_nan:
            selected = np.sum(pd.isnull(bins.values), axis=1) == 0
            bins = bins.iloc[selected,:]

        return bins
    

    def correct_bins(self,
                    bins: IntervalFrame,
                    normal_bins: IntervalFrame = None,
                    n_jobs: int = 1,
                    remove_nan: bool = True,
                    verbose: bool = False):
        """
        Correct bins
        """

        # Make sure bins have main chromosomes
        bins = bins.loc[self.genome["main_chromosomes"],:]

        # Match normals
        if normal_bins is not None:
            normal_bins = normal_bins.loc[self.genome["main_chromosomes"],:]
            bins = bins.exact_match(normal_bins)
            normal_bins = normal_bins.exact_match(bins)

        # Initialize bias records
        bias_record = self.genome.calculate_bias(bins.index)

        # Correct bins
        bins, bias_record = match_bins(bins,
                                       bias_record)
        bins, bias_record = filter_bins(bins,
                                        bias_record,
                                         blacklist_cutoff = 0.1,
                                         mapping_cutoff = 0.9,
                                         keep_sex_chroms = self.keep_sex_chroms)
        
        # Make sure bins is of correct dtype
        bins.df = bins.df.astype(float)

        # Correct genome biases
        if verbose: print("Correcting genome biases", flush=True)
        if n_jobs == 1:
            for sample in bins.df.columns:
                #bins.df.loc[:,sample] = binned_bias_correct_counts(bins.df.loc[:,sample].values,
                #                                                                bias_record)
                bins.df.loc[:,sample] = gc_mappability_bin_correct(bins.df.loc[:,sample].values,
                                                                bias_record)
        else:
            #corrected_values = Parallel(n_jobs=n_jobs)(delayed(binned_bias_correct_counts)(bins.values[:,i], bias_record) for i in range(bins.shape[1]))
            corrected_values = Parallel(n_jobs=n_jobs)(delayed(gc_mappability_bin_correct)(bins.values[:,i], bias_record) for i in range(bins.shape[1]))
            bins.df[:] = np.array(corrected_values).T
        
        # Correct normal_bins
        if normal_bins is not None:
            if verbose: print("Correcting normal biases", flush=True)
            # Match bins
            normal_bins = normal_bins.exact_match(bins)
            bins = bins.exact_match(normal_bins)

            # Regress out bias from normals
            if n_jobs == 1:
                for normal_sample in normal_bins.df.columns:
                    #normal_bins.df.loc[:,normal_sample] = binned_bias_correct_counts(normal_bins.df.loc[:,normal_sample].values,
                    #                                                                bias_record)
                    normal_bins.df.loc[:,normal_sample] = gc_mappability_bin_correct(normal_bins.df.loc[:,normal_sample].values,
                                                                                    bias_record)
            else:
                #corrected_values = Parallel(n_jobs=n_jobs)(delayed(binned_bias_correct_counts)(normal_bins.values[:,i], bias_record) for i in range(normal_bins.shape[1]))
                corrected_values = Parallel(n_jobs=n_jobs)(delayed(gc_mappability_bin_correct)(normal_bins.values[:,i], bias_record) for i in range(normal_bins.shape[1]))
                normal_bins.df[:] = np.array(corrected_values).T

            # Regress out normals
            if n_jobs == 1:
                for sample in bins.df.columns:
                    bins.df.loc[:,sample] = regress_normal(bins.df.loc[:,sample].values, normal_bins.df.values)
                    #bins.df.loc[:,sample] = binned_bias_correct_counts(bins.df.loc[:,sample].values,
                    #                                                                normal_bins)
            else:
                corrected_values = Parallel(n_jobs=n_jobs)(delayed(regress_normal)(bins.values[:,i], normal_bins.df) for i in range(bins.shape[1]))
                bins.df[:] = np.array(corrected_values).T

        # Remove bins with no signal
        bins.df.fillna(0, inplace=True)

        # Calculate ratios
        autosomes = self.genome["autosomes"]
        median = np.median(bins.loc[autosomes,:].df.values, axis=0)
        bins.df.iloc[:,:] = np.log2(bins.df.values / median)
        bins.df[np.isinf(bins.df.values)] = 0

        # Smooth ratios
        if verbose: print("Smoothing bins", flush=True)
        for sample in bins.df.columns:
            if verbose: print("   " +  sample, flush=True)
            for chrom in bins.index.unique_labels:
                if self.outlier_smooth:
                    bins.loc[chrom, sample] = smoothCNV(np.repeat(chrom, bins.index.label_counts[chrom]),
                                                                bins.loc[chrom, sample].values)
        
        # Remove nans
        if remove_nan:
            if verbose: print("Removing NaN's", flush=True)
            bins = bins.iloc[np.sum(pd.isnull(bins.values), axis=1) == 0,:]

        # Correct chr19 if necessary
        if "chr19" in bins.index.unique_labels:
            for sample in bins.df.columns:
                chr19 = np.median(bins.df.loc[bins.index.labels=="chr19",sample].values)
                if chr19 <= -0.05:
                    chr19_correction = min([0.1, chr19 * -1])
                    bins.df.loc[bins.index.labels=="chr19",sample] = bins.df.loc[bins.index.labels=="chr19",sample].values + chr19_correction

        return bins
    
    
    def predict_purity(self,
                       hmm_bins: IntervalFrame,
                       normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                       ploidy: List[int] = [2, 3, 4],
                       estimatePloidy: bool = True,
                       scStates: List[int] = [1, 3],
                       record: bool = True,
                       merge: bool = True,
                       merge_MAD: float = 1.4826,
                       verbose: bool = False):
        """
        """

        if verbose: print("Segmenting for HMM...", flush=True)
        hmm_states = {}
        for sample in hmm_bins.df.columns:
            if verbose: print("   " + sample, flush=True)
            # Segment bins
            hmm_cnv_segments = hmm_bins.segment(sample, method=self.method, cutoff=self.bcp_cutoff)
            # Annotate segments
            hmm_cnv_segments.annotate(hmm_bins, sample, method="median")
            if merge:
                hmm_cnv_segments = merge_segments(hmm_cnv_segments, hmm_bins.loc[:,[sample]], column = sample, merge_MAD = merge_MAD)
            #hmm_cnv_segments = merge_segments(hmm_cnv_segments, hmm_bins, column = sample)

            sample_bins = hmm_bins.loc[:,[sample]]
            sample_bins.df.columns = ["ratios"]

            # Classify calls be HMM
            hmm_states[sample] = train_hmm(sample_bins,
                                            normal = normal,
                                            ploidy = ploidy,
                                            gender = None,
                                            estimatePloidy = estimatePloidy,
                                            minSegmentBins = self.minSegmentBins,
                                            maxCN = self.maxCN,
                                            scStates = scStates,
                                            verbose = self.verbose)

            # Assign attributes
            if record:
                # Check all CNVs are not NEUT
                processed_cnvs = hmm_classify(hmm_cnv_segments, hmm_states[sample])
                processed_cnvs = validate_distributions(processed_cnvs, hmm_bins, column = sample)
                if np.sum(processed_cnvs.df.loc[:,"Corrected_Call"].values == "NEUT") == processed_cnvs.df.shape[0]:
                    self.pf.add_anno("purity", sample, 0.0)
                    self.pf.add_anno("ploidy", sample, hmm_states[sample]["phi"])
                    self.pf.add_anno("clonal", sample, 0.0)
                else:
                    self.pf.add_anno("purity", sample, 1 - hmm_states[sample]["n"])
                    self.pf.add_anno("ploidy", sample, hmm_states[sample]["phi"])
                    self.pf.add_anno("clonal", sample, hmm_states[sample]["Frac_genome_subclonal"])

        return hmm_states

    
    def predict_cnvs(self,
                    data: Fragments | IntervalFrame | Dict[str, Fragments],
                    normal_data: Fragments | IntervalFrame | Dict[str, Fragments] = None,
                    prebinned: bool = False,
                    merge: bool = True,
                    merge_MAD: float = 1.4826,
                    verbose: bool = True):
        """
        Fit HMM for copy number calling
        """

        # Get chrom lengths
        if isinstance(data, Fragments):
            chrom_ranges = data.frags.label_ranges
        elif isinstance(data, dict):
            key = list(data.keys())[0]
            chrom_ranges = data[key].frags.label_ranges
        else:
            chrom_ranges = data.index.label_ranges
        chrom_lengths = {chrom:chrom_ranges[chrom][1] for chrom in chrom_ranges}
        self.pf.uns["chrom_lengths"] = chrom_lengths
        
        # Bin data
        if prebinned == False:
            if verbose: print("Binning for HMMs...", flush=True)
            if self.hmm_binsize != self.cnv_binsize:
                hmm_bins = self.bin_data(data = data,
                                    use_hmm_bins = True)
            if verbose: print("Binning for CNVs...", flush=True)
            bins = self.bin_data(data = data,
                                use_hmm_bins = False)

        # Correct bins
        if normal_data is not None:

            if self.hmm_binsize != self.cnv_binsize:
                if prebinned == False:
                    if verbose: print("Binning normals for HMM...", flush=True)
                    if isinstance(normal_data, IntervalFrame):
                        normal_hmm_bins = IntervalFrame(intervals=hmm_bins.index)
                        normal_hmm_bins = annotate_bins(normal_data, normal_hmm_bins, method="nansum")
                    else:
                        normal_hmm_bins = self.bin_data(data = normal_data,
                                                use_hmm_bins = True)
                    
                hmm_bins = self.correct_bins(bins = hmm_bins,
                                        normal_bins = normal_hmm_bins,
                                        n_jobs = 1)
            if prebinned == False:
                if verbose: print("Binning normals for CNVs...", flush=True)
                if isinstance(normal_data, IntervalFrame):
                    normal_bins = IntervalFrame(intervals=bins.index)
                    normal_bins = annotate_bins(normal_data, normal_bins, method="nansum")
                else:
                    normal_bins = self.bin_data(data = normal_data,
                                            use_hmm_bins = False)
            bins = self.correct_bins(bins = bins,
                                    normal_bins = normal_bins,
                                    n_jobs = 1)
        else:
            if self.hmm_binsize != self.cnv_binsize:
                if verbose: print("Correcting for HMMs...", flush=True)
                hmm_bins = self.correct_bins(bins = hmm_bins,
                                            n_jobs = 1)
            if verbose: print("Correcting for CNVs...", flush=True)
            bins = self.correct_bins(bins = bins,
                                        n_jobs = 1)
        
        # Add bins
        self.pf.add_intervals("cnv_bins", bins)

        # Call HMM
        if verbose: print("Predicting purity...", flush=True)
        if self.hmm_binsize != self.cnv_binsize:
            hmm_states = self.predict_purity(hmm_bins,
                                            normal = self.normal,
                                            ploidy = self.ploidy,
                                            estimatePloidy = self.estimatePloidy,
                                            scStates = self.scStates,
                                            record = True,
                                            merge = merge,
                                            merge_MAD = merge_MAD)
            if verbose: print("Classifying segments..", flush=True)
            hmm_states = self.predict_purity(bins,
                                            normal = self.normal,
                                            ploidy = [2],
                                            estimatePloidy = False,
                                            scStates = self.scStates,
                                            record = False,
                                            merge = merge,
                                            merge_MAD = merge_MAD)
        else:
            if verbose: print("Classifying segments..", flush=True)
            hmm_states = self.predict_purity(bins,
                                            normal = self.normal,
                                            ploidy = [2],
                                            estimatePloidy = False,
                                            scStates = self.scStates,
                                            record = True,
                                            merge = merge,
                                            merge_MAD = merge_MAD)


        for sample in bins.df.columns:
            # Classify segments
            cnv_segments = bins.segment(sample, method=self.method, cutoff=self.bcp_cutoff)
            cnv_segments.annotate(bins, sample, method="median")
            if merge:
                cnv_segments = merge_segments(cnv_segments, bins, column = sample, merge_MAD = merge_MAD)
            processed_cnvs = hmm_classify(cnv_segments, hmm_states[sample])
            processed_cnvs = validate_distributions(processed_cnvs, bins, column = sample)
            #processed_cnvs.drop_columns(["copy_number","event","subclone_status","logR_Copy_Number"])

            #sample_bins = hmm_bins.loc[:,[sample]]
            #self.pf.add_intervals("cnv_bins", sample_bins)
            self.pf.add_obs_intervals(sample, "cnv_segments", processed_cnvs)


    def calculate_zscore(self,
                         verbose: bool = False):
        """
        Calculate zscores for CNV bins and segments
        """

        # Iterate over samples
        for sample in self.pf.obs:
            if verbose: print("Calculating zscores for " + sample, flush=True)

            # Get neutral bins
            neut = self.pf.obs_intervals[sample]["cnv_segments"].iloc[self.pf.obs_intervals[sample]["cnv_segments"].df.loc[:,"Corrected_Call"].values == "NEUT",:]
            n = neut.index.nhits_from_LabeledIntervalArray(self.pf.intervals["cnv_bins"].index)
            bins = self.pf.intervals["cnv_bins"].loc[:,[sample]]
            neut_bins = bins.iloc[n==1,:]

            # Calculate zscores
            mean = neut_bins.df.loc[:,sample].values.mean()
            std = neut_bins.df.loc[:,sample].values.std(ddof=1)
            zscores = (bins.df.loc[:,sample].values - mean) / std

            # Summarize zscores
            bins.df.loc[:,sample] = zscores
            segs = self.pf.obs_intervals[sample]["cnv_segments"]
            segs.annotate(bins, sample, method="median", column_name="zscore_median")
            
            self.pf.add_obs_intervals(sample, "cnv_segments", segs)
            self.pf.add_intervals("cnv_zscore_bins", bins)

    
    def rebin(self,
              bin_size: int = 100000,
              pf: ProjectFrame = None) -> ProjectFrame:
        """
        Re-size CNV bins and segments using old bins and segments
        """

        # Iterate over samples
        for sample in self.pf.obs:
            # Create new bins
            new_bins = LabeledIntervalArray.create_bin(self.genome["chrom_sizes"],
                                                    bin_size = bin_size)
            new_bins = IntervalFrame(intervals=new_bins)

            # Annotate new bins
            bins = self.pf.intervals["cnv_bins"].loc[:,[sample]]
            new_bins.annotate(bins, sample, method="mean")
            new_bins = new_bins.iloc[~pd.isnull(new_bins.df.loc[:,"mean"].values),:]

            # Annotate segs
            segs = self.pf.obs_intervals[sample]["cnv_segments"]
            segs.annotate(bins, sample, method="median")

            # Append to projectframe
            if pf is None:
                pf = ProjectFrame()
            pf.add_obs_intervals(sample, "cnv_segments", segs)
            pf.add_intervals("cnv_bins", new_bins)

        return pf
    

    def segment_to_bins(self,
                        bin_size: int = 100000) -> IntervalFrame:
        """
        Use CNV segments to generate bins
        """

        # Iterate over samples
        for sample in self.pf.obs:

            # Create new bins
            new_bins = LabeledIntervalArray.create_bin(self.genome["chrom_sizes"],
                                                    bin_size = bin_size)
            new_bins = IntervalFrame(intervals=new_bins)

            # Annotate new bins
            segs = self.pf.obs_intervals[sample]["cnv_segments"]
            new_bins.annotate(segs, "median", method="mean")
            new_bins = new_bins.iloc[~pd.isnull(new_bins.df.loc[:,"mean"].values),:]

        return new_bins



