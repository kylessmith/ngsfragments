from .fragments import fragments
from .peak_calling.CallPeaks import call_wps_peaks, normalize_signal
import numpy as np
import pandas as pd
import os
from pathlib import Path
import gzip

# Set data directory
data_dir = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
data_dir = os.path.join(data_dir, "data")


def wps_peaks(frags, wps_dict, shift=0, merge_distance=5, min_length=50, max_length=150):
    """
    """

    peaks = {}
    key_ranges = frags.frags.label_ranges

    for key in wps_dict:
        if len(wps_dict[key]) > 0:
            peaks[key] = call_wps_peaks(wps_dict[key].values, key_ranges[key][0], str(key), shift=shift, merge_distance=merge_distance,
                                        min_length=min_length, max_length=max_length)

    return peaks


def normalize_wps(wps_dict, window=1000, smooth=True, smooth_window=21, polyorder=2, n_threads=1, method="mean", verbose=False):
    """
    """

    for chrom in wps_dict:
        if verbose: print(chrom)
        if method == "median":
            wps_dict[chrom].iloc[:] = normalize_signal(wps_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=False)
        elif method == "mean":
            wps_dict[chrom].iloc[:] = normalize_signal(wps_dict[chrom].values, window, smooth, smooth_window, polyorder, n_threads, use_mean=True)


def window_scores(scores, bed_fn="hg19_TSS.bed.gz", upstream=1000, downstream=1000, stranded=False, verbose=False):
    """
    Calculate scores for each postion within a window
    """

    # Determine bed file
    if Path(bed_fn).is_file() == False:
        bed_fn = os.path.join(data_dir, bed_fn)
        if Path(bed_fn).is_file():
            if verbose: print("Found file:", bed_fn)
        else:
            raise FileExistsError("Could not find bed file.")

    # Initialize window values
    window_length = upstream + downstream
    window_total_values = []
    genes = []

    # Determine ranges for each chromosome
    ranges = {chrom:(scores[chrom].index.values[0], scores[chrom].index.values[-1]) for chrom in scores}

    # Check is gzipped
    open_func = open
    is_gzip = False
    if bed_fn.endswith(".gz"):
        is_gzip = True
        open_func = gzip.open
    
    # Iterate over BED records
    for i, line in enumerate(open_func(bed_fn, "r")):
        if is_gzip:
            line = line.decode()
        field = line.strip().split("\t") # find fields
        chrom = field[0]
        gene = field[-1]
        center = int((int(field[1]) + int(field[2])) / 2)
        
        # Determine if chrom is in values
        try:
            scores[chrom]
        except KeyError:
            continue

        # Check if region captured in scores
        if int(field[2]) < ranges[chrom][0] or int(field[1]) > ranges[chrom][1]:
            continue
        
        # Index window scores
        start_shift = 0
        end_shift = window_length
        window_values = np.zeros(window_length)
        if stranded and field[3] == "-":
            #start_index = max((center - downstream) - ranges[chrom[0]], 0)
            start_index = (center - downstream) - ranges[chrom][0]
            if start_index < 0:
                start_shift = 0 - start_index
                start_index = 0

            #end_index = min((center + upstream) - ranges[chrom][0], len(scores[chrom]))
            end_index = (center + upstream) - ranges[chrom][0]
            if end_index > len(scores[chrom]):
                end_shift = len(scores[chrom]) - end_index
                end_index = len(scores[chrom])
            
            window_values[start_shift:end_shift] = scores[chrom].values[start_index:end_index]
            window_values = window_values[::-1]
        
        else:
            #start_index = max((center - upstream) - ranges[chrom][0], 0)
            start_index = (center - upstream) - ranges[chrom][0]
            if start_index < 0:
                start_shift = 0 - start_index
                start_index = 0

            #end_index = min((center + downstream) - ranges[chrom][0], len(scores[chrom]))
            end_index = (center + downstream) - ranges[chrom][0]
            if end_index > len(scores[chrom]):
                end_shift = len(scores[chrom]) - end_index
                end_index = len(scores[chrom])

            #print(chrom, field[1], field[2], start_index, end_index, start_shift, end_shift)
            window_values[start_shift:end_shift] = scores[chrom].values[start_index:end_index]
        
        window_total_values.append(window_values)
        if stranded and len(field) >= 5:
            genes.append(gene)
        elif stranded == False and len(field) >= 4:
            genes.append(gene)

    # Convert to numpy array
    window_total_values = pd.DataFrame(np.array(window_total_values))
    if len(genes) > 0:
        window_total_values.index = genes
        
    return window_total_values


def nfr(values, center_x=1000, scale=False, order=35):
    """
    Calculate Nucleosome Free Region score

    Params
    ------
        values
            numpy.ndarray
        center_x
            int
        smooth
            bool
        scale
            bool
        order
            int

    Returns
    -------
        nfr_score
            float
        pos_score
            float
    """
    
    # Scale values between 0 and 1
    if scale:
        values = (values - np.min(values)) / (np.max(values) - np.min(values))

    # Calculate NFR
    nf = np.sum(values[center_x - 50:center_x + 50])
    n1 = np.sum(values[center_x - 200:center_x - 50])
    n2 = np.sum(values[center_x + 50:center_x + 200])
    nfr_score = np.log2(nf) - np.log2((n1 + n2) / 2)
    
    # Find position score
    control_mean = (np.sum(values[:100]) + np.sum(values[-100:])) / 200
    pos_score = np.log2(values[center_x] / control_mean)		
    
    return nfr_score, pos_score


def nfr_windows(window_values, center_x=1000, scale=False, order=35):
    """
    """
    nfr_scores = np.zeros(window_values.shape[0])
    for i in range(window_values.shape[0]):
        nfr_score, pos_score = nfr(window_values.values[i,:], center_x=center_x, scale=scale, order=order)
        nfr_scores[i] = nfr_score

    return nfr_scores


def get_genes(bed_fn):
    """
    """

    # Initialize genes
    genes = []

    # Iterate over BED records
    for i, line in enumerate(open(bed_fn, "r")):
        field = line.strip().split("\t") # find fields
        chrom = field[0]

        genes.append(field[-1])

    return np.array(genes)


def bin_coverage(frags, bin_size=100000, min_length=None, max_length=None):
    """
    """

    bins = frags.bin_nhits(bin_size=bin_size, min_length=min_length, max_length=max_length)

    return bins


def predict_sex(bin_coverage, fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = True):
    """
    """

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

    if chrXMedian < chrXMedianForMale:
        if useChrY and (chrYCov < fracReadsInChrYForMale):
            gender = "female"
        else:
            gender = "male"
    else:
        gender = "female"

    return {"gender":gender, "chrYCovRatio":chrYCov, "chrXMedian":chrXMedian}


def bin_array(scores, bin_size=100000, method="mean"):
    """
    """
    
    # Determine function
    if method == "mean":
        func = np.mean
    elif method == "median":
        func = np.median
    elif method == "sum":
        func = np.sum
    else:
        raise AttributeError("Unspecified method (mean, median, or sum)")

    # Determine bins
    if len(scores) == 0:
        return pd.Series([])

    start_bin = int(scores.index.values[0] / bin_size) * bin_size
    end_bin = int(scores.index.values[-1] / bin_size) * bin_size

    # Bin
    bin_results = np.zeros(int(((end_bin - start_bin) / bin_size) + 1))
    for i, b in enumerate(range(start_bin, end_bin+1, bin_size)):
        bin_results[i] = func(scores.loc[b:b+bin_size].values)

    # Create pandas.Series
    bin_results = pd.Series(bin_results,
                            index=np.arange(start_bin, end_bin+1, bin_size))

    return bin_results


def interval_annotate(intervals, scores_dict, method = "mean"):
    """
    """

    # Determine function
    if method == "mean":
        func = np.mean
    elif method == "median":
        func = np.median
    elif method == "sum":
        func = np.sum
    else:
        raise AttributeError("Unspecified method (mean, median, or sum)")

    # Initialize values
    values = np.zeros(len(intervals))

    # Iterate over intervals
    for i, interval in enumerate(intervals):
        start = scores_dict[interval.label].index[0]
        end = scores_dict[interval.label].index[-1]
        if interval.end < start or interval.start > end:
            continue
        values[i] = func(scores_dict[interval.label].values[interval.start - start:interval.end - start])

    return values


def write_bw(values_dict, out_fn, chrom=None):
    """
    """
    import pyBigWig

    # Determine chromosomes to query
    if chrom is not None:
        chroms = [chrom]
    else:
        chroms = self.chroms
    
    # Open BigWig file
    out_bw = pyBigWig.open(out_fn, "w")
    
    # Write header
    header = [(key, values_dict[key].index.values[-1]) for key in values_dict]
    out_bw.addHeader(header)
    
    # Write values to BigWig
    for c in chroms:
        if self.verbose: print(c)
        out_bw.addEntries(c, values_dict[c].index.values[0], values=values_dict[c].values.astype(float), span=1, step=1, validate=False)
    
    # Close file
    out_bw.close()


def write_h5(values_dict, out_fn, chrom=None):
    """
    """
    import h5py

    # Determine chromosomes to query
    if chrom is not None:
        chroms = [chrom]
    else:
        chroms = self.chroms
    
    # Open BigWig file
    out_h5 = h5py.File(out_fn, "w")

    # Write values to h5
    for c in chroms:
        if self.verbose: print(c)
        # Check values are not None
        if values_dict[c] is None:
            continue
        chrom_group = out_h5.create_group(c)
        chrom_group["values"] = values_dict[c].values
        chrom_group["position"] = values_dict[c].index.values
    
    # Close file
    out_h5.close()
    
    
def write_bedgraph(out_fn, interval_dict=None, chrom=None):
    """
    Write fragments as bedgraph format
    """
    import io

    # Determine chromosomes to query
    if chrom is not None:
        chroms = [chrom]
    else:
        chroms = self.chroms

    # Determine intervals to write
    if interval_dict is None:
        interval_dict = self.frags

    # Determine if out is a file or string
    if isinstance(out_fn, str):
        out_bed = open(out_fn, "w")
    elif isinstance(out_fn, io.TextIOWrapper):
        out_bed = out_fn
    
    # Iterate over fragments and write to file
    for c in chroms:
        for interval in interval_dict[c]:
            out_bed.write("%s\t%d\t%d\t%f\n)" % (c, interval.start, interval.end, interval.value))

    # If was string, close file
    if isinstance(out_fn, str):
        out_bed.close()