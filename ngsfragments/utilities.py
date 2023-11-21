import numpy as np
import pandas as pd
import gzip

#from numba import jit

# Set data directory
#data_dir = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
#data_dir = os.path.join(data_dir, "data")


def window_scores(scores, feature="tss", upstream=1000, downstream=1000, stranded=False, verbose=False):
    """
    Calculate scores for each postion within a window
    """

    # Determine bed file
    #bed_fn = get_data_file(bed_fn)

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
    with open_func(bed_fn, "r") as f:
        for i, line in enumerate(f):
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
    elif method == "var":
        func = np.var
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


