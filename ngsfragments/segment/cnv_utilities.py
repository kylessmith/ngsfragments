import numpy as np
import pandas as pd
from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
from projectframe import ProjectFrame
import hmmCNV
from typing import Any, Dict, List


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
                    column: str = "ratios"):
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
        
        segment_mad1 = 1.4826 * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
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
                segment_mad1 = 1.4826 * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
                segment_bounds1 = (segment_median1 - segment_mad1,
                                segment_median1 + segment_mad1)

            else:
                merged.add(previous_start, previous_end, chrom)
                merged_values.append(segment_median1)
                previous_start = chrom_segs.index[i].start
                previous_end = chrom_segs.index[i].end
                segment_median1 = float(chrom_segs.df.loc[:,"median"].values[i])
                segment_mad1 = 1.4826 * calculate_MAD(ratios.intersect(previous_start, previous_end, chrom).values)
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
                   column: str = "ratios"):
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
        merged = iterative_merge(merged, bins, column=column)
        current_length = merged.shape[0]

    return merged


def hmm_classify(segments: IntervalFrame,
                 hmm_states: Any,
                 includeHOMD: bool = False):
    """
    Classify segments using HMM states

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        hmm_states : Any
            HMM states
        includeHOMD : bool
            Include HOMD state

    Returns
    -------
        copy_number : IntervalFrame
            Copy number states
    """

    # Calculate probability for each segment value
    names = np.array(["HOMD","HETD","NEUT","GAIN","AMP","HLAMP"] + ["HLAMP"+str(i) for i in range(2,25)])
    values = segments.df.loc[:,"median"].values
    new_states = np.zeros(len(values), dtype=int)
    probs = np.zeros(len(values))
    # Iterate over states
    for state in range(len(hmm_states["mus"])):
        p = hmmCNV.em_utilities.tdistPDF(values[np.newaxis].T, [hmm_states["mus"][state]], [hmm_states["lambdas"][state]], hmm_states["nu"])
        new_states[p > probs] = state
        probs[p > probs] = p[p > probs]

    # Assign values
    copy_number = hmm_states["states"].values[new_states]
    segments.df.loc[:,"copy_number"] = copy_number
    segments.df.loc[:,"event"] = names[copy_number]
    segments.df.loc[:,"subclone_status"] = hmm_states["subclone"].values[new_states]

    # Create inputs for correctIntegerCN
    cn = hmm_states["cn"]
    purity = hmm_states["n"]
    ploidy = hmm_states["phi"]
    cellPrev = hmm_states["sp"]
    maxCN = len(hmm_states["mus"])

    # Correcte Integer CN
    correctedResults = hmmCNV.hmm_utilities.correctIntegerCN(cn = cn.copy(),
                                                        segs = segments.copy(), 
                                                        purity = 1 - purity,
                                                        ploidy = ploidy,
                                                        cellPrev = 1 - cellPrev, 
                                                        maxCNtoCorrect_autosomes = maxCN,
                                                        maxCNtoCorrect_X = maxCN,
                                                        minPurityToCorrect = 0.03, 
                                                        gender = None,
                                                        correctHOMD = includeHOMD)
    
    return correctedResults["segs"]


def validate_calls(segments: IntervalFrame,
                   bins: IntervalFrame,
                   n_mads: int = 1,
                   column: str = "ratios"):
    """
    Validate NEUT calls are in distribution

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        n_mads : int
            Number of MADs to use

    Returns
    -------
        segments : IntervalFrame
            IntervalFrame of segments
    """

    # Determine which bins are neutral
    is_neutral = segments.df.loc[:,"Corrected_Call"].values == "NEUT"

    # Keep correcting until no changes
    past_n_neutral = np.sum(is_neutral)
    n_neutral = 0
    while n_neutral != past_n_neutral:
        # Determine which bins are neutral
        is_neutral = segments.df.loc[:,"Corrected_Call"].values == "NEUT"
        medians = segments.df.loc[is_neutral,"median"].values

        # Update past_n_neutral
        past_n_neutral = np.sum(segments.loc[:,"Corrected_Call"].values == "NEUT")

        # Define a neutral distribution
        neutral_distribution = bins.overlap(segments.iloc[is_neutral,:]).df.loc[:,column].values
        neutral_median = np.median(neutral_distribution)
        neutral_mad = calculate_MAD(neutral_distribution)
        
        # Determine pvalues
        new_state = np.ones(is_neutral.sum())
        for i in range(is_neutral.sum()):
            valid = True
            if medians[i] > (neutral_median + (neutral_mad * n_mads)):
                valid = False
            elif medians[i] < (neutral_median - (neutral_mad * n_mads)):
                valid = False
            # Determine new state
            if not valid:
                if medians[i] < 0:
                    new_state[i] = 0
                else:
                    new_state[i] = 2

        # Correct neutral callings
        states = np.ones(segments.shape[0], dtype=int)
        states[is_neutral] = new_state
        
        # Assign new states
        segments.df.loc[states==0,"Corrected_Call"] = "HETD"
        segments.df.loc[states==0,"Corrected_Copy_Number"] = 0
        segments.df.loc[states==2,"Corrected_Call"] = "GAIN"
        segments.df.loc[states==2,"Corrected_Copy_Number"] = 2
        
        # Correct number of neutrals
        n_neutral = np.sum(segments.df.loc[:,"Corrected_Call"].values == "NEUT")

    return segments


def correct_neutral(segments: IntervalFrame,
                    bins: IntervalFrame,
                    n_mads: int = 1,
                    column: str = "ratios"):
    """
    Correct neutral calls

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        n_mads : int
            Number of MADs to use
    
    Returns
    -------
        segments : IntervalFrame
            IntervalFrame of segments
    """

    # Determine which bins are neutral
    is_neutral = segments.df.loc[:,"Corrected_Call"].values == "NEUT"

    # Keep correcting until no changes
    past_n_neutral = is_neutral.sum()
    n_neutral = 0
    while n_neutral != past_n_neutral:
        # Determine which bins are neutral
        is_neutral = segments.df.loc[:,"Corrected_Call"].values == "NEUT"

        # Update past_n_neutral
        past_n_neutral = is_neutral.sum()

        # Define a neutral distribution
        neutral_distribution = bins.overlap(segments.iloc[is_neutral,:]).df.loc[:,column].values
        neutral_median = np.median(neutral_distribution)
        neutral_mad = calculate_MAD(neutral_distribution)

        # Define non-neutral segments
        medians = segments.df.loc[~is_neutral,"median"].values
        
        # Determine pvalues
        valid = np.ones(len(medians), dtype=bool)
        for i in range(np.sum(~is_neutral)):
            if medians[i] > (neutral_median + (neutral_mad * n_mads)):
                valid[i] = False
            elif medians[i] < (neutral_median - (neutral_mad * n_mads)):
                valid[i] = False

        # Correct neutral callings
        new_neutral = np.zeros(segments.shape[0], dtype=bool)
        new_neutral[~is_neutral] = valid
        segments.df.loc[new_neutral,"Corrected_Call"] = "NEUT"
        segments.df.loc[new_neutral,"Corrected_Copy_Number"] = 2
        
        # Correct number of neutrals
        n_neutral = np.sum(segments.df.loc[:,"Corrected_Call"].values == "NEUT")

    return segments


def correct_direction(segments: IntervalFrame) -> IntervalFrame:
    """
    """

    # Determine direction
    for i in np.arange(segments.shape[0]):
        index = segments.df.index.values[i]
        if segments.df.loc[index,"Corrected_Call"] == "NEUT":
            continue
        elif segments.df.loc[index,"Corrected_Call"] == "HETD" and segments.df.loc[index,"median"] > 0:
            segments.df.loc[index,"Corrected_Call"] = "GAIN"
            segments.df.loc[index,"Corrected_Copy_Number"] == 3
        elif segments.df.loc[index,"Corrected_Call"] == "GAIN" and segments.df.loc[index,"median"] < 0:
            segments.df.loc[index,"Corrected_Call"] = "HETD"
            segments.df.loc[index,"Corrected_Copy_Number"] == 1

    return segments


def median_variance(segments: IntervalFrame,
                    bins: IntervalFrame,
                    key: str):
    """
    Calculate median segment variance

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        key : str
            Key to use for variance calculation
    
    Returns
    -------
        median_variance : float
            Median segment variance
    """

    # Calculate variance
    segments.annotate(bins, key, "var")
    median_variance = np.median(segments.df.loc[:,"var"].values)

    return median_variance


def train_hmm(bins: IntervalFrame,
              normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
              ploidy: List[float] = [2],
              gender: str | None = None,
              estimatePloidy: bool = False,
              minSegmentBins: int = 25,
              maxCN: int = 7,
              scStates: List[int] = [1, 3],
              verbose: int = False,
              **kwargs: Any):
    """
    Train HMM

    Parameters
    ----------
        bins : IntervalFrame
            IntervalFrame of bins
        normal : List[float]
            List of normal ratios
        ploidy : List[float]
            List of ploidies
        gender : str | None
            Gender of sample
        estimatePloidy : bool
            Whether to estimate ploidy
        minSegmentBins : int
            Minimum number of bins in a segment
        maxCN : int
            Maximum copy number
        verbose : int
            Verbosity level
        **kwargs : Any
            Additional arguments to pass to hmmCNV

    Returns
    -------
        hmm_loglik : Dict[str, float]
            DataFrame of log likelihoods
    """

    # Run HMMcopy
    hmm_loglik, hmm_results = hmmCNV.hmmCNV(bins,
                                            normal=normal,
                                            ploidy=ploidy,
                                            gender=gender,
                                            estimatePloidy=estimatePloidy,
                                            minSegmentBins=minSegmentBins,
                                            maxCN=maxCN,
                                            scStates=scStates,
                                            verbose=verbose,
                                            **kwargs)
    
    # Find optimal segments
    optimal_result = np.argmax(hmm_loglik.loc[:,"loglik"].values)
    #hmm_segments = hmm_results[optimal_result]["results"]["segs"]
    hmm_states = {}
    hmm_states["mus"] = hmm_results[optimal_result]["results"]["mus"][-1]
    hmm_states["lambdas"] = hmm_results[optimal_result]["results"]["lambdas"][-1]
    hmm_states["nu"] = hmm_results[optimal_result]["results"]["param"]["nu"]
    hmm_states["states"] = hmm_results[optimal_result]["results"]["param"]["jointCNstates"]["Sample_1"]
    hmm_states["subclone"] = hmm_results[optimal_result]["results"]["param"]["jointSCstatus"]["Sample_1"]
    hmm_states["n"] = hmm_results[optimal_result]["results"]["n"][0,-1]
    hmm_states["phi"] = hmm_results[optimal_result]["results"]["phi"][0,-1]
    hmm_states["sp"] = hmm_results[optimal_result]["results"]["sp"][0,-1]
    hmm_states["Frac_genome_subclonal"] = float(hmm_loglik.iloc[optimal_result].loc["Frac_genome_subclonal"])
    hmm_states["cn"] = hmm_results[optimal_result]["cna"]

    return hmm_states


def validate_distributions(segments: IntervalFrame,
                           bins: IntervalFrame,
                           n_mads: float = 1.4826,
                           column: str = "ratios"):
    """
    Validate segment distributions

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        bins : IntervalFrame
            IntervalFrame of bins
        n_mads : int
            Number of MADs to use

    Returns
    -------
        None
    """

    # Valicate segment classifications
    segments = correct_neutral(segments, bins, n_mads=n_mads, column = column)
    segments = validate_calls(segments, bins, n_mads=n_mads, column = column)
    segments = correct_direction(segments)

    return segments


def write_seg_file(pf: ProjectFrame,
                   out_fn: str,
                   sample: str):
    """
    Write segments to SEG file

    Parameters
    ----------
        segments : IntervalFrame
            IntervalFrame of segments
        out_fn : str
            Output file name
        sample : str
            Sample name

    Returns
    -------
        None
    """
    
    # Open output file
    out = open(out_fn, "w")
    fmt = "{sample}\t{chrom}\t{start}\t{end}\t{n_bins}\t{log2_ratio_median}\n"
    out.write(fmt.replace('}','').replace('{',''))

    # Get segments
    segs = pf.obs_intervals[sample]["cnv_segments"]
    chroms = segs.index.labels
    starts = segs.starts
    ends = segs.ends

    # Iterate over segment overlaps
    for i, overlap in enumerate(segs.iter_intersect(pf.intervals["cnv_bins"])):
        n_bins = overlap.shape[0]
        log2_ratio_median = segs.df.iloc[i,0]
        chrom = chroms[i]
        start = starts[i]
        end = ends[i]
        out.write(fmt.format(**locals()))
            
    # Close seg file
    out.close()