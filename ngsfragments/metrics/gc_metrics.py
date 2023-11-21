import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray

# Local imports
from ..fragments import Fragments


def gc_bias(intervals: Fragments | LabeledIntervalArray,
            genome_version: str = "hg19"):
    """
    Calculate GC bias

    Parameters
    ----------
        intervals : Fragments | LabeledIntervalArray
            Intervals to calculate GC bias for
        genome_version : str
            Genome version to use
        
    Returns
    -------
        gc_content : np.ndarray
            GC content of each interval
    """

    # Get intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags
    
    # Get genome
    if genome_version == "hg19":
        try:
            import hg19genome as genome
        except ImportError:
            print("hg19genome not installed. Please install hg19genome: 'pip install hg19genome'")
    elif genome_version == "hg38":
        try:
            import hg38genome as genome
        except ImportError:
            print("hg38genome not installed. Please install hg38genome: 'pip install hg38genome'")
    else:
        raise NotImplementedError("Only hg19 genome is currently supported")

    # Get bin bias
    gc_content = genome.kmers.kmer_reader.gc_percent(intervals)
    lengths = intervals.ends - intervals.starts

    return gc_content, lengths


def fragment_gc_bias(gc_content: np.ndarray,
                     lengths: np.ndarray):
    """
    Calculate GC bias per fragment

    Parameters
    ----------
        gc_content : np.ndarray
            GC content of each interval
        lengths : np.ndarray
            Length of each interval

    Returns
    -------
        record : pd.DataFrame
            GC bias per fragment
    """

    bins = np.arange(0,1,0.1)
    record = pd.DataFrame(np.zeros((len(np.arange(0, np.max(lengths), 100)),
                                    len(bins))),
                          index = np.arange(0, np.max(lengths), 100),
                          columns = bins)
    for i, n in enumerate(range(0, np.max(lengths), 100)):
        bin_nums = np.digitize(gc_content[np.logical_and(lengths > n, lengths < n + 100)], bins) - 1
        hist_vals = np.bincount(bin_nums)
        hist_vals = hist_vals / np.max(hist_vals)

        record.iloc[i,:len(hist_vals)] = hist_vals
    
    return record


def simulate_gc_bias(intervals,
                     n_simulations: int = 10):
    """
    Simulate GC bias

    Parameters
    ----------
        intervals : Fragments | LabeledIntervalArray
            Intervals to calculate GC bias for
        n_simulations : int
            Number of simulations to run
        
    Returns
    -------
        gc_content : np.ndarray
            GC content of each interval
    """

    # Get intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Initialize record
    bins = np.arange(0,1,0.1)
    lengths = intervals.ends - intervals.starts
    record = pd.DataFrame(np.zeros((len(np.arange(0, np.max(lengths), 100)),
                                    len(bins))),
                        index = np.arange(0, np.max(lengths), 100),
                        columns = bins)

    # Simulate GC bias
    for n_sim in range(n_simulations):
        print(n_sim)
        sim = intervals.simulate()
        # Calculate GC bias
        gc_content, lengths = gc_bias(sim)

        for i, n in enumerate(range(0, np.max(lengths), 100)):
            bin_nums = np.digitize(gc_content[np.logical_and(lengths > n, lengths < n + 100)], bins) - 1
            hist_vals = np.bincount(bin_nums)
            hist_vals = hist_vals / np.max(hist_vals)

            record.values[i,:len(hist_vals)] += hist_vals
        
    record.values = record.values / n_simulations
    
    return record

