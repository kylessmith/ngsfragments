from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
from math import e
import numpy as np
import pandas as pd
import glob
import os

# Local imports
from ..fragments import Fragments
from ..io.read_sam import from_sam
from ..segment.cnv import call_cnvs
from ..correct.correct_intervals import calculate_interval_bias


def score_interval(bins: IntervalFrame,
                    start: int,
                    end: int,
                    chrom: str):
    """
    Score an interval by weighted sum of the number of hits in the interval

    Parameters
    ----------
        bins : IntervalFrame
            IntervalFrame containing the number of hits per bin
        start : int
            Start of interval
        end : int
            End of interval
        chrom : str
            Chromosome of interval

    Returns
    -------
        score : float
            Score of interval
    """

    length = end - start
    nhits = bins.intersect(start, end, chrom).df.loc[:,"counts"].values
    #weight = (1 * (1/length)) + e**-1
    #weight = 1 + e**-1
    weight = (1 + (e**-1)) * (1/length)
    score = np.sum(nhits * weight)

    return score


def score_upstream(bins: IntervalFrame,
                    start: int,
                    end: int,
                    chrom: str,
                    upstream: int = 100000):
    """
    Score upstream of an interval by weighted sum of the number of hits in the interval

    Parameters
    ----------
        bins : IntervalFrame
            IntervalFrame containing the number of hits per bin
        start : int
            Start of interval
        end : int
            End of interval
        chrom : str
            Chromosome of interval
        upstream : int
            Number of basepairs upstream to score

    Returns
    -------
        score : float
            Score of interval
    """

    length = end - start
    end = start
    start = max(start - upstream, 0)
    nhits = bins.intersect(start, end, chrom)
    
    if nhits.shape[0] > 0:
        distance = end - nhits.ends
        #weight = ((e**(-abs(distance) / 5000)) * (1/length)) + (e**-1)
        #weight = (e**(-abs(distance) / 5000)) + (e**-1)
        weight = ((e**(-abs(distance) / 5000)) + (e**-1)) * (1/length)
        score = np.sum(nhits.df.loc[:,"counts"].values * weight)
    else:
        score = 0

    return score


def score_downstream(bins: IntervalFrame,
                    start: int,
                    end: int,
                    chrom: str,
                    downstream: int = 100000):
    """
    Score downstream of an interval by weighted sum of the number of hits in the interval

    Parameters
    ----------
        bins : IntervalFrame
            IntervalFrame containing the number of hits per bin
        start : int
            Start of interval
        end : int
            End of interval
        chrom : str
            Chromosome of interval
        downstream : int
            Number of basepairs downstream to score

    Returns
    -------
        score : float
            Score of interval
    """

    length = end - start
    start = end
    end = end + downstream
    nhits = bins.intersect(start, end, chrom)
    
    if nhits.shape[0] > 0:
        distance = end - nhits.starts
        #weight = ((e**(-abs(distance) / 5000)) * (1/length)) + (e**-1)
        #weight = (e**(-abs(distance) / 5000)) + (e**-1)
        weight = ((e**(-abs(distance) / 5000)) + (e**-1)) * (1/length)
        score = np.sum(nhits.df.loc[:,"counts"].values * weight)
    else:
        score = 0

    return score


def gene_activity(intervals: Fragments | LabeledIntervalArray,
                  genome_version: str = "hg19",
                  feature: str = "gene",
                  min_length: int = 120,
                  max_length: int = 220,
                  verbose: bool = False):
    """
    Determine the activity of genes

    Parameters
    ----------
        frags : Fragments
            Fragments object containing the fragments
        genome : str
            Genome to use
        feature : str
            Feature to use
        min_length : int
            Minimum length of fragments to use
        max_length : int
            Maximum length of fragments to use
        verbose : bool
            Print progress
    
    Returns
    -------
        scores : pd.Series
            Series containing the scores of the genes
    """

    # Get sample name
    if isinstance(intervals, Fragments):
        path = os.path.normpath(intervals.sam_file)
        file_name = path.split(os.sep)[-1]
        sample_name = file_name.split(".bam")[0]
    else:
        sample_name = "sample"

    # Get windows
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)

    # Get genes
    if feature == "gene":
        genes = genome.get_intervals("gene", upstream = 5000, filter_column = "gene_type", filter_selection="protein_coding")
    elif feature == "promoter":
        genes = genome.get_intervals("tss", upstream = 5000, downstream=1000, filter_column = "gene_type", filter_selection="protein_coding")
    else:
        raise NotImplementedError("Only gene and promoter are supported currently")

    # Convert intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Determine chromosomes
    chromosomes = intervals.unique_labels

    # Iterate over chromosomes
    scores = pd.Series(np.zeros(genes.shape[0]),
                       index = genes.loc[:,"gene_name"].values)
    scores.name = sample_name

    j = 0
    for chrom in chromosomes:
        if verbose: print(chrom)
        # Get genes on chromosome
        chrom_genes = genes.loc[chrom,:]
        # Determine nhits per bin
        bins = intervals.get(chrom).bin_nhits(500, min_length, max_length)

        # Create interval frame
        bins_iframe = IntervalFrame(bins[0])
        bins_iframe.df.loc[:,"counts"] = bins[1]

        for i in range(chrom_genes.shape[0]):
            score = score_interval(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            score += score_upstream(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            score += score_downstream(bins_iframe, chrom_genes.index[i].start, chrom_genes.index[i].end, chrom)
            scores[j] = score
            j += 1

    return scores


def correct_gene_activity(frags: Fragments,
                          scores: pd.Series,
                          correct_cnv: bool = False,
                          genome_version: str = "hg19",
                          feature: str = "gene",
                          verbose: bool = False):
    """
    Correct the activity of genes

    Parameters
    ----------
        scores : pd.Series
            Series containing the scores of the genes
        genome : str
            Genome to use
        feature : str
            Feature to use
        verbose : bool
            Print progress
    
    Returns
    -------
        scores : pd.Series
            Series containing the scores of the genes
    """
    # Get sample name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]
    sample_name = file_name.split(".bam")[0]

    # Get windows
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)

    # Get genes
    if feature == "gene":
        genes = genome.get_intervals("gene", upstream = 5000, filter_column = "gene_type", filter_selection="protein_coding")
    elif feature == "promoter":
        genes = genome.get_intervals("tss", upstream = 5000, downstream=1000, filter_column = "gene_type", filter_selection="protein_coding")
    else:
        raise NotImplementedError("Only gene and promoter are supported currently")

    # Remove unnecessary columns
    genes.drop_columns(["Strand","gene_type","level","Source"])

    # Filter scores
    lookup = pd.Series(np.arange(genes.shape[0]),
                       index = genes.loc[:,"gene_name"].values)
    gene_names = lookup.index.intersection(scores.index.values)
    genes = genes.iloc[lookup.loc[gene_names].values,:]
    scores = scores.loc[gene_names]

    # Calculate CNVs
    if correct_cnv:
        cnv_bins, cnv_segments = call_cnvs(frags,
                                            None,
                                            bin_size = 100000,
                                            genome = genome_version,
                                            method = "bcp_online_both",
                                            outlier_smooth = True,
                                            gauss_smooth = False,
                                            verbose = verbose)
    else:
        cnv_bins = None

    # Correct scores
    genes.loc[:,"score"] = scores.values
    genes = calculate_interval_bias(genes,
                                    column = "score",
                                    cnv_bins = cnv_bins,
                                    genome_version = genome_version)

    # Return genes
    corrected_scores = pd.Series(genes.loc[:,"corrected_values"].values,
                                 index = genes.loc[:,"gene_name"].values)
    corrected_scores.name = sample_name

    return corrected_scores


def multi_gene_activity(directory: str,
                        genome_version: str = "hg19",
                        feature: str = "gene",
                        min_length: int = 121,
                        max_length: int = 375,
                        nthreads: int = 1,
                        length_norm: bool = True,
                        verbose: bool = False):
    """
    Determine the activity of genes for multiple samples

    Parameters
    ----------
        directory : str
            Directory containing the bam files
        genome : str
            Genome to use
        feature : str
            Feature to use
        min_length : int
            Minimum length of fragments to use
        max_length : int
            Maximum length of fragments to use
        nthreads : int
            Number of threads to use
        length_norm : bool
            Whether to normalize by feature length
        verbose : bool
            Print progress
    
    Returns
    -------
        scores : pd.DataFrame
            DataFrame containing the scores of the genes
    """

    # Get files
    files = glob.glob(directory+"*.bam")

    # Iterate over files
    if verbose: print(files[0])
    frags = from_sam(files[0], nthreads=nthreads)
    scores = gene_activity(frags, feature=feature, genome_version=genome_version, min_length=min_length, max_length=max_length)
    scores = correct_gene_activity(frags,
                                   scores = scores,
                                   correct_cnv = False,
                                   genome_version = genome_version,
                                   feature = feature,
                                   verbose = verbose)
    scores = scores.to_frame().T
    scores.index = [files[0]]
    for bam in files[1:]:
        if verbose: print(bam)
        frags = from_sam(bam, nthreads=nthreads)
        bam_scores = gene_activity(frags, feature=feature, genome_version=genome_version, min_length=min_length, max_length=max_length)
        bam_scores = correct_gene_activity(frags,
                                    scores = bam_scores,
                                    correct_cnv = False,
                                    genome_version = genome_version,
                                    feature = feature,
                                    verbose = verbose)
        scores.loc[bam,:] = bam_scores

    return scores