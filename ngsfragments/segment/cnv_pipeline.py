import ngsfragments as ngs
import os
import numpy as np
import pandas as pd
from typing import List
import glob
from projectframe import ProjectFrame
np.seterr(all="ignore")

# Local imports
from ..fragments import Fragments
from .cnv import call_cnvs, sex_prediction
from .cnv_utilities import train_hmm, merge_segments, hmm_classify, validate_distributions
from ..io import from_sam


def log_fragments(pf: ProjectFrame,
                  frags: Fragments) -> ProjectFrame:
    """
    Log Fragments object to ProjectFrame

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragments object

    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Determine name
    file_name = frags.sam_file
    sample_name = os.path.split(file_name)[-1].split(".bam")[0]

    # Check is file_name exists
    #if sample_name not in pf.obs:
    pf.add_obs(sample_name)

    # Record
    if "length_dist" not in pf.values.keys or pd.isnull(pf.values["length_dist"].loc[sample_name,:].values).all():
        length_dist = frags.length_dist()
        length_dist.name = sample_name
        pf.add_values("length_dist", length_dist)

    if "n_fragments" not in pf.anno.columns or pd.isnull(pf.anno.loc[sample_name, "n_fragments"]):
        pf.add_anno("n_fragments", sample_name, len(frags.frags))
    
    if "chrom_lengths" not in pf.uns.keys:
        pf.uns["chrom_lengths"] = frags.genome["chrom_sizes"]
        pf.uns["chroms"] = frags.chroms

    if "chroms" not in pf.uns.keys:
        pf.uns["chroms"] = frags.chroms

    return pf


def call_cnv_pipeline(pf: ProjectFrame,
                      data: List[str] | str | Fragments,
                      cnv_binsize: int = 100000,
                      hmm_binsize: int = 1000000,
                      genome_version: str = "hg19",
                      nthreads: int = 1,
                      method: str = "online_both",
                      outlier_smooth: bool =True,
                      gauss_smooth: bool = False,
                      bcp_cutoff: float = 0.3,
                      hazard: int = 100,
                      shuffles: int = 5000,
                      p: float = 0.00005,
                      normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                      ploidy: List[int] = [2],
                      estimatePloidy: bool = False,
                      minSegmentBins: int = 25,
                      maxCN: int = 5,
                      use_normal: bool = False,
                      keep_sex_chroms: bool = False,
                      scStates: List[int] = [1, 3],
                      verbose: bool = False,
                      **kwargs) -> ProjectFrame:
    """
    Call CNVs from fragments

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        data : str | Fragments
            Fragments object or SAM file
        cnv_binsize : int
            Bin size for CNV calling
        hmm_binsize : int
            Bin size for HMM
        genome_version : str
            Genome version
        nthreads : int
            Number of threads
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        hazard : int
            Hazard
        shuffles : int
            Number of shuffles
        p : float
            P-value
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
        **kwargs : dict
            Additional arguments
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Log genome version
    pf.params["ref_genome"] = genome_version

    # Read bam file
    if isinstance(data, str):
        data = from_sam(data, nthreads=nthreads, verbose=verbose)
    
    # Process a single sample
    if isinstance(data, Fragments):
        pf = call_cnvs_single(pf,
                            data,
                            cnv_binsize,
                            hmm_binsize,
                            method,
                            outlier_smooth,
                            gauss_smooth,
                            bcp_cutoff,
                            normal,
                            ploidy,
                            estimatePloidy,
                            minSegmentBins,
                            maxCN,
                            hazard,
                            shuffles,
                            p,
                            use_normal,
                            keep_sex_chroms,
                            scStates,
                            verbose)

    if isinstance(data, list):
        pf = call_cnvs_multiple(pf,
                            data,
                            cnv_binsize,
                            hmm_binsize,
                            method,
                            outlier_smooth,
                            gauss_smooth,
                            bcp_cutoff,
                            normal,
                            ploidy,
                            estimatePloidy,
                            minSegmentBins,
                            maxCN,
                            hazard,
                            shuffles,
                            p,
                            nthreads,
                            use_normal,
                            verbose)

    return pf


def predict_cnvs(pf: ProjectFrame,
                    frags: ngs.Fragments,
                    bin_size: int = 100000,
                    method: str = "online_both",
                    outlier_smooth: bool =True,
                    gauss_smooth: bool = False,
                    bcp_cutoff: float = 0.3,
                    hazard: int = 100,
                    shuffles: int = 5000,
                    p: float = 0.00005,
                    use_normal: bool = False,
                    keep_sex_chroms: bool = False,
                    verbose: bool = False,
                    **kwargs):
    """
    Call Copy Number Variants (CNVs)

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragment
        bin_size : int
            Bin size
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        hazard : int
            Hazard
        shuffles : int
            Number of shuffles
        p : float
            P-value
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
        **kwargs : dict
            Additional arguments

    Returns
    -------
        pf : ProjectFrame
            ProjectFrame

    """

    # Find file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]
    sample_name = file_name.split(".bam")[0]

    # Log fragments
    pf = log_fragments(pf, frags)

    # Log bin size
    pf.params["cnv_binsize"] = bin_size

    # Call Copy Number Variations
    cnv_bins, cnv_segs = call_cnvs(data = frags,
                                    bin_size = bin_size,
                                    genome_version = str(pf.params["ref_genome"]),
                                    method = method,
                                    outlier_smooth = outlier_smooth,
                                    gauss_smooth = gauss_smooth,
                                    verbose = verbose,
                                    cutoff = bcp_cutoff,
                                    hazard = hazard,
                                    shuffles = shuffles,
                                    p = p,
                                    use_normal = use_normal,
                                    keep_sex_chroms = keep_sex_chroms)
    
    # Append bins
    cnv_bins = cnv_bins.loc[:,["ratios"]]
    cnv_bins.df.columns = [sample_name]
    pf.add_obs_intervals(sample_name, "cnv_segments", cnv_segs)
    pf.add_intervals("cnv_bins", cnv_bins)

    return pf


def predict_purity(pf: ProjectFrame,
                   frags: ngs.Fragments,
                   file_name: str = None,
                   bin_size: int = 1000000,
                   method: str = "online_both",
                   outlier_smooth: bool = True,
                   gauss_smooth: bool = False,
                   bcp_cutoff: float = 0.3,
                   normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                   ploidy: List[int] = [1, 2, 3],
                   estimatePloidy: bool = False,
                   minSegmentBins: int = 25,
                   maxCN: int = 5,
                   use_normal: bool = False,
                   scStates: List[int] = [1, 3],
                   verbose: bool = False,
                   **kwargs):
    """
    Predict tumor purity

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragment
        file_name : str
            File name
        bin_size : int
            Bin size
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
        **kwargs : dict
            Additional arguments
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Find file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]
    sample_name = file_name.split(".bam")[0]

    # Check if file was previously annotated
    pf = log_fragments(pf, frags)

    # Determine if remove chr19
    #if pf.params["ref_genome"] == "hg19" or pf.params["ref_genome"] == "hg38":
    #    remove_chroms.append("chr19")

    # Calculate bins
    hmm_binsize = bin_size
    bins, segments = call_cnvs(data = frags,
                                bin_size = bin_size,
                                genome_version = str(pf.params["ref_genome"]),
                                method = method,
                                outlier_smooth = outlier_smooth,
                                gauss_smooth = gauss_smooth,
                                verbose = verbose,
                                cutoff = bcp_cutoff,
                                use_normal = use_normal)

    # Define gender
    try:
        gender = pf.obs_values[sample_name]["gender"]
    except KeyError:
        gender = None
    
    # Classify calls be HMM
    hmm_states = train_hmm(bins,
                            normal = normal,
                            ploidy = ploidy,
                            gender = None,
                            estimatePloidy = estimatePloidy,
                            minSegmentBins = minSegmentBins,
                            maxCN = maxCN,
                            scStates = scStates,
                            verbose = verbose,
                            **kwargs)

    # Assign attributes
    pf.add_anno("purity", sample_name, 1 - hmm_states["n"])
    pf.add_anno("ploidy", sample_name, hmm_states["phi"])
    pf.add_anno("clonal", sample_name, hmm_states["Frac_genome_subclonal"])

    return pf


def calculate_hmm_states(pf: ProjectFrame,
                        frags: ngs.Fragments,
                        bin_size: int = 1000000,
                        method: str = "online_both",
                        outlier_smooth: bool = True,
                        gauss_smooth: bool = False,
                        bcp_cutoff: float = 0.3,
                        normal: List[float] = [0.1, 0.5, 0.9],
                        ploidy: List[int] = [2],
                        estimatePloidy: bool = False,
                        minSegmentBins: int = 25,
                        maxCN: int = 7,
                        use_normal: bool = False,
                        scStates: List[int] = [1, 3],
                        verbose: bool = False,
                        **kwargs):
    """
    Train Hidden Markov Model (HMM)

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragment
        bin_size : int
            Bin size
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
        **kwargs : dict
            Additional arguments
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]
    sample_name = file_name.split(".bam")[0]

    # Check if file was previously annotated
    pf = log_fragments(pf, frags)

    # Determine bin size
    pf.params["hmm_binsize"] = bin_size

    # Calculate bins
    bins, segments = call_cnvs(data = frags,
                                bin_size = bin_size,
                                genome_version = str(pf.params["ref_genome"]),
                                method = method,
                                outlier_smooth = outlier_smooth,
                                gauss_smooth = gauss_smooth,
                                verbose = verbose,
                                cutoff = bcp_cutoff,
                                use_normal = use_normal)

    # Define gender
    try:
        gender = pf.obs_values[sample_name]["gender"]
    except KeyError:
        gender = None

    # Train HMM
    hmm_states = train_hmm(bins,
                            normal = normal,
                            ploidy = ploidy,
                            gender = gender,
                            estimatePloidy = estimatePloidy,
                            minSegmentBins = minSegmentBins,
                            maxCN = maxCN,
                            verbose = verbose,
                            scStates = scStates,
                            **kwargs)

    pf.uns[sample_name] = {"hmm_states": hmm_states}

    return pf


def process_cnvs(pf: ProjectFrame,
                 frags: ngs.Fragments = None,
                 file_name: str = None,
                 merge: bool = True):
    """
    Process CNVs

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragment
        file_name : str
            File name
        merge : bool
            Merge segments
        
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Find file_name
    if frags is None and file_name is None:
        raise NameError("frags or file_name must be provided.")
    elif frags is not None:
        path = os.path.normpath(frags.sam_file)
        file_name = path.split(os.sep)[-1]
        sample_name = file_name.split(".bam")[0]

        # Check if file was previously annotated
        pf = log_fragments(pf, frags)
    
    # Check if hmm has been run
    try:
        pf.uns[sample_name]["hmm_states"]
    except KeyError:
        raise AttributeError("Must run .train_hmm() before .process_cnvs()")

    # Extract sample bins
    sample = pf.intervals["cnv_bins"].loc[:,[sample_name]]
    sample.df.columns = ["ratios"]

    # Merge segments
    if merge:
        merged = merge_segments(pf.obs_intervals[sample_name]["cnv_segments"], sample)
        pf.obs_intervals[sample_name]["cnv_segments"] = merged

    # Process CNVs
    processed_cnvs = hmm_classify(pf.obs_intervals[sample_name]["cnv_segments"],
                                  pf.uns[sample_name]["hmm_states"])
    validate_distributions(processed_cnvs, sample)
    pf.add_obs_intervals(sample_name, "cnv_segments", processed_cnvs)

    return pf


def call_cnvs_single(pf: ProjectFrame,
                     frags: ngs.Fragments,
                     cnv_binsize: int = 100000,
                     hmm_binsize: int = 1000000,
                     method: str = "online_both",
                     outlier_smooth: str = True,
                     gauss_smooth: bool = False,
                     bcp_cutoff: float = 0.3,
                     normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                     ploidy: List[int] = [2],
                     estimatePloidy: bool = False,
                     minSegmentBins: int = 25,
                     maxCN: int = 5,
                     hazard: int = 100,
                     shuffles: int = 5000,
                     p: float = 0.00005,
                     use_normal: bool = False,
                     keep_sex_chroms: bool = False,
                     scStates: List[int] = [1, 3],
                     verbose: bool = False):
    """
    Call CNVs in a single sample

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        frags : Fragments
            Fragment
        cnv_binsize : int
            CNV bin size
        hmm_binsize : int
            HMM bin size
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        hazard : int
            Hazard
        shuffles : int
            Shuffles
        p : float
            P value
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]
    sample_name = file_name.split(".bam")[0]

    # Check if file was previously annotated
    pf = log_fragments(pf, frags)

    # Check fragments have enough chromosomes
    label_counts = frags.frags.label_counts
    if sum([label_counts[c] > 10 for c in label_counts]) < 10:
        print("Not enough chromosomes in", sample_name)
        return pf
    
    # Predict sex
    gender = sex_prediction(frags, bin_size=cnv_binsize)
    pf.add_anno("predicted_gender", sample_name, gender["gender"])

    # Predict purity and ploidy
    pf = predict_purity(pf,
                        frags,
                        bin_size = hmm_binsize,
                        method = method,
                        outlier_smooth = outlier_smooth,
                        gauss_smooth = gauss_smooth,
                        bcp_cutoff = bcp_cutoff,
                        normal = normal,
                        ploidy = [1, 2, 3],
                        estimatePloidy = True,
                        minSegmentBins = minSegmentBins,
                        maxCN = maxCN,
                        verbose = verbose,
                        use_normal = use_normal,
                        scStates = scStates)

    # Train HMM
    pf = calculate_hmm_states(pf,
                            frags,
                            bin_size = hmm_binsize,
                            method = method,
                            outlier_smooth = outlier_smooth,
                            gauss_smooth = gauss_smooth,
                            bcp_cutoff = bcp_cutoff,
                            normal = normal,
                            ploidy = ploidy,
                            estimatePloidy = estimatePloidy,
                            minSegmentBins = minSegmentBins,
                            maxCN = maxCN,
                            use_normal = use_normal,
                            scStates = scStates,
                            verbose = verbose)

    # Call Copy Number Variations
    pf = predict_cnvs(pf,
              frags,
              bin_size = cnv_binsize,
              method = method,
              outlier_smooth = outlier_smooth,
              gauss_smooth = gauss_smooth,
              verbose = verbose,
              cutoff = bcp_cutoff,
              hazard = hazard,
              shuffles = shuffles,
              p = p,
              use_normal = use_normal,
              keep_sex_chroms = keep_sex_chroms)
    pf = process_cnvs(pf,
                    frags)

    # Calculate segment variance
    pf.obs_intervals[sample_name]["cnv_segments"].annotate(pf.intervals["cnv_bins"], sample_name, "var")

    return pf


def call_cnvs_multiple(pf: ProjectFrame,
                     directory: str,
                     cnv_binsize: int = 100000,
                     hmm_binsize: int = 1000000,
                     method: str = "online_both",
                     outlier_smooth: str = True,
                     gauss_smooth: bool = False,
                     bcp_cutoff: float = 0.3,
                     normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                     ploidy: List[int] = [2],
                     estimatePloidy: bool = False,
                     minSegmentBins: int = 25,
                     maxCN: int = 5,
                     hazard: int = 100,
                     shuffles: int = 5000,
                     p: float = 0.00005,
                     nthreads: int = 1,
                     use_normal: bool = False,
                     verbose: bool = False):
    """
    Call CNVs for all samples in a directory

    Parameters
    ----------
        pf : ProjectFrame
            ProjectFrame
        directory : str
            Directory
        cnv_binsize : int
            CNV bin size
        hmm_binsize : int
            HMM bin size
        method : str
            Method for CNV calling
        outlier_smooth : bool
            Smooth outliers
        gauss_smooth : bool
            Smooth Gaussian
        bcp_cutoff : float
            BCP cutoff
        normal : list
            Normal values
        ploidy : list
            Ploidy values
        estimatePloidy : bool
            Estimate ploidy
        minSegmentBins : int
            Minimum segment bins
        maxCN : int
            Maximum copy number
        hazard : int
            Hazard
        shuffles : int
            Shuffles
        p : float
            P value
        nthreads : int
            Number of threads
        wgbs : bool
            Whole genome bisulfite sequencing
        verbose : bool
            Verbose
    
    Returns
    -------
        pf : ProjectFrame
            ProjectFrame
    """

    for file in glob.glob(directory + "*.bam"):
        frags = from_sam(file, nthreads=nthreads, verbose=verbose)
        pf = call_cnvs_single(pf,
                     frags,
                     cnv_binsize,
                     hmm_binsize,
                     method,
                     outlier_smooth,
                     gauss_smooth,
                     bcp_cutoff,
                     normal,
                     ploidy,
                     estimatePloidy,
                     minSegmentBins,
                     maxCN,
                     hazard,
                     shuffles,
                     p,
                     use_normal,
                     verbose)

    return pf