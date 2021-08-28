import numpy as np
import math
import matplotlib.pyplot as plt
import h5py
import pysam
from collections import Counter
import pyBigWig
import urllib
import gzip
import os
from ailist import AIList


def read_bed(bed_fn, genome_style = "UCSC"):
    """
    Read bed file
    
    Params
    ---------
        bed_fn
            str (Name of the bed file to read)
        genome_style
            str (Naming style: 'UCSC'='chr1')
    
    Returns
    ----------
        bed_record
            dict (Dictionary of AILists per chromosome of bed regions)
    """
    
    # Determine if bed is url
    if "http" in bed_fn:
        # Read url bed
        response = urllib.request.urlopen(bed_fn)
        if ".gz" in url:
            gunzip_response = gzip.GzipFile(fileobj=response)
            bed_file = gunzip_response.readlines()
        else:
            bed_file = response.readlines()
    else:
        bed_file = open(bed_fn, "r")

    # Initialize record
    bed_record = {}
    previous_chrom = None

    # Iterate over bed
    for line in bed_file:
        # Determine fields
        if isinstance(line, bytes):
            fields = line.strip().split(b'\t')
        else:
            fields = line.strip().split('\t')
            
        # Determine naming style
        chrom = fields[0]
        if genome_style == "UCSC":
            if "chr" not in chrom:
                chrom = "chr" + chrom
        start = int(fields[1])
        end = int(fields[2])

        # Check if new chromosome
        if chrom != previous_chrom:
            bed_record[chrom] = AIList()
            previous_chrom = chrom

        # Add fields to record
        bed_record[chrom].add(start, end)

    return bed_record


def get_bias_files(fasta_fn, mappability_fn="wgEncodeCrgMapabilityAlign50mer",
                   blacklist_fn="hg19", genome_style = "UCSC", verbose=False):
    """
    Organize bias files
    
    Params
    ---------
        fasta_fn
            str (Name of FASTA file)
        mappability_fn
            str (Name of mappability file)
        blacklist_fn
            str (Name of mappability file)
        genome_style
            str (Naming style: 'UCSC'='chr1')
        verbose
            bool (Flag for print statements)
    
    Returns
    ----------
        bias_files
            dict (Names of bias files)
    """

    # Set data directory
    data_dir = os.path.split(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0])[0]
    data_dir = os.path.join(data_dir, "data")	

    # Read Fasta
    if verbose: print("   FASTA:", fasta_fn)
    if fasta_fn is not None:
        fasta = pysam.Fastafile(fasta_fn)
    elif fasta_fn is not None:
        fasta = None

    # Read mappability
    if mappability_fn == "wgEncodeCrgMapabilityAlign50mer":
        if verbose: print("   MAPPABILITY:", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig")
        mappability = pyBigWig.open("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig")
    elif mappability_fn is not None:
        if verbose: print("   MAPPABILITY:", "mappability")
        mappability = pyBigWig.open(mappability_fn)
    elif mappability_fn is not None:
        if verbose: print("   MAPPABILITY:", "None")
        mappability = None

    # Read blacklist
    if ".bed" not in blacklist_fn:
        blacklist_fn = os.path.join(data_dir, blacklist_fn + "_QDNAseq_blacklisted.bed")
    if blacklist_fn is not None:
        if verbose: print("   BLACKLIST:", blacklist_fn)
        blacklist = read_bed(blacklist_fn, genome_style)
    elif blacklist_fn is not None:
        if verbose: print("   BLACKLIST:", "None")
        blacklist = None

    # Record file handles
    bias_files = {"fasta":fasta,
                  "map":mappability,
                  "black":blacklist}

    return bias_files

    
def calculate_bias(bias_files, bin_size=100000, verbose=False):
    """
    Calculate bias per bin
    
    Params
    ---------
        bias_files
            dict (Names of bias files)
        bin_size
            int (Number of positions in bin)
        verbose
            bool (Flag for print statements)
    
    Returns
    ----------
        bias_record
            dict (Bias bins for given chromosome)
    """
    
    # Determine chromosomes
    chromosomes = list(zip(bias_files["fasta"].references,
                           bias_files["fasta"].lengths))

    # Initialize bias records
    bias_record = {}
    bias_record["gc"] = {}
    bias_record["repeat"] = {}
    bias_record["n"] = {}
    bias_record["mappability"] = {}
    bias_record["blacklist"] = {}
    
    # Iterate over chromosomes
    for chrom, chrom_length in chromosomes:
        if verbose: print(chrom)
        # Determine number of bins
        n_bins = max(1, int(chrom_length / bin_size))
        if verbose: print("   n_bins:", n_bins)

        # Initialize bias records
        bias_record["gc"][chrom] = np.zeros(n_bins)
        bias_record["repeat"][chrom] = np.zeros(n_bins)
        bias_record["n"][chrom] = np.zeros(n_bins)
        bias_record["mappability"][chrom] = np.zeros(n_bins)
        bias_record["blacklist"][chrom] = np.zeros(n_bins)

        # Iterate over bins
        for i in range(0, n_bins):
            # Define bin
            start = i * bin_size
            end = start + bin_size

            # Iterate over fasta
            sequence = bias_files["fasta"].fetch(reference=chrom, start=start, end=end)
            counts = Counter(sequence)
            # Record content
            bias_record["n"][chrom][i] = counts["N"] / bin_size
            if bias_record["n"][chrom][i] == 1.0:
                bias_record["gc"][chrom][i] = 0.0
                bias_record["repeat"][chrom][i] = 0.0
            else:
                bias_record["gc"][chrom][i] = (counts["G"] + counts["g"] + counts["C"] + counts["c"]) / (bin_size - counts["N"])
                bias_record["repeat"][chrom][i] = (counts["a"] + counts["g"] + counts["t"] + counts["c"]) / (bin_size - counts["N"])

            # Record mappability content
            if bias_files["map"] is not None:
                try:
                    bias_record["mappability"][chrom][i] = np.nansum(bias_files["map"].values(chrom, start, end)) / bin_size
                except RuntimeError:
                    # Try adding the chr
                    try:
                        bias_record["mappability"][chrom][i] = np.nansum(bias_files["map"].values("chr"+chrom, start, end)) / bin_size
                    except RuntimeError:
                        bias_record["mappability"][chrom][i] = 0
            
            # Record blacklist content
            if bias_files["black"] is not None:
                try:
                    bias_record["blacklist"][chrom][i] = np.sum(bias_files["black"][chrom].interval_coverage(start, end).values > 0) / bin_size
                except KeyError:
                    # Try adding the chr
                    try:
                        bias_record["blacklist"][chrom][i] = np.sum(bias_files["black"]["chr"+chrom].interval_coverage(start, end).values > 0) / bin_size
                    except KeyError:
                        bias_record["blacklist"][chrom][i] = 1

    # Close bias files
    bias_files["fasta"].close()
    bias_files["map"].close()
    for chrom in bias_files["black"]:
        bias_files["black"][chrom].close()

    return bias_record


def write_bias_record(bias_record, out_h5_fn="genome_bias.h5", bin_size=100000):
    """
    Write the bias record to h5 file
    
    Params
    ---------
        bias_record
            dict (Bias bins for given chromosome)
        out_h5_fn
            str (Name of output h5 file [default: 'genome_bias.h5'])
        bin_size
            int (Number of positions in bin)
    
    Returns
    ----------
        None
    """

    # Open output h5 file
    out_h5 = h5py.File(out_h5_fn, "w")

    # Iterate over chromosomes
    chroms = np.array([])
    starts = np.array([])
    for chrom in bias_record["gc"].keys():
        n_bins = len(bias_record["gc"][chrom])
        chroms = np.append(chroms, np.repeat(chrom, n_bins))
        starts = np.append(starts, np.arange(0, n_bins * bin_size, bin_size))
    out_h5["chroms"] = chroms.astype(bytes)
    out_h5["starts"] = starts
    out_h5["ends"] = starts + bin_size

    # Iterate over biases
    for bias in bias_record:
        bias_array = np.array([])
        for chrom in bias_record["gc"].keys():
            bias_array = np.append(bias_array, bias_record[bias][chrom])
        out_h5[bias] = bias_array

    # Close output h5 file
    out_h5.close()