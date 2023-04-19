import os
import pysam

from ..fragments import Fragments
from .parse_sam import read_fragments


def get_chromosomes(sam_file_fn: str):
	"""
	Read genome file
	"""

	# Open SAM file
	sam_file = pysam.Samfile(sam_file_fn)

	# Read chromosomes in SAM file
	chroms = sam_file.references

	# Close SAM file
	sam_file.close()

	return chroms


def from_sam(sam_fn: str = None,
                chrom: str = None,
                min_size: int = 1,
                max_size: int = 1000,
                paired: bool = True,
                qcfail: bool = False,
                mapq_cutoff: int = 25,
                verbose: bool = False,
                nthreads: int = 1,
                proportion: float = 1.0,
                genome: str = "hg19",
                n_frags: int = None):
    """
    Initialize fragments class

    Params
    ------
        sam_fn : str
            SAM file
        chrom : str
            Chromosome
        sam_file : str
            SAM file
        min_size : int
            Minimum fragment size
        max_size : int
            Maximum fragment size
        paired : bool
            Whether to use paired-end reads
        qcfail : bool
            Whether to use QC failed reads
        mapq_cutoff : int
            Minimum mapping quality
        verbose : bool
            Whether to print progress
        n_jobs : int
            Number of threads
        proportion : float
            Proportion of fragments to use
        genome : str
            Genome version
        n_frags : int
            Number of fragments to use
    
    Returns
    -------
        fragments : :class:`~fragments.fragments`
            :class:`~fragments.fragments`
    """

    # Assign bam file
    path = os.path.normpath(sam_fn)
    file_name = path.split(os.sep)[-1]
    sam_file = file_name

    # Read genome
    chromosomes = get_chromosomes(sam_fn)
    add_chr = True
    if chromosomes[0].startswith("chr"):
        add_chr = False

    # Add fragment intervals
    if verbose: print("Reading")
    if chrom is None:
        frags = read_fragments(sam_fn, min_size, max_size, paired, qcfail, mapq_cutoff, proportion, nthreads=nthreads, add_chr=add_chr)
    else:
        frags = read_fragments(sam_fn, min_size, max_size, paired, qcfail, mapq_cutoff, proportion, chrom, 1, genome[chrom])

    # Downsample
    if n_frags is not None:
        proportion = int(frags.size * proportion)
    if proportion != 1:
        frags = frags.downsample(proportion)

    # Build
    fragments = Fragments(frags,
                            sam_file=sam_file,
                            genome=genome)

    return fragments