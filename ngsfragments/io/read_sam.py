import os
import pysam

from ..fragments import Fragments
from .parse_sam import read_fragments


def get_chromosomes(sam_file_fn: str):
	"""
	Read genome file
    
    Parameters
    ----------
        sam_file_fn : str
            SAM file
    
    Returns
    -------
        chroms : list
            List of chromosomes
	"""

	# Open SAM file
	sam_file = pysam.Samfile(sam_file_fn)

	# Read chromosomes in SAM file
	chroms = sam_file.references

	# Close SAM file
	sam_file.close()

	return chroms


def from_sam(sam_fn: str = None,
                min_size: int = 1,
                max_size: int = 1000,
                paired: bool = True,
                qcfail: bool = False,
                mapq_cutoff: int = 25,
                verbose: bool = False,
                nthreads: int = 1,
                proportion: float = 1.0,
                genome_version: str = "hg19",
                n_frags: int = None):
    """
    Initialize fragments class

    Parameters
    ----------
        sam_fn : str
            SAM file
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
    frags = read_fragments(sam_fn, min_size, max_size, paired, qcfail, mapq_cutoff, proportion, nthreads=nthreads, add_chr=add_chr)

    # Downsample
    if n_frags is not None:
        proportion = int(frags.size * proportion)
    if proportion != 1:
        frags = frags.downsample(proportion)

    # Build
    fragments = Fragments(frags,
                            sam_file=sam_file,
                            genome_version=genome_version)

    return fragments


def methyl_length_match(filename: str,
                        genome_version: str = "hg38",
                        n_jobs: int = 1):
    """
    """

    from multiprocessing import Pool
    import os
    from .parse_sam.ReadSam import methyl_length_decompose, merge_bams

    # Get genome file
    if genome_version == "hg19":
        from hg19genome import Hg19Genome as Genome
    elif genome_version == "hg38":
        from hg38genome import Hg38Genome as Genome
    else:
        raise ValueError("Genome version not recognized")
    genome = Genome()
    chroms = genome.main_chromosomes
    chroms = chroms[chroms != "chrM"]

    # Assign parameters
    prefix = filename.replace(".bam", "")
    parameters = [(filename, str(c), prefix, genome_version) for c in chroms]

    # Run
    print("Decomposing...")
    with Pool(processes=n_jobs) as pool:
        results = pool.starmap(methyl_length_decompose, parameters)
    print("Done")

    # Merge
    print("Merging...")
    inputs1 = [prefix + "_" + chromosome + ".1.bam" for chromosome in chroms]
    inputs2 = [prefix + "_" + chromosome + ".2.bam" for chromosome in chroms]
    output1 = filename.replace(".bam", ".1.bam")
    output2 = filename.replace(".bam", ".2.bam")
    merge_bams(inputs1, output1)
    merge_bams(inputs2, output2)
    print("Done")

    # Remove intermediate files
    for f in inputs1 + inputs2:
        os.remove(f)

    return