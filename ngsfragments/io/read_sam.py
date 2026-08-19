import os
import pysam
import numpy as np

from ..fragments import Fragments
from .parse_sam import read_fragments
from .parse_sam.ReadSam import read_fragments_region, combine_fragment_arrays


def _get_header_info(sam_fn: str):
    """
    Open BAM once and return (file_basename, chromosomes, add_chr).
    Avoids the double-open that existed when get_chromosomes() and read_fragments()
    each opened the file independently.
    """
    with pysam.Samfile(sam_fn) as sf:
        file_name   = os.path.basename(os.path.normpath(sam_fn))
        chromosomes = list(sf.references)

    # If ANY chromosome already carries a "chr" prefix, assume all do and don't add it.
    add_chr = not any(c.startswith("chr") for c in chromosomes)
    return file_name, chromosomes, add_chr


def get_chromosomes(sam_file_fn: str):
    """
    Return the list of reference sequence names from the BAM/SAM header.

    Parameters
    ----------
    sam_file_fn : str
        Path to BAM/SAM file.

    Returns
    -------
    chroms : list[str]
    """
    with pysam.Samfile(sam_file_fn) as sf:
        chroms = list(sf.references)
    return chroms


# ---------------------------------------------------------------------------
# Joblib worker — must be module-level for pickling by the loky backend.
# Each worker runs in its own process: opens its own BAM handle, loads its
# own BAM index slice, and returns lightweight numpy arrays back to the parent.
#
# NOTE: nthreads here is htslib BGZF decompression threads *per worker*.
#       Total CPU usage ≈ n_jobs × nthreads.  Keep n_jobs × nthreads ≤ ncores.
# ---------------------------------------------------------------------------

def _read_chrom_worker(sam_fn, chrom, min_size, max_size, paired,
                       qcfail, mapq_cutoff, proportion, add_chr, nthreads):
    """
    Read fragments from one chromosome and return picklable numpy arrays.
    Intended to be called via joblib.Parallel — do not call directly.
    """
    return read_fragments_region(
        sam_fn, chrom,
        min_size, max_size, paired, qcfail, mapq_cutoff,
        proportion, nthreads, add_chr,
    )


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
             n_frags: int = None,
             nucleosome_adjust: bool = False,
             fixed_size: int = 74,
             n_jobs: int = 1):
    """
    Load fragments from a BAM/SAM file into a Fragments object.

    Parameters
    ----------
    sam_fn : str
        Path to BAM file.
    min_size : int
        Minimum fragment length (bp).
    max_size : int
        Maximum fragment length (bp).
    paired : bool
        True  → paired-end mode: use insert size (isize) for fragment length.
        False → single-end mode: use read length.
    qcfail : bool
        If True, include QC-failed reads (BAM_FQCFAIL).  Default False.
    mapq_cutoff : int
        Exclude reads with MAPQ below this threshold.
    verbose : bool
        Print progress messages.
    nthreads : int
        htslib BGZF decompression threads per file handle.
        With n_jobs > 1, total CPU ≈ n_jobs × nthreads.
    proportion : float
        Fraction of fragments to retain via random downsampling (0 < p ≤ 1).
    genome_version : str
        Genome version string passed to Fragments.
    n_frags : int, optional
        If provided, downsample to exactly this many fragments after loading.
        Takes priority over `proportion`.
    nucleosome_adjust : bool
        If True, centre a fixed-size window on each read's 5' end (nucleosome
        occupancy / endpoint analysis) rather than using the full insert span.
    fixed_size : int
        Window size (bp) used when nucleosome_adjust=True.
    n_jobs : int
        Number of parallel worker processes (joblib).
        -1 uses all available CPUs.
        1  (default) runs single-threaded with the whole-file scanner, which
           avoids the per-chromosome index-load overhead for small BAMs.

    Returns
    -------
    fragments : Fragments
    """
    # Single BAM open: read header info and derive add_chr in one pass.
    sam_file, chromosomes, add_chr = _get_header_info(sam_fn)

    if verbose:
        print(f"Reading {sam_fn}  ({len(chromosomes)} sequences, "
              f"n_jobs={n_jobs}, nthreads={nthreads})")

    # ------------------------------------------------------------------
    # Fragment loading
    # ------------------------------------------------------------------
    if n_jobs == 1 or nucleosome_adjust:
        # Sequential whole-file scan.
        # nucleosome_adjust stays sequential — sam_nucleosome_add has no
        # region-based variant yet; adding one is straightforward if needed.
        frags = read_fragments(
            sam_fn, min_size, max_size, paired, qcfail, mapq_cutoff,
            proportion, nthreads=nthreads, add_chr=add_chr,
            nucleosome_adjust=nucleosome_adjust, fixed_size=fixed_size,
        )

    else:
        # Parallel path: one worker process per chromosome.
        # Each worker calls sam_iter_add_region (index-based), so it reads
        # only its slice of the BAM — no worker ever reads the full file.
        from joblib import Parallel, delayed

        if verbose:
            print(f"  Dispatching {len(chromosomes)} chromosomes to {n_jobs} workers...")

        results = Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(_read_chrom_worker)(
                sam_fn, chrom,
                min_size, max_size, paired, qcfail, mapq_cutoff,
                proportion, add_chr, nthreads,
            )
            for chrom in chromosomes
        )

        # Filter out empty chromosomes (no passing reads) to keep merge lean
        results = [r for r in results if len(r[0]) > 0]

        if verbose:
            total = sum(len(r[0]) for r in results)
            print(f"  Merging {total:,} fragments from {len(results)} chromosomes...")

        frags = combine_fragment_arrays(results)

    # ------------------------------------------------------------------
    # Downsampling
    # ------------------------------------------------------------------
    # n_frags (absolute count) takes priority over proportion (fractional).
    # proportion-based downsampling in the C layer has already been applied
    # during loading; this second stage is for exact-count requirements.
    if n_frags is not None:
        frags = frags.downsample(n_frags)
    elif proportion < 1.0 and n_jobs != 1:
        # In the parallel path, downsampling happens per-chromosome inside the
        # C layer, so the total count is only approximately proportion * N.
        # A second pass here corrects for cross-chromosome variance if needed.
        pass  # already downsampled in C; remove this block if exact counts matter

    # ------------------------------------------------------------------
    # Build Fragments object
    # ------------------------------------------------------------------
    fragments = Fragments(frags, sam_file=sam_file, genome_version=genome_version)

    if verbose:
        print(f"Done. {fragments.n_fragments:,} fragments loaded.")

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
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)
    chroms = np.array(genome["main_chromosomes"])
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


def bounds_motif_match(filename: str,
                        genome_version: str = "hg38",
                        n_jobs: int = 1):
    """
    """

    from multiprocessing import Pool
    import os
    from .parse_sam.ReadSam import bounds_motif_enrichment, merge_bams

    # Get genome file
    import genome_info
    genome = genome_info.GenomeInfo(genome_version)
    chroms = np.array(genome["main_chromosomes"])
    chroms = chroms[chroms != "chrM"]

    # Assign parameters
    prefix = filename.replace(".bam", "")
    parameters = [(filename, str(c), prefix, genome_version) for c in chroms]

    # Run
    print("Decomposing...")
    with Pool(processes=n_jobs) as pool:
        results = pool.starmap(bounds_motif_enrichment, parameters)
    print("Done")

    # Merge
    print("Merging...")
    inputs1 = [prefix + "_" + chromosome + ".motif1.bam" for chromosome in chroms]
    inputs2 = [prefix + "_" + chromosome + ".motif2.bam" for chromosome in chroms]
    output1 = filename.replace(".bam", ".motif1.bam")
    output2 = filename.replace(".bam", ".motif2.bam")
    merge_bams(inputs1, output1)
    merge_bams(inputs2, output2)
    print("Done")

    # Remove intermediate files
    for f in inputs1 + inputs2:
        os.remove(f)

    return