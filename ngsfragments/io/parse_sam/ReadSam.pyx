#cython: embedsignature=True
#cython: profile=False
#cython: language_level=3

cimport cython
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.stdlib cimport malloc, free, rand, RAND_MAX, srand
from libc.stdio cimport printf
import math
from intervalframe import IntervalFrame
import pysam
import pandas as pd
cimport pysam.libcalignedsegment as pysam_aligned
cimport pysam.libchtslib as htslib
from pysam.libcalignedsegment cimport AlignedSegment
from pysam import AlignedSegment as PyAlignedSegment

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add

#import pysam
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, hts_itr_next, hts_get_bgzfp
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion, IteratorRowAllRefs

import numpy as np
cimport numpy as np
np.import_array()
import time


cdef np.ndarray pointer_to_numpy_array_int(void *ptr, np.npy_intp size):
	"""
	Convert c pointer to numpy array.
	The memory will be freed as soon as the ndarray is deallocated.

	Parameters
	----------
		ptr : void
			Pointer to be given to numpy
		size : np.npy_intp
			Size of the array

	Returns
	-------
		arr : numpy.ndarray
			Numpy array from given pointer

	"""

	# Import functions for numpy C header
	cdef extern from "numpy/arrayobject.h":
		void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

	# Create shape of ndarray
	cdef np.npy_intp dims[1]
	dims[0] = size

	# Create ndarray from C pointer
	cdef np.ndarray arr = np.PyArray_SimpleNewFromData(1, &dims[0], np.NPY_INT32, ptr)

	# Hand control of data freeing to numpy
	PyArray_ENABLEFLAGS(arr, np.NPY_ARRAY_OWNDATA)
	#np.PyArray_UpdateFlags(arr, arr.flags.num | np.NPY_OWNDATA)

	return arr


cdef labeled_aiarray_t *sam_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion, int nthreads, int add_chr):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# RNG seeds itself on first use per thread via fast_rand_float() — no srand() needed

	# Initialize
	sam_iter_add(samfile_name, cintervals, min_size, max_size, paired, qcfail, mapq_cutoff,
					  proportion, nthreads, add_chr)

	return cintervals


cdef labeled_aiarray_t *sam_nucleosome_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 			float proportion, int nthreads, int add_chr, int fixed_size):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# RNG seeds itself on first use per thread via fast_rand_float() — no srand() needed

	# Initialize
	sam_nucleosome_add(samfile_name, cintervals, min_size, max_size, paired, fixed_size, qcfail, mapq_cutoff,
					  proportion, nthreads, add_chr)

	return cintervals


def read_fragments(str samfile_name,
					int min_size,
					int max_size,
					int paired,
					int qcfail,
					int mapq_cutoff,
					float proportion,
					nthreads = 1,
					add_chr = False,
					nucleosome_adjust = False,
					fixed_size = 74):
	"""
	"""
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef str chrom_pos_str
	cdef bytes chrom_pos_bytes
	cdef const char* chrom_pos

	cdef labeled_aiarray_t *cintervals
	if nucleosome_adjust:
		cintervals = sam_nucleosome_read(name, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion, nthreads, int(add_chr), fixed_size)
	else:
		cintervals = sam_read(name, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion, nthreads, int(add_chr))

	print("Done reading", flush=True)
	cdef LabeledIntervalArray fragments = LabeledIntervalArray()
	fragments.set_list(cintervals)

	return fragments


cdef AlignedSegment bam1_to_aligned_segment(bam1_t *bam_record, bam_hdr_t *header):
	"""
	Convert a bam1_t pointer to an AlignedSegment object.

	Parameters:
	-----------
	bam_record : htslib.bam1_t *
	Pointer to the BAM record
	header : htslib.bam_hdr_t * (optional)
	Pointer to the BAM header (needed for some operations)

	Returns:
	--------
	AlignedSegment
	The constructed AlignedSegment object
	"""
	cdef AlignedSegment aligned_segment = AlignedSegment()

	# Copy the bam1_t data into the AlignedSegment
	# First, ensure the internal bam1_t is allocated
	if aligned_segment._delegate == NULL:
		aligned_segment._delegate = htslib.bam_init1()

	# Copy the data from the source bam1_t to the AlignedSegment's internal bam1_t
	htslib.bam_copy1(aligned_segment._delegate, bam_record)

	# If header is provided, associate it
	#if header != NULL:
	#	aligned_segment.header = header

	return aligned_segment


cdef class MethylFragment(object):
	"""
	"""

	def __cinit__(self):
		self.name = ""
		self.start = 0
		self.end = 0
		self.length = 0
		self.read1_start = 0
		self.read1_end = 0
		self.read2_start = 0
		self.read2_end = 0
		self.methyl = np.zeros(0, dtype = np.int8)
		self.pos = np.zeros(0, dtype = np.int64)
		self.read1_pos = np.zeros(0, dtype = np.int64)
		self.read1_methyl = np.zeros(0, dtype = np.int8)
		self.read2_pos = np.zeros(0, dtype = np.int64)
		self.read2_methyl = np.zeros(0, dtype = np.int8)

		self.read1 = None
		self.read2 = None


	def __init__(self):
		"""
		"""
		pass

	def __repr__(self):
		"""
		"""
		return "MethylFragment(start = {}, end = {}, length = {}, methyl = {})".format(self.start, self.end, self.length, self.methyl)

	
	cdef void set_frag(MethylFragment self, methyl_read_t *frag, methyl_read_t *read1, methyl_read_t *read2):
		"""
		"""

		#cdef bytes bchrom = frag.chrom
		#self.chrom = bchrom.decode()
		self.name = frag.name.decode()
		self.start = frag.start
		self.end = frag.end
		self.length = frag.length

		self.read1 = bam1_to_aligned_segment(read1.read, read1.header)
		self.read2 = bam1_to_aligned_segment(read2.read, read2.header)

		self.read1_start = read1.start
		self.read1_end = read1.end
		self.read2_start = read2.start
		self.read2_end = read2.end
		self.read1_strand = read1.strand
		self.read2_strand = read2.strand

		#self.cpgs.set_list(frag.cpgs)
		self.pos = np.zeros(frag.ncpgs, dtype = np.int64)
		self.methyl = np.zeros(frag.ncpgs, dtype = np.int8)
		self.read1_pos = np.zeros(read1.ncpgs, dtype = np.int64)
		self.read1_methyl = np.zeros(read1.ncpgs, dtype = np.int8)
		self.read2_pos = np.zeros(read2.ncpgs, dtype = np.int64)
		self.read2_methyl = np.zeros(read2.ncpgs, dtype = np.int8)
		cdef int i
		for i in range(frag.ncpgs):
			self.pos[i] = frag.pos[i]
			self.methyl[i] = frag.methyl[i]
		for i in range(read1.ncpgs):
			self.read1_pos[i] = read1.pos[i]
			self.read1_methyl[i] = read1.methyl[i]
		for i in range(read2.ncpgs):
			self.read2_pos[i] = read2.pos[i]
			self.read2_methyl[i] = read2.methyl[i]
		return


def read_fragments_region(str samfile_name,
						str chromosome,
						int min_size,
						int max_size,
						int paired,
						int qcfail,
						int mapq_cutoff,
						float proportion,
						nthreads = 1,
						add_chr = False):
	"""
	Read fragments from a single chromosome using the BAM index.

	Returns (starts, ends, chromosome) as numpy arrays — picklable, so safe to
	call from a joblib worker process and return to the parent for merging.

	Parameters
	----------
	samfile_name : str
		Path to BAM file (must have a .bai or .csi index).
	chromosome : str
		Chromosome name as it appears in the BAM header (e.g. "chr1" or "1").
	min_size, max_size : int
		Fragment length filter bounds.
	paired : int
		1 = paired-end (use isize), 0 = single-end (use read length).
	qcfail : int
		0 = exclude QC-failed reads (default), 1 = include them.
	mapq_cutoff : int
		Minimum mapping quality.
	proportion : float
		Downsampling fraction in (0, 1].
	nthreads : int
		htslib BGZF decompression threads.
	add_chr : bool
		Prepend "chr" to chromosome names if not already present.
	"""
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef bytes bchrom = chromosome.encode()
	cdef char *c_chrom = bchrom

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()

	sam_iter_add_region(name, cintervals, c_chrom,
						min_size, max_size, paired, qcfail, mapq_cutoff,
						proportion, nthreads, int(add_chr))

	cdef LabeledIntervalArray la = LabeledIntervalArray()
	la.set_list(cintervals)

	# Empty chromosomes: .starts/.ends raise IndexError on a zero-size buffer.
	# Return empty arrays — the parent process filters these out before merging.
	if la.size == 0:
		return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32), chromosome

	la.construct()  # build augmented index — required before .starts/.ends are accessible

	# Return numpy arrays + chromosome name — all picklable for joblib inter-process return
	return la.starts.copy(), la.ends.copy(), chromosome


def combine_fragment_arrays(list results):
	"""
	Merge per-chromosome (starts, ends, chrom) tuples into a single LabeledIntervalArray.

	Called in the main process after collecting joblib worker results.
	The chromosome name is encoded once per chunk — O(n_chroms) string ops total
	rather than O(n_frags).

	Parameters
	----------
	results : list of (starts: np.ndarray, ends: np.ndarray, chrom: str)
		Output of read_fragments_region calls collected from parallel workers.

	Returns
	-------
	LabeledIntervalArray
	"""
	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	cdef np.ndarray[np.int32_t, ndim=1] starts
	cdef np.ndarray[np.int32_t, ndim=1] ends
	cdef bytes clabel
	cdef char *c_clabel
	cdef int i, n

	for starts_arr, ends_arr, chrom in results:
		starts  = np.asarray(starts_arr, dtype=np.int32)
		ends    = np.asarray(ends_arr,   dtype=np.int32)
		# Encode chromosome name once per chunk, not once per interval
		clabel  = chrom.encode() if isinstance(chrom, str) else chrom
		c_clabel = clabel
		n = starts.shape[0]
		for i in range(n):
			labeled_aiarray_add(cintervals, starts[i], ends[i], c_clabel)

	cdef LabeledIntervalArray fragments = LabeledIntervalArray()
	fragments.set_list(cintervals)
	return fragments


def read_methyl(str filename,
				str chrom,
				int start = -1,
				int end = -1,
				str genome_version = "hg19",
				str ref_file = None,
				int min_size = 1,
				int max_size = 1000,
				int mapq_cutoff = 13,
				float proportion = 1.0,
				int nthreads = 1):
	"""
	Iterate over methylated read pairs from a BAM file.

	Parameters
	----------
	filename : str
		Path to the BAM file.
	chrom : str
		Chromosome name.
	start, end : int, optional
		Region bounds (default: whole chromosome).
	genome_version : str, optional
		Genome version passed to genome_info when ref_file is not provided.
	ref_file : str, optional
		Explicit path to a reference file (.2bit or FASTA/.fa/.fa.gz).
		When provided, genome_version is ignored.  A FASTA must be indexed
		(samtools faidx); the index will be created automatically if absent.
	"""

	# RNG seeds itself on first use per thread via fast_rand_float() — no srand() needed

	# Resolve reference path: explicit file takes priority over genome_info
	if ref_file is None:
		import genome_info
		genome = genome_info.GenomeInfo(genome_version)
		ref_file = genome.seq_file

	# Open BAM header then close the file handle immediately
	with pysam.AlignmentFile(filename, "rb") as _bam:
		header = _bam.header

	cdef bytes bfname = filename.encode()
	cdef char *fname = bfname
	cdef bytes bref = ref_file.encode()
	cdef char *c_ref = bref
	cdef bytes bchrom = chrom.encode()
	cdef char *c_chrom = bchrom

	cdef methyl_read_iterator_t *iterator = methyl_read_iterator_init(fname,
																c_ref,
																c_chrom,
																start,
																end,
																min_size,
																max_size,
																0,
																mapq_cutoff,
																proportion,
																nthreads)
	if iterator == NULL:
		raise ValueError(
			f"Could not initialise methylation iterator for chromosome '{chrom}'. "
			f"Failed to open reference file '{ref_file}', or the chromosome was not "
			"found under any known naming convention. Check that the reference path "
			"is correct and the genome build matches the BAM."
		)

	cdef MethylFragment frag
	cdef int i
	cdef int n = 0
	cdef int status
	try:
		while True:
			status = methyl_read_iterator_next(iterator)
			if status < 0:
				raise RuntimeError(
					f"Failed to fetch reference sequence window for chromosome '{chrom}'. "
					"The 2bit file may be truncated or corrupted."
				)
			if status < 1:
				break
			n += 1
			frag = MethylFragment()
			frag.set_frag(iterator.methyl_pair, iterator.read1, iterator.read2)
			frag.read1.header = header
			frag.read2.header = header
			yield frag
	finally:
		methyl_read_iterator_destroy(iterator)

	return


def fetch_cpgs(str chrom,
				str genome_version = "hg19"):
	"""
	"""

	# RNG seeds itself on first use per thread via fast_rand_float() — no srand() needed

	# Get genome file
	import genome_info
	genome = genome_info.GenomeInfo(genome_version)
	ref_filename = genome.seq_file

	cdef bytes bname = ref_filename.encode()
	cdef char *name = bname
	cdef bytes bchrom = chrom.encode()
	cdef char *c_chrom = bchrom

	cdef reference_cpgs_t *cpgs = fetch_reference_cpgs(name, c_chrom)
	cdef np.ndarray pos = pointer_to_numpy_array(cpgs.pos, cpgs.size)
	cdef np.ndarray strand = pointer_to_numpy_array_int(cpgs.strand, cpgs.size)
	
	values = pd.Series(strand, index=pos)

	return values


def methyl_length_decompose(str bam_file_path,
						str chromosome,
                        str prefix = "",
						str genome_version = "hg38",
                        int min_size1 = 1,
                        int max_size1 = 150,
                        int min_size2 = 151,
                        int max_size2 = 1000,
                        int min_distance = 10,
                        int tolerance = 10,
						int max_iter = 10,
                        int qcfail = 0,
                        int mapq_cutoff = 13,
                        float proportion = 1.0,
                        int nthreads = 1):
	"""
	"""

	# RNG seeds itself on first use per thread via fast_rand_float() — no srand() needed

	# Get genome file
	import genome_info
	genome = genome_info.GenomeInfo(genome_version)
	ref_filename = genome.seq_file

	if prefix == "":
		prefix = bam_file_path.replace(".bam", "")
	
	output_bam_file_path1 = prefix + "_" + chromosome + ".1.bam"
	output_bam_file_path2 = prefix + "_" + chromosome + ".2.bam"

	cdef bytes ref_2bit_bytes = ref_filename.encode()
	#cdef char *name = bname
	cdef bytes chromosome_bytes = chromosome.encode()
	#cdef char *c_chrom = bchrom
	cdef bytes bam_file_path_bytes = bam_file_path.encode()
	cdef bytes output_bam_file_path1_bytes = output_bam_file_path1.encode()
	cdef bytes output_bam_file_path2_bytes = output_bam_file_path2.encode()

	split_methyl_bam(bam_file_path_bytes,
                     output_bam_file_path1_bytes,
                     output_bam_file_path2_bytes,
                     ref_2bit_bytes,
                     chromosome_bytes,
                     min_size1,
                     max_size1,
                     min_size2,
                     max_size2,
                     min_distance,
                     tolerance,
					 max_iter,
                     qcfail,
                     mapq_cutoff,
                     proportion,
					 nthreads)

	return


def merge_bams(inputs, output):
	"""
	"""
	import pysam

	# Open output BAM file
	input_bam = pysam.AlignmentFile(inputs[0], "rb")
	output_bam = pysam.AlignmentFile(output, "wb", template=input_bam)
	input_bam.close()

	for input1 in inputs:
		# Open input BAM files
		input_bam = pysam.AlignmentFile(input1, "rb")

		# Merge BAM files
		for read in input_bam:
			output_bam.write(read)

		# Close BAM files
		input_bam.close()

	output_bam.close()

	return

