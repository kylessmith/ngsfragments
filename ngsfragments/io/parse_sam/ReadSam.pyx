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
	# Set random generator seed
	srand(time.time())

	# Initialize
	sam_iter_add(samfile_name, cintervals, min_size, max_size, paired, qcfail, mapq_cutoff,
					  proportion, nthreads, add_chr)

	return cintervals


cdef labeled_aiarray_t *sam_nucleosome_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 			float proportion, int nthreads, int add_chr, int fixed_size):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# Set random generator seed
	srand(time.time())

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


cdef class MethylFragment(object):
	"""
	"""

	def __cinit__(self):
		#self.chrom = ""
		self.start = 0
		self.end = 0
		self.length = 0
		self.methyl = np.zeros(0, dtype = np.int8)
		self.pos = np.zeros(0, dtype = np.int32)


	def __init__(self):
		"""
		"""
		pass

	def __repr__(self):
		"""
		"""
		return "MethylFragment(start = {}, end = {}, length = {}, methyl = {})".format(self.start, self.end, self.length, self.methyl)

	
	cdef void set_frag(MethylFragment self, methyl_read_t *frag):
		"""
		"""

		#cdef bytes bchrom = frag.chrom
		#self.chrom = bchrom.decode()
		self.start = frag.start
		self.end = frag.end
		self.length = frag.length
		#self.cpgs.set_list(frag.cpgs)
		self.pos = np.zeros(frag.ncpgs, dtype = np.int32)
		self.methyl = np.zeros(frag.ncpgs, dtype = np.int8)
		cdef int i
		for i in range(frag.ncpgs):
			self.pos[i] = frag.pos[i]
			self.methyl[i] = frag.methyl[i]
		
		return


def read_methyl(str filename,
				str chrom,
				str genome_version = "hg19",
				int min_size = 1,
				int max_size = 1000,
				int mapq_cutoff = 13,
				float proportion = 1.0,
				int nthreads = 1):
	"""
	"""

	# Set random generator seed
	srand(time.time())

	# Get genome file
	import genome_info
	genome = genome_info.GenomeInfo(genome_version)
	ref_filename = genome.seq_file

	cdef bytes bfname = filename.encode()
	cdef char *fname = bfname
	cdef bytes bname = ref_filename.encode()
	cdef char *name = bname
	cdef bytes bchrom = chrom.encode()
	cdef char *c_chrom = bchrom

	#iterate_over_methyl(fname, bname, c_chrom, -1000, 1000, 1, 0, 13, 1.0, 2)

	cdef methyl_read_iterator_t *iterator = methyl_read_iterator_init(fname,
																bname,
																bchrom,
																min_size,
																max_size,
																0,
																mapq_cutoff,
																proportion,
																nthreads)
	
	cdef MethylFragment frag
	
	cdef int i
	cdef int n = 0
	while methyl_read_iterator_next(iterator) >= 1:
		n += 1
		frag = MethylFragment()
		frag.set_frag(iterator.methyl_pair)
		yield frag

	methyl_read_iterator_destroy(iterator)

	return


def fetch_cpgs(str chrom,
				str genome_version = "hg19"):
	"""
	"""

	# Set random generator seed
	srand(time.time())

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

	# Set random generator seed
	srand(time.time())

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

def bounds_motif_enrichment(str bam_file_path,
							str chromosome,
							str prefix = "",
							str genome_version = "hg38",
							int n_bases = 21,
							int min_size = 1,
							int max_size = 1000,
							int qcfail = 0,
							int mapq_cutoff = 13,
							float proportion = 1.0,
							int nthreads = 1):
	"""
	"""

	# Set random generator seed
	srand(time.time())

	# Get genome file
	import genome_info
	genome = genome_info.GenomeInfo(genome_version)
	ref_filename = genome.seq_file

	if prefix == "":
		prefix = bam_file_path.replace(".bam", "")
	
	output_bam_file_path1 = prefix + "_" + chromosome + ".motif1.bam"
	output_bam_file_path2 = prefix + "_" + chromosome + ".motif2.bam"

	cdef bytes ref_2bit_bytes = ref_filename.encode()
	cdef bytes chromosome_bytes = chromosome.encode()
	cdef bytes bam_file_path_bytes = bam_file_path.encode()
	cdef bytes output_bam_file_path1_bytes = output_bam_file_path1.encode()
	cdef bytes output_bam_file_path2_bytes = output_bam_file_path2.encode()

	bounds_motif_split(bam_file_path_bytes,
						ref_2bit_bytes,
						chromosome_bytes,
						output_bam_file_path1_bytes,
						output_bam_file_path2_bytes,
						n_bases,
						min_size,
						max_size,
						1,
						qcfail,
						mapq_cutoff,
						proportion,
						nthreads,
						0)

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


def test_memory(str bam_file_path,
                        str output_bam_file_path1,
                        str output_bam_file_path2,
                        str chromosome,
						str genome_version = "hg19",
                        int min_size1 = 1,
                        int max_size1 = 150,
                        int min_size2 = 151,
                        int max_size2 = 1000,
                        int min_distance = 10,
                        int tolerance = 1,
                        int qcfail = 0,
                        int mapq_cutoff = 13,
                        float proportion = 1.0,
                        int nthreads = 1):
	"""
	"""

	# Set random generator seed
	srand(time.time())

	# Get genome file
	import genome_info
	genome = genome_info.GenomeInfo(genome_version)
	ref_filename = genome.seq_file

	cdef bytes ref_2bit_bytes = ref_filename.encode()
	#cdef char *name = bname
	cdef bytes chromosome_bytes = chromosome.encode()
	#cdef char *c_chrom = bchrom
	cdef bytes bam_file_path_bytes = bam_file_path.encode()
	cdef bytes output_bam_file_path1_bytes = output_bam_file_path1.encode()
	cdef bytes output_bam_file_path2_bytes = output_bam_file_path2.encode()

	test_memory_leak(bam_file_path_bytes,
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
                     qcfail,
                     mapq_cutoff,
                     proportion,
					 nthreads)

	return