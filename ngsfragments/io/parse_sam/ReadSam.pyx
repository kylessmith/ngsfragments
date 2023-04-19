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

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add

#import pysam
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, hts_itr_next, hts_get_bgzfp
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion, IteratorRowAllRefs

import numpy as np
cimport numpy as np

import time



cdef labeled_aiarray_t *sam_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion, int nthreads, int add_chr):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# Set random generator seed
	srand(time.time())

	# Initialize
	sam_iter_add(samfile_name, cintervals, min_size, max_size, paired, qcfail, mapq_cutoff,
					  proportion, nthreads, add_chr)

	return cintervals


cdef labeled_aiarray_t *chrom_read(char *samfile_name, const char *chrom_pos, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# Set random generator seed
	srand(time.time())

	# Initialize
	cdef int result = chrom_iter_add(samfile_name, chrom_pos, cintervals, min_size, max_size, paired, qcfail, mapq_cutoff,
					  				 proportion)

	return cintervals


cdef void chrom_bin_read(char *samfile_name, const char *chrom_pos, long[::1] nhits, int bin_size, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
						 float proportion):
	
	# Set random generator seed
	srand(time.time())

	# Initialize
	cdef int result = chrom_bin_iter_add(samfile_name, chrom_pos, bin_size, &nhits[0], min_size, max_size, paired, qcfail, mapq_cutoff,
					  					proportion)

	return


cdef void chrom_bin_read_float(char *samfile_name, const char *chrom_pos, long *nhits, int bin_size, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
						 float proportion):
	
	# Set random generator seed
	srand(time.time())

	# Initialize
	cdef int result = chrom_bin_iter_add(samfile_name, chrom_pos, bin_size, &nhits[0], min_size, max_size, paired, qcfail, mapq_cutoff,
					  					proportion)

	return


def read_fragments(str samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
					float proportion, chrom = None,  start = None, end = None, nthreads = 1, add_chr = False):
	"""
	"""
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef str chrom_pos_str
	cdef bytes chrom_pos_bytes
	cdef const char* chrom_pos

	cdef labeled_aiarray_t *cintervals
	if chrom is None:
		cintervals = sam_read(name, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion, nthreads, int(add_chr))
	else:
		chrom_pos_str = chrom + ":" + str(start) + "-" + str(end)
		chrom_pos_bytes = chrom_pos_str.encode()
		chrom_pos = chrom_pos_bytes
		cintervals = chrom_read(name, chrom_pos, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion)

	print("Done reading", flush=True)
	cdef LabeledIntervalArray fragments = LabeledIntervalArray()
	fragments.set_list(cintervals)

	return fragments
	

def read_chrom_bin_fragments(str samfile_name, int bin_size, str chrom, int chrom_length, int min_size=1, int max_size=1000, int paired=1, int qcfail=0, int mapq_cutoff=25,
							 float proportion=1.0):
	"""
	"""
	
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef str chrom_pos_str
	cdef bytes chrom_pos_bytes
	cdef const char* chrom_pos

	# Determine bins
	cdef int n_bins = math.ceil(chrom_length / bin_size)
	cdef np.ndarray nhits = np.zeros(n_bins, dtype=np.int_)
	cdef long[::1] nhits_mem = nhits

	chrom_pos_str = chrom + ":0-" + str(chrom_length)
	chrom_pos_bytes = chrom_pos_str.encode()
	chrom_pos = chrom_pos_bytes
	chrom_bin_read(name, chrom_pos, nhits_mem, bin_size, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion)

	cdef LabeledIntervalArray bins = LabeledIntervalArray()
	bins.create_bin({chrom:chrom_length}, bin_size)

	iframe_bins = IntervalFrame(bins)
	iframe_bins.df.loc[:,"counts"] = nhits

	return iframe_bins


def read_bin_fragments(str samfile_name, int bin_size=100000, int min_size=1, int max_size=1000, int paired=1, int qcfail=0, int mapq_cutoff=25,
					   float proportion=1.0):
	"""
	"""

	# Read chomosome lengths
	sam_file = pysam.Samfile(samfile_name)
	genome = {}
	genome.update(list(zip(sam_file.references, sam_file.lengths)))
	sam_file.close()

	# Create bins
	cdef LabeledIntervalArray bins = LabeledIntervalArray()
	bins.create_bin(genome, bin_size)

	# Create nhits
	cdef long nbins = bins.size
	cdef np.ndarray nhits = np.zeros(nbins, dtype=int)
	cdef long[::1] nhits_mem = nhits

	# Statis typing
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef str chrom_pos_str
	cdef bytes chrom_pos_bytes
	cdef const char* chrom_pos

	# Iterate over chroms
	cdef int chrom_bins
	cdef str chrom
	cdef bytes bchrom
	cdef long shift = 0
	for chrom in genome:
		bchrom = chrom.encode()
		chrom_bins = math.ceil(genome[chrom] / bin_size)

		chrom_pos_str = chrom + ":0-" + str(genome[chrom])
		chrom_pos_bytes = chrom_pos_str.encode()
		chrom_pos = chrom_pos_bytes
		chrom_bin_read_float(name, chrom_pos, &nhits_mem[shift], bin_size, min_size, max_size, paired, qcfail, mapq_cutoff,
							proportion)
		
		shift += chrom_bins

	# Create IntervalFrame
	iframe = IntervalFrame(bins)
	iframe.df.loc[:,"counts"] = nhits

	return iframe