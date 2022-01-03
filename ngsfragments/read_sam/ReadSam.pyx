#cython: embedsignature=True
#cython: profile=False
#cython: language_level=3

cimport cython
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.stdlib cimport malloc, free, rand, RAND_MAX, srand
from libc.stdio cimport printf

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add

#import pysam
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, hts_itr_next, hts_get_bgzfp
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion, IteratorRowAllRefs

import numpy as np
cimport numpy as np

import time



cdef labeled_aiarray_t *sam_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion):

	cdef labeled_aiarray_t *cintervals = labeled_aiarray_init()
	# Set random generator seed
	srand(time.time())

	# Initialize
	sam_iter_add(samfile_name, cintervals, min_size, max_size, paired, qcfail, mapq_cutoff,
					  proportion)

	return cintervals



def read_fragments(str samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
					float proportion):
	"""
	"""
	cdef bytes bname = samfile_name.encode()
	cdef char *name = bname
	cdef labeled_aiarray_t *cintervals = sam_read(name, min_size, max_size, paired, qcfail, mapq_cutoff,
				  								  proportion)

	print("Done reading")
	cdef LabeledIntervalArray fragments = LabeledIntervalArray()
	fragments.set_list(cintervals)

	return fragments
	