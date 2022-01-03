# cython: language_level=3

cimport cython
import numpy as np
cimport numpy as np

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add

from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t, uint16_t, int32_t
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, htsFile, hts_itr_t, hts_idx_t
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion


cdef extern from "read_intervals.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "htslib/hts.h":
	# C is include here so that it doesn't need to be compiled externally
	
	ctypedef int64_t hts_pos_t

cdef extern from "htslib/sam.h":
	# C is include here so that it doesn't need to be compiled externally

	ctypedef struct bam1_core_t:
		hts_pos_t pos
		int32_t tid
		uint16_t bin # NB: invalid on 64-bit pos
		uint8_t qual
		uint8_t l_extranul
		uint16_t flag
		uint16_t l_qname
		uint32_t n_cigar
		int32_t l_qseq
		int32_t mtid
		hts_pos_t mpos
		hts_pos_t isize

	ctypedef struct bam1_t:
		bam1_core_t core
		uint64_t id
		uint8_t *data
		int l_data
		uint32_t m_data
		#uint32_t mempolicy:2, :30 # Reserved

cdef extern from "read_intervals.h":
	# C is include here so that it doesn't need to be compiled externally

	int check_read(bam1_t *aln, int min_size, int max_size, int paired, int qcfail,
			   		int mapq_cutoff, float proportion) nogil

	void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
					int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
					float proportion) nogil


cdef labeled_aiarray_t *sam_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion)



