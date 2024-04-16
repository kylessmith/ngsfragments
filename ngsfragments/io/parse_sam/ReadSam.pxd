 # cython: language_level=3

cimport cython
import numpy as np
cimport numpy as np

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add
from ailist.array_query_core cimport pointer_to_numpy_array
#from hg19genome.kmers.kmer_reader cimport kmer_count_t

from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t, uint16_t, int32_t, int8_t, int16_t
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, htsFile, hts_itr_t, hts_idx_t
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion


cdef extern from "read_intervals.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "read_iterator.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "methyl_fragment_iter.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "reference_methyl.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "read_name_store.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "methyl_size_split.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "methyl_record.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "merge_bams.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "bounds_motif.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "kmers/interval_kmer.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "kmers/2bit.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "kmers/2bit.h":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "kmers/interval_kmer.h":
	# C is include here so that it doesn't need to be compiled externally
	#pass

	ctypedef struct kmer_t:
		char *name
		int count

	ctypedef struct kmer_count_t:
		int max_kmers
		int n_kmers
		kmer_t *kmers
		void *kmer_lookup

	int fetch_kmer(kmer_count_t *kc, char *seq) nogil

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

	ctypedef struct methyl_read_t:
		char *name
		int start
		int end
		int length
		long *pos
		int8_t *methyl
		uint8_t *qual
		int ncpgs
		int size
		int max_size

	ctypedef struct methyl_read_iterator_t:
		const char *chrom
		methyl_read_t *methyl_pair

	ctypedef struct reference_cpgs_t:
		long *pos
		int *strand
		int size
		int max_size
	
	#==================================================================================================
	# read_intervals.c
	#--------------------------------------------------------------------------------------------------

	# Check that read passed QC
	int check_read(bam1_t *aln, int min_size, int max_size, int paired, int qcfail,
			   		int mapq_cutoff, float proportion) nogil

	# Add reads from sam file to interval list
	void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
					int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
					float proportion, int nthreads, int add_chr) nogil

	# Add reads from sam file to inteval list and adjust for nucleosome occupancy
	void sam_nucleosome_add(char *samfile_name,
							labeled_aiarray_t *intervals,
							int min_size,
							int max_size,
							int paired,
							int fixed_size,
							int qcfail,
							int mapq_cutoff,
							float proportion,
							int nthreads,
							int add_chr) nogil

	#==================================================================================================
	# methyl_fragment_iter.c
	#--------------------------------------------------------------------------------------------------

	# Initialize methyl_read_iterator_t struct
	methyl_read_iterator_t *methyl_read_iterator_init(const char *bam_file_path,
														char *ref_2bit,
														const char *chromosome,
														int min_size,
														int max_size,
														int qcfail,
														int mapq_cutoff,
														float proportion,
														int nthreads) nogil

	# Free memory allocated for methyl_read_iterator_t struct
	void methyl_read_iterator_destroy(methyl_read_iterator_t *iter) nogil

	# Iterate over BAM file
	int methyl_read_iterator_next(methyl_read_iterator_t *iter) nogil


	#==================================================================================================
	# reference_methyl.c
	#--------------------------------------------------------------------------------------------------

	# Initialize reference CpGs struct
	reference_cpgs_t *reference_cpgs_init() nogil

	# Free memory allocated for reference CpGs struct
	void reference_cpgs_destroy(reference_cpgs_t *cpgs) nogil

	# Append CpG to reference CpGs struc
	void reference_cpgs_append(reference_cpgs_t *cpgs, long pos, int strand) nogil

	# Fetch reference CpGs
	reference_cpgs_t *fetch_reference_cpgs(char *ref_2bit, char *chrom) nogil


	#==================================================================================================
	# methyl_size_split.c
	#--------------------------------------------------------------------------------------------------
	
	# Split reads by length/profile and write to output BAM files
	void split_methyl_bam(const char *bam_file_path,
							const char *output_bam_file_path1,
							const char *output_bam_file_path2,
							char *ref_2bit,
							const char *chromosome,
							int min_size1,
							int max_size1,
							int min_size2,
							int max_size2,
							int min_distance,
							int tolerance,
							int max_iter,
							int qcfail,
							int mapq_cutoff,
							float proportion,
                        	int nthreads) nogil

	
	void test_memory_leak(const char *bam_file_path,
                        const char *output_bam_file_path1,
                        const char *output_bam_file_path2,
                        char *ref_2bit,
                        const char *chromosome,
                        int min_size1,
                        int max_size1,
                        int min_size2,
                        int max_size2,
                        int min_distance,
                        int tolerance,
                        int qcfail,
                        int mapq_cutoff,
                        float proportion,
                        int nthreads) nogil

	#==================================================================================================
	# bounds_motif.c
	#--------------------------------------------------------------------------------------------------

	void bounds_motif_split(char *samfile_name,
                                            char *fname,
											char *chromosome,
                                            const char *output_bam_file_path1,
                                            const char *output_bam_file_path2,
                                            int n_bases,
                                            int min_size,
                                            int max_size,
                                            int paired,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads,
                                            int add_chr) nogil

	#==================================================================================================
	# merge_bams.c
	#--------------------------------------------------------------------------------------------------

	#void merge_bams(const char *input1_bam_path, const char *input2_bam_path, const char *output_bam_path) nogil


	#==================================================================================================


cdef labeled_aiarray_t *sam_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 float proportion, int nthreads, int add_chr)

cdef labeled_aiarray_t *sam_nucleosome_read(char *samfile_name, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
								 			float proportion, int nthreads, int add_chr, int fixed_size)

cdef np.ndarray pointer_to_numpy_array_int(void *ptr, np.npy_intp size)


cdef class MethylFragment(object):

	# Attributes
	cdef public str chrom
	cdef public int start
	cdef public int end
	cdef public int length
	cdef public np.ndarray pos
	cdef public np.ndarray methyl

	# Methods
	cdef void set_frag(MethylFragment self, methyl_read_t *frag)
