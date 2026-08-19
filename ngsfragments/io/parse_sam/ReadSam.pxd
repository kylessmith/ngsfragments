 # cython: language_level=3

cimport cython
import numpy as np
cimport numpy as np
np.import_array()

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init, labeled_aiarray_add
from ailist.array_query_core cimport pointer_to_numpy_array
#from hg19genome.kmers.kmer_reader cimport kmer_count_t

from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t, uint16_t, int32_t, int8_t, int16_t
#from pysam.libchtslib cimport bam1_t, bam_pileup1_t, htsFile, hts_itr_t, hts_idx_t
#from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion
cimport pysam.libcalignedsegment as pysam_aligned
cimport pysam.libchtslib as htslib
from pysam.libcalignedsegment cimport AlignedSegment
from pysam import AlignedSegment as PyAlignedSegment


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

	#ctypedef struct sam_hrecs_t sam_hrecs_t

	ctypedef struct sam_hdr_t:
		int32_t n_targets
		int32_t ignore_sam_err
		size_t l_text
		uint32_t *target_len
		const int8_t *cigar_tab
		char **target_name
		char *text
		void *sdict
		#sam_hrecs_t *hrecs
		uint32_t ref_count

	ctypedef sam_hdr_t bam_hdr_t


cdef extern from "htslib/faidx.h":
	ctypedef struct faidx_t:
		pass

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
		int strand
		bam1_t *read
		bam_hdr_t *header

	ctypedef struct methyl_read_iterator_t:
		char *chrom
		faidx_t *fai
		methyl_read_t *methyl_pair
		methyl_read_t *read1
		methyl_read_t *read2

	ctypedef struct reference_cpgs_t:
		long *pos
		int *strand
		int size
		int max_size
	
	#==================================================================================================
	# read_intervals.c
	#--------------------------------------------------------------------------------------------------

	# Check that read passed QC (static inline in read_intervals.h)
	# reject_mask is precomputed via build_reject_mask(qcfail) — avoids per-read recompute
	int check_read(bam1_t *aln, uint32_t reject_mask, int min_size, int max_size,
	               int paired, int mapq_cutoff, float proportion) nogil

	# Add reads from sam file to interval list
	void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
					int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
					float proportion, int nthreads, int add_chr) nogil

	# Add reads from one chromosome using BAM index (parallel-friendly, index-based)
	void sam_iter_add_region(char *samfile_name, labeled_aiarray_t *intervals,
					const char *chromosome,
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
														const char *ref_file,
														const char *chromosome,
														int start,
														int end,
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

cdef AlignedSegment bam1_to_aligned_segment(bam1_t *bam_record, bam_hdr_t *header)

cdef class MethylFragment(object):

	# Attributes
	cdef public str name
	cdef public str chrom
	cdef public int start
	cdef public int end
	cdef public int length
	cdef public AlignedSegment read1
	cdef public AlignedSegment read2

	cdef public int read1_start
	cdef public int read1_end
	cdef public int read1_strand
	cdef public int read2_start
	cdef public int read2_end
	cdef public int read2_strand
	cdef public np.ndarray pos
	cdef public np.ndarray methyl
	cdef public np.ndarray read1_pos
	cdef public np.ndarray read1_methyl
	cdef public np.ndarray read2_pos
	cdef public np.ndarray read2_methyl

	# Methods
	cdef void set_frag(MethylFragment self, methyl_read_t *frag, methyl_read_t *read1, methyl_read_t *read2)
