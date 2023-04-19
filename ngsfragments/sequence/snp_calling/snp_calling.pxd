from pysam.libchtslib cimport bam1_t, bam_pileup1_t, bam_mplp_t
from pysam.libcalignedsegment cimport PileupColumn
from pysam.libcalignmentfile cimport IteratorColumn
from pysam.libcsamfile cimport Samfile, \
    pysam_bam_get_seq, pysam_bam_get_qual
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from cpython cimport PyBytes_FromStringAndSize

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t
from ailist.AIList_core cimport ailist_t



cdef struct _CountAllele:
	int fwd, rev, pres, count, mapqual
	double qual


cdef inline object _get_seq_base(bam1_t *src, uint32_t k)
cdef inline _init_count(_CountAllele* c)
cdef inline _incr_count(_CountAllele* c, bint is_reverse, double baseq, int mapq)
cdef int process_chrom(bam_mplp_t bam_record, ailist_t *ail) nogil