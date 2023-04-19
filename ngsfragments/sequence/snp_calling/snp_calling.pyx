#cython: language_level=3

from cython.parallel import prange

from pysam.libchtslib cimport bam1_t, bam_pileup1_t, bam_mplp_t, bam_mplp_auto
from pysam.libcalignedsegment cimport PileupColumn
from pysam.libcalignmentfile cimport IteratorColumn
from pysam.libcsamfile cimport Samfile, \
    pysam_bam_get_seq, pysam_bam_get_qual
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from cpython cimport PyBytes_FromStringAndSize
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, label_t
from ailist.AIList_core cimport ailist_t

## These are bits set in the flag.
## have to put these definitions here, in csamtools.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
DEF BAM_FPAIRED       =1
## @abstract the read is mapped in a proper pair */
DEF BAM_FPROPER_PAIR  =2
## @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
DEF BAM_FUNMAP        =4
## @abstract the mate is unmapped */
DEF BAM_FMUNMAP       =8
## @abstract the read is mapped to the reverse strand */
DEF BAM_FREVERSE      =16
## @abstract the mate is mapped to the reverse strand */
DEF BAM_FMREVERSE     =32
## @abstract this is read1 */
DEF BAM_FREAD1        =64
## @abstract this is read2 */
DEF BAM_FREAD2       =128
## @abstract not primary alignment */
DEF BAM_FSECONDARY   =256
## @abstract QC failure */
DEF BAM_FQCFAIL      =512
## @abstract optical or PCR duplicate */
DEF BAM_FDUP        =1024


cdef char* CODE2CIGAR= "MIDNSHP=X"

cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"


cdef inline object _get_seq_base(bam1_t *src, uint32_t k):
    cdef uint8_t * p
    cdef char * s

    if not src.core.l_qseq:
        return None

    seq = PyBytes_FromStringAndSize(NULL, 1)
    s   = <char*>seq
    p   = pysam_bam_get_seq(src)

    s[0] = bam_nt16_rev_table[p[k//2] >> 4 * (1 - k%2) & 0xf]

    return seq


cdef inline _init_count(_CountAllele* c):
    c.qual = c.fwd = c.rev = c.pres = c.count = c.mapqual = 0


cdef inline _incr_count(_CountAllele* c, bint is_reverse, double baseq, int mapq):
    c.count+=1
    c.pres=1
    c.qual+=baseq
    c.mapqual+=mapq
    if is_reverse:
        c.rev += 1
    else:
        c.fwd += 1


cdef int process_chrom(bam_mplp_t bam_record, ailist_t *ail) nogil:

    cdef int tid = 0
    cdef int pos = 0
    cdef int n_plp = 0
    cdef const bam_pileup1_t * plp = NULL
    cdef int ret = bam_mplp_auto(bam_record,
                                 &tid,
                                 &pos,
                                 &n_plp,
                                 &plp)
    printf("ret: %d\n", ret)

    return 0


def simple_seq( LabeledIntervalArray intervals,
                str samfile_name,
                int min_reads=10,
                int min_support=4, 
                double min_af=0.05,
                double min_pvalue=0.0001,
                double min_f_r=0.05,
                int min_qual=25,
                double min_se=0.001,
                uint64_t min_mapq=55, 
                region_chrom=None,
                region_start=None,
                region_end=None,
                int min_baseq=13):

    #statically type variables
    cdef int current_pos = 0
    cdef set seen_chroms = set()
    cdef tuple germline_pos
    cdef list germline_observed = []
    cdef list possible_snps = []
    cdef list coverage = []
    cdef list total_positions = []
    cdef bint in_germline
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef uint32_t flag
    cdef bint is_reverse, is_proper_pair, germline
    cdef int i  # loop index
    cdef int n  # total number of reads in column
    cdef double baseq
    cdef double total_baseq = 0.0
    cdef PileupColumn col
    cdef double alt_allele_freq
    cdef str chrom
    cdef bytes alnbase
    cdef int pos
    cdef _CountAllele A, T, G, C
    cdef int total_alleles
    cdef double reverse, forward

    cdef Samfile samfile
    cdef IteratorColumn it
    cdef int n_chromosomes = len(intervals.unique_labels)
    cdef bam_mplp_t *bam_record = <bam_mplp_t *>malloc(n_chromosomes * sizeof(bam_mplp_t))
    i = 0
    for chrom in intervals.unique_labels:
        samfile = Samfile(samfile_name)
        it = samfile.pileup(reference=chrom,
                            start=None,
                            end=None,
                            truncate=False,
                            max_depth=80000,
                            ignore_orphans=True,
                            min_base_quality=min_qual,
                            min_mapping_quality=min_mapq,
                            ignore_overlaps=True,
                            redo_baq=True)
        bam_record[i] = it.pileup_iter
        i += 1

    cdef labeled_aiarray_t *laia = intervals.laia
    cdef label_t *p
    cdef Py_ssize_t t
    for t in prange(n_chromosomes, nogil=True):
        p = &laia.labels[t]
        process_chrom(bam_record[t], p.ail)

    return