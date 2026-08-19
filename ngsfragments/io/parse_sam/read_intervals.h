//=====================================================================================
// Common structs, parameters, functions
// Original: https://github.com/databio/aiarray/tree/master/src
// by Kyle S. Smith
//-------------------------------------------------------------------------------------
#ifndef __READ_INTERVALS_H__
#define __READ_INTERVALS_H__
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "kmers/interval_kmer.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"

//=====================================================================================
// Thread-local xorshift32 RNG — faster than rand(), no global lock, no float division
//-------------------------------------------------------------------------------------

static __thread uint32_t _rng_state = 0;

static inline float fast_rand_float(void)
{
    /* Self-seeds from wall-clock time XOR stack address on first use per thread. */
    if (_rng_state == 0)
        _rng_state = (uint32_t)time(NULL) ^ (uint32_t)(uintptr_t)&_rng_state;
    _rng_state ^= _rng_state << 13;
    _rng_state ^= _rng_state >> 17;
    _rng_state ^= _rng_state << 5;
    /* Multiply by 1/2^32 instead of dividing — single instruction on modern CPUs */
    return _rng_state * (1.0f / 4294967296.0f);
}

//=====================================================================================
// Flag rejection mask helpers
// BAM_FSECONDARY | BAM_FSUPPLEMENTARY were previously unfiltered — chimeric/multi-map
// reads with garbage isize values could corrupt fragment intervals.
//-------------------------------------------------------------------------------------

#define BASE_REJECT_MASK (BAM_FUNMAP | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)

static inline uint32_t build_reject_mask(int qcfail)
{
    uint32_t mask = BASE_REJECT_MASK;
    if (!qcfail) mask |= BAM_FQCFAIL;  /* exclude QC-fails unless caller wants them */
    return mask;
}

//=====================================================================================
// Define macros
//-------------------------------------------------------------------------------------

// Append to methylation read
/* NOTE: the write at the end is GUARDED by (read->size < read->max_size). On a
 * realloc failure the growth block `break`s out of the do{}while(0); previously
 * execution then fell through to the unconditional write, which (a) wrote one
 * element past a full buffer, and (b) if only the first of the three reallocs had
 * succeeded, indexed three arrays of MISMATCHED capacity — a heap overflow. The
 * guard makes a failed grow drop the element instead of corrupting the heap.
 * read->size is only incremented when the element is actually stored, so it stays
 * an accurate count (see methyl_read_append, which mirrors ncpgs to size). */
#define push_methyl(read, pos, methyl, qual) do {                                           \
    if (read->size == read->max_size) {                                                     \
        int _new_max = read->max_size ? read->max_size << 1 : 2;                            \
        long *new_pos = (long*)realloc(read->pos, sizeof(long) * _new_max);                 \
        if (new_pos == NULL) break;                                                         \
        read->pos = new_pos;                                                                \
        int8_t *new_methyl = (int8_t*)realloc(read->methyl, sizeof(int8_t) * _new_max);    \
        if (new_methyl == NULL) break;                                                      \
        read->methyl = new_methyl;                                                          \
        uint8_t *new_qual = (uint8_t*)realloc(read->qual, sizeof(uint8_t) * _new_max);     \
        if (new_qual == NULL) break;                                                        \
        read->qual = new_qual;                                                              \
        read->max_size = _new_max; /* only commit after all three reallocs succeed */       \
    }                                                                                       \
    if (read->size < read->max_size) { /* only write if capacity actually exists */         \
        read->pos[read->size] = pos;                                                        \
        read->methyl[read->size] = (methyl == 1) ? 1 : 0;                                  \
        read->qual[read->size++] = qual;                                                    \
    }                                                                                       \
} while (0)

// Append to cpg record
#define push_cpg(cpgs, pos, strand) do {                                                    \
    if (cpgs->size == cpgs->max_size) {                                                     \
        int _new_max = cpgs->max_size ? cpgs->max_size << 1 : 2;                            \
        long *new_pos = (long*)realloc(cpgs->pos, sizeof(long) * _new_max);                 \
        if (new_pos == NULL) break;                                                         \
        cpgs->pos = new_pos;                                                                \
        int *new_strand = (int*)realloc(cpgs->strand, sizeof(int) * _new_max);              \
        if (new_strand == NULL) break;                                                      \
        cpgs->strand = new_strand;                                                          \
        cpgs->max_size = _new_max; /* only commit after both reallocs succeed */            \
    }                                                                                       \
    if (cpgs->size < cpgs->max_size) { /* only write if capacity actually exists */         \
        cpgs->pos[cpgs->size] = pos;                                                        \
        cpgs->strand[cpgs->size++] = strand;                                                \
    }                                                                                       \
} while (0)

KHASH_SET_INIT_STR(read_name_set);

//=====================================================================================
// Define structs and functions
//-------------------------------------------------------------------------------------

typedef struct {
    samFile *fp;
    bam_hdr_t *header;
    hts_itr_t *iter;
    bam1_t *aln;
    int min_size;
    int max_size;
    int paired;
    uint32_t reject_mask;  /* precomputed from qcfail at init time — avoids per-read recompute */
    int mapq_cutoff;
    float proportion;
} read_iter_t;

typedef struct {
    char *name;
    int start;
    int end;
    int length;
    long *pos;
    int8_t *methyl;
    uint8_t *qual;
    int ncpgs;
    int size;
    int max_size;
    int strand;
    bam1_t *read;
    bam_hdr_t *header;
} methyl_read_t;

typedef struct {
    read_iter_t *read_iter;
    void *methyl_hash;
    TwoBit   *tb;          /* non-NULL when using a .2bit reference  */
    faidx_t  *fai;         /* non-NULL when using a FASTA reference  */
    char *chrom;
    int ref_start;
    int ref_end;
    char *ref_seq;
    int seq_len;
    int chrom_length;
    int use_xm;            /* 1 = use XM-tag path (no ref); 0 = reference-based path */
    methyl_read_t *methyl_pair;
    methyl_read_t *read1;
    methyl_read_t *read2;
} methyl_read_iterator_t;

typedef struct {
    long *pos;
    int *strand;
    int size;
    int max_size;
} reference_cpgs_t;

typedef struct {
    void *lookup;
    long *pos;
    int size;
    int max_size;
} int_index_t;

typedef struct {
    int_index_t *index;
    int16_t *methyl;
    int16_t *unmethyl;
} methyl_record_t;

typedef struct {
    methyl_record_t *record1;
    methyl_record_t *record2;
} methyl_record_pair_t;

typedef struct {
    khash_t(read_name_set) *set1;
    khash_t(read_name_set) *set2;
} read_name_sets_t;


//==================================================================================================
// read_intervals.c
//--------------------------------------------------------------------------------------------------

/* check_read is static inline so both read_intervals.c and read_iterator.c get an inlined copy
 * with no function-call overhead on the hot per-read path.
 *
 * reject_mask: precomputed with build_reject_mask(qcfail) once per file/iterator.
 *   Combines BAM_FUNMAP | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY and optionally
 *   BAM_FQCFAIL into a single bitwise test — one branch instead of four.
 */
static inline int check_read(bam1_t *aln,
                              uint32_t reject_mask,
                              int min_size,
                              int max_size,
                              int paired,
                              int mapq_cutoff,
                              float proportion)
{
    if (aln->core.qual < mapq_cutoff)
        return 0;

    /* Single combined flag test — branch-predictor-friendly on clean BAMs */
    if (aln->core.flag & reject_mask)
        return 0;

    int tlen;
    if (paired == 1)
    {
        if (!(aln->core.flag & BAM_FPROPER_PAIR))
            return 0;
        if (aln->core.flag & BAM_FMUNMAP)
            return 0;
        /* Fragment length is |isize| regardless of which mate this is.
         * NOTE: do NOT filter to BAM_FREAD1 here. check_read is shared with the
         * name-paired methyl iterator (read_iter_next), which must see BOTH mates
         * to form a pair. Per-fragment de-duplication (one interval per pair) is
         * the interval builder's job and is done there via the isize>0 test. */
        tlen = (int)llabs(aln->core.isize);
    }
    else
    {
        tlen = aln->core.l_qseq;
    }

    if (tlen < min_size || tlen > max_size)
        return 0;

    /* Downsample: fast_rand_float() uses thread-local xorshift32, no global lock */
    if (proportion < 1.0f && fast_rand_float() >= proportion)
        return 0;

    return 1;
}

// Add reads from sam file to interval list (whole-file sequential scan)
void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
                  int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
                  float proportion, int nthreads, int add_chr);

// Add reads from one chromosome using the BAM index (parallel-friendly)
void sam_iter_add_region(char *samfile_name, labeled_aiarray_t *intervals,
                         const char *chromosome,
                         int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
                         float proportion, int nthreads, int add_chr);

// Add reads from sam file to interval list and adjust for nucleosome occupancy
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
                        int add_chr);


//==================================================================================================
// read_iterator.c
//--------------------------------------------------------------------------------------------------

// Initialize read iterator for specific genomic region
read_iter_t *read_iter_init(const char *bam_file_path,
                            const char *chromosome,
                            int start_pos,
                            int end_pos,
                            int min_size,
                            int max_size,
                            int paired,
                            int qcfail,
                            int mapq_cutoff,
                            float proportion,
                            int nthreads);

// Initialize read iterator for entire chromosome (convenience function)
read_iter_t *read_iter_init_chromosome(const char *bam_file_path,
                                       const char *chromosome,
                                       int min_size,
                                       int max_size,
                                       int paired,
                                       int qcfail,
                                       int mapq_cutoff,
                                       float proportion,
                                       int nthreads);

// Free memory allocated for read iterator
void read_iter_destroy(read_iter_t *read_iter);

// Iterate over BAM file
int read_iter_next(read_iter_t *read_iter);


//==================================================================================================
// methyl_fragment_iter.c
//--------------------------------------------------------------------------------------------------

// Initialize methyl_read_t struct
methyl_read_t *methyl_read_init(bam1_t *aln, bam_hdr_t *header);

// Free memory allocated for methyl_read_t struct
void methyl_read_destroy(methyl_read_t *read);

//Check if position is CpG
static inline int isCpG(char *seq, int pos, int seqlen);

// Get strand of read
int getStrand(bam1_t *b);

// Process read using reference sequence
void processRead(bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, methyl_read_t *read, int minPhred);

// Process read using Bismark XM tag (reference-free fallback); returns 1 on success, 0 if tag absent
static int processReadXM(bam1_t *b, methyl_read_t *read, int minPhred);

// Process methylated read pair
methyl_read_t *methyl_pair_process(methyl_read_t *read1, methyl_read_t *read2);

// Initialize methyl_read_iterator_t struct
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
                                                    int nthreads);

// Free memory allocated for methyl_read_iterator_t struct
void methyl_read_iterator_destroy(methyl_read_iterator_t *iter);

// Iterate over BAM file
int methyl_read_iterator_next(methyl_read_iterator_t *iter);


//==================================================================================================
// reference_methyl.c
//--------------------------------------------------------------------------------------------------

// Initialize reference CpGs struct
reference_cpgs_t *reference_cpgs_init();

// Free memory allocated for reference CpGs struct
void reference_cpgs_destroy(reference_cpgs_t *cpgs);

// Append CpG to reference CpGs struc
void reference_cpgs_append(reference_cpgs_t *cpgs, long pos, int strand);

// Fetch reference CpGs
reference_cpgs_t *fetch_reference_cpgs(char *ref_2bit, char *chrom);

// Fetch reference methyl record
methyl_record_t *fetch_reference_methyl_record(char *ref_2bit, char *chrom);


//==================================================================================================
// methyl_record.c
//--------------------------------------------------------------------------------------------------

// Initialize methyl record
methyl_record_t *methyl_record_init(long *pos, int size);

// Free memory allocated for methyl record
void methyl_record_destroy(methyl_record_t *series);

void methyl_record_pair_transfer_null(methyl_record_pair_t *pair);

// Add methylated positions to methyl record
void methyl_record_add(methyl_record_t *series, long *pos, int8_t *methyl, int size);

// Get methylated positions from methyl record
int16_t *methyl_record_get(methyl_record_t *series, long pos);

// Initialize methyl record pair
methyl_record_pair_t *methyl_record_pair_init(methyl_record_t *record1,
                                              methyl_record_t *record2);

// Free memory allocated for methyl record pair
void methyl_record_pair_destroy(methyl_record_pair_t *pair);

// Assign read to methyl record pair
int assign_methyl_read(methyl_record_pair_t *pair, methyl_read_t *read);

// Compare methyl records using Euclidean distance
double compare_methyl_records(methyl_record_pair_t *pair);

int methyl_record_pair_write(methyl_record_pair_t *pair, char *file_fn);


//==================================================================================================
// read_name_store.c
//--------------------------------------------------------------------------------------------------

// Initialize a read name sets struct
read_name_sets_t *init_read_name_sets();

// Destroy a read name sets struct
void destroy_read_name_sets(read_name_sets_t *read_name_sets);


//==================================================================================================
// methyl_size_split.c
//--------------------------------------------------------------------------------------------------

// Split reads by length and quantify methylation
methyl_record_pair_t *methyl_length_split(const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size1,
                                            int max_size1,
                                            int min_size2,
                                            int max_size2,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads);

// Split reads based on similarity to other methylation profiles
methyl_record_pair_t *methyl_profile_split(methyl_record_pair_t *pair,
                                            const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size,
                                            int max_size,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads);

// Split reads based on similarity to other methylation profiles and return read names
read_name_sets_t *methyl_profile_split_names(methyl_record_pair_t *pair,
                                            const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size,
                                            int max_size,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads);

// Write reads to output BAM files based on read name
void write_split_reads(const char *bam_file_path,
                        const char *output_bam_file_path1,
                        const char *output_bam_file_path2,
                        read_name_sets_t *read_names,
                        char *chromosome,
                        int min_size,
                        int max_size,
                        int qcfail,
                        int mapq_cutoff,
                        int nthreads);

// Split reads by length/profile and write to output BAM files
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
                        int nthreads);



//==================================================================================================
// merge_bams.c
//--------------------------------------------------------------------------------------------------

//void merge_bams(const char *input1_bam_path, const char *input2_bam_path, const char *output_bam_path);

//==================================================================================================

#endif
