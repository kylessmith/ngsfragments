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
#include "kmers/interval_kmer.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"

//=====================================================================================
// Define macros
//-------------------------------------------------------------------------------------

// Append to mathylation read
#define push_methyl(read, pos, methyl, qual) do {											\
    if (read->size == read->max_size){													    \
		read->max_size = read->max_size? read->max_size<<1 : 2;							    \
		read->pos = (long*)realloc(read->pos, sizeof(long) * read->max_size);		        \
        read->methyl = (int8_t*)realloc(read->methyl, sizeof(int8_t) * read->max_size);	    \
        read->qual = (uint8_t*)realloc(read->qual, sizeof(uint8_t) * read->max_size);	    \
	}																				        \
	if (methyl == 1) {																        \
        read->pos[read->size] = pos;												        \
        read->methyl[read->size] = 1;												        \
        read->qual[read->size++] = qual;												    \
    } else {																		        \
        read->pos[read->size] = pos;												        \
        read->methyl[read->size] = 0;												        \
        read->qual[read->size++] = qual;												    \
    }																				        \
} while (0)

// Append to cpg record
#define push_cpg(cpgs, pos, strand) do {											        \
    if (cpgs->size == cpgs->max_size){													    \
		cpgs->max_size = cpgs->max_size? cpgs->max_size<<1 : 2;							    \
		cpgs->pos = (long*)realloc(cpgs->pos, sizeof(long) * cpgs->max_size);		        \
        cpgs->strand = (int*)realloc(cpgs->strand, sizeof(int) * cpgs->max_size);	        \
	}																				        \
    cpgs->pos[cpgs->size] = pos;												            \
    cpgs->strand[cpgs->size++] = strand;												    \
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
    int qcfail;
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
} methyl_read_t;

typedef struct {
    read_iter_t *read_iter;
    void *methyl_hash;
    TwoBit *tb;
    const char *chrom;
    int ref_start;
    int ref_end;
    char *ref_seq;
    int seq_len;
    int chrom_length;
    methyl_read_t *methyl_pair;
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

// Check that read passed QC
int check_read(bam1_t *aln, int min_size, int max_size, int paired, int qcfail,
		   		int mapq_cutoff, float proportion);

// Add reads from sam file to interval list
void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
                  int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
				  float proportion, int nthreads, int add_chr);

// Add reads from sam file to inteval list and adjust for nucleosome occupancy
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

// Initialize read iterator
read_iter_t *read_iter_init(const char *bam_file_path,
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
methyl_read_t *methyl_read_init(bam1_t *aln);

// Free memory allocated for methyl_read_t struct
void methyl_read_destroy(methyl_read_t *read);

//Check if position is CpG
inline int isCpG(char *seq, int pos, int seqlen);

// Get strand of read
int getStrand(bam1_t *b);

// Process read
void processRead(bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, methyl_read_t *read);

// Process methylated read pair
methyl_read_t *methyl_pair_process(methyl_read_t *read1, methyl_read_t *read2);

// Initialize methyl_read_iterator_t struct
methyl_read_iterator_t *methyl_read_iterator_init(const char *bam_file_path,
                                                    char *ref_2bit,
                                                    const char *chromosome,
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
                        int nthreads);


//==================================================================================================
// bounds_motif.c
//--------------------------------------------------------------------------------------------------

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
                                            int add_chr);

//==================================================================================================
// merge_bams.c
//--------------------------------------------------------------------------------------------------

//void merge_bams(const char *input1_bam_path, const char *input2_bam_path, const char *output_bam_path);

//==================================================================================================
// .c
//--------------------------------------------------------------------------------------------------


//==================================================================================================

#endif
