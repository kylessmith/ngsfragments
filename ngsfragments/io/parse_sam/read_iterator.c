//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


read_iter_t *read_iter_init(const char *bam_file_path,
                            const char *chromosome,
                            int min_size,
                            int max_size,
                            int paired,
                            int qcfail,
                            int mapq_cutoff,
                            float proportion,
                            int nthreads)
{   /* Initialize read iterator */

    // Initialize variables
    read_iter_t *read_iter = malloc(sizeof(read_iter_t));
    read_iter->fp = sam_open(bam_file_path, "r");
    if (read_iter->fp == NULL)
    {
        fprintf(stderr, "Error opening BAM file: %s\n", bam_file_path);
        exit(1);
    }
    // Set number of threads
    if (nthreads > 1)
    {
        hts_set_threads(read_iter->fp, nthreads);
    }

    // Read header
    read_iter->header = sam_hdr_read(read_iter->fp);
    if (read_iter->header == NULL)
    {
        fprintf(stderr, "Failed to read header from BAM file: %s\n", bam_file_path);
        sam_close(read_iter->fp);
        exit(1);
    }

    // Load the index
    hts_idx_t *idx = sam_index_load(read_iter->fp, bam_file_path);
    if (idx == NULL)
    {
        fprintf(stderr, "Failed to load index for BAM file: %s\n", bam_file_path);
        bam_hdr_destroy(read_iter->header);
        sam_close(read_iter->fp);
        exit(1);
    }

    // Fetch the region
    read_iter->iter = sam_itr_querys(idx, read_iter->header, chromosome);
    if (read_iter->iter == NULL)
    {
        fprintf(stderr, "Failed to fetch region: %s\n", chromosome);
        hts_idx_destroy(idx);
        bam_hdr_destroy(read_iter->header);
        sam_close(read_iter->fp);
        exit(1);
    }

    // Close index
    hts_idx_destroy(idx);

    read_iter->aln = bam_init1();

    // Set variables
    read_iter->min_size = min_size;
    read_iter->max_size = max_size;
    read_iter->paired = paired;
    read_iter->qcfail = qcfail;
    read_iter->mapq_cutoff = mapq_cutoff;
    read_iter->proportion = proportion;

    return read_iter;
}


void read_iter_destroy(read_iter_t *read_iter)
{   /* Free memory allocated for read iterator */

    // Close BAM file
    bam_destroy1(read_iter->aln);
    bam_hdr_destroy(read_iter->header);
    hts_itr_destroy(read_iter->iter);
    sam_close(read_iter->fp);
    free(read_iter);
}


int read_iter_next(read_iter_t *read_iter)
{   /* Iterate over BAM file */

    // Iterate over BAM file
    while (sam_itr_next(read_iter->fp, read_iter->iter, read_iter->aln) >= 0)
    {
        // Check read
        int passing = check_read(read_iter->aln,
                                 read_iter->min_size,
                                 read_iter->max_size,
                                 read_iter->paired,
                                 read_iter->qcfail,
                                 read_iter->mapq_cutoff,
                                 read_iter->proportion);
        // Check if passing, else continue iterating
        if (passing == 0)
        {
            continue;
        } else {
            return 1;
        }
    }

    return 0;
}