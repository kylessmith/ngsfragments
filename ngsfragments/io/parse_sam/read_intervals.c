//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


int check_read(bam1_t *aln, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff, float proportion)
{   /* Check that read passed QC */
    
    // Initialize variables
    uint32_t flag;
    int is_proper_pair, is_duplicate;
    int is_qcfail, mate_is_unmapped;
    int is_unmapped;
    int mapq, tlen, start;
    float r;
    float max = (float)RAND_MAX;

    // Check mapping quality
    if (aln->core.qual < mapq_cutoff)
    {
        return 0;
    }

    // Set flag
    flag = aln->core.flag;

    // Check qcfail
    //is_qcfail = (int)(flag & BAM_FQCFAIL);
    if (aln->core.flag & BAM_FQCFAIL)
    {
        return 0;
    }

    // Check mapping status
    //is_unmapped = (int)(flag & BAM_FUNMAP);
    if (aln->core.flag & BAM_FUNMAP)
    {
        return 0;
    }

    // Check duplicate status
    //is_duplicate = (int)(flag & BAM_FDUP);
    if (aln->core.flag & BAM_FDUP)
    {
        return 0;
    }

    // Check paired status
    if (paired == 1)
    {
        //is_proper_pair = (int)(flag & BAM_FPROPER_PAIR);
        //if (aln->core.flag & BAM_FPROPER_PAIR)
        //{
        //    return 0;
        //}
        //mate_is_unmapped = (int)(flag & BAM_FMUNMAP);
        if (aln->core.flag & BAM_FMUNMAP)
        {
            return 0;
        }
    }

    // Record insert/fragment length
    start = aln->core.pos;
    if (paired == 1)
    {
        tlen = aln->core.isize; // insert size
    }
    else {
        tlen = aln->core.l_qseq; // length of read
    }
    
    // Insert read into interval list
    if (aln->core.flag & BAM_FPROPER_PAIR)
    {
        if (tlen >= min_size && tlen <= max_size)
        {
            // Randomly downsample
            if (proportion < 1.0)
            {
                //r = rand() / RAND_MAX;
                float random = (float)rand();
                r = (random / max);
                if (r < proportion)
                {
                    return 1;
                }
            } else {
                return 1;
            }
        }
    }

    return 0;
}


void sam_iter_add(char *samfile_name,
                  labeled_aiarray_t *intervals,
                  int min_size,
                  int max_size,
                  int paired,
                  int qcfail,
                  int mapq_cutoff,
				  float proportion,
                  int nthreads,
                  int add_chr)
{   /* Add reads from sam file to interval list */
    
    // Open sam files
    samFile *fp_in = hts_open(samfile_name, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();
	int start;
	int tlen;
	int end;
	int passing;

    // Set number of threads
    if (nthreads > 1)
    {
        hts_set_threads(fp_in, nthreads);
    }

    // Iterate over bam reads
    while (sam_read1(fp_in, sam_hdr, aln) >= 0)
    {
        // Check read

        //printf("chrom: %s, chrom0: %s\n", chrom, chrom0);

        passing = check_read(aln, min_size, max_size, paired, qcfail, mapq_cutoff, proportion);

        //printf("passing: %d\n", passing);
	    
		// Add read
		if (passing == 1){
			start = aln->core.pos;
		    if (paired == 1)
		    {
		        tlen = aln->core.isize; // insert size
		    }
		    else {
		        tlen = aln->core.l_qseq; // length of read
		    }

            const char *chrom = sam_hdr->target_name[aln->core.tid];
            if (add_chr == 0)
            {
			    labeled_aiarray_add(intervals, start, start+tlen, chrom);
            }
            else {
                //char *tmp_chrom = (char*)chrom;
                char new_chrom[100] = "chr";
                strcat(new_chrom, chrom);
                labeled_aiarray_add(intervals, start, start+tlen, new_chrom);
            }
		}
    }

    // Clean up
    //printf("cleaning up\n");
    bam_destroy1(aln);
    sam_close(fp_in);
    sam_hdr_destroy(sam_hdr);

    return;
}


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
                        int add_chr)
{   /* Add reads from sam file to inteval list and adjust for nucleosome occupancy */
    
    // Open sam files
    samFile *fp_in = hts_open(samfile_name, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();
	int start;
	int tlen;
	int end;
	int passing;

    // Set number of threads
    if (nthreads > 1)
    {
        hts_set_threads(fp_in, nthreads);
    }

    // Iterate over bam reads
    while (sam_read1(fp_in, sam_hdr, aln) >= 0)
    {
        // Check read

        //printf("chrom: %s, chrom0: %s\n", chrom, chrom0);

        passing = check_read(aln, min_size, max_size, paired, qcfail, mapq_cutoff, proportion);

        //printf("passing: %d\n", passing);
	    
		// Add read
		if (passing == 1){
			start = aln->core.pos;
		    if (paired == 1)
		    {
		        tlen = aln->core.isize; // insert size
		    }
		    else {
		        tlen = aln->core.l_qseq; // length of read
		    }

            const char *chrom = sam_hdr->target_name[aln->core.tid];
            //int strand = bam_is_rev(b);
            if (add_chr == 0)
            {
                if (bam_is_rev(aln))
                {
                    int midpoint = start;
                    //midpoint = midpoint - (uint32_t)(tlen / 2);
                    int fixed_start = midpoint - (fixed_size / 2);
                    int fixed_end = fixed_start + fixed_size;
                    labeled_aiarray_add(intervals, fixed_start, fixed_end, chrom);
                }
                else
                {
                    int midpoint = start + tlen;
                    //midpoint = midpoint + (uint32_t)(tlen / 2);
                    int fixed_start = midpoint - (fixed_size / 2);
                    int fixed_end = fixed_start + fixed_size;
                    labeled_aiarray_add(intervals, fixed_start, fixed_end, chrom);
                }
            }
            else {
                //char *tmp_chrom = (char*)chrom;
                char new_chrom[100] = "chr";
                strcat(new_chrom, chrom);
                if (bam_is_rev(aln))
                {
                    int midpoint = start;
                    //midpoint = midpoint + (uint32_t)(tlen / 2);
                    int fixed_start = midpoint - (fixed_size / 2);
                    int fixed_end = fixed_start + fixed_size;
                    labeled_aiarray_add(intervals, fixed_start, fixed_end, new_chrom);
                }
                else
                {
                    int midpoint = start + tlen;
                    //midpoint = midpoint - (uint32_t)(tlen / 2);
                    int fixed_start = midpoint - (fixed_size / 2);
                    int fixed_end = fixed_start + fixed_size;
                    labeled_aiarray_add(intervals, fixed_start, fixed_end, new_chrom);
                }
            }
		}
    }

    // Clean up
    //printf("cleaning up\n");
    bam_destroy1(aln);
    sam_close(fp_in);
    sam_hdr_destroy(sam_hdr);

    return;
}