//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include "read_intervals.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "src/labeled_aiarray/labeled_augmented_array.h"


//-----------------------------------------------------------------------------


int check_read(bam1_t *aln, int min_size, int max_size, int paired, int qcfail, int mapq_cutoff, float proportion)
{   /* Filter read and add to augmented interval list */
    
    // Initialize variables
    uint32_t flag;
    int is_proper_pair, is_duplicate;
    int is_qcfail, mate_is_unmapped;
    int is_unmapped;
    int mapq, tlen, start;
    float r;

    // Check mapping quality
    if (aln->core.qual < mapq_cutoff)
    {
        return 0;
    }

    // Set flag
    flag = aln->core.flag;

    // Check qcfail
    is_qcfail = (int)(flag & BAM_FQCFAIL);
    if (is_qcfail != qcfail)
    {
        return 0;
    }

    // Check mapping status
    is_unmapped = (int)(flag & BAM_FUNMAP);
    if (is_unmapped == 1)
    {
        return 0;
    }

    // Check duplicate status
    is_duplicate = (int)(flag & BAM_FDUP);
    if (is_duplicate == 1)
    {
        return 0;
    }

    // Check paired status
    if (paired == 1)
    {
        is_proper_pair = (int)(flag & BAM_FPROPER_PAIR);
        if (is_proper_pair == 0)
        {
            return 0;
        }
        mate_is_unmapped = (int)(flag & BAM_FMUNMAP);
        if (mate_is_unmapped == 1)
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
    if (tlen >= min_size && tlen <= max_size)
    {
        // Randomly downsample
        if (proportion < 1.0)
        {
            r = rand() / RAND_MAX;
            if (r < proportion)
            {
                return 1;
				//labeled_aiarray_add(intervals, start, start+tlen, chrom_id);
            }
        } else {
            return 1;
        }
            
    }

    return 0;
}


void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
                            int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
							float proportion, int nthreads, int trim_name)
{
    
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

            char chrom0[100];
            const char *chrom;
            const char *sam_chrom = sam_hdr->target_name[aln->core.tid];
            strcpy(chrom0, sam_chrom);
            
            // Trim chrom name by 2 characters if requested
            //printf("chrom0: %s\n", chrom0);
            if (trim_name == 1)
            {
                int len = strlen(chrom0);
                //printf("len: %d\n", len);
                if (len > 2)
                {
                    chrom0[len-2] = '\0';
                    //strcpy(chrom, chrom0);
                    chrom = chrom0;
                }
            } else {
                //strcpy(chrom, chrom0);
                chrom = chrom0;
            }

            //printf("adding\n");
			labeled_aiarray_add(intervals, start, start+tlen, chrom);
            //printf("added\n");
		}
    }

    // Clean up
    //printf("cleaning up\n");
    bam_destroy1(aln);
    sam_close(fp_in);
    sam_hdr_destroy(sam_hdr);

    return;
}


int chrom_iter_add(char *samfile_name, const char *chrom_pos, labeled_aiarray_t *intervals,
                            int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
							float proportion)
{
    // Open Sam file
    bam1_t *aln = bam_init1();
    samFile *fp_in = sam_open(samfile_name, "r");
    if (fp_in == NULL)
    {
        return -1;
    }
    bam_hdr_t *header = sam_hdr_read(fp_in);
    if (header == NULL)
    {
        return -1;
    }

    // initiate an iterator
    hts_itr_t *iter;
    // initiate an BAM index
    hts_idx_t *idx;

    // second parameter is same as BAM file path
    idx = sam_index_load(fp_in, samfile_name);
    if (idx == NULL)
    {
        return -1;
    }

    //  third parameter is like chr1:1000-2000, locate to the start of the region
    iter  = sam_itr_querys(idx, header, chrom_pos);

    // Iterate over sam
    int start;
	int tlen;
	int end;
	int passing;
    while ( sam_itr_next(fp_in, iter, aln) >= 0)
    {
        // Check read
		char *chrom = header->target_name[aln->core.tid];
        passing = check_read(aln, min_size, max_size, paired, qcfail, mapq_cutoff, proportion);
	    
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
			labeled_aiarray_add(intervals, start, start+tlen, chrom);
		}
    }

    // Clean up
    bam_destroy1(aln);
    sam_close(fp_in);
    hts_itr_destroy(iter);
    bam_hdr_destroy(header);

    return 0;
}


int chrom_bin_iter_add(char *samfile_name, const char *chrom_pos, int bin_size,
                        long nhits[], int min_size, int max_size, int paired,
                        int qcfail, int mapq_cutoff, float proportion)
{
    // Open Sam file
    bam1_t *aln = bam_init1();
    samFile *fp_in = sam_open(samfile_name, "r");
    if (fp_in == NULL)
    {
        return -1;
    }
    bam_hdr_t *header = sam_hdr_read(fp_in);
    if (header == NULL)
    {
        return -1;
    }

    // initiate an iterator
    hts_itr_t *iter;
    // initiate an BAM index
    hts_idx_t *idx;

    // second parameter is same as BAM file path
    idx = sam_index_load(fp_in, samfile_name);
    if (idx == NULL)
    {
        return -1;
    }

    //  third parameter is like chr1:1000-2000, locate to the start of the region
    iter  = sam_itr_querys(idx, header, chrom_pos);

    //printf("here1\n");

    // Iterate over sam
    int start;
	int tlen;
	int end;
	int passing;
    while ( sam_itr_next(fp_in, iter, aln) >= 0)
    {
        // Check read
		char *chrom = header->target_name[aln->core.tid];
        passing = check_read(aln, min_size, max_size, paired, qcfail, mapq_cutoff, proportion);
	    
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

            int start_bin = start / bin_size;
            double length = (double)(tlen);
            int n_bins = ceil(((double)(start % bin_size) / bin_size) + (length / bin_size));
            int n;
            for (n = 0; n < n_bins; n++)
            {
                int bin = start_bin + n;
                nhits[bin] = nhits[bin] + 1;
            }
		}
    }

    // Clean up
    bam_destroy1(aln);
    sam_close(fp_in);
    hts_itr_destroy(iter);
    bam_hdr_destroy(header);

    return 0;
}