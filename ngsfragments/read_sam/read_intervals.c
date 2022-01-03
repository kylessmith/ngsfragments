//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include "read_intervals.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "src/labeled_augmented_array.h"


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
							float proportion)
{
    
    // Open sam files
    samFile *fp_in = hts_open(samfile_name, "r");
    bam_hdr_t *bam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();
	int start;
	int tlen;
	int end;
	int passing;

    // Iterate over bam reads
    while (sam_read1(fp_in, bam_hdr, aln) > 0)
    {
        // Check read
		char *chrom = bam_hdr->target_name[aln->core.tid];
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

    return;
}
