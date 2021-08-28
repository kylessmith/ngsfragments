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


void read_add(labeled_aiarray_t *intervals, bam1_t *aln, char *chrom_id,
                      int min_size, int max_size, int paired, int qcfail, int mapq_cutoff, float proportion)
{   /* Filter read and add to augmented interval list */
    
    // Initialize variables
    uint32_t flag;
    int is_proper_pair, is_duplicate;
    int is_qcfail, mate_is_unmapped;
    int is_unmapped;
    int mapq, tlen, start;
    float r;
    int passing = 1;

    // Check mapping quality
    if (aln->core.qual < mapq_cutoff)
    {
        passing = 0;
    }

    // Set flag
    flag = aln->core.flag;

    // Check qcfail
    is_qcfail = (int)(flag & BAM_FQCFAIL);
    if (passing == 1 && is_qcfail != qcfail)
    {
        passing = 0;
    }

    // Check mapping status
    is_unmapped = (int)(flag & BAM_FUNMAP);
    if (passing == 1 && is_unmapped == 1)
    {
        passing = 0;
    }

    // Check duplicate status
    is_duplicate = (int)(flag & BAM_FDUP);
    if (passing == 1 && is_duplicate == 1)
    {
        passing = 0;
    }

    // Check paired status
    if (passing == 1 && paired == 1)
    {
        is_proper_pair = (int)(flag & BAM_FPROPER_PAIR);
        if (is_proper_pair == 0)
        {
            passing = 0;
        }
        mate_is_unmapped = (int)(flag & BAM_FMUNMAP);
        if (mate_is_unmapped == 1)
        {
            passing = 0;
        }
    }

    // Record insert/fragment length
    if (passing == 1)
    {
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
                    labeled_aiarray_add(intervals, start, start+tlen, chrom_id);
                }
            } else {
                labeled_aiarray_add(intervals, start, start+tlen, chrom_id);
            }
                
        }
    }

    return;
}


void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
                            int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
							float proportion)
{
    
    // Open sam files
    samFile *fp_in = hts_open(samfile_name, "r");
    bam_hdr_t *bam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();

    // Iterate over bam reads
    while (sam_read1(fp_in, bam_hdr, aln) > 0)
    {
        char *chrom = bam_hdr->target_name[aln->core.tid];
        read_add(intervals, aln, chrom, min_size, max_size, paired, qcfail, mapq_cutoff, proportion);
    }

    // Clean up
    bam_destroy1(aln);
    sam_close(fp_in);

    return;
}