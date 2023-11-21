//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


reference_cpgs_t *reference_cpgs_init()
{   /* Initialize reference CpGs struct */

    // Initialize struct
    reference_cpgs_t *cpgs = malloc(sizeof(reference_cpgs_t));

    // Initialize arrays for storing CpGs
    cpgs->pos = (long *)malloc(sizeof(long) * 64);
    cpgs->strand = (int *)malloc(sizeof(int) * 64);
    cpgs->size = 0;
    cpgs->max_size = 64;

    return cpgs;
}


void reference_cpgs_destroy(reference_cpgs_t *cpgs)
{   /* Free memory allocated for reference CpGs struct */

    free(cpgs->pos);
    free(cpgs->strand);
    free(cpgs);
}


void reference_cpgs_append(reference_cpgs_t *cpgs, long pos, int strand)
{   /* Append CpG to reference CpGs struct */

    push_cpg(cpgs, pos, strand);
}


reference_cpgs_t *fetch_reference_cpgs(char *ref_2bit, char *chrom)
{   /* Fetch reference CpGs */

    // Initialize
    reference_cpgs_t *cpgs = reference_cpgs_init();
    TwoBit *ref_file = twobitOpen(ref_2bit, 0);
    int chrom_length = twobitChromLen(ref_file, chrom);
    int position = 0;
    int start = 0;
    int end = 1000;
    char *sequence = twobitSequence(ref_file, chrom, start, end + 1);
    int seq_len = strlen(sequence);

    // Check if chromosome is in the 2bit file
    if (chrom_in(chrom, ref_file->cl->chrom, ref_file->hdr->nChroms) == 0)
    {
        return 0;
    }

    
    // Iterate over chromosome
    int i = 0;
    while (position < chrom_length - 1)
    {
        // Increment position
        position++;

        // Check if current position is in sequence
        if (position >= end)
        {
            // Check sequence bounds
            start = end;
            if ((end + 1000) > chrom_length - 1)
            {
                end = chrom_length - 1;
            } else {
                end = end + 1000;
            }

            // Fetch sequence
            if ((end - start) > 0)
            {
                free(sequence);
                sequence = twobitSequence(ref_file, chrom, start, end + 1);
                seq_len = strlen(sequence);
            }

            // Reset position
            i = 0;
        }

        // Check if CpG
        int direction = isCpG(sequence, i, seq_len);
        if (direction != 0)
        {
            //printf("%s\t%d\n", chrom, position);
            reference_cpgs_append(cpgs, position, direction);
        }

        i++;
    }

    twobitClose(ref_file);
    free(sequence);

    return cpgs;
}


methyl_record_t *fetch_reference_methyl_record(char *ref_2bit, char *chrom)
{   /* Fetch reference methyl record */

    // Fetch CpGs
    //printf("  Fetching reference CpGs...\n");
    reference_cpgs_t *cpgs = fetch_reference_cpgs(ref_2bit, chrom);
    //printf("  done\n");

    // Initialize methyl record
    //printf("  Initializing methyl record...\n");
    //printf("   CpGs: %d\n", cpgs->size);
    methyl_record_t *series = methyl_record_init(cpgs->pos, cpgs->size);
    //printf("  done\n");
    free(cpgs->strand);

    return series;
}