#include "read_intervals.h"

//-----------------------------------------------------------------------------

static const int khStrMeth = 56;
KHASH_MAP_INIT_STR(khStrMeth, methyl_read_t*);
typedef khash_t(khStrMeth) methhash_t;

//-----------------------------------------------------------------------------

methyl_read_t *methyl_read_init(bam1_t *aln)
{   /* Initialize methyl_read_t struct */

    // Initialize struct
    methyl_read_t *read = malloc(sizeof(methyl_read_t));
    if (read == NULL) {
        return NULL;
    }

    // Set variables
    const char *read_name = bam_get_qname(aln);
    read->name = strdup(read_name);
    //read->name = read_name;
    read->start = aln->core.pos;
    read->end = aln->core.pos + aln->core.l_qseq + 1;
    read->length = aln->core.isize;
    read->pos = (long *)malloc(sizeof(long) * 64);
    if (read->pos == NULL) {
        free(read->name);
        free(read);
        return NULL;
    }
    read->methyl = (int8_t *)malloc(sizeof(int8_t) * 64);
    if (read->methyl == NULL) {
        free(read->name);
        free(read->pos);
        free(read);
        return NULL;
    }
    read->qual = (uint8_t *)malloc(sizeof(uint8_t) * 64);
    if (read->qual == NULL) {
        free(read->name);
        free(read->pos);
        free(read->methyl);
        free(read);
        return NULL;
    }
    read->ncpgs = 0;
    read->size = 0;
    read->max_size = 64;

    return read;
}


void methyl_read_destroy(methyl_read_t *read)
{   /* Free memory allocated for methyl_read_t struct */

    free(read->name);
    free(read->pos);
    free(read->methyl);
    free(read->qual);
    free(read);
}


void methyl_read_append(methyl_read_t *read, long pos, int8_t methyl, uint8_t qual)
{
    // Append
    push_methyl(read, pos, methyl, qual);
    read->ncpgs++;

    return;
}


inline int isCpG(char *seq, int pos, int seqlen)
{   /* Check if position is CpG */
    if (pos >= seqlen)
    {
        return 0;
    }

    if (*(seq+pos) == 'C' || *(seq+pos) == 'c')
    {
        if (pos+1 == seqlen)
        {
            return 0;
        }
        if (*(seq+pos+1) == 'G' || *(seq+pos+1) == 'g')
        {
            return 1;
        }

        return 0;

    } else if (*(seq+pos) == 'G' || *(seq+pos) == 'g')
    {
        if (pos == 0)
        {
            return 0;
        }
        if (*(seq+pos-1) == 'C' || *(seq+pos-1) == 'c')
        {
            return -1;
        }

        return 0;
    }

    return 0;
}


int getStrand(bam1_t *b)
{   /* Get strand of read */

    char *XG = (char *) bam_aux_get(b, "XG");
    //Only bismark uses the XG tag like this. Some other aligners use it for other purposes...
    if(XG != NULL && *(XG+1) != 'C' && *(XG+1) != 'G') XG = NULL;
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 1; //OT
        }
    } else {
        if(*(XG+1) == 'C') { //OT or CTOT, due to C->T converted genome
            if((b->core.flag & 0x51) == 0x41) return 1; //Read#1 forward == OT
            else if((b->core.flag & 0x51) == 0x51) return 3; //Read #1 reverse == CTOT
            else if((b->core.flag & 0x91) == 0x81) return 3; //Read #2 forward == CTOT
            else if((b->core.flag & 0x91) == 0x91) return 1; //Read #2 reverse == OT
            else if(b->core.flag & 0x10) return 3; //Single-end reverse == CTOT
            else return 1; //Single-end forward == OT
        } else {
            if((b->core.flag & 0x51) == 0x41) return 4; //Read#1 forward == CTOB
            else if((b->core.flag & 0x51) == 0x51) return 2; //Read #1 reverse == OB
            else if((b->core.flag & 0x91) == 0x81) return 2; //Read #2 forward == OB
            else if((b->core.flag & 0x91) == 0x91) return 4; //Read #2 reverse == CTOB
            else if(b->core.flag & 0x10) return 2; //Single-end reverse == OB
            else return 4; //Single-end forward == CTOB
        }
    }
}


void processRead(bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, methyl_read_t *read)
{   /* Process read */

    uint32_t readPosition = 0;
    uint32_t mappedPosition = b->core.pos;
    int cigarOPNumber = 0;
    int cigarOPOffset = 0;
    uint32_t *CIGAR = bam_get_cigar(b);
    uint8_t *readSeq = bam_get_seq(b);
    uint8_t *readQual = bam_get_qual(b);
    int strand = getStrand(b);
    int cigarOPType;
    int direction;
    int base;

    // Parameters
    int minPhred = 5;

    while(readPosition < b->core.l_qseq && cigarOPNumber < b->core.n_cigar)
    {
        if(cigarOPOffset >= bam_cigar_oplen(CIGAR[cigarOPNumber]))
        {
            cigarOPOffset = 0;
            cigarOPNumber++;
        }
        cigarOPType = bam_cigar_type(CIGAR[cigarOPNumber]);
        if(cigarOPType & 2) { //not ISHPB
            if(cigarOPType & 1) { //M=X
                // Skip poor base calls
                if(readQual[readPosition] < minPhred) {
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                }

                direction = isCpG(seq, mappedPosition - sequenceStart, seqLen);
                if(direction)
                {
                    base = bam_seqi(readSeq, readPosition);  // Filtering by quality goes here
                    if(direction == 1 && (strand & 1) == 1) { // C & OT/CTOT
                        // methylated
                        if(base == 2) //C
                        {
                            methyl_read_append(read, mappedPosition, 1, readQual[readPosition]);
                        }
                        // unmethylated
                        else if(base == 8) //T
                        {
                            methyl_read_append(read, mappedPosition, 0, readQual[readPosition]);
                        }
                    } else if(direction == -1 && (strand & 1) == 0) { // G & OB/CTOB
                        // methylated
                        if(base == 4)  //G
                        {
                            methyl_read_append(read, mappedPosition, 1, readQual[readPosition]);
                        }
                        // unmethylated
                        else if(base == 1) //A
                        {
                            methyl_read_append(read, mappedPosition, 0, readQual[readPosition]);
                        }
                    }
                }
                mappedPosition++;
                readPosition++;
                cigarOPOffset++;
            } else { //DN
                mappedPosition += bam_cigar_oplen(CIGAR[cigarOPNumber++]);
                cigarOPOffset = 0;
                continue;
            }
        } else if(cigarOPType & 1) { // IS
            readPosition += bam_cigar_oplen(CIGAR[cigarOPNumber++]);
            cigarOPOffset = 0;
            continue;
        } else { // HPB Note that B is not handled properly, but it doesn't currently exist in the wild
            cigarOPOffset = 0;
            cigarOPNumber++;
            continue;
        }
    }

    return;
}


methyl_read_t *methyl_pair_process(methyl_read_t *read1, methyl_read_t *read2)
{   /* Process methylated read pair */

    // Initialize struct
    methyl_read_t *read = malloc(sizeof(methyl_read_t));

    // Set variables
    read->name = strdup(read1->name);
    //read->name = read1->name;
    read->start = read1->start;
    read->end = read2->end;
    read->length = read1->length;
    read->pos = (long *)malloc(sizeof(long) * 64);
    read->methyl = (int8_t *)malloc(sizeof(int8_t) * 64);
    read->qual = (uint8_t *)malloc(sizeof(uint8_t) * 64);
    read->ncpgs = 0;
    read->size = 0;
    read->max_size = 64;
    int i;

    // No CpGs in read 2
    if (read1->ncpgs > 0 && read2->ncpgs == 0)
    {
        for (i = 0; i < read1->ncpgs; i++)
        {
            methyl_read_append(read, read1->pos[i], read1->methyl[i], read1->qual[i]);
        }
    }
    else if (read1->ncpgs == 0 && read2->ncpgs > 0) // No CpGs in read 1
    {
        for (i = 0; i < read2->ncpgs; i++)
        {
            methyl_read_append(read, read2->pos[i], read2->methyl[i], read2->qual[i]);
        }
    } else {
        //int pos;
        //int qual;
        int read1_i = 0;
        int read2_i = 0;


        // Iterate over reads
        //printf("Iterating over reads\n");
        while (read1_i < read1->ncpgs || read2_i < read2->ncpgs)
        {
            //printf("   %d, %d, %d, %d\n", read1_i, read2_i, read1->ncpgs, read2->ncpgs);
            //if (read1_i > 50 || read2_i > 50) {exit(1);}

            // Check if read 1 is finished
            if (read1_i >= read1->ncpgs)
            {
                methyl_read_append(read, read2->pos[read2_i], read2->methyl[read2_i], read2->qual[read2_i]);
                read2_i++;
                continue;
            }

            // Check if read 2 is finished
            if (read2_i >= read2->ncpgs)
            {
                methyl_read_append(read, read1->pos[read1_i], read1->methyl[read1_i], read1->qual[read1_i]);
                read1_i++;
                continue;
            }

            // Check if read 1 position is before read 2 position
            if (read1->pos[read1_i] < read2->pos[read2_i])
            {
                methyl_read_append(read, read1->pos[read1_i], read1->methyl[read1_i], read1->qual[read1_i]);
                read1_i++;
            } else if (read1->pos[read1_i] == read2->pos[read2_i]) // Read 1 and read 2 positions are the same
            {
                // Take best quality read
                if (read1->qual[read1_i] >= read2->qual[read2_i])
                {
                    methyl_read_append(read, read1->pos[read1_i], read1->methyl[read1_i], read1->qual[read1_i]);
                    read1_i++;
                    read2_i++;
                } else {
                    methyl_read_append(read, read2->pos[read2_i], read2->methyl[read2_i], read2->qual[read2_i]);
                    read1_i++;
                    read2_i++;
                }
            } else { // Read 2 position is before read 1 position
                methyl_read_append(read, read2->pos[read2_i], read2->methyl[read2_i], read2->qual[read2_i]);
                read2_i++;
            }
            
        }
    }

    return read;
}


methyl_read_iterator_t *methyl_read_iterator_init(const char *bam_file_path,
                                                    char *ref_2bit,
                                                    const char *chromosome,
                                                    int min_size,
                                                    int max_size,
                                                    int qcfail,
                                                    int mapq_cutoff,
                                                    float proportion,
                                                    int nthreads)
{   /* Initialize methyl_read_iterator_t struct */

    // Initialize struct
    methyl_read_iterator_t *iter = malloc(sizeof(methyl_read_iterator_t));

    // Create read iterator
    iter->read_iter = read_iter_init(bam_file_path,
                                            chromosome,
                                            min_size,
                                            max_size,
                                            1,
                                            qcfail,
                                            mapq_cutoff,
                                            proportion,
                                            nthreads);

    // Initialize the hash table for storing read_pairs
    iter->methyl_hash = kh_init(khStrMeth);
    
    // Get reference sequence
    iter->tb = twobitOpen(ref_2bit, 0);
    iter->chrom = strdup(chromosome);
    iter->ref_start = 0;
    iter->ref_end = 1000;
    iter->ref_seq = twobitSequence(iter->tb, iter->chrom, iter->ref_start, iter->ref_end);
    if (iter->ref_seq == NULL)
    {
        fprintf(stderr, "Failed to fetch reference sequence: %s\n", iter->chrom);
        exit(1);
    }
    iter->seq_len = strlen(iter->ref_seq);
    iter->chrom_length = twobitChromLen(iter->tb, iter->chrom);
    iter->methyl_pair = NULL;

    return iter;
}


void methyl_read_iterator_destroy(methyl_read_iterator_t *iter)
{   /* Free memory allocated for methyl_read_iterator_t struct */

    // Destroy read iterator
    read_iter_destroy(iter->read_iter);

    // Destroy hash table
    methhash_t *h = (methhash_t*)iter->methyl_hash;
    khiter_t k;
    for (k = 0; k < kh_end(h); ++k)
    {
        if (kh_exist(h, k))
        {
            methyl_read_t *mread = kh_value(h, k);
            methyl_read_destroy(mread);
            //free((char*)kh_key(h, k));
        }
    }
    kh_destroy(khStrMeth, iter->methyl_hash);

    // Destroy reference sequence
    free(iter->ref_seq);
    free(iter->chrom);
    twobitClose(iter->tb);
    if (iter->methyl_pair != NULL)
    {
        methyl_read_destroy(iter->methyl_pair);
    }

    // Destroy struct
    free(iter);

    return;
}


int methyl_read_iterator_next(methyl_read_iterator_t *iter)
{   /* Iterate over BAM file */

    methhash_t *h = (methhash_t*)iter->methyl_hash;
    khiter_t k;

    // Iterator over reads
    while (read_iter_next(iter->read_iter) >= 1)
    {
        // Get reference sequence
        while (iter->read_iter->aln->core.pos > iter->ref_end && iter->ref_end < iter->chrom_length)
        {
            iter->ref_start = iter->ref_end;
            iter->ref_end = iter->ref_end + 1000;
            if (iter->ref_end > iter->chrom_length)
            {
                iter->ref_end = iter->chrom_length;
            }
            free(iter->ref_seq);
            iter->ref_seq = twobitSequence(iter->tb, iter->chrom, iter->ref_start, iter->ref_end);
            if (iter->ref_seq == NULL)
            {
                fprintf(stderr, "Failed to fetch reference sequence: %s\n", iter->chrom);
                exit(1);
            }
            iter->seq_len = strlen(iter->ref_seq);
        }

        // Initialize read
        methyl_read_t *read = methyl_read_init(iter->read_iter->aln);
        processRead(iter->read_iter->aln, iter->ref_seq, iter->ref_start, iter->seq_len, read);

        // Check if read name is in hash table
        k = kh_get(khStrMeth, h, read->name);
        if (k == kh_end(h))
        {   
            // Add read name to hash table
            int ret;
            k = kh_put(khStrMeth, h, read->name, &ret);
            kh_value(h, k) = read;

        } else {
            methyl_read_t *read2 = kh_value(h, k);

            // Free order read
            //methyl_read_t *tmp_read = iter->methyl_pair;
            if (iter->methyl_pair != NULL)
            {
                methyl_read_destroy(iter->methyl_pair);
            }
            // Determine read order
            if (read->start < read2->start)
            {
                iter->methyl_pair = methyl_pair_process(read, read2);
            } else {
                iter->methyl_pair = methyl_pair_process(read2, read);
            }
            kh_del(khStrMeth, h, k);
            methyl_read_destroy(read);
            methyl_read_destroy(read2);

            return 1;
        }

        // TMP
        //methyl_read_destroy(read);
    }

    return 0;
}