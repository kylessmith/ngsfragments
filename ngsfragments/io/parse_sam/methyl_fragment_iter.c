#include "read_intervals.h"

//-----------------------------------------------------------------------------

#define REF_WINDOW_SIZE 1000

KHASH_MAP_INIT_STR(khStrMeth, methyl_read_t*);
typedef khash_t(khStrMeth) methhash_t;

//-----------------------------------------------------------------------------

methyl_read_t *methyl_read_init(bam1_t *aln, bam_hdr_t *header)
{   /* Initialize methyl_read_t struct */

    // Initialize struct
    methyl_read_t *read = malloc(sizeof(methyl_read_t));
    if (read == NULL) {
        return NULL;
    }

    // Set variables
    read->header = header;
    read->read = bam_init1();
    bam_copy1(read->read, aln);
    const char *read_name = bam_get_qname(aln);
    read->name = strdup(read_name);
    read->start = aln->core.pos;
    read->end = bam_endpos(aln);           /* CIGAR-aware; l_qseq ignores deletions */
    read->length = llabs(aln->core.isize);  /* isize is hts_pos_t (long long); llabs avoids truncation */
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

    if (read == NULL)
    {
        return;
    }

    free(read->name);
    free(read->pos);
    free(read->methyl);
    free(read->qual);
    
    // Only free BAM data if it exists (combined reads have NULL here)
    if (read->read != NULL) {
        bam_destroy1(read->read);
    }

    free(read);
}


void methyl_read_append(methyl_read_t *read, long pos, int8_t methyl, uint8_t qual)
{
    // Append. push_methyl bumps read->size ONLY when the element is actually stored
    // (a failed realloc now drops it instead of corrupting the heap — see the macro).
    // Previously ncpgs was incremented unconditionally, so a dropped append left
    // ncpgs > size; every consumer loops on ncpgs, so it would read pos[size..ncpgs-1]
    // past the valid region. Mirroring ncpgs to size makes desync impossible.
    push_methyl(read, pos, methyl, qual);
    read->ncpgs = read->size;

    return;
}


static inline int isCpG(char *seq, int pos, int seqlen)
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


void processRead(bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, methyl_read_t *read, int minPhred)
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

    // Set strand in read
    read->strand = strand;

    // getStrand returns 0 when the strand cannot be determined (e.g. a paired read
    // with an unexpected flag combination and no usable XG tag). We cannot place
    // CpG calls without knowing the strand — the direction==-1 / (strand&1)==0
    // branch below would otherwise MIS-process such a read as OB and shift every
    // call by -1. Skip extraction; the read keeps 0 CpGs and still pairs normally
    // (its mate's calls, if on a determinable strand, are used via methyl_pair_process).
    if (strand == 0)
        return;

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
                // Skip poor base calls; advance counters and move on
                if(readQual[readPosition] < minPhred) {
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                } else {
                    direction = isCpG(seq, mappedPosition - sequenceStart, seqLen);
                    if(direction)
                    {
                        base = bam_seqi(readSeq, readPosition);
                        if(direction == 1 && (strand & 1) == 1) { // C on OT/CTOT strand
                            // methylated
                            if(base == 2) //C
                            {
                                if (strand == 1 || strand == 3) // OT or CTOT: C is at the CpG position
                                {
                                    methyl_read_append(read, mappedPosition, 1, readQual[readPosition]);
                                }
                                else if (strand == 2 || strand == 4) // OB or CTOB: shift to C position
                                {
                                    methyl_read_append(read, mappedPosition - 1, 1, readQual[readPosition]);
                                }
                            }
                            // unmethylated
                            else if(base == 8) //T
                            {
                                if (strand == 1 || strand == 3) // OT or CTOT
                                {
                                    methyl_read_append(read, mappedPosition, 0, readQual[readPosition]);
                                }
                                else if (strand == 2 || strand == 4) // OB or CTOB: shift to C position
                                {
                                    methyl_read_append(read, mappedPosition - 1, 0, readQual[readPosition]);
                                }
                            }
                        } else if(direction == -1 && (strand & 1) == 0) { // G on OB/CTOB strand
                            // Outer gate guarantees strand is 2 or 4 (OB/CTOB).
                            // G is always one base right of the C; shift by -1 to report at CpG (C) position.
                            // methylated
                            if(base == 4) //G
                            {
                                methyl_read_append(read, mappedPosition - 1, 1, readQual[readPosition]);
                            }
                            // unmethylated
                            else if(base == 1) //A
                            {
                                methyl_read_append(read, mappedPosition - 1, 0, readQual[readPosition]);
                            }
                        }
                    }
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                } // end else (quality filter passed)
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


/*
 * processReadXM — reference-free methylation extraction using the Bismark XM tag.
 *
 * The XM auxiliary tag encodes the methylation context of every query base in
 * alignment order (gaps in the read skipped, deletions represented by dots):
 *
 *   Z / z  — methylated / unmethylated CpG
 *   X / x  — methylated / unmethylated CHG
 *   H / h  — methylated / unmethylated CHH
 *   U / u  — methylated / unmethylated unknown context
 *   .       — non-cytosine or deletion in reference
 *
 * Only CpG calls (Z/z) are recorded.  Position is reported at the C of the CpG
 * on the forward strand, matching the coordinate convention used by processRead:
 *   OT/CTOT reads (C observed): position = mappedPosition
 *   OB/CTOB reads (G observed): position = mappedPosition - 1
 *
 * Returns 1 on success, 0 if the XM tag is absent (caller should fall back to
 * the reference-based path or skip the read).
 */
static int processReadXM(bam1_t *b, methyl_read_t *read, int minPhred)
{
    uint8_t *xm_raw = bam_aux_get(b, "XM");
    if (xm_raw == NULL) return 0;               /* tag absent — can't use this path */
    const char *XM = bam_aux2Z(xm_raw);         /* pointer into the BAM data block  */
    if (XM == NULL) return 0;

    int strand = getStrand(b);
    read->strand = strand;

    // Undeterminable strand (see processRead): the OB/CTOB coordinate shift (-1)
    // is applied via the `else` branch below, so a strand-0 read would have every
    // CpG mis-shifted. Return success with 0 CpGs recorded — consistent with the
    // reference-based path, and the read still pairs normally.
    if (strand == 0)
        return 1;

    uint8_t *readQual = bam_get_qual(b);
    uint32_t *CIGAR   = bam_get_cigar(b);

    uint32_t readPosition   = 0;
    uint32_t mappedPosition = b->core.pos;
    int      xmPosition     = 0;           /* index into XM string (query bases only) */
    int      cigarOPNumber  = 0;
    int      cigarOPOffset  = 0;
    int      cigarOPType;

    while (readPosition < (uint32_t)b->core.l_qseq && cigarOPNumber < b->core.n_cigar)
    {
        if (cigarOPOffset >= (int)bam_cigar_oplen(CIGAR[cigarOPNumber]))
        {
            cigarOPOffset = 0;
            cigarOPNumber++;
        }
        cigarOPType = bam_cigar_type(CIGAR[cigarOPNumber]);

        if (cigarOPType & 2) { /* consumes reference */
            if (cigarOPType & 1) { /* also consumes query (M=X) */
                char ctx = XM[xmPosition];
                if (readQual[readPosition] >= minPhred && (ctx == 'Z' || ctx == 'z'))
                {
                    int8_t  methyl = (ctx == 'Z') ? 1 : 0;
                    uint32_t cpg_pos;
                    if (strand == 1 || strand == 3)       /* OT / CTOT: C at this pos */
                        cpg_pos = mappedPosition;
                    else                                   /* OB / CTOB: G, C is at -1 */
                        cpg_pos = mappedPosition - 1;
                    methyl_read_append(read, (long)cpg_pos, methyl, readQual[readPosition]);
                }
                mappedPosition++;
                readPosition++;
                xmPosition++;
                cigarOPOffset++;
            } else { /* D/N: consumes reference only; XM uses '.' for ref deletions */
                uint32_t oplen = bam_cigar_oplen(CIGAR[cigarOPNumber++]);
                mappedPosition += oplen;
                xmPosition     += oplen;   /* XM dots cover deleted reference bases */
                cigarOPOffset   = 0;
                continue;
            }
        } else if (cigarOPType & 1) { /* I/S: consumes query only */
            uint32_t oplen = bam_cigar_oplen(CIGAR[cigarOPNumber++]);
            readPosition += oplen;
            xmPosition   += oplen;
            cigarOPOffset = 0;
            continue;
        } else { /* H/P/B */
            cigarOPOffset = 0;
            cigarOPNumber++;
            continue;
        }
    }

    return 1;
}



methyl_read_t *methyl_pair_process(methyl_read_t *read1, methyl_read_t *read2)
{   /* Process methylated read pair */

    // Initialize struct
    methyl_read_t *read = malloc(sizeof(methyl_read_t));
    if (read == NULL)
    {
        return NULL;
    }

    // Set variables
    read->name = strdup(read1->name);
    if (read->name == NULL)
    {
        free(read);
        return NULL;
    }

    read->start = read1->start;
    read->end = read2->end;
    read->length = read1->length;
    read->strand = read1->strand;

    // Initialize BAM-related fields to NULL (combined reads don't need these)
    read->read = NULL;
    read->header = NULL;

    read->pos = (long *)malloc(sizeof(long) * 64);
    if (read->pos == NULL)
    {
        free(read->name);
        free(read);
        return NULL;
    }

    read->methyl = (int8_t *)malloc(sizeof(int8_t) * 64);
    if (read->methyl == NULL)
    {
        free(read->name);
        free(read->pos);
        free(read);
        return NULL;
    }
    read->qual = (uint8_t *)malloc(sizeof(uint8_t) * 64);
    if (read->qual == NULL)
    {
        free(read->name);
        free(read->pos);
        free(read->methyl);
        free(read);
        return NULL;
    }

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
        int read1_i = 0;
        int read2_i = 0;

        // Iterate over reads, merging by position (stable sort; ties go to higher quality)
        while (read1_i < read1->ncpgs || read2_i < read2->ncpgs)
        {
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


/*
 * Reference backend abstraction
 * ─────────────────────────────
 * Callers use ref_open / ref_chrom_length / ref_fetch_seq / ref_close and
 * never touch TwoBit or faidx_t directly.  Exactly one of *tb_out / *fai_out
 * will be non-NULL after a successful ref_open.
 */

/* Detect file type by extension and open the appropriate handle.
 * Returns 1 for .2bit, 2 for FASTA, 0 on failure.
 * Recognised FASTA extensions: .fa  .fasta  .fa.gz  .fasta.gz */
static int ref_open(const char *path, TwoBit **tb_out, faidx_t **fai_out)
{
    *tb_out  = NULL;
    *fai_out = NULL;

    if (path == NULL) return 0;

    /* Check for .2bit suffix */
    const char *dot = strrchr(path, '.');
    if (dot != NULL && strcmp(dot, ".2bit") == 0) {
        *tb_out = twobitOpen(path, 0);
        return (*tb_out != NULL) ? 1 : 0;
    }

    /* Everything else treated as FASTA; fai_load will create a .fai index
     * automatically if one is not already present. */
    *fai_out = fai_load(path);
    return (*fai_out != NULL) ? 2 : 0;
}

/* Return chromosome length, or ≤0 if not found.
 * faidx_seq_len returns -1 for unknown sequences. */
static long long ref_chrom_length(TwoBit *tb, faidx_t *fai, const char *chrom)
{
    if (tb  != NULL) return (long long)twobitChromLen(tb,  (char *)chrom);
    if (fai != NULL) return (long long)faidx_seq_len(fai,           chrom);
    return 0LL;
}

/* Fetch the reference sequence for the window [start, end).
 * twobitSequence uses 0-based half-open [start, end).
 * faidx_fetch_seq uses 0-based inclusive [start, end-1].
 * The returned pointer must be freed by the caller. */
static char *ref_fetch_seq(TwoBit *tb, faidx_t *fai,
                           const char *chrom, int start, int end)
{
    if (tb  != NULL) return twobitSequence(tb, (char *)chrom, start, end);
    if (fai != NULL) {
        int len = 0;
        return faidx_fetch_seq(fai, chrom, start, end - 1, &len);
    }
    return NULL;
}

/* Close whichever handle is non-NULL. */
static void ref_close(TwoBit *tb, faidx_t *fai)
{
    if (tb  != NULL) twobitClose(tb);
    if (fai != NULL) fai_destroy(fai);
}


methyl_read_iterator_t *methyl_read_iterator_init(const char *bam_file_path,
                                                    const char *ref_file,
                                                    const char *chromosome,
                                                    int start_pos,
                                                    int end_pos,
                                                    int min_size,
                                                    int max_size,
                                                    int qcfail,
                                                    int mapq_cutoff,
                                                    float proportion,
                                                    int nthreads)
{   /* Initialize methyl_read_iterator_t struct */

    // Initialize struct
    methyl_read_iterator_t *iter = malloc(sizeof(methyl_read_iterator_t));
    if (iter == NULL) {
        fprintf(stderr, "Failed to allocate methyl_read_iterator_t\n");
        return NULL;
    }

    // Create read iterator
    iter->read_iter = read_iter_init(bam_file_path,
                                            chromosome,
                                            start_pos,
                                            end_pos,
                                            min_size,
                                            max_size,
                                            1,
                                            qcfail,
                                            mapq_cutoff,
                                            proportion,
                                            nthreads);

    // Initialize the hash table for storing read_pairs
    iter->methyl_hash = kh_init(khStrMeth);

    // Open reference file (auto-detected as .2bit or FASTA by extension)
    int ref_type = ref_open(ref_file, &iter->tb, &iter->fai);
    if (ref_type == 0) {
        fprintf(stderr, "Failed to open reference file: %s\n",
                ref_file ? ref_file : "(null)");
        kh_destroy(khStrMeth, iter->methyl_hash);
        read_iter_destroy(iter->read_iter);
        free(iter);
        return NULL;
    }

    iter->chrom = strdup(chromosome);

    // Validate that the chromosome exists in the reference before fetching sequence.
    iter->chrom_length = (int)ref_chrom_length(iter->tb, iter->fai, iter->chrom);
    if (iter->chrom_length <= 0)
    {
        // Common failure: BAM uses "chr1" but reference uses "1", or vice versa.
        // Try toggling the "chr" prefix before giving up.
        char *alt_chrom = NULL;
        if (strncmp(chromosome, "chr", 3) == 0) {
            alt_chrom = strdup(chromosome + 3);          /* "chr1" -> "1"   */
        } else {
            alt_chrom = malloc(strlen(chromosome) + 4);  /* "1"   -> "chr1" */
            if (alt_chrom != NULL) sprintf(alt_chrom, "chr%s", chromosome);
        }

        long long alt_len = (alt_chrom != NULL)
            ? ref_chrom_length(iter->tb, iter->fai, alt_chrom)
            : 0LL;

        if (alt_len > 0) {
            fprintf(stderr,
                    "Warning: chromosome '%s' not found in reference; "
                    "using '%s' instead\n", chromosome, alt_chrom);
            free(iter->chrom);
            iter->chrom = alt_chrom;
            iter->chrom_length = (int)alt_len;
        } else {
            free(alt_chrom);
            fprintf(stderr,
                    "Warning: chromosome '%s' not found in reference under any known "
                    "naming convention; falling back to Bismark XM tag (reference-free)\n",
                    chromosome);
            iter->ref_seq = NULL;
            iter->seq_len = 0;
            iter->use_xm  = 1;
            iter->methyl_pair = NULL;
            iter->read1 = NULL;
            iter->read2 = NULL;
            return iter;
        }
    }

    // Seed the reference window at the query start position so region queries
    // don't have to slide forward from position 0 one window at a time.
    // For whole-chromosome queries (start_pos <= 0) begin at the chromosome start.
    iter->ref_start = (start_pos > 0) ? start_pos : 0;
    iter->ref_end   = iter->ref_start + REF_WINDOW_SIZE;
    if (iter->ref_end > iter->chrom_length)
        iter->ref_end = iter->chrom_length;
    iter->ref_seq   = ref_fetch_seq(iter->tb, iter->fai, iter->chrom,
                                    iter->ref_start, iter->ref_end);
    if (iter->ref_seq == NULL)
    {
        fprintf(stderr, "Failed to fetch reference sequence for '%s'\n", iter->chrom);
        free(iter->chrom);
        ref_close(iter->tb, iter->fai);
        kh_destroy(khStrMeth, iter->methyl_hash);
        read_iter_destroy(iter->read_iter);
        free(iter);
        return NULL;
    }
    iter->seq_len = strlen(iter->ref_seq);
    iter->use_xm = 0;
    iter->methyl_pair = NULL;
    iter->read1 = NULL;
    iter->read2 = NULL;

    return iter;
}


void methyl_read_iterator_destroy(methyl_read_iterator_t *iter)
{   /* Free memory allocated for methyl_read_iterator_t struct */

    // Destroy read iterator
    read_iter_destroy(iter->read_iter);

    // Destroy hash table — free any unpaired reads and their key copies
    methhash_t *h = (methhash_t*)iter->methyl_hash;
    khiter_t k;
    for (k = 0; k < kh_end(h); ++k)
    {
        if (kh_exist(h, k))
        {
            free((char*)kh_key(h, k));       /* free the strdup'd key */
            methyl_read_t *mread = kh_value(h, k);
            methyl_read_destroy(mread);
        }
    }
    kh_destroy(khStrMeth, iter->methyl_hash);

    // Destroy reference
    free(iter->ref_seq);
    free(iter->chrom);
    ref_close(iter->tb, iter->fai);
    if (iter->methyl_pair != NULL)
    {
        methyl_read_destroy(iter->methyl_pair);
    }
    // Destroy read1 and read2 if they exist
    if (iter->read1 != NULL) {
        methyl_read_destroy(iter->read1);
    }
    if (iter->read2 != NULL) {
        methyl_read_destroy(iter->read2);
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
        // Advance reference sequence window (reference-based path only).
        // Jump directly to the window containing the read rather than sliding
        // one step at a time, which would be O(distance/window_size) fetches.
        if (!iter->use_xm) {
            hts_pos_t read_pos = iter->read_iter->aln->core.pos;
            if (read_pos >= iter->ref_end && iter->ref_end < iter->chrom_length)
            {
                iter->ref_start = (int)(read_pos / REF_WINDOW_SIZE) * REF_WINDOW_SIZE;
                iter->ref_end   = iter->ref_start + REF_WINDOW_SIZE;
                if (iter->ref_end > iter->chrom_length)
                    iter->ref_end = iter->chrom_length;
                free(iter->ref_seq);
                iter->ref_seq = ref_fetch_seq(iter->tb, iter->fai, iter->chrom,
                                              iter->ref_start, iter->ref_end);
                if (iter->ref_seq == NULL)
                {
                    fprintf(stderr, "Failed to fetch reference window [%d, %d) for '%s'\n",
                            iter->ref_start, iter->ref_end, iter->chrom);
                    return -1;
                }
                iter->seq_len = strlen(iter->ref_seq);
            }
        } /* end !use_xm */

        // Initialize read and extract methylation calls
        methyl_read_t *read = methyl_read_init(iter->read_iter->aln, iter->read_iter->header);
        if (iter->use_xm) {
            if (!processReadXM(iter->read_iter->aln, read, 20)) {
                /* XM tag absent on this read — skip it silently */
                methyl_read_destroy(read);
                continue;
            }
        } else {
            processRead(iter->read_iter->aln, iter->ref_seq, iter->ref_start, iter->seq_len, read, 20);
        }

        // Check if read name is in hash table
        k = kh_get(khStrMeth, h, read->name);
        if (k == kh_end(h))
        {   
            // Insert a copy of the name as the key so it outlives the read struct
            int ret;
            k = kh_put(khStrMeth, h, strdup(read->name), &ret);
            kh_value(h, k) = read;

        } else {
            methyl_read_t *read2 = kh_value(h, k);

            // Destroy previous pair and individual reads before creating new ones
            if (iter->methyl_pair != NULL)
            {
                methyl_read_destroy(iter->methyl_pair);
            }
            if (iter->read1 != NULL) {
                methyl_read_destroy(iter->read1);
            }
            if (iter->read2 != NULL) {
                methyl_read_destroy(iter->read2);
            }

            // Determine read order
            if (read->start < read2->start)
            {
                iter->methyl_pair = methyl_pair_process(read, read2);
                iter->read1 = read;
                iter->read2 = read2;
            } else {
                iter->methyl_pair = methyl_pair_process(read2, read);
                iter->read1 = read2;
                iter->read2 = read;
            }
            // Free the key copy that was strdup'd on insertion
            free((char*)kh_key(h, k));
            kh_del(khStrMeth, h, k);

            return 1;
        }
    }

    return 0;
}