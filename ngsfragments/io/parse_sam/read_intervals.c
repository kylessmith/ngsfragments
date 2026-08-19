//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

// check_read is now static inline in read_intervals.h — inlined at every call site.

//-----------------------------------------------------------------------------
// Chromosome name table helpers
// Build once per file open: maps tid -> "chr"-prefixed (or raw) name.
// Eliminates strcat + 100-byte stack init + branch from the per-read hot loop.
//-----------------------------------------------------------------------------

static char **build_chrom_table(sam_hdr_t *sam_hdr, int add_chr)
{
    int n = sam_hdr->n_targets;
    char **table = (char **)malloc(n * sizeof(char *));
    if (table == NULL) return NULL;

    for (int i = 0; i < n; i++)
    {
        const char *raw = sam_hdr->target_name[i];
        if (add_chr)
        {
            size_t len = strlen(raw);
            table[i] = (char *)malloc(len + 4);
            if (table[i]) { table[i][0]='c'; table[i][1]='h'; table[i][2]='r';
                            memcpy(table[i]+3, raw, len+1); }
        }
        else
        {
            table[i] = strdup(raw);
        }
    }
    return table;
}

static void free_chrom_table(char **table, int n)
{
    for (int i = 0; i < n; i++) free(table[i]);
    free(table);
}

//-----------------------------------------------------------------------------


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
{   /* Add reads from sam file to interval list (whole-file sequential scan) */

    samFile *fp_in = hts_open(samfile_name, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();

    if (nthreads > 1)
        hts_set_threads(fp_in, nthreads);

    // Precompute: reject mask (one bitwise op replaces 4+ flag checks per read)
    uint32_t reject_mask = build_reject_mask(qcfail);

    // Precompute: chrom name table (eliminates strcat + branch from hot loop)
    char **chrom_table = build_chrom_table(sam_hdr, add_chr);

    while (sam_read1(fp_in, sam_hdr, aln) >= 0)
    {
        if (!check_read(aln, reject_mask, min_size, max_size, paired, mapq_cutoff, proportion))
            continue;

        /* Paired: emit one interval per fragment. The upstream mate carries the
         * positive isize and sits at the fragment's leftmost coordinate, so it
         * alone defines [pos, pos + isize). Skipping isize<=0 de-duplicates the
         * pair without depending on which mate is R1/R2. */
        int start = aln->core.pos;
        int tlen;
        if (paired == 1)
        {
            if (aln->core.isize <= 0)
                continue;                 /* downstream mate (or unset isize) — already counted */
            tlen = (int)aln->core.isize;
        }
        else
        {
            tlen = aln->core.l_qseq;
        }
        labeled_aiarray_add(intervals, start, start + tlen, chrom_table[aln->core.tid]);
    }

    free_chrom_table(chrom_table, sam_hdr->n_targets);
    bam_destroy1(aln);
    sam_close(fp_in);
    sam_hdr_destroy(sam_hdr);
}


void sam_iter_add_region(char *samfile_name,
                         labeled_aiarray_t *intervals,
                         const char *chromosome,
                         int min_size,
                         int max_size,
                         int paired,
                         int qcfail,
                         int mapq_cutoff,
                         float proportion,
                         int nthreads,
                         int add_chr)
{   /* Add reads from one chromosome using the BAM index.
     * Each call is independent: safe to run per-chromosome in parallel processes.
     * Workers open their own file handle — no shared state. */

    samFile *fp_in = hts_open(samfile_name, "r");
    if (fp_in == NULL)
    {
        fprintf(stderr, "sam_iter_add_region: failed to open %s\n", samfile_name);
        return;
    }

    if (nthreads > 1)
        hts_set_threads(fp_in, nthreads);

    sam_hdr_t *sam_hdr = sam_hdr_read(fp_in);
    if (sam_hdr == NULL)
    {
        fprintf(stderr, "sam_iter_add_region: failed to read header from %s\n", samfile_name);
        sam_close(fp_in);
        return;
    }

    hts_idx_t *idx = sam_index_load(fp_in, samfile_name);
    if (idx == NULL)
    {
        fprintf(stderr, "sam_iter_add_region: failed to load index for %s\n", samfile_name);
        sam_hdr_destroy(sam_hdr);
        sam_close(fp_in);
        return;
    }

    hts_itr_t *iter = sam_itr_querys(idx, sam_hdr, chromosome);
    if (iter == NULL)
    {
        fprintf(stderr, "sam_iter_add_region: chromosome '%s' not found in %s\n",
                chromosome, samfile_name);
        hts_idx_destroy(idx);
        sam_hdr_destroy(sam_hdr);
        sam_close(fp_in);
        return;
    }

    uint32_t reject_mask = build_reject_mask(qcfail);
    char **chrom_table   = build_chrom_table(sam_hdr, add_chr);

    bam1_t *aln = bam_init1();
    while (sam_itr_next(fp_in, iter, aln) >= 0)
    {
        if (!check_read(aln, reject_mask, min_size, max_size, paired, mapq_cutoff, proportion))
            continue;

        /* Paired: emit one interval per fragment. The upstream mate carries the
         * positive isize and sits at the fragment's leftmost coordinate, so it
         * alone defines [pos, pos + isize). Skipping isize<=0 de-duplicates the
         * pair without depending on which mate is R1/R2. */
        int start = aln->core.pos;
        int tlen;
        if (paired == 1)
        {
            if (aln->core.isize <= 0)
                continue;                 /* downstream mate (or unset isize) — already counted */
            tlen = (int)aln->core.isize;
        }
        else
        {
            tlen = aln->core.l_qseq;
        }
        labeled_aiarray_add(intervals, start, start + tlen, chrom_table[aln->core.tid]);
    }

    free_chrom_table(chrom_table, sam_hdr->n_targets);
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(sam_hdr);
    sam_close(fp_in);
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
{   /* Add reads and centre a fixed-size window on the 5' end of each read strand. */

    samFile *fp_in = hts_open(samfile_name, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();

    if (nthreads > 1)
        hts_set_threads(fp_in, nthreads);

    uint32_t reject_mask = build_reject_mask(qcfail);
    char **chrom_table   = build_chrom_table(sam_hdr, add_chr);

    int half = fixed_size / 2;

    while (sam_read1(fp_in, sam_hdr, aln) >= 0)
    {
        if (!check_read(aln, reject_mask, min_size, max_size, paired, mapq_cutoff, proportion))
            continue;

        int start = aln->core.pos;
        /* One window per fragment: keep only the upstream mate (isize>0) when paired. */
        if (paired == 1 && aln->core.isize <= 0)
            continue;
        int tlen  = (paired == 1) ? (int)aln->core.isize : aln->core.l_qseq;

        /* Centre a fixed-size window on the 5' end of each read strand.
         * For the reverse-strand read the 5' end is at aln->core.pos (leftmost coord).
         * For the forward-strand read the 5' end is at pos + tlen (right fragment endpoint). */
        int five_prime  = bam_is_rev(aln) ? start : start + tlen;
        int fixed_start = five_prime - half;
        int fixed_end   = fixed_start + fixed_size;

        labeled_aiarray_add(intervals, fixed_start, fixed_end, chrom_table[aln->core.tid]);
    }

    free_chrom_table(chrom_table, sam_hdr->n_targets);
    bam_destroy1(aln);
    sam_close(fp_in);
    sam_hdr_destroy(sam_hdr);
}
