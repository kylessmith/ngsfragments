#include <stdio.h>
#include <htslib/sam.h>


typedef struct {
    char *base_calls;
    int base_calls_count;
} BaseCalls;

typedef struct {
    samFile *in;
    bam_hdr_t *header;
    bam_plp_t pileup;
    const bam_pileup1_t *pl;
    int tid, pos, n;
    BaseCalls base_calls;
} PileupIterator;


static int record_base_calls(const bam_pileup1_t *pl, int n, void *data)
{
    BaseCalls *base_calls = (BaseCalls *)data;
    base_calls->base_calls = malloc(n * sizeof(char));
    base_calls->base_calls_count = 0;

    int i;
    for (i = 0; i < n; i++)
    {
        const bam_pileup1_t *p = pl + i;
        if (!p->is_del && !p->is_refskip)
        {
            int base = bam_seqi(bam_get_seq(p->b), p->qpos);
            int qual = bam_get_qual(p->b)[p->qpos];
            if (qual >= 30) // Quality threshold
            {
                if (bam_is_rev(p->b)) // Check if the read is reverse complemented
                {
                    char *base_str = seq_nt16_str[base];
                    base_calls->base_calls[base_calls->base_calls_count++] = tolower(base_str[0]);
                }
                else
                {
                    base_calls->base_calls[base_calls->base_calls_count++] = seq_nt16_str[base];
                }
            }
        }
    }
    return 0;
}


PileupIterator *pileup_iterator_init(const char *bam_file)
{
    PileupIterator *iter = malloc(sizeof(PileupIterator));
    iter->in = sam_open(bam_file, "r");
    iter->header = sam_hdr_read(iter->in);
    iter->pileup = bam_plp_init(record_base_calls, &iter->base_calls);
    bam_plp_set_maxcnt(iter->pileup, 8000);
    return iter;
}

void pileup_iterator_destroy(PileupIterator *iter)
{
    free(iter->base_calls.base_calls);
    bam_plp_destroy(iter->pileup);
    bam_hdr_destroy(iter->header);
    sam_close(iter->in);
    free(iter);
}

int pileup_iterator_next(PileupIterator *iter)
{
    iter->pl = bam_plp_auto(iter->pileup, &iter->tid, &iter->pos, &iter->n);
    return iter->pl != 0;
}