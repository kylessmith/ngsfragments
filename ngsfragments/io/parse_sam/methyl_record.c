//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include <math.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------

static const int khLongLong = 7;
KHASH_MAP_INIT_INT64(khLongLong, long);
typedef khash_t(khLongLong) longhash_t;

//-----------------------------------------------------------------------------


int_index_t *int_index_init(long *pos, int size)
{
    // Initialize struct
    int_index_t *index = (int_index_t*)malloc(sizeof(int_index_t));
    if (index == NULL)
    {
        fprintf(stderr, "Error initializing int_index_t\n");
        exit(1);
    }

    // Initialize variables
    index->lookup = kh_init(khLongLong);
    index->pos = pos;
    index->size = size;
    index->max_size = size;

    // Add values to lookup table
    longhash_t *h = (longhash_t*)index->lookup;
    int absent;
    khiter_t k;
    long i;
    for (i = 0; i < size; i++)
    {
        k = kh_put(khLongLong, h, pos[i], &absent);
        kh_value(h, k) = i;
    }

    return index;
}


void int_index_destroy(int_index_t *index)
{   /* Free memory allocated for index */

    khiter_t k;
    longhash_t *h = (longhash_t*)index->lookup;
    kh_destroy(khLongLong, h);
    free(index->pos);
    free(index);
}


long int_index_get(int_index_t *index, long pos)
{
    longhash_t *h = (longhash_t*)index->lookup;
    khiter_t k = kh_get(khLongLong, h, pos);
    if (k == kh_end(h))
    {
        return -1;
    }

    long i = kh_value(h, k);
    
    return i;
}


methyl_record_t *methyl_record_init(long *pos, int size)
{   /* Initialize methyl record */

    // Initialize struct
    methyl_record_t *series = malloc(sizeof(methyl_record_t));
    if (series == NULL)
    {
        fprintf(stderr, "Error initializing methyl_record_t\n");
        exit(1);
    }

    // Initialize variables
    series->index = int_index_init(pos, size);

    // Initialize arrays
    series->methyl = (int16_t*)calloc(size, sizeof(int16_t));
    if (series->methyl == NULL)
    {
        fprintf(stderr, "Error initializing methyl_record_t\n");
        exit(1);
    }
    series->unmethyl = (int16_t*)calloc(size, sizeof(int16_t));
    if (series->unmethyl == NULL)
    {
        fprintf(stderr, "Error initializing methyl_record_t\n");
        exit(1);
    }

    return series;
}


void methyl_record_destroy(methyl_record_t *series)
{   /* Free memory allocated for methyl record */

    int_index_destroy(series->index);
    free(series->methyl);
    free(series->unmethyl);
    free(series);
}


void methyl_record_add(methyl_record_t *series, long *pos, int8_t *methyl, int size)
{   /* Add methylated positions to methyl record */

    int_index_t *index = series->index;

    // Iterate over positions
    int i;
    for (i = 0; i < size; i++)
    {
        long pos_i = pos[i];
        long j = int_index_get(index, pos_i);
        if (j == -1)
        {
            fprintf(stderr, "Error: position not found in series\n");
            //printf("%ld\n", pos_i);
            //exit(1);
            continue;
        }

        // Add methylated positions
        if (methyl[i] == 1)
        {
            series->methyl[j] += 1;
        }
        else if (methyl[i] == 0)
        {
            series->unmethyl[j] += 1;
            //printf("%ld\t%d\t%d\n", pos_i, methyl[i], series->unmethyl[j]);
        }
    }

    return;
}


int16_t *methyl_record_get(methyl_record_t *series, long pos)
{   /* Get methylated positions from methyl record */

    int_index_t *index = series->index;
    long i = int_index_get(index, pos);
    if (i == -1)
    {
        return NULL;
    }

    int16_t *values = malloc(sizeof(int16_t) * 2);
    values[0] = series->methyl[i];
    values[1] = series->unmethyl[i];

    return values;
}


methyl_record_pair_t *methyl_record_pair_init(methyl_record_t *record1,
                                              methyl_record_t *record2)
{   /* Initialize methyl record pair */

    // Initialize struct
    methyl_record_pair_t *pair = malloc(sizeof(methyl_record_pair_t));
    if (pair == NULL)
    {
        fprintf(stderr, "Error initializing methyl_record_pair_t\n");
        exit(1);
    }

    // Initialize variables
    pair->record1 = record1;
    pair->record2 = record2;

    return pair;
}


void methyl_record_pair_destroy(methyl_record_pair_t *pair)
{   /* Free memory allocated for methyl record pair */

    methyl_record_destroy(pair->record1);
    methyl_record_destroy(pair->record2);
    free(pair);
}


void methyl_record_pair_transfer_null(methyl_record_pair_t *pair)
{   /* Transfer null methyl record */

    int i;
    for (i = 0; i < pair->record1->index->size; i++)
    {
        if (pair->record1->methyl[i] == 0 && pair->record1->unmethyl[i] == 0)
        {
            pair->record1->methyl[i] = pair->record2->methyl[i];
            pair->record1->unmethyl[i] = pair->record2->unmethyl[i];
        }
        else if (pair->record2->methyl[i] == 0 && pair->record2->unmethyl[i] == 0)
        {
            pair->record2->methyl[i] = pair->record1->methyl[i];
            pair->record2->unmethyl[i] = pair->record1->unmethyl[i];
        }
    }

    return;
}


/* Jeffreys prior pseudocount for the Beta-Bernoulli posterior methylation rate.
 * beta_hat = (methyl + A) / (methyl + unmethyl + 2A) is the posterior mean under a
 * Beta(A, A) prior. A = 0.5 shrinks low-depth sites toward 0.5 so a single 1/0 read
 * (depth 1) no longer produces a hard 1.0/0.0 that dominates the distance. */
#define METHYL_BETA_PRIOR 0.1

int assign_methyl_read(methyl_record_pair_t *pair, methyl_read_t *read)
{   /* Assign a read to whichever profile its per-CpG methylation better matches.
     *
     * Distance to profile p is a COVERAGE-WEIGHTED mean absolute deviation between
     * the profile's Beta-Bernoulli posterior methylation rate and this read's
     * per-CpG call:
     *
     *     dist_p = sum_i [ depth_pi * |beta_pi - r_i| ]  /  sum_i depth_pi
     *
     * where depth_pi = methyl+unmethyl in profile p at CpG i, r_i in {0,1} is the
     * read's call, and beta_pi is the pseudocount-shrunk posterior rate. High-depth
     * (confident) CpGs dominate; depth-0 CpGs contribute nothing to either the
     * numerator or the denominator and so drop out automatically.
     *
     * Fixes vs. the previous version:
     *  - No division by zero. Previously sum/n with n==0 produced NaN, and
     *    `sum1 <= NaN` is false, so every read with no coverage in profile 2 was
     *    silently sent to profile 2 — the side with NO evidence. Now the read goes
     *    to the side that HAS evidence (or a deterministic default if neither does).
     *  - Coverage weighting + shrinkage: depth-1 sites no longer count the same as
     *    depth-100 sites, which is what matters in the low-coverage cfDNA regime.
     *  - No per-CpG malloc/free: reads record arrays directly via the shared index
     *    instead of methyl_record_get()'s 2-int16 heap allocation (previously
     *    4 mallocs per CpG per read, every EM iteration).
     *
     * Return contract is unchanged: 0 -> record1, 1 -> record2, ties -> record1. */

    int_index_t *idx1 = pair->record1->index;
    int_index_t *idx2 = pair->record2->index;

    double num1 = 0.0, den1 = 0.0;   /* weighted deviation sum / total weight, profile 1 */
    double num2 = 0.0, den2 = 0.0;   /* ... profile 2 */

    int i;
    for (i = 0; i < read->ncpgs; i++)
    {
        long pos = read->pos[i];
        double r = (double)read->methyl[i];   /* 0 or 1 */

        /* Profile 1 */
        long j1 = int_index_get(idx1, pos);
        if (j1 != -1)
        {
            double m = (double)pair->record1->methyl[j1];
            double u = (double)pair->record1->unmethyl[j1];
            double depth = m + u;
            if (depth > 0.0)
            {
                double beta = (m + METHYL_BETA_PRIOR) / (depth + 2.0 * METHYL_BETA_PRIOR);
                num1 += depth * fabs(beta - r);
                den1 += depth;
            }
        }

        /* Profile 2 */
        long j2 = int_index_get(idx2, pos);
        if (j2 != -1)
        {
            double m = (double)pair->record2->methyl[j2];
            double u = (double)pair->record2->unmethyl[j2];
            double depth = m + u;
            if (depth > 0.0)
            {
                double beta = (m + METHYL_BETA_PRIOR) / (depth + 2.0 * METHYL_BETA_PRIOR);
                num2 += depth * fabs(beta - r);
                den2 += depth;
            }
        }
    }

    /* Assignment with explicit handling of the zero-evidence cases. */
    if (den1 == 0.0 && den2 == 0.0)
        return 0;              /* no coverage in EITHER profile at this read's CpGs:
                                 * genuinely uninformative -> deterministic default
                                 * (record1), consistent with tie-breaking below. */
    if (den2 == 0.0)
        return 0;              /* only profile 1 has evidence -> record1 */
    if (den1 == 0.0)
        return 1;              /* only profile 2 has evidence -> record2 */

    double dist1 = num1 / den1;
    double dist2 = num2 / den2;

    return (dist1 <= dist2) ? 0 : 1;   /* ties -> record1, matching original */
}


double compare_methyl_records(methyl_record_pair_t *pair)
{   /* Compare methyl records using Euclidean distance */

    // Calculate distance
    double sum = 0;

    // Iterate over positions
    int size = pair->record1->index->size;
    int i;
    for (i = 0; i < size; i++)
    {
        if (pair->record1->methyl[i] == 0 && pair->record1->unmethyl[i] == 0)
        {
            continue;
        }
        if (pair->record2->methyl[i] == 0 && pair->record2->unmethyl[i] == 0)
        {
            continue;
        }
        double beta1 = (double)pair->record1->methyl[i] / ((double)pair->record1->unmethyl[i] + (double)pair->record1->methyl[i]);
        double beta2 = (double)pair->record2->methyl[i] / ((double)pair->record2->unmethyl[i] + (double)pair->record2->methyl[i]);
        // Binarize
        //beta1 = (beta1 >= 0.5) ? 1 : 0;
        //beta2 = (beta2 >= 0.5) ? 1 : 0;
        // Euclidean distance
        double diff = beta1 - beta2;
        sum += diff * diff;
        // Manhattan distance
        //sum += fabs(beta1 - beta2);
    }
    sum = sqrt(sum);

    return sum;
}


int methyl_record_pair_write(methyl_record_pair_t *pair, char *file_fn)
{   /* Write methyl record pair to file */

    // Open file
    FILE *fp = fopen(file_fn, "w");

    // Iterate over positions
    int size = pair->record1->index->size;
    int i;
    for (i = 0; i < size; i++)
    {
        //double beta1 = (double)pair->record1->methyl[i] / ((double)pair->record1->unmethyl[i] + (double)pair->record1->methyl[i]);
        //double beta2 = (double)pair->record2->methyl[i] / ((double)pair->record2->unmethyl[i] + (double)pair->record2->methyl[i]);
        // Binarize
        //beta1 = (beta1 >= 0.5) ? 1 : 0;
        //beta2 = (beta2 >= 0.5) ? 1 : 0;
        // Write to file
        fprintf(fp, "%d\t%d\t%d\t%d\n", (int)pair->record1->methyl[i], (int)pair->record1->unmethyl[i], (int)pair->record2->methyl[i], (int)pair->record2->unmethyl[i]);
    }

    // Close file
    fclose(fp);

    return 0;
}