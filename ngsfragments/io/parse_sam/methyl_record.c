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


int assign_methyl_read(methyl_record_pair_t *pair, methyl_read_t *read)
{   /* Assign read to methyl record pair */

    // Calculate distances
    int n1 = 0;
    int n2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    int i;
    for (i = 0; i < read->ncpgs; i++)
    {
        long pos = read->pos[i];
        int8_t methyl = read->methyl[i];
        int16_t *values1 = methyl_record_get(pair->record1, pos);
        int16_t *values2 = methyl_record_get(pair->record2, pos);
        //printf("%ld\t%d\t%d\t%d\t%d\n", pos, values1[0], values1[1], values2[0], values2[1]);
        /*if (values1[0] == 0 && values1[1] == 0)
        {
            continue;
        }
        if (values2[0] == 0 && values2[1] == 0)
        {
            continue;
        }*/

        // Calculate beta values
        double beta1;
        if (values1[0] != 0 || values1[1] != 0)
        {
            beta1 = (double)values1[0] / ((double)values1[0] + (double)values1[1]);
            //beta1 = (beta1 >= 0.5) ? 1 : 0;
            sum1 += fabs(beta1 - (double)methyl);
            n1++;
        }

        double beta2;
        if (values2[0] != 0 || values2[1] != 0)
        {
            beta2 = (double)values2[0] / ((double)values2[0] + (double)values2[1]);
            //beta2 = (beta2 >= 0.5) ? 1 : 0;
            sum2 += fabs(beta2 - (double)methyl);
            n2++;
        }

        //double beta1 = (double)values1[0] / ((double)values1[0] + (double)values1[1]);
        //double beta2 = (double)values2[0] / ((double)values2[0] + (double)values2[1]);
        // Binarize
        //beta1 = (beta1 >= 0.5) ? 1 : 0;
        //beta2 = (beta2 >= 0.5) ? 1 : 0;
        //sum1 += fabs(beta1 - (double)methyl);
        //sum2 += fabs(beta2 - (double)methyl);
        //printf("   %f\t%f\t%d\n", beta1, beta2, methyl);
        // Euclidean distance
        //double diff1 = beta1 - (double)methyl;
        //double diff2 = beta2 - (double)methyl;
        //sum1 += diff1 * diff1;
        //sum2 += diff2 * diff2;
        // Manhattan distance
        //sum1 += fabs(beta1 - (double)methyl);
        //sum2 += fabs(beta2 - (double)methyl);

        // Free memory
        free(values1);
        free(values2);
    }

    //sum1 = sqrt(sum1);
    //sum2 = sqrt(sum2);
    sum1 = sum1 / (double)n1;
    sum2 = sum2 / (double)n2;

    // Assign read to record
    if (sum1 <= sum2)
    //if (sum1 <= sum2 && read->ncpgs > 0)
    {
        return 0;
    }
    else
    {
       return 1;
    }    
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