//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


read_name_sets_t *init_read_name_sets()
{   /* Initialize a read name sets struct */
    
    read_name_sets_t *read_name_sets = malloc(sizeof(read_name_sets_t));
    if (read_name_sets == NULL)
    {
        fprintf(stderr, "Error allocating memory for read name sets\n");
        exit(1);
    }
    read_name_sets->set1 = kh_init(read_name_set);
    read_name_sets->set2 = kh_init(read_name_set);
    return read_name_sets;
}


void destroy_read_name_sets(read_name_sets_t *read_name_sets)
{   /* Destroy a read name sets struct */

    khint_t k;
    for (k = 0; k < kh_end(read_name_sets->set1); ++k)
        if (kh_exist(read_name_sets->set1, k))
            free((char*)kh_key(read_name_sets->set1, k));
    kh_destroy(read_name_set, read_name_sets->set1);

    for (k = 0; k < kh_end(read_name_sets->set2); ++k)
        if (kh_exist(read_name_sets->set2, k))
            free((char*)kh_key(read_name_sets->set2, k));
    kh_destroy(read_name_set, read_name_sets->set2);

    free(read_name_sets);
}