//=====================================================================================
// Common structs, parameters, functions
// Original: https://github.com/databio/aiarray/tree/master/src
// by Kyle S. Smith
//-------------------------------------------------------------------------------------
#ifndef __READ_INTERVALS_H__
#define __READ_INTERVALS_H__
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "src/labeled_augmented_array.h"


void read_add(labeled_aiarray_t *intervals, bam1_t *aln, char *chrom_id,
              int min_size, int max_size, int paired, int qcfail, int mapq_cutoff, float proportion);

void sam_iter_add(char *samfile_name, labeled_aiarray_t *intervals,
                  int min_size, int max_size, int paired, int qcfail, int mapq_cutoff,
				  float proportion);

#endif
