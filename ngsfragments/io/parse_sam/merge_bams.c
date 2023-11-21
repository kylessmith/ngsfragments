//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


void merge_bams(const char *input1_bam_path, const char *input2_bam_path, const char *output_bam_path)
{
    // Open input BAM files
    samFile *input1_bam = sam_open(input1_bam_path, "r");
    if (input1_bam == NULL) {
        fprintf(stderr, "Error opening input BAM file: %s\n", input1_bam_path);
        exit(1);
    }
    samFile *input2_bam = sam_open(input2_bam_path, "r");
    if (input2_bam == NULL) {
        fprintf(stderr, "Error opening input BAM file: %s\n", input2_bam_path);
        sam_close(input1_bam);
        exit(1);
    }

    // Open output BAM file
    samFile *output_bam = sam_open(output_bam_path, "wb");
    if (output_bam == NULL) {
        fprintf(stderr, "Error opening output BAM file: %s\n", output_bam_path);
        sam_close(input1_bam);
        sam_close(input2_bam);
        exit(1);
    }

    // Merge BAM files
    bam_hdr_t *input1_header = sam_hdr_read(input1_bam);
    bam_hdr_t *input2_header = sam_hdr_read(input2_bam);
    bam_hdr_t *output_header = bam_hdr_dup(input1_header);
    /*if (bam_hdr_combine(output_header, input2_header) != 0) {
        fprintf(stderr, "Error combining BAM headers\n");
        bam_hdr_destroy(input1_header);
        bam_hdr_destroy(input2_header);
        bam_hdr_destroy(output_header);
        sam_close(input1_bam);
        sam_close(input2_bam);
        sam_close(output_bam);
        exit(1);
    }*/

    // Write header to output BAM file
    if (sam_hdr_write(output_bam, output_header) < 0)
    {
        fprintf(stderr, "Error writing BAM header to output BAM file: %s\n", output_bam_path);
        exit(1);
    }

    bam1_t *read = bam_init1();
    while (sam_read1(input1_bam, input1_header, read) >= 0) {
        sam_write1(output_bam, output_header, read);
    }
    while (sam_read1(input2_bam, input2_header, read) >= 0) {
        sam_write1(output_bam, output_header, read);
    }
    bam_destroy1(read);

    // Close BAM files
    bam_hdr_destroy(input1_header);
    bam_hdr_destroy(input2_header);
    bam_hdr_destroy(output_header);
    sam_close(input1_bam);
    sam_close(input2_bam);
    sam_close(output_bam);

    return;
}