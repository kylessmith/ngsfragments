//=============================================================================
// Store BAM fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <string.h>
#include "read_intervals.h"

//-----------------------------------------------------------------------------


methyl_record_pair_t *methyl_length_split(const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size1,
                                            int max_size1,
                                            int min_size2,
                                            int max_size2,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads)
{   /* Split reads by length and quantify methylation */
    
    // Initialize reference CpGs
    //printf("Initializing reference methyl record\n");
    methyl_record_t *record1 = fetch_reference_methyl_record(ref_2bit, chromosome);
    methyl_record_t *record2 = fetch_reference_methyl_record(ref_2bit, chromosome);
    //printf("done\n");

    // Initialize read iterator
    int min_size = max_size2 * -1;
    int max_size = max_size2;
    methyl_read_iterator_t *methyl_iter = methyl_read_iterator_init(bam_file_path,
                                                                    ref_2bit,
                                                                    chromosome,
                                                                    min_size,
                                                                    max_size,
                                                                    qcfail,
                                                                    mapq_cutoff,
                                                                    proportion,
                                                                    nthreads);

    // Iterate through reads
    methyl_read_t *read;
    while (methyl_read_iterator_next(methyl_iter) != 0)
    {
        read = methyl_iter->methyl_pair;
        /*if (read->length >= 1 && read->length <= 135)
        {
            methyl_record_add(record1, read->pos, read->methyl, read->size);
        }
        else if (read->length >= 240 && read->length <= 324)
        {
            methyl_record_add(record1, read->pos, read->methyl, read->size);
        }
        else if (read->length >= 136 && read->length <= 239)
        {
            methyl_record_add(record2, read->pos, read->methyl, read->size);
        }
        else if (read->length >= 325 && read->length <= 1000)
        {
            methyl_record_add(record2, read->pos, read->methyl, read->size);
        }*/

        if (read->length >= min_size1 && read->length <= max_size1)
        {
            methyl_record_add(record1, read->pos, read->methyl, read->size);
        }
        else if (read->length >= min_size2 && read->length <= max_size2)
        {
            methyl_record_add(record2, read->pos, read->methyl, read->size);
        }
    }

    //printf("Created methyl records\n");
    //int i;
    //for (i = 0; i < 10000; i++)
    //{
    //    printf("%d\t%d\t%d\t%d\t%d\n", i, record1->methyl[i], record1->unmethyl[i], record2->methyl[i], record2->unmethyl[i]);
    //}

    // Clean up
    methyl_read_iterator_destroy(methyl_iter);

    // Create pair
    methyl_record_pair_t *pair = methyl_record_pair_init(record1, record2);

    return pair;
}


methyl_record_pair_t *methyl_profile_split(methyl_record_pair_t *pair,
                                            const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size,
                                            int max_size,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads)
{   /* Split reads based on similarity to other methylation profiles */
    
    //printf("Initializing reference methyl record\n");
    //methyl_record_t *record1 = fetch_reference_methyl_record(ref_2bit, chromosome);
    //methyl_record_t *record2 = fetch_reference_methyl_record(ref_2bit, chromosome);

    // Initialize read iterator
    //printf("Initializing read iterator\n");
    //printf("   min_size: %d\n", min_size);
    //printf("   max_size: %d\n", max_size);
    methyl_read_iterator_t *methyl_iter = methyl_read_iterator_init(bam_file_path,
                                                                    ref_2bit,
                                                                    chromosome,
                                                                    min_size,
                                                                    max_size,
                                                                    qcfail,
                                                                    mapq_cutoff,
                                                                    proportion,
                                                                    nthreads);

    // Iterate through reads
    methyl_read_t *read;
    while (methyl_read_iterator_next(methyl_iter) != 0)
    {
        read = methyl_iter->methyl_pair;
        int assignment = assign_methyl_read(pair, read);
        if (assignment == 0)
        {
            //methyl_record_add(record1, read->pos, read->methyl, read->size);
            methyl_record_add(pair->record1, read->pos, read->methyl, read->size);
        }
        else
        {
            methyl_record_add(pair->record2, read->pos, read->methyl, read->size);
        }
    }

    // Clean up
    methyl_read_iterator_destroy(methyl_iter);
    //methyl_record_pair_destroy(pair);

    // Create pair
    //methyl_record_pair_t *new_pair = methyl_record_pair_init(record1, record2);

    return pair;
}


methyl_record_pair_t *methyl_profile_split_new(methyl_record_pair_t *pair,
                                                const char *bam_file_path,
                                                char *ref_2bit,
                                                const char *chromosome,
                                                int min_size,
                                                int max_size,
                                                int qcfail,
                                                int mapq_cutoff,
                                                float proportion,
                                                int nthreads)
{   /* Split reads based on similarity to other methylation profiles */
    
    //printf("Initializing reference methyl record\n");
    methyl_record_t *record1 = fetch_reference_methyl_record(ref_2bit, chromosome);
    methyl_record_t *record2 = fetch_reference_methyl_record(ref_2bit, chromosome);

    // Initialize read iterator
    //printf("Initializing read iterator\n");
    //printf("   min_size: %d\n", min_size);
    //printf("   max_size: %d\n", max_size);
    methyl_read_iterator_t *methyl_iter = methyl_read_iterator_init(bam_file_path,
                                                                    ref_2bit,
                                                                    chromosome,
                                                                    min_size,
                                                                    max_size,
                                                                    qcfail,
                                                                    mapq_cutoff,
                                                                    proportion,
                                                                    nthreads);

    // Iterate through reads
    methyl_read_t *read;
    while (methyl_read_iterator_next(methyl_iter) != 0)
    {
        read = methyl_iter->methyl_pair;
        int assignment = assign_methyl_read(pair, read);
        if (assignment == 0)
        {
            methyl_record_add(record1, read->pos, read->methyl, read->size);
            //methyl_record_add(pair->record1, read->pos, read->methyl, read->size);
        }
        else
        {
            //methyl_record_add(pair->record2, read->pos, read->methyl, read->size);
            methyl_record_add(record2, read->pos, read->methyl, read->size);
        }
    }

    // Clean up
    methyl_read_iterator_destroy(methyl_iter);
    methyl_record_pair_destroy(pair);

    // Create pair
    methyl_record_pair_t *new_pair = methyl_record_pair_init(record1, record2);

    return new_pair;
}


read_name_sets_t *methyl_profile_split_names(methyl_record_pair_t *pair,
                                            const char *bam_file_path,
                                            char *ref_2bit,
                                            const char *chromosome,
                                            int min_size,
                                            int max_size,
                                            int qcfail,
                                            int mapq_cutoff,
                                            float proportion,
                                            int nthreads)
{   /* Split reads based on similarity to other methylation profiles and return read names */
    
    // Initialize read iterator
    methyl_read_iterator_t *methyl_iter = methyl_read_iterator_init(bam_file_path,
                                                                    ref_2bit,
                                                                    chromosome,
                                                                    min_size,
                                                                    max_size,
                                                                    qcfail,
                                                                    mapq_cutoff,
                                                                    proportion,
                                                                    nthreads);

    // Initialize hash table for storing read names
    read_name_sets_t *read_names = init_read_name_sets();
    
    // Iterate through reads
    methyl_read_t *read;
    while (methyl_read_iterator_next(methyl_iter) != 0)
    {
        read = methyl_iter->methyl_pair;
        //char *read_name = strdup(read->name);
        char *read_name = strdup(read->name);
        int assignment = assign_methyl_read(pair, read);
        if (assignment == 0)
        {
            int absent;
            khint_t k = kh_put(read_name_set, read_names->set1, read_name, &absent);
            //if (absent)
            //{
            //    fprintf(stderr, "Error storing read name: %s\n", read_name);
            //    exit(1);
            //}
        }
        else
        {
            int absent;
            khint_t k = kh_put(read_name_set, read_names->set2, read_name, &absent);
            //if (absent)
            //{
            //    fprintf(stderr, "Error storing read name: %s\n", read_name);
            //    exit(1);
            //}
        }
    }
    //printf("found %d in set1\n", kh_end(read_names->set1));
    //printf("found %d in set2\n", kh_end(read_names->set2));

    // Clean up
    methyl_read_iterator_destroy(methyl_iter);

    return read_names;
}


void write_split_reads(const char *bam_file_path,
                        const char *output_bam_file_path1,
                        const char *output_bam_file_path2,
                        read_name_sets_t *read_names,
                        char *chromosome,
                        int min_size,
                        int max_size,
                        int qcfail,
                        int mapq_cutoff,
                        int nthreads)
{   /* Write reads to output BAM files based on read name */

    // Open input BAM file
    //samFile *bam_file = sam_open(bam_file_path, "r");
    //if (bam_file == NULL)
    //{
    //    fprintf(stderr, "Error opening BAM file: %s\n", bam_file_path);
    //    exit(1);
    //}

    // Create read iterator
    read_iter_t *read_iter = read_iter_init(bam_file_path,
                                            chromosome,
                                            min_size,
                                            max_size,
                                            1,
                                            qcfail,
                                            mapq_cutoff,
                                            1.0,
                                            nthreads);

    // Open output BAM file
    samFile *output_bam_file1 = sam_open(output_bam_file_path1, "wb");
    if (output_bam_file1 == NULL)
    {
        fprintf(stderr, "Error opening output BAM file: %s\n", output_bam_file_path1);
        exit(1);
    }
    samFile *output_bam_file2 = sam_open(output_bam_file_path2, "wb");
    if (output_bam_file2 == NULL)
    {
        fprintf(stderr, "Error opening output BAM file: %s\n", output_bam_file_path2);
        exit(1);
    }

    // Copy header from input BAM file to output BAM file
    //bam_hdr_t *bam_header = sam_hdr_read(bam_file);
    if (sam_hdr_write(output_bam_file1, read_iter->header) < 0)
    {
        fprintf(stderr, "Error writing BAM header to output BAM file: %s\n", output_bam_file_path1);
        exit(1);
    }
    if (sam_hdr_write(output_bam_file2, read_iter->header) < 0)
    {
        fprintf(stderr, "Error writing BAM header to output BAM file: %s\n", output_bam_file_path2);
        exit(1);
    }

    // Initialize BAM alignment
    //bam1_t *bam_alignment = bam_init1();

    // Iterate through input BAM file and write matching reads to output BAM file
    //printf("Iterating over BAM file\n");
    while (read_iter_next(read_iter) >= 1)
    {
        char *read_name = bam_get_qname(read_iter->aln);

        khint_t k = kh_get(read_name_set, read_names->set1, read_name);
        if (k != kh_end(read_names->set1))
        {
            if (sam_write1(output_bam_file1, read_iter->header, read_iter->aln) < 0)
            {
                fprintf(stderr, "Error writing BAM alignment to output BAM file: %s\n", output_bam_file_path1);
                exit(1);
            }
        } else {
            k = kh_get(read_name_set, read_names->set2, read_name);
            if (k != kh_end(read_names->set2))
            {
                if (sam_write1(output_bam_file2, read_iter->header, read_iter->aln) < 0)
                {
                    fprintf(stderr, "Error writing BAM alignment to output BAM file: %s\n", output_bam_file_path2);
                    exit(1);
                }
            }
        }
    }

    // Clean up
    //bam_destroy1(bam_alignment);
    //bam_hdr_destroy(bam_header);
    //sam_close(bam_file);
    read_iter_destroy(read_iter);
    sam_close(output_bam_file1);
    sam_close(output_bam_file2);

    return;
}


void split_methyl_bam(const char *bam_file_path,
                        const char *output_bam_file_path1,
                        const char *output_bam_file_path2,
                        char *ref_2bit,
                        const char *chromosome,
                        int min_size1,
                        int max_size1,
                        int min_size2,
                        int max_size2,
                        int min_distance,
                        int tolerance,
                        int max_iter,
                        int qcfail,
                        int mapq_cutoff,
                        float proportion,
                        int nthreads)
{   /* Split reads by length/profile and write to output BAM files */
    
    printf("Splitting reads: %s\n", chromosome);

    // Split reads by length
    int min_size = max_size2 * -1;
    int max_size = max_size2;
	//printf("Spliting on lengths\n");
    printf("   splitting lengths %s\n", chromosome);
    methyl_record_pair_t *pair = methyl_length_split(bam_file_path,
                                                        ref_2bit,
                                                        chromosome,
                                                        min_size1,
                                                        max_size1,
                                                        min_size2,
                                                        max_size2,
                                                        qcfail,
                                                        mapq_cutoff,
                                                        proportion,
                                                        nthreads);

    // Determine similarity between methylation profiles
	//printf("dones\n");
	//printf("Checking distances\n");
    double current_distance = compare_methyl_records(pair);
    double previous_distance = current_distance + tolerance + 1;
	//printf("Done\n");

    //printf("Current distance: %f\n", current_distance);
    //printf("Previous distance: %f\n", previous_distance);
    //printf("min_distance: %d\n", min_distance);
    //methyl_record_pair_write(pair, "test.txt");
    //exit(1);

    // Split reads by similarity to other methylation profiles
    //printf("Spliting on profiles\n");
    printf("   splitting profiles %s\n", chromosome);
    //methyl_record_pair_transfer_null(pair);
    int n = 0;
    while (fabs(previous_distance - current_distance) > tolerance && n < max_iter)
    {
        //printf("   splitting\n");
        previous_distance = current_distance;
        pair = methyl_profile_split_new(pair,
                                    bam_file_path,
                                    ref_2bit,
                                    chromosome,
                                    min_size,
                                    max_size,
                                    qcfail,
                                    mapq_cutoff,
                                    proportion,
                                    nthreads);
        //methyl_record_pair_transfer_null(pair);
        current_distance = compare_methyl_records(pair);
        //printf("   Current distance: %f\n", current_distance);
        //printf("   Previous distance: %f\n", previous_distance);
        n++;
    }
    //printf("done\n");

    //printf("Done with all splitting\n");
    //int k;
    //for (k = 0; k < 600; k++)
    //{
    //    printf("%d\t%d\t%d\t%d\t%d\n", k, pair->record1->methyl[k], pair->record1->unmethyl[k], pair->record2->methyl[k], pair->record2->unmethyl[k]);
    //}

    // Split reads by similarity to other methylation profiles and return read names
    //printf("Spliting on profiles and returning read names\n");
    
    printf("   splitting names %s\n", chromosome);
    read_name_sets_t *read_names = methyl_profile_split_names(pair,
                                                                bam_file_path,
                                                                ref_2bit,
                                                                chromosome,
                                                                min_size,
                                                                max_size,
                                                                qcfail,
                                                                mapq_cutoff,
                                                                proportion,
                                                                nthreads);
    //printf("done\n");

    // Write reads to output BAM files based on read name
    //printf("Writing reads to output BAM files\n");
    printf("    writing %s\n", chromosome);
    write_split_reads(bam_file_path,
                        output_bam_file_path1,
                        output_bam_file_path2,
                        read_names,
                        chromosome,
                        min_size,
                        max_size,
                        qcfail,
                        mapq_cutoff,
                        nthreads);
    //printf("done\n");

    // Clean up
    destroy_read_name_sets(read_names);
    methyl_record_pair_destroy(pair);

    printf("DONE %s\n", chromosome);

    return;
}


void test_memory_leak(const char *bam_file_path,
                        const char *output_bam_file_path1,
                        const char *output_bam_file_path2,
                        char *ref_2bit,
                        const char *chromosome,
                        int min_size1,
                        int max_size1,
                        int min_size2,
                        int max_size2,
                        int min_distance,
                        int tolerance,
                        int qcfail,
                        int mapq_cutoff,
                        float proportion,
                        int nthreads)
{
    // Iterate over methyl fragments
    int min_size = max_size2 * -1;
    int max_size = max_size2;
    methyl_read_iterator_t *iter = methyl_read_iterator_init(bam_file_path,
                                                            ref_2bit,
                                                            chromosome,
                                                            min_size,
                                                            max_size,
                                                            qcfail,
                                                            mapq_cutoff,
                                                            proportion,
                                                            nthreads);

    while (methyl_read_iterator_next(iter) != 0)
    {
        
        //printf("%s\n", iter->methyl_pair->name);
    }

    // Clean up
    methyl_read_iterator_destroy(iter);

    return;
}