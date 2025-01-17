//=============================================================================
// Store fragment information
// by Kyle S. Smith
//
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "read_fragments.h"

#define BUFFER_SIZE 4096

//-----------------------------------------------------------------------------


void process_line(char *line)
{
    // Split the line by tabs and process each field
    char *token = strtok(line, "\t");
    while (token != NULL)
    {
        printf("%s\n", token); // Replace this with your processing logic
        token = strtok(NULL, "\t");
    }
}

void iterate_gzipped_file(const char *filename)
{
    gzFile file = gzopen(filename, "rb");
    if (!file)
    {
        perror("Could not open gzipped file");
        return;
    }

    char buffer[BUFFER_SIZE];
    while (gzgets(file, buffer, sizeof(buffer)) != NULL)
    {
        // Skip lines starting with '#'
        if (buffer[0] == '#')
        {
            continue;
        }
        
        // Remove the newline character at the end of the line, if present
        buffer[strcspn(buffer, "\n")] = 0;

        // Process the line
        process_line(buffer);
    }

    gzclose(file);
}

