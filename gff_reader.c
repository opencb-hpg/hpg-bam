
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "commons.h"
#include "gff_reader.h"
#include "log.h"

unsigned int gff_file_read(char* filename, gff_line_t* gff_lines_p) {
    unsigned int count = 0;
    FILE* file = fopen(filename, "r");

    if (file == NULL) {
        LOG_ERROR("Error opening gff file");
    }

    char* buffer = (char*) calloc(1, MAX_GFF_FILE_LINE_LENGTH);
    const char delimiters[] = "\t";
    char* token = NULL;
    int pos;

    while (fgets(buffer, MAX_GFF_FILE_LINE_LENGTH, file)) {
        if ((buffer[0] == '#') || (buffer[0] != 'c') || (buffer[1] != 'h') || (buffer[2] != 'r')) continue;

        pos = 0;
        token = strtok(buffer, delimiters);

        while (token != NULL) {
            switch (pos) {
                case 0:
                    gff_lines_p[count].seqname = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].seqname, token);
                    break;
                case 1:
                    gff_lines_p[count].source = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].source, token);
                    break;
                case 2:
                    gff_lines_p[count].feature = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].feature, token);
                    break;
                case 3:
                    sscanf(token, "%i", &(gff_lines_p[count].start));
                    break;
                case 4:
                    sscanf(token, "%i", &(gff_lines_p[count].end));
                    break;
                case 5:
                    gff_lines_p[count].score = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].score, token);
                    break;
                case 6:
                    gff_lines_p[count].strand = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].strand, token);
                    break;
                case 7:
                    gff_lines_p[count].frame = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].frame, token);
                    break;
                case 8:
                    gff_lines_p[count].group = (char*) calloc(strlen(token), sizeof(char));
                    strcpy(gff_lines_p[count].group, token);
                    break;
                default:
                    break;
            }

            token = strtok(NULL, delimiters);
            pos++;
        }

        count++;
    }

    fclose(file);
    free(buffer);

    return count;
}
