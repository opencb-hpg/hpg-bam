
#ifndef CONVERT_C
#define CONVERT_C

#include <stdlib.h>
#include <stdio.h>

#include "convert.h"

/* ******************************************************
 *       	Function implementations      		*
 * ******************************************************/

void convert_bam_to_sam(char* bam_input, char* sam_input) {
    int read_bytes;
    bam1_t* bam_p = bam_init1();
    char* bam_string;

    LOG_DEBUG("CONVERT-START: bam to sam\n");

    //open BAM file for read
    if (time_flag) {
        start_timer(t1_convert);
    }
    bam_file_t* bam_file_p =  bam_fopen_mode(bam_input, NULL, "r");

    //open SAM file for write, SAM file is a text file!!!
    FILE* sam_fd = fopen(sam_input, "w");

    if (sam_fd == NULL) {
        char log_message[200];
        sprintf(log_message, "Error opening file '%.150s' in mode 'r' !!!!!\n", sam_input);
        LOG_FATAL(log_message);
    }

    //header for BAM file has been done in the opening
    bam_header_t* bam_header_p = bam_file_p->bam_header_p;

    //write header text to SAM file
    fprintf(sam_fd, "%s", bam_header_p->text);

    //write string alignments to SAM file
    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        bam_string = bam_format1(bam_header_p, bam_p);
        fprintf(sam_fd, "%s\n", bam_string);
        free(bam_string); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    //close BAM and SAM files, free bam alignment and bam file object
    bam_fclose(bam_file_p);
    fclose(sam_fd);
    bam_destroy1(bam_p);
    if (time_flag) {
        stop_timer(t1_convert, t2_convert, convert_time);
    }

    //number_of_batchs = 1, convention value for statistics (not real batch)
    number_of_batchs = 1;

    LOG_DEBUG("CONVERT-START: bam to sam\n");
}

void convert_sam_to_bam(char* sam_input, char* bam_input) {
    bam1_t* bam_p = bam_init1();

    LOG_DEBUG("CONVERT-START: sam to bam\n");

    //open SAM file for read
    if (time_flag) {
        start_timer(t1_convert);
    }
    tamFile sam_fd = sam_open(sam_input);

    //open BAM file for write
    bam_file_t* bam_file_p =  bam_fopen_mode(bam_input, NULL, "w");

    //read header from SAM file
    bam_header_t* bam_header_p = sam_header_read(sam_fd);

    //write header to BAM file
    bam_header_write(bam_file_p->bam_fd, bam_header_p);

    //write alignments to BAM file
    while (sam_read1(sam_fd, bam_header_p, bam_p) > 0) {
        bam_write1(bam_file_p->bam_fd, bam_p);
        num_alignments++;
    }

    //close BAM and SAM files, free bam alignment and bam file object
    bam_fclose(bam_file_p);
    sam_close(sam_fd);
    bam_header_destroy(bam_header_p);
    bam_destroy1(bam_p);
    if (time_flag) {
        stop_timer(t1_convert, t2_convert, convert_time);
    }

    //number_of_batchs = 1, convention value for statistics (not real batch)
    number_of_batchs = 1;
}

#endif  /* CONVERT_C */
