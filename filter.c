
#include <stdlib.h>
#include <stdio.h>

#include "filter.h"

/* **************************************************************
 *       	Private functions declarations      		*
 * *************************************************************/

int write_bam_by_mismatches_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int max_mismatches);
int write_bam_by_chromosome_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, short int chromosome);
int write_bam_by_length_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int min_length);
int write_bam_by_quality_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int min_quality, int max_quality);
int write_bam_by_distance_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int max_distance);
int write_bam_by_criteria_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, short int chromosome, int min_length, int min_quality, int max_distance);
int num_errors_from_cigar_(uint32_t* cigar, int num_cigar_operations);

/* ******************************************************
 *       	Function implementations      		*
 * ******************************************************/

void filter_bam_by_chromosome(char* bam_input, char* output_directory, short int chromosome) {
    char* bam_output;

    //open BAM file for read
    if (time_flag) {
        start_timer(t1_filter);
    }

    bam_file_t* bam_file_p =  bam_fopen_mode(bam_input, NULL, "r");

    //header for BAM file has been done in the opening
    bam_header_t* bam_header_p = bam_file_p->bam_header_p;

    //open BAM file for write!!!
    char in_shortname[MAX_FULL_PATH_LENGTH];
    char* str_chr_suffix = "chr";

    get_filename_from_path(bam_input, in_shortname);
    bam_output = (char*) calloc((strlen(output_directory) + strlen(in_shortname) + strlen(str_chr_suffix) + 10), sizeof(char));

    sprintf(bam_output, "%s/%s.%s%hi", output_directory, in_shortname, str_chr_suffix, chromosome);
    
    char log_message[200];
    sprintf(log_message, "FILTER-START: filter bam by chromosome %i in output file %s\n", chromosome, bam_output);
    LOG_DEBUG(log_message);
    
    bam_file_t* bam_slice_file_p = bam_fopen_mode(bam_output, bam_header_p, "w");

    if (bam_slice_file_p == NULL) {
        char log_message[200];
        sprintf(log_message, "Error opening file '%.150s' in mode 'r' !!!!!\n", bam_output);
        LOG_FATAL(log_message);
    }

    //write header text to SAM file
    bam_header_write(bam_slice_file_p->bam_fd, bam_header_p);

    //chromosome minus 1 for matching
    chromosome--;
    
    //write bam alignments that meet the filter conditions
    write_bam_by_chromosome_(bam_file_p, bam_slice_file_p, chromosome);    

    //close BAM and SAM files, free bam alignment and bam file object
    bam_fclose(bam_file_p);
    bam_slice_file_p->bam_header_p = NULL; //it has been freed in the former line
    bam_fclose(bam_slice_file_p);
    free(bam_output);
    if (time_flag) {
        stop_timer(t1_filter, t2_filter, filter_time);
    }

    //number_of_batchs = 1, convention value for statistics (not real batch)
    number_of_batchs = 1;
    
    LOG_DEBUG("FILTER-END: filter bam by chromosome ended\n");
}

void filter_bam_by_criteria(char* bam_input, char* output_directory, int max_mismatches, short int chromosome, int min_length, int min_quality, int max_quality, int max_distance) {
    char* bam_output;
    int num_of_filtering_criteria = 0;
    
    //count the number of criteria
    if (max_mismatches != -1) {
        num_of_filtering_criteria++;
    }    
    if (chromosome != -1) {
        num_of_filtering_criteria++;
    }
    if (min_length != -1) {
        num_of_filtering_criteria++;
    }
    if ((min_quality != -1) || (max_quality != -1)) {
        num_of_filtering_criteria++;
    }
    if (max_distance != -1) {
        num_of_filtering_criteria++;
    }    
    
    //open BAM file for read
    if (time_flag) {
        start_timer(t1_filter);
    }

    bam_file_t* bam_file_p =  bam_fopen_mode(bam_input, NULL, "r");

    //header for BAM file has been done in the opening
    bam_header_t* bam_header_p = bam_file_p->bam_header_p;

    //open BAM file for write!!!
    char in_shortname[MAX_FULL_PATH_LENGTH];
    char* str_suffix1;
    
    get_filename_from_path(bam_input, in_shortname);
    bam_output = (char*) calloc((strlen(output_directory) + strlen(in_shortname) + 25), sizeof(char));
    
    if ((num_of_filtering_criteria == 1) && (max_mismatches != -1)) {
        sprintf(bam_output, "%s/%s.%s%i", output_directory, in_shortname, "maxmismatches", max_mismatches);
    } else if ((num_of_filtering_criteria == 1) && (chromosome != -1)) {
        sprintf(bam_output, "%s/%s.%s%hi", output_directory, in_shortname, "chr", chromosome);
    } else if ((num_of_filtering_criteria == 1) && (min_length != -1)) {
        sprintf(bam_output, "%s/%s.%s%i", output_directory, in_shortname, "minlength", min_length);
    } else if ((num_of_filtering_criteria == 1) && (min_quality != -1) && (max_quality != -1)) {
        sprintf(bam_output, "%s/%s.%s[%i,%i]", output_directory, in_shortname, "qualityrange", min_quality, max_quality);
    } else if ((num_of_filtering_criteria == 1) && (min_quality != -1)) {
        sprintf(bam_output, "%s/%s.%s%i", output_directory, in_shortname, "minquality", min_quality);
    } else if ((num_of_filtering_criteria == 1) && (max_quality != -1)) {
        sprintf(bam_output, "%s/%s.%s%i", output_directory, in_shortname, "maxquality", max_quality);
    } else if ((num_of_filtering_criteria == 1) && (max_distance != -1)) {
        sprintf(bam_output, "%s/%s.%s%i", output_directory, in_shortname, "maxdistance", max_distance);
    } else if (num_of_filtering_criteria > 1) {        
        sprintf(bam_output, "%s/%s.%s", output_directory, in_shortname, "filtered");
    }
    
    char log_message[200];
    sprintf(log_message, "FILTER-START: filter bam in output file %s\n", bam_output);
    LOG_DEBUG(log_message);
    
    bam_file_t* bam_slice_file_p = bam_fopen_mode(bam_output, bam_header_p, "w");

    if (bam_slice_file_p == NULL) {
        char log_message[200];
        sprintf(log_message, "Error opening file '%.150s' in mode 'r' !!!!!\n", bam_output);
        LOG_FATAL(log_message);
    }

    //write header text to SAM file
    bam_header_write(bam_slice_file_p->bam_fd, bam_header_p);
    
    //if exists a filter by chromosome it must be decreased by 1
    if (chromosome != -1) {
        chromosome--;  //chromosome minus 1 for matching
    }
    
    //write bam alignments that meet the filter conditions
    if ((num_of_filtering_criteria == 1) && (max_mismatches != -1)) {
        write_bam_by_mismatches_(bam_file_p, bam_slice_file_p, max_mismatches);
    } else if ((num_of_filtering_criteria == 1) && (chromosome != -1)) {
        write_bam_by_chromosome_(bam_file_p, bam_slice_file_p, chromosome);    
    } else if ((num_of_filtering_criteria == 1) && (min_length != -1)) {
        write_bam_by_length_(bam_file_p, bam_slice_file_p, min_length);    
    } else if ((num_of_filtering_criteria == 1) && ((min_quality != -1) || (max_quality != -1))) {
        write_bam_by_quality_(bam_file_p, bam_slice_file_p, min_quality, max_quality);    
    } else if ((num_of_filtering_criteria == 1) && (max_distance != -1)) {
        write_bam_by_distance_(bam_file_p, bam_slice_file_p, max_distance);    
    } else if (num_of_filtering_criteria > 1) {        
        write_bam_by_criteria_(bam_file_p, bam_slice_file_p, chromosome, min_length, min_quality, max_distance);    
    }

    //close BAM and SAM files, free bam alignment and bam file object
    bam_fclose(bam_file_p);
    bam_slice_file_p->bam_header_p = NULL; //it has been freed in the former line
    bam_fclose(bam_slice_file_p);
    free(bam_output);
    if (time_flag) {
        stop_timer(t1_filter, t2_filter, filter_time);
    }

    //number_of_batchs = 1, convention value for statistics (not real batch)
    number_of_batchs = 1;
    
    LOG_DEBUG("FILTER-END: filter bam by chromosome ended\n");
}

/* **************************************************************
 *       	Private functions implementations      		*
 * *************************************************************/

int write_bam_by_mismatches_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int max_mismatches) {
    LOG_DEBUG("FILTER: by number of mismatches\n");
    
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0;

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        if (num_errors_from_cigar_(bam1_cigar(bam_p), bam_p->core.n_cigar) <= max_mismatches) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);

    return num_write_alignments;  
}

int write_bam_by_chromosome_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, short int chromosome) {
    LOG_DEBUG("FILTER: by chromosome\n");
    
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0;

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        if (bam_p->core.tid == chromosome) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);

    return num_write_alignments;
}

int write_bam_by_length_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int min_length) {
    LOG_DEBUG("FILTER: by length\n");
    
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0;

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        if (bam_cigar2qlen(&(bam_p->core), bam1_cigar(bam_p)) >= min_length) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);

    return num_write_alignments;
}

int write_bam_by_quality_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int min_quality, int max_quality) {
    LOG_DEBUG("FILTER: by quality\n");
   
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0;
    
    // if max_quality is not informed max value is assumed 
    if (max_quality == -1) {
        max_quality = 256;  // max value of the quality
    }

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        if ((bam_p->core.qual >= min_quality) && (bam_p->core.qual <= max_quality)) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);
    printf("num_alignments: %i, num_write_alignments: %i\n", num_alignments, num_write_alignments);

    return num_write_alignments;
} 

int write_bam_by_distance_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, int max_distance) {
    LOG_DEBUG("FILTER: by distance between paired ends\n");
    
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0;

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        //if ((bam_p->core.flag & BAM_FPROPER_PAIR) && (abs(bam_p->core.pos - bam_p->core.mpos) <= max_distance)) {
        if ((bam_p->core.flag & BAM_FPROPER_PAIR) && (abs(bam_p->core.isize) <= max_distance)) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);

    return num_write_alignments;
}
        
int write_bam_by_criteria_(bam_file_t* bam_file_p, bam_file_t* bam_slice_file_p, short int chromosome, int min_length, int min_quality, int max_distance) {
    LOG_DEBUG("FILTER: by several criteria\n");
  
    bam1_t* bam_p = bam_init1();
    int read_bytes, num_write_alignments = 0, criteria_failed = 0;

    while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_p)) > 0) {
        if ((chromosome != -1) && (bam_p->core.tid != chromosome)) {
            criteria_failed++;
        }
        if ((min_length != -1) && (bam_cigar2qlen(&(bam_p->core), bam1_cigar(bam_p)) < min_length)) {
            criteria_failed++;
        }
        if ((min_quality != -1) && (bam_p->core.qual < min_quality)) {
            criteria_failed++;
        }
        if ((max_distance != -1) && (abs(bam_p->core.pos - bam_p->core.mpos) > max_distance)) {
            criteria_failed++;
        }
      
        if (!criteria_failed) {
            bam_write1(bam_slice_file_p->bam_fd, bam_p);
            num_write_alignments++;
        }
        
        //bam_destroy1(bam_p); // it was allocated by the sam-tools, we must free it !!
        num_alignments++;
    }

    bam_destroy1(bam_p);

    return num_write_alignments;
}

int num_errors_from_cigar_(uint32_t* cigar_p, int num_cigar_operations) {
    uint32_t cigar_int;
    int num_errors = 0;
    
    for (int i = 0; i < num_cigar_operations; i++) {
        cigar_int = cigar_p[i];

        switch (cigar_int&BAM_CIGAR_MASK) {
            case BAM_CINS:  //Insertion
                 num_errors += (cigar_int >> BAM_CIGAR_SHIFT);
                 break;
            case BAM_CDEL:  //Deletion
                 num_errors += (cigar_int >> BAM_CIGAR_SHIFT);
                 break;
            case BAM_CDIFF:  //Mismatch
                 num_errors += (cigar_int >> BAM_CIGAR_SHIFT);
                 break;
        }
    }    

    return num_errors;
}
