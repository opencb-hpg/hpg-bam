
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "bam.h"
#include "bam_coverage.h"
#include "gff_data.h"

#define LINES_PER_COVERAGE_PRINT_BUFFER 	 7
#define COVERAGE_LINE_LENGTH   			21

/* ******************************************************
 *      		Private functions     		*
 * *****************************************************/

void bam_coverage_compute_alignment_in_region_(int chromosome, int start_coordinate, uint32_t* cigar_data_p, int cigar_index, int num_cigar_operations, gff_region_t* region_p, bam_chromosome_coverage_t* bam_chromosome_coverage_p);
void bam_coverage_compute_alignment_(int chromosome, int start_coordinate, uint32_t* cigar_data_p, int cigar_index, int num_cigar_operations, bam_chromosome_coverage_t* bam_chromosome_coverage_p);
void bam_coverage_counter_incr_between_(int interval_start, int interval_end, bam_coverage_counter_t* counter_p);
char* num_chromosome_to_string_(int chromosome);
int nt_coordinate_to_string_(char* str_coordinate, int coordinate);
int nt_coverage_to_string_(char* str_coverage, int coverage);

/* ******************************************************
 *      		Global variables     		*
 * *****************************************************/

//char* str_chromosomes[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25"};
//int strlen_chromosomes[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
char* str_chromosomes[] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
                            "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", 
                            "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", 
                            "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", 
                            "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", 
                            "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", 
                            "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", 
                            "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", 
                            "91", "92", "93", "94", "95", "96", "97", "98", "99", "100" };
int strlen_chromosomes[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                             2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

/* ******************************************************
 *      	Function implementations    		*
 * *****************************************************/

/* ******************************************************
 *      	Coverage counters functions     	*
 * *****************************************************/

void bam_chromosome_coverage_init(bam_chromosome_coverage_t* bam_chromosome_coverage_p) {
    bam_chromosome_coverage_p->bam_coverage_counter_p = (bam_coverage_counter_t**) calloc((MAX_NTS_PER_CHROMOSOME / NTS_PER_COUNTER), sizeof(bam_coverage_counter_t*));
}

void bam_chromosome_coverage_clear(bam_chromosome_coverage_t* bam_chromosome_coverage_p) {
    for (int i = 0; i < (MAX_NTS_PER_CHROMOSOME / NTS_PER_COUNTER) ; i++) {
        free(bam_chromosome_coverage_p->bam_coverage_counter_p[i]);
    }
    free(bam_chromosome_coverage_p->bam_coverage_counter_p);
}

void bam_coverage_counter_mark_to_print(bam_chromosome_coverage_t* bam_chromosome_coverage_p, int all) {
    bam_coverage_counter_t** counter_p;
    int count_not_null = 0;

    for (int i = (num_of_chromosomes - 1); i >= 0; i--) {  
        counter_p = bam_chromosome_coverage_p[i].bam_coverage_counter_p;

        for (int j = ((MAX_NTS_PER_CHROMOSOME / NTS_PER_COUNTER) - 1); j >= 0; j--) {
            if (counter_p[j] != NULL) {
                count_not_null++;

                if ((count_not_null > 1) || (all)) {
                    counter_p[j]->print = 1;
                }
            }
        }
    }
}

char coverage_lines[1000000];
char* p = coverage_lines;
int chunk = 750000;

char* str_coverage_matrix[65536];
short int strlen_coverage_matrix[65536];

void bam_coverage_counter_print(bam_chromosome_coverage_t* bam_chromosome_coverage_p, char* output_directory, char* input_filename) {
    char coverage_filename[MAX_FULL_PATH_LENGTH];
    char in_shortname[MAX_FULL_PATH_LENGTH];

    get_filename_from_path(input_filename, in_shortname);
    sprintf(coverage_filename, "%s/%s%s", output_directory, in_shortname, COVERAGE_FILE_SUFFIX);

    char log_message[200];
    sprintf(log_message, "COVERAGE FILENAME: %.170s\n", coverage_filename);
    LOG_DEBUG(log_message);

    FILE* fd = (coverage_filename == NULL ? stdout : fopen(coverage_filename, "a"));

    int counter_start_position, strlen_chromosome, strlen_coordinate;    
    bam_coverage_counter_t** bam_coverage_counter_p;
    unsigned short int* coverage_counter_p;

    char* str_chromosome;
    char* str_coordinate = (char*) calloc(9, sizeof(char));
    char* str_coverage = (char*) calloc(5, sizeof(char));

    if (time_flag) {
        stop_timer(t1_cpu, t2_cpu, cpu_time);
    }
    if (time_flag) {
        start_timer(t1_write);
    }

    for (int k = 0; k < num_of_chromosomes; k++) {
        bam_coverage_counter_p = bam_chromosome_coverage_p[k].bam_coverage_counter_p;

        for (int i = 0; i < (MAX_NTS_PER_CHROMOSOME / NTS_PER_COUNTER); i++) {
            if ((bam_coverage_counter_p[i] != NULL) && (bam_coverage_counter_p[i]->print)) {
                // printf referenced to chromosome 1 (considered in hash definition) and position 1
                str_chromosome = str_chromosomes[k];
                strlen_chromosome = strlen_chromosomes[k];

                counter_start_position = i * NTS_PER_COUNTER + 1;
                coverage_counter_p = bam_coverage_counter_p[i]->coverage_counter;

                for (int j = 0; j < NTS_PER_COUNTER; j++) {
                    // condition for printing coverage, could be > 0 or not
                    if (coverage_counter_p[j] > 0) {
                        strlen_coordinate = nt_coordinate_to_string_(str_coordinate, counter_start_position + j);

                        memcpy(p, str_chromosome, strlen_chromosome);
                        p += strlen_chromosome;
                        memcpy(p, "\t", 1);
                        p += 1;
                        memcpy(p, str_coordinate, strlen_coordinate);
                        p += strlen_coordinate;
                        memcpy(p, "\t", 1);
                        p += 1;
                        memcpy(p, str_coverage_matrix[coverage_counter_p[j]], strlen_coverage_matrix[coverage_counter_p[j]]);
                        p += strlen_coverage_matrix[coverage_counter_p[j]];
                        memcpy(p, "\n", 1);
                        p += 1;

                        if ((p - coverage_lines) > chunk) {
                            *p = '\0';
                            fwrite(coverage_lines, 1, p - coverage_lines, fd);
                            p = coverage_lines;
                        }

                    } // end of if coverage_counter_p

                    // condition for computing mean coverage, compute always if > 0
                    if (coverage_counter_p[j] > 0) {
                        nts_with_coverage++;
                        mean_coverage += coverage_counter_p[j];
                    }
                }

                free(bam_coverage_counter_p[i]);
                bam_coverage_counter_p[i] = NULL;
            }
        }
    }

    if (p != coverage_lines) {
        *p = '\0';
        fwrite(coverage_lines, 1, p - coverage_lines, fd);
        p = coverage_lines;
    }

    if (time_flag) {
        stop_timer(t1_write, t2_write, write_time);
    }
    if (time_flag) {
        start_timer(t1_cpu);
    }

    free(str_coordinate);
    free(str_coverage);
    fclose(fd);
}

void bam_coverage_counter_print_block(bam_chromosome_coverage_t* bam_chromosome_coverage_p, char* output_directory, char* input_filename) {
    char coverage_filename[MAX_FULL_PATH_LENGTH];
    char in_shortname[MAX_FULL_PATH_LENGTH];

    get_filename_from_path(input_filename, in_shortname);
    sprintf(coverage_filename, "%s/%s%s", output_directory, in_shortname, COVERAGE_FILE_SUFFIX);

    char log_message[200];
    sprintf(log_message, "COVERAGE FILENAME: %.170s\n", coverage_filename);
    LOG_DEBUG(log_message);

    FILE* fd = (coverage_filename == NULL ? stdout : fopen(coverage_filename, "a"));

    int chromosome;
    bam_coverage_counter_t** bam_coverage_counter_p;
    unsigned short int* coverage_counter_p;
    int counter_start_position;

    char* coverage_line = (char*) calloc(1, COVERAGE_LINE_LENGTH * sizeof(char));
    char* coverage_lines_block = (char*) calloc(1, COVERAGE_LINE_LENGTH * LINES_PER_COVERAGE_PRINT_BUFFER * sizeof(char));

    if (time_flag) {
        start_timer(t1_write);
    }
    
    for (int k = 0; k < num_of_chromosomes; k++) {
        bam_coverage_counter_p = bam_chromosome_coverage_p[k].bam_coverage_counter_p;

        for (int i = 0; i < (MAX_NTS_PER_CHROMOSOME / NTS_PER_COUNTER); i++) {
            if ((bam_coverage_counter_p[i] != NULL) && (bam_coverage_counter_p[i]->print)) {
                // printf referenced to chromosome 1 and position 1
                chromosome = k + 1;
                counter_start_position = i * NTS_PER_COUNTER + 1;
                coverage_counter_p = bam_coverage_counter_p[i]->coverage_counter;

                for (int j = 0; j < NTS_PER_COUNTER; j++) {

                    // condition for printing coverage, could be > 0 or not
                    if (coverage_counter_p[j] > 0) {
                        sprintf(coverage_lines_block + strlen(coverage_lines_block), "%i\t%i\t%i\n", chromosome, (counter_start_position + j), coverage_counter_p[j]);
                    }

                    if (j % LINES_PER_COVERAGE_PRINT_BUFFER) {
                        fprintf(fd, "%s", coverage_lines_block);
                        memset(coverage_lines_block, 0, COVERAGE_LINE_LENGTH * LINES_PER_COVERAGE_PRINT_BUFFER);
                    }

                    // condition for computing mean coverage, compute always if > 0
                    if (coverage_counter_p[j] > 0) {
                        nts_with_coverage++;
                        mean_coverage += coverage_counter_p[j];
                    }
                }

                fprintf(fd, "%s", coverage_lines_block);

                free(bam_coverage_counter_p[i]);
                bam_coverage_counter_p[i] = NULL;
            }
        }
    }
    if (time_flag) {
        stop_timer(t1_write, t2_write, write_time);
    }

    free(coverage_line);
    free(coverage_lines_block);

    fclose(fd);
}

void bam_coverage_counter_delete_file(char* output_directory, char* input_filename) {
    char coverage_filename[MAX_FULL_PATH_LENGTH];
    char in_shortname[MAX_FULL_PATH_LENGTH];

    get_filename_from_path(input_filename, in_shortname);
    sprintf(coverage_filename, "%s/%s%s", output_directory, in_shortname, COVERAGE_FILE_SUFFIX);
    remove(coverage_filename);
}

/* **************************************************************
 *      	Coverage computation function     		*
 * *************************************************************/

void bam_coverage_compute(bam_data_batch_t* batch_p, bam_chromosome_coverage_t* bam_chromosome_coverage_p, gff_data_t* gff_data_p, char* output_directory, char* input_filename, int cpu_num_threads) {
    int i, chromosome, alignment_start, alignment_max_end, cigar_index, num_cigar_operations, num_regions;
    int regions[100];
    uint32_t* cigar_data_p;

#pragma omp parallel for num_threads(cpu_num_threads) shared(batch_p, bam_chromosome_coverage_p, gff_data_p) private(i, chromosome, alignment_start, alignment_max_end, cigar_index, num_cigar_operations, num_regions, regions, cigar_data_p)
    for (int i = 0; i < batch_p->num_alignments; i++) {
        chromosome = batch_p->core_data_p[i].chromosome;
        alignment_start = batch_p->core_data_p[i].start_coordinate;
        alignment_max_end = batch_p->core_data_p[i].start_coordinate + (2 * batch_p->core_data_p[i].alignment_length - 1);
        cigar_data_p = batch_p->cigar_data_p;
        cigar_index = batch_p->core_data_p[i].cigar_index;
        num_cigar_operations = batch_p->core_data_p[i+1].cigar_index - batch_p->core_data_p[i].cigar_index;

        if (gff_data_p != NULL) {
            num_regions = gff_data_alignment_in_region(gff_data_p, chromosome, alignment_start, alignment_max_end, regions);

            for (int j = 0; j < num_regions; j++) {
                bam_coverage_compute_alignment_in_region_(chromosome, alignment_start, cigar_data_p, cigar_index, num_cigar_operations, &gff_data_p->gff_regions_p[regions[j]], bam_chromosome_coverage_p);
            }
        } else {
            bam_coverage_compute_alignment_(chromosome, alignment_start, cigar_data_p, cigar_index, num_cigar_operations, bam_chromosome_coverage_p);
        }
    }
#pragma omp barrier

    bam_coverage_counter_mark_to_print(bam_chromosome_coverage_p, 0);
    bam_coverage_counter_print(bam_chromosome_coverage_p, output_directory, input_filename);
}

void str_coverage_matrix_init() {
    for (int i = 0; i < 10; i++) {
        strlen_coverage_matrix[i] = 1;
    }

    for (int i = 10; i < 100; i++) {
        strlen_coverage_matrix[i] = 2;
    }

    for (int i = 100; i < 1000; i++) {
        strlen_coverage_matrix[i] = 3;
    }
    
    for (int i = 1000; i < 10000; i++) {
        strlen_coverage_matrix[i] = 4;
    }

    for (int i = 10000; i < 65536; i++) {
        strlen_coverage_matrix[i] = 5;
    }

    for (int i = 0; i < 65536; i++) {
        str_coverage_matrix[i] = (char*) calloc(strlen_coverage_matrix[i], sizeof(char));
        nt_coverage_to_string_(str_coverage_matrix[i], i);
    }
}

/* **************************************************************
 *      	Private function implementations    		*
 * *************************************************************/

void bam_coverage_compute_alignment_in_region_(int chromosome, int start_coordinate, uint32_t* cigar_data_p, int cigar_index, int num_cigar_operations, gff_region_t* region_p, bam_chromosome_coverage_t* bam_chromosome_coverage_p) {
    int region_start, region_end;
    int interval_count_start, interval_count_end, cigar, cigar_num_nts, cigar_operation;
    int fragment_start_coordinate, fragment_end_coordinate, counter_offset_start, counter_offset_end;

    region_start = region_p->start;
    region_end = region_p->end;

    fragment_start_coordinate = start_coordinate;

#pragma omp critical
    for (int i = 0; i < num_cigar_operations; i++) {
        cigar = cigar_data_p[cigar_index++];
        cigar_operation = (cigar & BAM_CIGAR_MASK);
        cigar_num_nts = cigar >> BAM_CIGAR_SHIFT;

        //TO DO: review with complete datasets and adjust conditions

        if ((cigar_operation == BAM_CMATCH) || (cigar_operation == BAM_CEQUAL) || (cigar_operation == BAM_CDIFF)) {
            fragment_end_coordinate = fragment_start_coordinate + cigar_num_nts;
        } else if ((cigar_operation == BAM_CREF_SKIP) || (cigar_operation == BAM_CPAD) || (cigar_operation == BAM_CDEL)) {
            fragment_start_coordinate += cigar_num_nts;
            continue;
        } else if ((cigar_operation == BAM_CINS) || (cigar_operation == BAM_CSOFT_CLIP)) {
            continue;
        }

        if (((fragment_start_coordinate >= region_start) && (fragment_end_coordinate <= region_end)) ||
                ((fragment_end_coordinate >= region_start) || (fragment_start_coordinate <= region_end)) ||
                ((fragment_start_coordinate < region_start) && (fragment_end_coordinate > region_end))) {

            interval_count_start = max(fragment_start_coordinate, region_start);
            interval_count_end = min(fragment_end_coordinate, region_end + 1);

            counter_offset_start = interval_count_start / NTS_PER_COUNTER;
            counter_offset_end = interval_count_end / NTS_PER_COUNTER;

            if (counter_offset_start == counter_offset_end) {
                if (bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] == NULL) bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] = (bam_coverage_counter_t*) calloc(1, sizeof(bam_coverage_counter_t));
                bam_coverage_counter_incr_between_((interval_count_start % NTS_PER_COUNTER), (interval_count_end % NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]);
            } else {
                if (bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end] == NULL) bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end] = (bam_coverage_counter_t*) calloc(1, sizeof(bam_coverage_counter_t));
                bam_coverage_counter_incr_between_((interval_count_start % NTS_PER_COUNTER), (NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]);
                bam_coverage_counter_incr_between_(0, (interval_count_end % NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end]);
            }
        }

        fragment_start_coordinate = fragment_end_coordinate;
    }
}

void bam_coverage_compute_alignment_(int chromosome, int start_coordinate, uint32_t* cigar_data_p, int cigar_index, int num_cigar_operations, bam_chromosome_coverage_t* bam_chromosome_coverage_p) {
    int cigar, cigar_num_nts, cigar_operation;
    int fragment_start_coordinate, fragment_end_coordinate, counter_offset_start, counter_offset_end;

    fragment_start_coordinate = start_coordinate;

    for (int i = 0; i < num_cigar_operations; i++) {
        cigar = cigar_data_p[cigar_index++];
        cigar_operation = (cigar & BAM_CIGAR_MASK);
        cigar_num_nts = cigar >> BAM_CIGAR_SHIFT;

        //TO DO: review with complete datasets and adjust conditions

        if ((cigar_operation == BAM_CMATCH) || (cigar_operation == BAM_CEQUAL) || (cigar_operation == BAM_CDIFF)) {
            fragment_end_coordinate = fragment_start_coordinate + cigar_num_nts;
        } else if ((cigar_operation == BAM_CREF_SKIP) || (cigar_operation == BAM_CPAD) || (cigar_operation == BAM_CDEL)) {
            fragment_start_coordinate += cigar_num_nts;
            continue;
        } else if ((cigar_operation == BAM_CINS) || (cigar_operation == BAM_CSOFT_CLIP)) {
            continue;
        }

        counter_offset_start = fragment_start_coordinate / NTS_PER_COUNTER;
        counter_offset_end = fragment_end_coordinate / NTS_PER_COUNTER;

        if (counter_offset_start == counter_offset_end) {

#pragma omp critical
            {
                if (bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] == NULL) {
                    bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] = (bam_coverage_counter_t*) calloc(1, sizeof(bam_coverage_counter_t));
                    pthread_mutex_init(&(bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]->lock), NULL);
                }
            }

            bam_coverage_counter_incr_between_((fragment_start_coordinate % NTS_PER_COUNTER), (fragment_end_coordinate % NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]);
        } else {
#pragma omp critical
            {
                if (bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] == NULL) {
                    bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start] = (bam_coverage_counter_t*) calloc(1, sizeof(bam_coverage_counter_t));
                    pthread_mutex_init(&(bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]->lock), NULL);
                }
      
                if (bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end] == NULL) {
                    bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end] = (bam_coverage_counter_t*) calloc(1, sizeof(bam_coverage_counter_t));
                    pthread_mutex_init(&(bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end]->lock), NULL);
                }
            }

            bam_coverage_counter_incr_between_((fragment_start_coordinate % NTS_PER_COUNTER), (NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_start]);
            bam_coverage_counter_incr_between_(0, (fragment_end_coordinate % NTS_PER_COUNTER), bam_chromosome_coverage_p[chromosome].bam_coverage_counter_p[counter_offset_end]);
        }

        fragment_start_coordinate = fragment_end_coordinate;
    }
}

void bam_coverage_counter_incr_between_(int interval_start, int interval_end, bam_coverage_counter_t* counter_p) {
    pthread_mutex_lock(&counter_p->lock);

    for (int i = interval_start; i < interval_end; i++) {
        counter_p->coverage_counter[i]++;
    }

    pthread_mutex_unlock(&counter_p->lock);
}

char* num_chromosome_to_string_(int chromosome) {
    return str_chromosomes[chromosome];
}

int nt_coordinate_to_string_(char* str_coordinate, int coordinate) {
    int i, quotient;

    if ((coordinate >= 10000000) && (coordinate < 100000000)) {
        i = 8;
    } else if (coordinate >= 100000000) {
        i = 9;
    } else if ((coordinate >= 1000000) && (coordinate < 10000000)) {
        i = 7;
    } else if ((coordinate >= 100000) && (coordinate < 1000000)) {
        i = 6;
    } else if ((coordinate >= 10000) && (coordinate < 100000)) {
        i = 5;
    } else if ((coordinate >= 1000) && (coordinate < 10000)) {
        i = 4;
    } else if ((coordinate >= 100) && (coordinate < 1000)) {
        i = 3;
    } else if ((coordinate >= 10) && (coordinate < 100)) {
        i = 2;
    } else if (coordinate < 10) {
        i = 1;
    }

    quotient = coordinate;

    for (int j = (i - 1); j >= 0; j--) {
        str_coordinate[j] = 48 + (quotient % 10);
        quotient /= 10;
    }

    return i;
}

int nt_coverage_to_string_(char* str_coverage, int coverage) {
    int i, quotient;

    if ((coverage >= 10) && (coverage < 100)) {
        i = 2;
    } else if ((coverage >= 100) && (coverage < 1000)) {
        i = 3;
    } else if (coverage < 10) {
        i = 1;
    } else if ((coverage >= 1000) && (coverage < 10000)) {
        i = 4;
    } else if (coverage >= 10000) {
        i = 5;
    }

    quotient = coverage;

    for (int j = (i - 1); j >= 0; j--) {
        str_coverage[j] = 48 + (quotient % 10);
        quotient /= 10;
    }

    return i;
}
