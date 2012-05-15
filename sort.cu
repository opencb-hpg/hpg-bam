
#ifndef SORT_CU
#define SORT_CU

#include <stdlib.h>
#include <stdio.h>

#ifdef THRUST-GPU
    #include <thrust/host_vector.h>
    #include <thrust/device_vector.h>
    #include <thrust/sort.h>
    #include <thrust/copy.h>
#endif

extern "C" {
    #include "aligner_dataset.h"
    #include "aligner_dataset_file.h"
    #include "bam.h"
    #include "bam_reader.h"
    #include "bam_writer.h"
    #include "chrom_alignments.h"
    #include "commons.h"
    #include "file_utils.h"
    #include "log.h"
    #include "sam.h"
    #include "sort.h"
    #include "system_utils.h"
}

#include "sort_thrust.h"

#define BLOCK_SIZE 16

/* **********************************************
 *    		Global variables  		*
 * *********************************************/

int bam_reader_alive = 1;

/* ******************************************************
 *    		Function implementations  		*
 * *****************************************************/

void sort_bam_file(size_t batch_size, char* input_filename, char* output_directory) {
    int num_aligments_read;
    char input_shortname[MAX_FULL_PATH_LENGTH];
    char output_filename[MAX_FULL_PATH_LENGTH];
    char** split_filename = (char**) calloc(NUM_OF_CHROMOSOMES, sizeof(char*));

    get_filename_from_path(input_filename, input_shortname);
    sprintf(output_filename, "%s/%s%s", output_directory, input_shortname, SORTED_FILE_SUFFIX);

    alignments_list_t* list_p = alignments_list_new(NUM_OF_CHROMOSOMES);

    // first phase: one-reader vs multi-writers
    // calling threads to read alignments from file and write back to file (one per chromosome)
    bam_reader_t* bam_first_reader_p = bam_reader_new(input_filename, batch_size, 0, list_p, CHROMOSOME_MODE, NO_SORT, ALL_CHROMOSOMES); //base_quality = 0
    bam_writer_t** bam_split_writer_p = (bam_writer_t**) malloc(NUM_OF_CHROMOSOMES * sizeof(bam_writer_t*));

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        split_filename[i] = (char*) malloc(MAX_FULL_PATH_LENGTH * sizeof(char));
        sprintf(split_filename[i], "%s/%s.%i", output_directory, input_shortname, i);
        bam_split_writer_p[i] = bam_writer_new(split_filename[i], list_p, bam_first_reader_p->bam_file_p->bam_header_p, CHROMOSOME_MODE, i);
    }

    bam_reader_start(bam_first_reader_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_writer_start(bam_split_writer_p[i]);
    }

    if (time_flag) {
        start_timer(t1_write);
    }

    num_alignments = bam_reader_join(bam_first_reader_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_writer_join(bam_split_writer_p[i]);
    }

    // second phase: one-reader vs one-writer
    //    one-reader is in charged of reading all segmented-files
    //¿set counters to 0?
    alignments_list_free(list_p);
    list_p = alignments_list_new(NUM_OF_CHROMOSOMES);

    bam_reader_t* bam_split_reader_p;

    bam_file_t* bam_file_p = bam_fopen(input_filename);
    //bam_writer_t* bam_sorted_writer_p = bam_writer_new(output_filename, list_p, bam_first_reader_p->bam_file_p->bam_header_p, SEQUENTIAL_MODE, ALL_CHROMOSOMES);
    bam_writer_t* bam_sorted_writer_p = bam_writer_new(output_filename, list_p, bam_file_p->bam_header_p, SEQUENTIAL_MODE, ALL_CHROMOSOMES);
    bam_fclose(bam_file_p);

    bam_reader_alive = 1;

    bam_writer_start(bam_sorted_writer_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_split_reader_p = bam_reader_new(split_filename[i], batch_size, 0, list_p, SEQUENTIAL_MODE, SORT_BY_POSITION, i);
        bam_reader_start(bam_split_reader_p);
        num_aligments_read = bam_reader_join(bam_split_reader_p);
        //bam_reader_free(bam_split_reader_p);
    }

    bam_reader_alive = 0;

    bam_writer_join(bam_sorted_writer_p);

    // delete bam files per chromosome
    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        remove(split_filename[i]);
        free(split_filename[i]);
    }

    alignments_list_free(list_p);
}

void sort_bam_file_by_id(size_t batch_size, char* input_filename, char* output_directory) {
    int num_aligments_read;
    char input_shortname[MAX_FULL_PATH_LENGTH];
    char output_filename[MAX_FULL_PATH_LENGTH];
    char** split_filename = (char**) calloc(NUM_OF_CHROMOSOMES, sizeof(char*));

    get_filename_from_path(input_filename, input_shortname);
    sprintf(output_filename, "%s/%s%s", output_directory, input_shortname, SORTED_FILE_SUFFIX);

    alignments_list_t* list_p = alignments_list_new(NUM_OF_CHROMOSOMES);

    // first phase: one-reader vs multi-writers
    // calling threads to read alignments from file and write back to file (one per chromosome)
    bam_reader_t* bam_first_reader_p = bam_reader_new(input_filename, batch_size, 0, list_p, CHROMOSOME_MODE, NO_SORT, ALL_CHROMOSOMES); //base_quality = 0
    bam_writer_t** bam_split_writer_p = (bam_writer_t**) malloc(NUM_OF_CHROMOSOMES * sizeof(bam_writer_t*));

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        split_filename[i] = (char*) malloc(MAX_FULL_PATH_LENGTH * sizeof(char));
        sprintf(split_filename[i], "%s/%s.%i", output_directory, input_shortname, i);
        bam_split_writer_p[i] = bam_writer_new(split_filename[i], list_p, bam_first_reader_p->bam_file_p->bam_header_p, CHROMOSOME_MODE, i);
    }

    bam_reader_start(bam_first_reader_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_writer_start(bam_split_writer_p[i]);
    }

    if (time_flag) {
        start_timer(t1_write);
    }

    num_alignments = bam_reader_join(bam_first_reader_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_writer_join(bam_split_writer_p[i]);
    }

    // second phase: one-reader vs one-writer
    //    one-reader is in charged of reading all segmented-files
    //¿set counters to 0?
    alignments_list_free(list_p);
    list_p = alignments_list_new(NUM_OF_CHROMOSOMES);

    bam_reader_t* bam_split_reader_p;

    bam_file_t* bam_file_p = bam_fopen(input_filename);
    //bam_writer_t* bam_sorted_writer_p = bam_writer_new(output_filename, list_p, bam_first_reader_p->bam_file_p->bam_header_p, SEQUENTIAL_MODE, ALL_CHROMOSOMES);
    bam_writer_t* bam_sorted_writer_p = bam_writer_new(output_filename, list_p, bam_file_p->bam_header_p, SEQUENTIAL_MODE, ALL_CHROMOSOMES);
    bam_fclose(bam_file_p);

    bam_reader_alive = 1;

    bam_writer_start(bam_sorted_writer_p);

    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        bam_split_reader_p = bam_reader_new(split_filename[i], batch_size, 0, list_p, SEQUENTIAL_MODE, SORT_BY_POSITION, i);
        bam_reader_start(bam_split_reader_p);
        num_aligments_read = bam_reader_join(bam_split_reader_p);
        //bam_reader_free(bam_split_reader_p);
    }

    bam_reader_alive = 0;
    bam_writer_join(bam_sorted_writer_p);

    // delete bam files per chromosome
    for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
        remove(split_filename[i]);
        free(split_filename[i]);
    }

    alignments_list_free(list_p);
}

/*
// merge implementation
void sort_bam_file(char* bam_input, char* output_directory) {  // test for performance measurements

 int N = 150000000;
 int X = 1;
//  int ordena[N];
//  int ordena_join[X*N];
 int* ordena = (int*) malloc(N * sizeof(int));
 memset(ordena, 0, N * sizeof(int));
 int* ordena_before_merge = (int*) malloc(N * X * sizeof(int));
 memset(ordena_before_merge, 0, N * X * sizeof(int));
 int* ordena_merged = (int*) malloc(N * X * sizeof(int));
 memset(ordena_merged, 0, N * X * sizeof(int));

 unsigned int iseed = (unsigned int)time(NULL);
 srand(iseed);

 if (time_flag) { start_timer(t1_write); } //time_write stores random generator time
 for (int i=0; i<N; i++) {
     ordena[i] = rand();
     //if (N - i < 20) { printf("ordena[%i]: %i\n", i, ordena[i]);}
 }
 if (time_flag) { stop_timer(t1_write, t2_write, write_time); }

 for (int j=0; j<X; j++) {
     //sort_vector(ordena, N);
     //sort_key_value(ordena, ordena, N);
 }

 //sort_key_value(ordena_join, ordena_join, N);

 int ordena_indices[X];
 memset(&ordena_indices, 0, X * sizeof(int));

 for (int i=0; i<X; i++) {
     for (int j=0; j<N; j++) {
  ordena_before_merge[i*N + j] = ordena[j];
     }
 }

//  for (int i=0; i<N*X; i++) {
//      printf("ordena_before_merge[%i]: %i\n", i, ordena_before_merge[i]);
//  }

 if (time_flag) { start_timer(t1_result); } //time_result stores merge time
 int next_pos;
 for (int k=0; k<N*X; k++) {

     for (int m=0; m<X; m++) {
  if (ordena_before_merge[m * N + ordena_indices[m]] < ordena_before_merge[next_pos * N + ordena_indices[next_pos]]) {
      next_pos = m;
  }
     }

     ordena_merged[k] = ordena_before_merge[next_pos * N + ordena_indices[next_pos]];
     ordena_indices[next_pos]++;
     if (ordena_indices[next_pos] == N) { ordena_indices[next_pos]--; }
     next_pos = 0;

 }
 if (time_flag) { stop_timer(t1_result, t2_result, result_time); }

//  for (int i=0; i<N*X; i++) {
//      printf("ordena_merged[%i]: %i\n", i, ordena_merged[i]);
//  }

}*/

void sort_dataset_by_id(char* dataset_input, char* output_directory) {
    unsigned int read_lines = 0;
    aligner_dataset_list_t* list_p = aligner_dataset_list_new(20000000);

    aligner_dataset_file_t* dataset_file_p = aligner_dataset_fopen(dataset_input);

    //fill dataset list in the reader
    read_lines = aligner_dataset_read_list(dataset_file_p, list_p, 0);

    //sort dataset
    aligner_dataset_list_sort_by_id(list_p);

    //write sorted file to disk and free resources
    char input_shortname[MAX_FULL_PATH_LENGTH];
    char output_filename[MAX_FULL_PATH_LENGTH];

    get_filename_from_path(dataset_input, input_shortname);
    sprintf(output_filename, "%s/%s%s", output_directory, input_shortname, ALIGNER_DATASET_SORTED_FILE_SUFFIX);

    char** split_filename = (char**) calloc(NUM_OF_CHROMOSOMES, sizeof(char*));
    aligner_dataset_list_write(list_p, output_filename);

    aligner_dataset_list_free(list_p);
    aligner_dataset_fclose(dataset_file_p);
}

#endif /* SORT_CU */
