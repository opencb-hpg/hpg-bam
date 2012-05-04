
#ifndef CHROM_ALIGNMENTS_H_
#define CHROM_ALIGNMENTS_H_

#include <stdio.h>

#include <pthread.h>

#ifdef THRUST-GPU
  #include <thrust/host_vector.h>
  #include <thrust/device_vector.h>
  #include <thrust/sort.h>
  #include <thrust/copy.h>
#endif

#include "bam.h"
#include "bam_file.h"
#include "commons.h"


//=====================================================
// structures
//=====================================================

// typedef struct chrom_alignments {
// 
//     int initial_allocation;
//     int allocated_alignment;
//     
//     int alignment_count;
//     pthread_mutex_t alignment_count_lock;
//     
//     int complete;
//     pthread_mutex_t complete_lock;
//     
//     int* bam_alignment_coordinates_p;
//     int* indices_p;
//     bam1_t** bam_alignments_p;
// } chrom_alignments_t;

typedef struct chrom_alignments {

    int initial_allocation;
    int allocated_alignment;
    
    int alignment_count;
    pthread_mutex_t alignment_count_lock;
    
    int complete;
    pthread_mutex_t complete_lock;

#ifdef THRUST-GPU
    thrust::host_vector<int> bam_alignment_coordinates_p;
    thrust::host_vector<int> indices_p;
#else
    int* bam_alignment_coordinates_p;
    int* indices_p;
#endif    

    bam1_t** bam_alignments_p;
} chrom_alignments_t;


typedef struct alignments_list {
    int num_chromosomes;
    chrom_alignments_t** chromosomes_p;
} alignments_list_t;


//=====================================================
// chromosome alignment functions
//=====================================================

chrom_alignments_t* chrom_alignments_new(int num_alignments);
chrom_alignments_t* chrom_alignments_realloc(chrom_alignments_t* chrom_alignments_p, int num_alignments);
void chrom_alignments_free(chrom_alignments_t* chrom_alignments_p);

bam1_t* chrom_alignments_get_alignment(chrom_alignments_t* chrom_alignments_p, int index);
int chrom_alignments_is_complete(chrom_alignments_t* chrom_alignments_p);
void chrom_alignments_set_complete(chrom_alignments_t* chrom_alignments_p, int complete);


//=====================================================
// alignment list functions
//=====================================================

alignments_list_t* alignments_list_new(int num_chromosomes);
void alignments_list_free(alignments_list_t* list_p);

void alignments_list_new_chrom_alignment(int chromosome, int num_alignments, alignments_list_t* list_p);
chrom_alignments_t* alignments_list_get_chrom_alignment(int chromosome, alignments_list_t* list_p);

void alignments_list_insert_batch(bam_batch_t* batch_p, alignments_list_t* list_p);


//bam_batch_t* get_bam_batch_form_alignments_list(size_t batch_size, alignments_list_t* list_p);


/*
chrom_alignments_t* bam_alignments_init(chrom_alignments_t* bam_chrom_alignments_p);
chrom_alignments_t chrom_alignments_init(chrom_alignments_t chrom_alignments, int num_alignments);
void bam_alignments_free(chrom_alignments_t* bam_chrom_alignments_p);
void chrom_alignments_free(chrom_alignments_t chrom_alignments, int num_alignments);
int* bam_read_alignment_count_init(int* bam_read_alignment_count_p);
*/
#endif /* BAM_ALIGNMENT_MATRIX_H_ */
