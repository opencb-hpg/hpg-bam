/*
 * qc_hash.h
 *
 *  Created on: Aug 5, 2011
 *      Author: victor
 */

#ifndef QC_HASH_H
#define QC_HASH_H

#include <pthread.h>

#include "bam_qc_report.h"
#include "commons.h"
#include "qc_hash_list.h"

//====================================================================================
//  macros for extracting tid and count for each paired end
//====================================================================================

#define QC_HASH_LENGTH		10000000

#define PAIRED_END1		1
#define PAIRED_END2		2
#define NO_PAIRED_END		0

#define LOWER_8_BITS_MASK	255
#define HIGHER_8_BITS_MASK	65280

#define GET_TID1( tid )		tid&LOWER_8_BITS_MASK
#define GET_TID2( tid )		(tid&HIGHER_8_BITS_MASK)>>8

#define SET_TID1( tid )		tid&LOWER_8_BITS_MASK
#define SET_TID2( tid )		(tid&HIGHER_8_BITS_MASK)>>8

#define GET_COUNT1( count )	count&LOWER_8_BITS_MASK
#define GET_COUNT2( count )	(count&HIGHER_8_BITS_MASK)>>8

#define COUNT_INCR1( count ) 	count++
#define COUNT_INCR2( count )  count += 256

#define MAX_MAPPING_COUNT_IN_HISTOGRAM	10


//====================================================================================
//  qc_hash.h
//
//  structures and methods for inserting, calculating and reading from qc hash table
//====================================================================================

typedef struct qc_hash {
  int length;
  int count;
  qc_hash_list_t* qc_hash_list_p;
  pthread_mutex_t lock;
} qc_hash_t;

typedef struct qc_mapping_counter {
  unsigned int num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 2];
  //unsigned int* num_mappings_histogram;
  unsigned long mean_paired_end_distance;
  pthread_mutex_t lock;
} qc_mapping_counter_t;

// ------------------------------------------------
// bam_qc_hash functions
// ------------------------------------------------

qc_hash_t* qc_hash_new(int length);
void qc_hash_free(qc_hash_t* qc_hash_p, int all);
void qc_hash_init(qc_hash_t* qc_hash_p);
void qc_hash_lock(qc_hash_t* qc_hash_p);
void qc_hash_unlock(qc_hash_t* qc_hash_p);

void qc_hash_insert_alignment(qc_hash_t* qc_hash_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end);
void qc_hash_print(qc_hash_t* qc_hash_p);

//void qc_hash_perform_calculations(qc_hash_t* qc_hash_p, unsigned int* num_mappings_histogram, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads);
void qc_hash_perform_calculations(qc_hash_t* qc_hash_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads);

void qc_mapping_counter_init(qc_mapping_counter_t* qc_mapping_counter_p);

#endif /* QC_HASH_H */