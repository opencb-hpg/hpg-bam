/*
 * qc_hash.h
 *
 *  Created on: Aug 5, 2011
 *      Author: victor
 */

#ifndef QC_HASH_LIST_H
#define QC_HASH_LIST_H

#include <pthread.h>

#include "commons.h"

#define ALLOCATED_MAPPINGS_PER_READ		8

//====================================================================================
//  qc_hash_list.h
//
//  structures and methods for inserting and reading from qc hash table lists
//====================================================================================
typedef struct read_position {
  short int tid;
  int position;
} read_position_t;

typedef struct qc_hash_list_item {
  read_position_t* pairend1_p;
  read_position_t* pairend2_p;
  short int allocated_pairends1;
  short int allocated_pairends2;
  short int num_pairends1;
  short int num_pairends2;
  //int pair_distance;
  char* id_seq;
  pthread_mutex_t lock;
  struct qc_hash_list_item* prev_p;
  struct qc_hash_list_item* next_p;
} qc_hash_list_item_t;

typedef struct qc_hash_list {
  int length;
  pthread_mutex_t lock;
  qc_hash_list_item_t* first_p;
  qc_hash_list_item_t* last_p;  
} qc_hash_list_t;


// ------------------------------------------------
// qc hash list functions
// ------------------------------------------------

qc_hash_list_item_t* qc_hash_list_item_new(char* id_seq, int tid, int start_coordinate, int seq_length, int paired_end);
void qc_hash_list_item_free(qc_hash_list_item_t* item_p, int all);
void qc_hash_item_free(qc_hash_list_item_t* qc_hash_list_item_p, int all);

void qc_hash_list_init(qc_hash_list_t* list_p);
//void qc_hash_list_init(qc_hash_list_t* qc_hash_list_p, int producers);
void qc_hash_list_lock(qc_hash_list_t* list_p);
void qc_hash_list_unlock(qc_hash_list_t* list_p);
void qc_hash_list_insert(qc_hash_list_item_t* qc_hash_list_item_p, qc_hash_list_t* qc_hash_list_p);
qc_hash_list_item_t* qc_hash_list_remove(qc_hash_list_t* qc_hash_list_p);
int qc_hash_list_length(qc_hash_list_t* qc_hash_list_p);
void qc_hash_list_print(qc_hash_list_t* qc_hash_list_p);

void qc_hash_list_items_free(qc_hash_list_t* list_p);


#endif /* QC_HASH_LIST_H */