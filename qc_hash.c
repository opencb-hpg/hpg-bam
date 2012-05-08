/*
 * qc_hash.c
 *
 *  Created on: Aug 5, 2011
 *      Author: victor
 */

#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "bam.h"
#include "bam_qc_report.h"
#include "qc_hash.h"
#include "GeneralHashFunctions.h"
#include "string_utils.h"

//====================================================================================
//  qc_hash.c
//
//  structures and methods for inserting, calculating and reading from qc hash table
//====================================================================================

// ------------------------------------------------
//  p r i v a t e  f u n c t i o n s
// ------------------------------------------------

qc_hash_list_item_t* qc_hash_find_id_seq_(qc_hash_t* qc_hash_p, char* id_seq, unsigned int* alignment_hash);
void qc_hash_insert_id_seq_(qc_hash_t* qc_hash_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end, unsigned int alignment_hash);
void qc_hash_update_id_seq_(qc_hash_list_item_t* item_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end);
void qc_hash_incr_count_(qc_hash_list_item_t* item_p, char* id_seq, short int paired_end);
int qc_hash_determine_paired_end_(char* id_seq);
//int qc_hash_list_item_perform_calculations_(qc_hash_list_item_t* item_p, unsigned int* num_mappings_histogram, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads);
int qc_hash_list_item_perform_calculations_(qc_hash_list_item_t* item_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads);

// ------------------------------------------------
// qc_hash_new functions
// ------------------------------------------------

qc_hash_t* qc_hash_new(int qc_hash_length) {
  
  qc_hash_t* qc_hash_p = (qc_hash_t*) calloc(1, sizeof(qc_hash_t));
  qc_hash_p->length = qc_hash_length;
  qc_hash_p->count = 0;
  qc_hash_p->qc_hash_list_p = (qc_hash_list_t*) calloc(qc_hash_length, sizeof(qc_hash_list_t));
  
  return qc_hash_p; 
}

//-----------------------------------------------------
// qc_hash_free
//
// free the hash table
//-----------------------------------------------------

void qc_hash_free(qc_hash_t* qc_hash_p, int all) {
  
  if (all) {
    
//     for (int i=0; i < qc_hash_p->length; i++) {
//       qc_hash_list_items_free(&qc_hash_p->qc_hash_list_p[i]);
//     }
        
    if (qc_hash_p->qc_hash_list_p != NULL) free(qc_hash_p->qc_hash_list_p);
  }
  
  free(qc_hash_p);  
}

//-----------------------------------------------------
// qc_hash_init
//
// init the hash mutex
//-----------------------------------------------------

void qc_hash_init(qc_hash_t* qc_hash_p) {
  //qc_hash_p->lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&(qc_hash_p->lock), NULL);
}

//-----------------------------------------------------
// qc_hash_lock
//
// lock the hash mutex
//-----------------------------------------------------

void qc_hash_lock(qc_hash_t* qc_hash_p) {
  pthread_mutex_lock(&(qc_hash_p->lock));
}

//-----------------------------------------------------
// qc_hash_unlock
//
// unlock the hash mutex
//-----------------------------------------------------

void qc_hash_unlock(qc_hash_t* qc_hash_p) {
  pthread_mutex_unlock(&(qc_hash_p->lock));
}

//-----------------------------------------------------
// qc_hash_update_count
//
// update the count for paired end 1 or 2
//-----------------------------------------------------

void qc_hash_insert_alignment(qc_hash_t* qc_hash_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end) {

  unsigned int alignment_hash;
  qc_hash_list_item_t* found_item_p = qc_hash_find_id_seq_(qc_hash_p, id_seq, &alignment_hash);


  if (found_item_p != NULL) {
    qc_hash_update_id_seq_(found_item_p, id_seq, tid, start_coordinate, seq_length, paired_end);
    //printf("UPDATE paired end: %i\n", paired_end);
  } else {
    qc_hash_insert_id_seq_(qc_hash_p, id_seq, tid, start_coordinate, seq_length, paired_end, alignment_hash);
    //printf("INSERT paired end: %i\n", paired_end);
  }  
}

//-----------------------------------------------------
// qc_hash_perform_calculations
//
//    calculate over the qc hash table to obtain:
//    - Mean distance between paired ends
//    - Histogram of mappings per reads
//-----------------------------------------------------

/*
void qc_hash_perform_calculations(qc_hash_t* qc_hash_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads) {
  
  int num_mappings = 0;
  
  qc_hash_list_t* list_p;
  qc_hash_list_item_t* item_p;
  
  //#pragma omp parallel for num_threads(cpu_num_threads) shared(qc_hash_p, counter_p, mean_paired_end_distance, max_distance_size)
  for (int i=0; i < qc_hash_p->length; i++) {
    list_p = &(qc_hash_p->qc_hash_list_p[i]);
    
    for (int j=0; j < list_p->length; j++) {
      //printf("list_length: %i\n", list_p->length);	//
      item_p = qc_hash_list_remove(list_p);
      num_mappings += qc_hash_list_item_perform_calculations_(item_p, counter_p, mean_paired_end_distance, max_distance_size,cpu_num_threads);
    }    
  }
  //#pragma omp barrier

  //num_mappings = 1;
  if (num_mappings > 0) *mean_paired_end_distance /= (unsigned long) num_mappings;
}*/


void qc_hash_perform_calculations(qc_hash_t* qc_hash_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads) {
  
  int num_mappings = 0;
  
  qc_hash_list_t* list_p;
  qc_hash_list_item_t* item_p;

  for (int i=0; i < qc_hash_p->length; i++) {
     
    for (int j=0, length = qc_hash_p->qc_hash_list_p[i].length; j != length; j++) {
      
      item_p = qc_hash_list_remove(&(qc_hash_p->qc_hash_list_p[i]));
      
      num_mappings += qc_hash_list_item_perform_calculations_(item_p, counter_p, mean_paired_end_distance, max_distance_size,cpu_num_threads);

      qc_hash_list_item_free(item_p, 1);
    }
  }

  //printf("mean_paired_end_distance: %li, num_mappings: %i\n", *mean_paired_end_distance, num_mappings);
  *mean_paired_end_distance /= (unsigned long) num_mappings;  
}


//-----------------------------------------------------
// qc_hash_print
//
// print the hash table
//-----------------------------------------------------

void qc_hash_print(qc_hash_t* qc_hash_p) {  
  for (int i=0; i < qc_hash_p->length; i++) {
    printf("hash position: %i\n", i);
    qc_hash_list_print(&(qc_hash_p->qc_hash_list_p[i]));
  }  
}

//-----------------------------------------------------
// qc_mapping_counter_init
//
// init the qc_mapping_counter
//-----------------------------------------------------

void qc_mapping_counter_init(qc_mapping_counter_t* qc_mapping_counter_p) {
  
  //qc_mapping_counter_p->num_mappings_histogram = {0};  
  //for (int i = 0; i < (MAX_MAPPING_COUNT_IN_HISTOGRAM + 2); i++) {
  //    qc_mapping_counter_p->num_mappings_histogram[i] = 0;
  //}
  
  memset(&qc_mapping_counter_p->num_mappings_histogram, 0, (MAX_MAPPING_COUNT_IN_HISTOGRAM + 2) * sizeof(int));
  //qc_mapping_counter_p->lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&(qc_mapping_counter_p->lock), NULL);
  
}

// ------------------------------------------------
//  p r i v a t e  f u n c t i o n s  i m p l
// ------------------------------------------------

// ------------------------------------------------
// qc_hash_find_id_seq_
// ------------------------------------------------

qc_hash_list_item_t* qc_hash_find_id_seq_(qc_hash_t* qc_hash_p, char* id_seq, unsigned int* alignment_hash) {
  
  int id_seq_length = strlen(id_seq);
  
  //unsigned int find_index = (RSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
  unsigned int find_index = (RSHash(id_seq, id_seq_length) % QC_HASH_LENGTH);
  *alignment_hash = find_index;
  
  //printf("id_seq: %s, find_index: %i\n", id_seq, find_index);
  
// --------------------------- T E S T -----------------  
//   unsigned int find_index = APHash(id_seq, strlen(id_seq) - 2);
//   unsigned int find_index = (RSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (JSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (PJWHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (ELFHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (BKDRHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (SDBMHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (DJBHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (DEKHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (BPHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//   unsigned int find_index = (FNVHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
// --------------------------- T E S T -----------------    
  
//   printf("find index: %u\n");
//   printf("--------------------->");
 
  qc_hash_list_t* list_p = &(qc_hash_p->qc_hash_list_p[find_index%QC_HASH_LENGTH]);
  qc_hash_list_lock(list_p);
  qc_hash_list_item_t* item_p = list_p->first_p;

  //for (int i=1; i < list_p->length; i++) {
  for (int i=0; i < list_p->length; i++) {
    //printf("item_p->id_seq: %s, id_seq: %s\n", item_p->id_seq, id_seq);
    
    //if (strcmp(item_p->id_seq, id_seq)) {
    if (strcmp(item_p->id_seq, id_seq) == 0) {
    
      //printf("SAME ID SEQ...\n");
      qc_hash_list_unlock(list_p);
      return item_p;
    }

    item_p = item_p->next_p;
  }
  
  qc_hash_list_unlock(list_p);
  return NULL;
}

// ------------------------------------------------
// qc_hash_insert_id_seq_
// ------------------------------------------------

void qc_hash_insert_id_seq_(qc_hash_t* qc_hash_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end, unsigned int alignment_hash) {
  
  int id_seq_length = strlen(id_seq);
  
  //unsigned int index = (RSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH); 
  //unsigned int index = (RSHash(id_seq, id_seq_length) % QC_HASH_LENGTH);
  unsigned int index = alignment_hash; 
  
// --------------------------- T E S T -----------------  
// if (time_flag) { start_timer(t1_write); }
//  for (int i = 0; i < 10; i++) {
//    index = (APHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);  
//    index = (RSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (JSHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (PJWHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (ELFHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (BKDRHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (SDBMHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (DJBHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (DEKHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (BPHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//    index = (FNVHash(id_seq, id_seq_length - 2) % QC_HASH_LENGTH);
//  }
// if (time_flag) { stop_timer(t1_write, t2_write, write_time); }
// --------------------------- T E S T -----------------    
  
  //int paired_end = qc_hash_determine_paired_end_(id_seq);
  //printf("paired_end______%i\n", paired_end);

  qc_hash_list_item_t* item_p = qc_hash_list_item_new(id_seq, tid, start_coordinate, seq_length, paired_end);
  qc_hash_list_insert(item_p, &(qc_hash_p->qc_hash_list_p[index]));  
}

// ------------------------------------------------
// qc_hash_update_id_seq_
// ------------------------------------------------

void qc_hash_update_id_seq_(qc_hash_list_item_t* item_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end) {
  
  //int paired_end = qc_hash_determine_paired_end_(id_seq);
  
  if (seq_length > 0) {
    if (paired_end == PAIRED_END1) {
      item_p->pairend1_p[item_p->num_pairends1].tid = tid;
      item_p->pairend1_p[item_p->num_pairends1].position = start_coordinate + seq_length;
      item_p->num_pairends1++;
    } else if (paired_end == PAIRED_END2) {
      item_p->pairend2_p[item_p->num_pairends2].tid = tid;
      item_p->pairend2_p[item_p->num_pairends2].position = start_coordinate;
      item_p->num_pairends2++;
    } else {
      //printf("id_seq do not match /1 or /2 termination, NO PAIRED END FOUND!!!!!\n");
    }
  }
  
  //reallocating memory for read positions in paired end 1
  //
  if (item_p->num_pairends1 == item_p->allocated_pairends1) {
    read_position_t* read_position_aux1_p = (read_position_t*) calloc(item_p->allocated_pairends1, sizeof(read_position_t));
    memcpy(read_position_aux1_p, item_p->pairend1_p, item_p->allocated_pairends1 * sizeof(read_position_t));
    free(item_p->pairend1_p);
    
    item_p->allocated_pairends1 += ALLOCATED_MAPPINGS_PER_READ;
    item_p->pairend1_p = (read_position_t*) calloc(item_p->allocated_pairends1, sizeof(read_position_t));
    memcpy(item_p->pairend1_p, read_position_aux1_p, (item_p->allocated_pairends1 - ALLOCATED_MAPPINGS_PER_READ) * sizeof(read_position_t));
    
    free(read_position_aux1_p);
  }
  
  //reallocating memory for read positions in paired end 2
  //
  if (item_p->num_pairends2 == item_p->allocated_pairends2) {
    read_position_t* read_position_aux2_p = (read_position_t*) calloc(item_p->allocated_pairends2, sizeof(read_position_t));
    memcpy(read_position_aux2_p, item_p->pairend2_p, item_p->allocated_pairends2 * sizeof(read_position_t));
    free(item_p->pairend2_p);
    
    item_p->allocated_pairends2 += ALLOCATED_MAPPINGS_PER_READ;
    item_p->pairend2_p = (read_position_t*) calloc(item_p->allocated_pairends2, sizeof(read_position_t));
    memcpy(item_p->pairend2_p, read_position_aux2_p, (item_p->allocated_pairends2 - ALLOCATED_MAPPINGS_PER_READ) * sizeof(read_position_t));
    
    free(read_position_aux2_p);
  }
  
}

// ------------------------------------------------
// qc_hash_incr_count_
// ------------------------------------------------

void qc_hash_incr_count_(qc_hash_list_item_t* item_p, char* id_seq, short int paired_end) {
  
  //int paired_end = qc_hash_determine_paired_end_(id_seq);
  
   if (paired_end == PAIRED_END1) {
    item_p->num_pairends1++;
  } else if (paired_end == PAIRED_END2) {
    item_p->num_pairends2++;
  } else {
    //printf("qc_hash_incr_count_: id_seq do not match /1 or /2 termination, NO PAIRED END FOUND!!!!!\n");
  }
}

// ------------------------------------------------
// qc_hash_determine_paired_end_
// ------------------------------------------------

int qc_hash_determine_paired_end_(char* id_seq) {
  
  int id_seq_length = strlen(id_seq);
  
  if (id_seq[id_seq_length-2] == '/') {
    if (id_seq[id_seq_length-1] == '1') {
      return PAIRED_END1;
    } else if (id_seq[id_seq_length-1] == '2') {
      return PAIRED_END2;
    } else {
      LOG_FATAL("id_seq do not match /1 or /2 termination, NO PAIRED END FOUND!!!!!\n");
      return NO_PAIRED_END;
    }    
  } else {
    return NO_PAIRED_END;
  }  
}

// ------------------------------------------------
// qc_hash_list_item_perform_calculations_
// ------------------------------------------------
/*
int qc_hash_list_item_perform_calculations_(qc_hash_list_item_t* item_p, unsigned int* num_mappings_histogram, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads) {
    
  int distance_size, num_mappings = 0;
  
  int num_pairends1 = item_p->num_pairends1;
  int num_pairends2 = item_p->num_pairends2;
  
  read_position_t* pairend1_p = item_p->pairend1_p;
  read_position_t* pairend2_p = item_p->pairend2_p;
    
  //printf("entering double for... num_pairends1: %i, num_pairends2: %i\n", num_pairends1, num_pairends2);
  
  #pragma omp parallel for num_threads(cpu_num_threads) shared(num_pairends1, num_pairends2, pairend1_p, pairend2_p, num_mappings_histogram, num_mappings, mean_paired_end_distance) private(distance_size)
  for (int i=0; i < num_pairends1; i++) {    
    for (int j=0; j < num_pairends2; j++) {
      
      if (pairend1_p[i].tid == pairend2_p[j].tid) {
	//printf("same chromosome...\n");
	distance_size = abs(pairend2_p[j].position - pairend1_p[i].position);  

	if (distance_size <= max_distance_size) {
	    #pragma omp critical 
	    {
	      num_mappings++;
	      *mean_paired_end_distance += distance_size;
	    }
	    //printf("num_mappings###: %i\n", num_mappings);
	    //printf("distance_size: %i\n", abs(distance_size));
	}
      } else {
	//printf("distinct chromosome...\n");
      }
    }
  }
  
  if (num_mappings <= MAX_MAPPING_COUNT_IN_HISTOGRAM) {
    num_mappings_histogram[num_mappings]++;
  } else {
    num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 1]++;    
  }
  
  return num_mappings;
}*/


int qc_hash_list_item_perform_calculations_(qc_hash_list_item_t* item_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads) {  
    
  int distance_size, num_mappings = 0;
  
  int num_pairends1 = item_p->num_pairends1;
  int num_pairends2 = item_p->num_pairends2;
  
  read_position_t* pairend1_p = item_p->pairend1_p;
  read_position_t* pairend2_p = item_p->pairend2_p;
  
  //printf("entering double for... num_pairends1: %i, num_pairends2: %i\n", num_pairends1, num_pairends2);
  
  #pragma omp parallel for num_threads(cpu_num_threads) shared(num_pairends1, num_pairends2, pairend1_p, pairend2_p, counter_p, num_mappings, mean_paired_end_distance) private(distance_size)
  for (int i=0; i < num_pairends1; i++) {    
    for (int j=0; j < num_pairends2; j++) {
      
      if (pairend1_p[i].tid == pairend2_p[j].tid) {

	distance_size = abs(pairend2_p[j].position - pairend1_p[i].position);  

	if (distance_size <= max_distance_size) {
	    #pragma omp critical 
	    {
	      num_mappings++;	      
	      *mean_paired_end_distance += distance_size;
	    }
	    //printf("num_mappings: %i\n", num_mappings);
	    //printf("distance_size: %i\n", abs(distance_size));
	}
      } else {
	//printf("distinct chromosome...\n");
      }
    }
  }
  
  if (num_mappings <= MAX_MAPPING_COUNT_IN_HISTOGRAM) {
    pthread_mutex_lock(&(counter_p->lock));
    counter_p->num_mappings_histogram[num_mappings]++;
    pthread_mutex_unlock(&(counter_p->lock));
  } else {
    pthread_mutex_lock(&(counter_p->lock));
    counter_p->num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 1]++;    
    pthread_mutex_unlock(&(counter_p->lock));
  }
  
  return num_mappings;
}



























