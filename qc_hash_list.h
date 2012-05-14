
#ifndef QC_HASH_LIST_H
#define QC_HASH_LIST_H

#include <pthread.h>

#include "commons.h"

#define ALLOCATED_MAPPINGS_PER_READ  8

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for storing reads positions
*
* Structure for storing reads positions
*/
typedef struct read_position {
    short int tid;			/**< Chromosome. */
    int position;			/**< Start position. */
} read_position_t;

/**
* @brief Item of the qc_hash_list
*
* Item of the qc_hash_list containing information for mapping histogram
*/
typedef struct qc_hash_list_item {
    read_position_t* pairend1_p;		/**< Chromosome. */
    read_position_t* pairend2_p;		/**< Chromosome. */
    short int allocated_pairends1;		/**< Chromosome. */
    short int allocated_pairends2;		/**< Chromosome. */
    short int num_pairends1;			/**< Chromosome. */
    short int num_pairends2;			/**< Chromosome. */
    //int pair_distance;
    char* id_seq;				/**< Chromosome. */
    pthread_mutex_t lock;			/**< Chromosome. */
    struct qc_hash_list_item* prev_p;		/**< Chromosome. */
    struct qc_hash_list_item* next_p;		/**< Chromosome. */
} qc_hash_list_item_t;

/**
* @brief List containing qc hash items
*
* Item of the qc_hash_list
*/
typedef struct qc_hash_list {
    int length;					/**< Chromosome. */
    pthread_mutex_t lock;			/**< Chromosome. */
    qc_hash_list_item_t* first_p;		/**< Chromosome. */
    qc_hash_list_item_t* last_p;		/**< Chromosome. */
} qc_hash_list_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

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
