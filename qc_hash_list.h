
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
    read_position_t* pairend1_p;	/**< Mapping positions for the pairend 1. */
    read_position_t* pairend2_p;	/**< Mapping positions for the pairend 2. */
    short int allocated_pairends1;	/**< Number of pairend 1 allocated read_position. */
    short int allocated_pairends2;	/**< Number of pairend 2 allocated read_position. */
    short int num_pairends1;		/**< Number of pairend 1 mappings. */
    short int num_pairends2;		/**< Number of pairend 2 mappings. */
    char* id_seq;			/**< Id. of the sequence. */
    pthread_mutex_t lock;		/**< Lock of the item. */
    struct qc_hash_list_item* prev_p;	/**< Pointer to the previous element in the list. */
    struct qc_hash_list_item* next_p;	/**< Pointer to the next element in the list. */
} qc_hash_list_item_t;

/**
* @brief List containing qc hash items
*
* Item of the qc_hash_list
*/
typedef struct qc_hash_list {
    int length;				/**< Length of the list. */
    pthread_mutex_t lock;		/**< Lock of the list. */
    qc_hash_list_item_t* first_p;	/**< Pointer to the first element of the list. */
    qc_hash_list_item_t* last_p;	/**< Pointer to the last element of the list. */
} qc_hash_list_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a new qc hash list item
*  @param id_seq id of the sequence in the item
*  @param tid chromosome
*  @param start_coordinate start coordinate
*  @param seq_length sequence length
*  @param paired_end flag indicated paired end
*  @return pointer to the created qc hash list item
*  
*  Creates and returns a new qc hash list item
*/
qc_hash_list_item_t* qc_hash_list_item_new(char* id_seq, int tid, int start_coordinate, int seq_length, int paired_end);

/**
*  @brief Frees a given qc hash list item
*  @param[in,out] item_p pointer to the item
*  @param all flag to free the inner mallocs of the item
*  @return void
*  
*  Frees a given qc hash list item
*/
void qc_hash_list_item_free(qc_hash_list_item_t* item_p, int all);

/**
*  @brief Inits a qc hash list
*  @param[in,out] list_p pointer to the qc_hash_list
*  @return void
*  
*  Inits a qc hash list
*/
void qc_hash_list_init(qc_hash_list_t* list_p);

/**
*  @brief Locks a qc hash list
*  @param[in,out] list_p pointer to the qc_hash_list
*  @return void
*  
*  Locks a qc hash list
*/
void qc_hash_list_lock(qc_hash_list_t* list_p);

/**
*  @brief Unlocks a qc hash list
*  @param[in,out] list_p pointer to the qc_hash_list
*  @return void
*  
*  Unlocks a qc hash list
*/
void qc_hash_list_unlock(qc_hash_list_t* list_p);

/**
*  @brief Inserts an item into the qc hash list
*  @param qc_hash_list_item_p pointer to the qc_hash_list_item to insert
*  @param qc_hash_list_p pointer to the qc_hash_list
*  @return void
*  
*  Inserts an item into the qc hash list
*/
void qc_hash_list_insert(qc_hash_list_item_t* qc_hash_list_item_p, qc_hash_list_t* qc_hash_list_p);


qc_hash_list_item_t* qc_hash_list_remove(qc_hash_list_t* qc_hash_list_p);

/**
*  @brief Gets the length of a qc_hash_list
*  @param qc_hash_list_p pointer to the qc_hash_list
*  @return void
*  
*  Gets the length of a qc_hash_list
*/
int qc_hash_list_length(qc_hash_list_t* qc_hash_list_p);

/**
*  @brief Prints a qc hash list
*  @param qc_hash_list_p pointer to the qc_hash_list
*  @return void
*  
*  Prints a qc hash list and all its items
*/
void qc_hash_list_print(qc_hash_list_t* qc_hash_list_p);

/**
*  @brief Frees ALL qc hash list items
*  @param[in,out] list_p pointer to the qc hash list
*  @return void
*  
*  Frees ALL qc hash list items
*/
void qc_hash_list_items_free(qc_hash_list_t* list_p);

#endif /* QC_HASH_LIST_H */
