
#ifndef QC_HASH_H
#define QC_HASH_H

#include <pthread.h>

#include "bam_qc_report.h"
#include "commons.h"
#include "qc_hash_list.h"

#define QC_HASH_LENGTH  	10000000

#define PAIRED_END1  		1
#define PAIRED_END2  		2
#define NO_PAIRED_END  		0

#define LOWER_8_BITS_MASK 255
#define HIGHER_8_BITS_MASK 65280

#define GET_TID1( tid )  tid&LOWER_8_BITS_MASK
#define GET_TID2( tid )  (tid&HIGHER_8_BITS_MASK)>>8

#define SET_TID1( tid )  tid&LOWER_8_BITS_MASK
#define SET_TID2( tid )  (tid&HIGHER_8_BITS_MASK)>>8

#define GET_COUNT1( count ) count&LOWER_8_BITS_MASK
#define GET_COUNT2( count ) (count&HIGHER_8_BITS_MASK)>>8

#define COUNT_INCR1( count )  count++
#define COUNT_INCR2( count )  count += 256

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for storing qc hash sublists
*
* Structure for storing qc hash sublists
*/
typedef struct qc_hash {
    int length;   			/**< Length of the structure. */
    int count;   			/**< Number of elements inserted. */
    qc_hash_list_t* qc_hash_list_p;   	/**< Pointer to the qc_hash_list. */
    pthread_mutex_t lock;   		/**< Lock for the structure. */
} qc_hash_t;

/**
* @brief Counter for mapping histogram and paired end distance
*
* Counter for mapping histogram and paired end distance
*/
typedef struct qc_mapping_counter {
    unsigned int num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 2];   	/**< Vector for mapping histogram. */
    unsigned long mean_paired_end_distance;   					/**< Mean distance between paired ends. */
    pthread_mutex_t lock;   							/**< Lock for the structure. */
} qc_mapping_counter_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a new qc hash
*  @param length length of the structure
*  @return pointer to the created qc hash
*  
*  Creates and returns a new qc hash
*/
qc_hash_t* qc_hash_new(int length);

/**
*  @brief Frees a given qc hash
*  @param[in,out] qc_hash_p pointer to the qc_hash
*  @param all flag to free the inner list of the qc hash
*  @return void
*  
*  Frees a given qc hash and optionally its inner list
*/
void qc_hash_free(qc_hash_t* qc_hash_p, int all);

/**
*  @brief Inits a qc hash 
*  @param[in,out] qc_hash_p pointer to the qc_hash
*  @return void
*  
*  Inits a qc hash
*/
void qc_hash_init(qc_hash_t* qc_hash_p);

/**
*  @brief Locks a qc hash 
*  @param[in,out] qc_hash_p pointer to the qc_hash
*  @return void
*  
*  Locks a qc hash
*/
void qc_hash_lock(qc_hash_t* qc_hash_p);

/**
*  @brief Unlocks a qc hash 
*  @param[in,out] qc_hash_p pointer to the qc_hash
*  @return void
*  
*  Unlocks a qc hash
*/
void qc_hash_unlock(qc_hash_t* qc_hash_p);

/**
*  @brief Inserts an aligment into the qc hash
*  @param qc_hash_p pointer to the qc_hash
*  @param id_seq id of the sequence
*  @param tid chromosome
*  @param start_coordinate start coordinate of the sequence
*  @param seq_length length of the inserted sequence
*  @param pairend_end flag of paired-end
*  @return void
*  
*  Inserts an aligment into the qc hash
*/
void qc_hash_insert_alignment(qc_hash_t* qc_hash_p, char* id_seq, int tid, int start_coordinate, int seq_length, short int paired_end);

/**
*  @brief Prints a qc hash 
*  @param qc_hash_p pointer to the qc_hash
*  @return void
*  
*  Prints a qc hash
*/
void qc_hash_print(qc_hash_t* qc_hash_p);

/**
*  @brief Perform calculations over the qc hash
*  @param qc_hash_p pointer to the qc_hash
*  @param[in,out] counter_p pointer to the qc_mapping_counter
*  @param[in,out] mean_paired_end_distance mean distance between paired ends
*  @param max_distance_size 
*  @param cpu_num_threads
*  @return void
*  
*  Perform two calculations over the qc hash:
*      - Mean distance between paired ends
*      - Histogram of mappings per reads  
*/
void qc_hash_perform_calculations(qc_hash_t* qc_hash_p, qc_mapping_counter_t* counter_p, unsigned long* mean_paired_end_distance, int max_distance_size, int cpu_num_threads);

/**
*  @brief Inits a qc mapping counter
*  @param[in,out] qc_mapping_counter_p pointer to the qc_mapping_counter
*  @return void
*  
*  Inits a qc mapping counter
*/
void qc_mapping_counter_init(qc_mapping_counter_t* qc_mapping_counter_p);

#endif /* QC_HASH_H */
