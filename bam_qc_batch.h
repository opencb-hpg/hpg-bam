
#ifndef BAM_QC_BATCH_H
#define BAM_QC_BATCH_H

#include <pthread.h>

#include "bam_file.h"
#include "commons.h"
#include "qc_hash.h"

#define X  		0
#define I  		1
#define D  		2
#define M  		3
#define EQUAL  		4
#define MISMATCHES 	5

#define COUNTERS_SIZE 	6

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for qc alignment data
*
* Structure containing qc data of an alignment
*/
typedef struct qc_alignment {
    unsigned int counters[COUNTERS_SIZE];  /**< Counters with qc data of an alignment. */
} qc_alignment_t;

/**
* @brief Batch of bam qc data
*
* Structure containing quality control data
*/
typedef struct bam_qc_batch {
    int id;    				/**< Id of the batch. */
    int num_alignments;   		/**< Number of alignments. */
    int num_blocks;   			/**< Number of GPU blocks (for some calculations). */
    bam1_t** alignments_p;  		/**< Pointers to the alignments. */
    qc_alignment_t* qc_alignment_p; 	/**< Pointer to qc alignment structures. */
    int* strand_counter_p;  		/**< Pointer to strand counter vector. */
    int* map_quality_p;   		/**< Pointer to map quality vector. */
    int* alignment_length_p;  		/**< Pointer to alignment length vector. */
    struct bam_qc_batch* prev_p;  	/**< Pointer to the previous bam_qc_batch. */
    struct bam_qc_batch* next_p;  	/**< Pointer to the next bam_qc_batch. */
} bam_qc_batch_t;

/**
* @brief List of bam_qc_batch elements
*
* List containing and linking bam_qc_batch elements
*/
typedef struct bam_qc_batch_list {
    int length;    			/**< Length of the list. */
    pthread_mutex_t lock;   		/**< Lock. */
    struct bam_qc_batch* first_p;  	/**< Pointer to the first item in the list. */
    struct bam_qc_batch* last_p;   	/**< Pointer to the last item in the list. */
} bam_qc_batch_list_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a bam qc batch
*  @param id id of the batch to create
*  @return void
*
*  Creates a new bam qc batch with the given id
*/
bam_qc_batch_t* bam_qc_batch_new(int id);

/**
*  @brief Frees a bam qc batch
*  @param[in,out] batch_p pointer to the bam qc batch to free
*  @param all flag to indicate if alignments are freed
*  @return void
*
*  Free a bam qc batch (all = 0 -> alignments not freed, all = 1 -> alignments freed)
*/
void bam_qc_batch_free(bam_qc_batch_t* batch_p, int all);

#endif /* BAM_QC_BATCH_H */
