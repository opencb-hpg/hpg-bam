
#ifndef CHROM_ALIGNMENTS_H
#define CHROM_ALIGNMENTS_H

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

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for single-chromosome alignments data
*
* Structure containing alignments data for a single chromosome
*/
typedef struct chrom_alignments {
    int initial_allocation;   			/**< Initial allocated alignments. */
    int allocated_alignment;   			/**< Actual allocated alignments (after resize if it is necessary). */

    int alignment_count;   			/**< Number of alignments. */
    pthread_mutex_t alignment_count_lock;   	/**< Lock for alignment_count variable. */

    int complete;   				/**< Flag of completed chromosome (in coverage calculation). */		
    pthread_mutex_t complete_lock;   		/**< Lock for complete variable. */

#ifdef THRUST-GPU
    thrust::host_vector<int> bam_alignment_coordinates_p;	/**< Vector of coordinates of the bam alignments. */
    thrust::host_vector<int> indices_p;				/**< Indices of the coordinate (for sorting). */
#else
    int* bam_alignment_coordinates_p;   	/**< Vector of coordinates of the bam alignments. */
    int* indices_p;   				/**< Indices of the coordinate (for sorting). */
#endif

    bam1_t** bam_alignments_p;   		/**< Pointers to bam1_t alignments. */
} chrom_alignments_t;

/**
* @brief Structure for storing chrom_alignments structures for each chromosome
*
* Structure for storing chrom_alignments structures for each chromosome
*/
typedef struct alignments_list {
    int num_chromosomes;			/**< Number of chromosomes. */
    chrom_alignments_t** chromosomes_p;		/**< Pointers to chrom_alignments structures. */
} alignments_list_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a new chrom_alignment
*  @param num_alignments number of alignments for memory allocation
*  @return pointer to chrom_alignment structure
*  
*  Creates a new chrom_alignment with the allocated num_alignments
*/
chrom_alignments_t* chrom_alignments_new(int num_alignments);

/**
*  @brief Reallocates space for more alignments in the chrom_alignment
*  @param chrom_alignments_p pointer to chrom_alignment structure
*  @param num_alignments new number of alignments for memory reallocation
*  @return pointer to chrom_alignment structure
*  
*  Reallocates space for more alignments in the given chrom_alignment
*/
chrom_alignments_t* chrom_alignments_realloc(chrom_alignments_t* chrom_alignments_p, int num_alignments);

/**
*  @brief Frees the given chrom_alignment
*  @param chrom_alignments_p pointer to chrom_alignment structure
*  @return void
*  
*  Frees the given chrom_alignment
*/
void chrom_alignments_free(chrom_alignments_t* chrom_alignments_p);

/**
*  @brief Gets an alignment from the chrom_aligment by index
*  @param chrom_alignments_p pointer to chrom_alignment structure
*  @param index index position to get from the chrom_alignment
*  @return pointer to bam1_t structure
*  
*  Reallocates space for more alignments in the given chrom_alignment
*/
bam1_t* chrom_alignments_get_alignment(chrom_alignments_t* chrom_alignments_p, int index);

/**
*  @brief Indicates if a chrom_alignment is completely filled
*  @param chrom_alignments_p pointer to chrom_alignment structure
*  @return 1 if complete, 0 if not complete
*  
*  Indicates if a chrom_alignment is completely filled for coverage calculations
*/
int chrom_alignments_is_complete(chrom_alignments_t* chrom_alignments_p);

/**
*  @brief Sets the completion flag in a chrom alignment
*  @param chrom_alignments_p pointer to chrom_alignment structure
*  @param complete completion flag (0 not completed, 1 completed)
*  @return void
*  
*  Sets the completion flag in a chrom alignment
*/
void chrom_alignments_set_complete(chrom_alignments_t* chrom_alignments_p, int complete);

/**
*  @brief Creates a new alignment list
*  @param num_chromosomes number of chromosomes of the list
*  @return pointer to the alignment list
*  
*  Creates a new alignment list
*/
alignments_list_t* alignments_list_new(int num_chromosomes);

/**
*  @brief Frees an alignment list
*  @param list_p pointer to the alignment list
*  @return void
*  
*  Frees an alignment list
*/
void alignments_list_free(alignments_list_t* list_p);

/**
*  @brief Creates a new chrom_alignment in the list
*  @param chromosome number of chromosome of the associated chrom_alignment
*  @param num_alignments number of alignments for memory allocation
*  @param list_p pointer to the alignment list
*  @return void
*  
*  Creates a new chrom_alignment in the list, if the chrom_alignment is yet created 
*  then does nothing
*/
void alignments_list_new_chrom_alignment(int chromosome, int num_alignments, alignments_list_t* list_p);

/**
*  @brief Gets a chrom_alignment from the list
*  @param chromosome number of chromosome to get its associated chrom_alignment
*  @param list_p pointer to the alignment list
*  @return chrom_alignment for the given chromosome
*  
*  Gets the chrom_alignment from the list for the given chromosome
*/
chrom_alignments_t* alignments_list_get_chrom_alignment(int chromosome, alignments_list_t* list_p);

/**
*  @brief Inserts the alignments of a given bam batch in a list of chrom_alignment
*  @param batch_p pointer to the bam batch
*  @param[in,out] list_p pointer to the alignment list
*  @return void
*  
*  Inserts the alignments of a given bam batch in a list of chrom_alignment
*/
void alignments_list_insert_batch(bam_batch_t* batch_p, alignments_list_t* list_p);

#endif /* BAM_ALIGNMENT_MATRIX_H */
