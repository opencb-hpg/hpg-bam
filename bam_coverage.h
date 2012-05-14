
#ifndef BAM_COVERAGE_H
#define BAM_COVERAGE_H

#include <pthread.h>

#include "bam_commons.h"
#include "gff_data.h"

#include "bam_data_batch_list.h"
#include "commons.h"

#define COVERAGE_FILE_SUFFIX  ".coverage"

#define MAX_NTS_PER_CHROMOSOME 250000000
#define NTS_PER_COUNTER   10000000


/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Coverage counter
*
* Coverage counter for nucleotide coverage store
*/
typedef struct bam_coverage_counter {
    short int chromosome;     /**< Chromosome. */
    short int print;     /**< Flag to print or not the counter (printed when completely filled). */
    unsigned short int coverage_counter[NTS_PER_COUNTER]; /**< Vector for storing coverage count by nucleotide. */
    pthread_mutex_t lock;     /**< Lock. */
} bam_coverage_counter_t;

/**
* @brief Container for bam_chromosome_counter structures
*
* Container that points to partial bam_chromosome_counter to complete
* the required positions to compute a complete chromosome count
*/
typedef struct bam_chromosome_coverage {
    bam_coverage_counter_t** bam_coverage_counter_p; /**< Pointers to partial bam_coverage_counter. */
} bam_chromosome_coverage_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
 *  @brief Inits a bam_chromosome_coverage structure
 *  @param bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to initialize
 *  @return void
 *
 *  Initializes bam_chromosome_coverage structure
 */
void bam_chromosome_coverage_init(bam_chromosome_coverage_t* bam_chromosome_coverage_p);

/**
 *  @brief Clears a bam_chromosome_coverage structure
 *  @param bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to clear
 *  @return void
 *
 *  Clears/frees a given bam_chromosome_coverage structure
 */
void bam_chromosome_coverage_clear(bam_chromosome_coverage_t* bam_chromosome_coverage_p);

/**
 *  @brief Marks those counters within a bam_chromosome_coverage structure to print
 *  @param bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to mark
 *  @param all flag to mark all counters pointed in the bam_chromosome_coverage structure
 *  @return void
 *
 *  Marks those counters within a bam_chromosome_coverage structure to print. The criteria to mark
 *  the counter is the completeness of the counter
 */
void bam_coverage_counter_mark_to_print(bam_chromosome_coverage_t* bam_chromosome_coverage_p, int all);

/**
 *  @brief Prints to disk the marked-to-print counters (line by line)
 *  @param bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to print
 *  @param output_directory directory where coverage file will be written
 *  @param input_filename file name to write the computed coverage
 *  @return void
 *
 *  Marks those counters within a bam_chromosome_coverage structure to print. The criteria to mark
 *  the counter is the completeness of the counter. The file is printed line by line.
 */
void bam_coverage_counter_print(bam_chromosome_coverage_t* bam_chromosome_coverage_p, char* output_directory, char* input_filename);

/**
 *  @brief Prints to disk the marked-to-print counters (block by block)
 *  @param bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to print
 *  @param output_directory directory where coverage file will be written
 *  @param input_filename file name to write the computed coverage
 *  @return void
 *
 *  Marks those counters within a bam_chromosome_coverage structure to print. The criteria to mark
 *  the counter is the completeness of the counter. The file is printed block by block. A block is
 *  a buffer of coverage lines.
 */
void bam_coverage_counter_print_block(bam_chromosome_coverage_t* bam_chromosome_coverage_p, char* output_directory, char* input_filename);

/**
 *  @brief Deletes a coverage file
 *  @param output_directory directory where coverage file will be deleted
 *  @param input_filename file name to delete
 *  @return void
 *
 *  Delete a coverage file. This operation is done previously to write a new coverage file
 */
void bam_coverage_counter_delete_file(char* output_directory, char* input_filename);

/**
 *  @brief Computes coverage of alignments in a batch
 *  @param batch_p pointer to the bam batch for computing the coverage
 *  @param[in,out] bam_chromosome_coverage_p pointer to bam_chromosome_coverage structure to store coverage values
 *  @param gff_data_p pointer to gff data with regions of interest in which coverage has to be calculated
 *  @param output_directory directory where coverage file will be written
 *  @param input_filename file name in which coverage is written
 *  @param cpu_num_threads number of cpu threads to compute coverage (OpenMP)
 *  @return void
 *
 *  Delete a coverage file. This operation is done previously to write a new coverage file
 */
void bam_coverage_compute(bam_data_batch_t* batch_p, bam_chromosome_coverage_t* bam_chromosome_coverage_p, gff_data_t* gff_data_p, char* output_directory, char* input_filename, int cpu_num_threads);

/**
 *  @brief Inits str_coverage_matrix and strlen_coverage_matrix global structures
 *  @return void
 *
 *  Initializes str_coverage_matrix and strlen_coverage_matrix global structures
 */
void str_coverage_matrix_init();

#endif /* BAM_COVERAGE_H */
