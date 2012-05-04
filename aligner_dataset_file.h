
#ifndef ALIGNER_DATASET_FILE_H
#define ALIGNER_DATASET_FILE_H

#include <stdio.h>

#include "aligner_dataset.h"
#include "bam_commons.h"
#include "commons.h"
#include "limits.h"

#define MAX_DATASET_LINE_LENGTH		512

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Aligner dataset file handler
* 
* Handler structure for managing aligner dataset files
*/
typedef struct aligner_dataset_file {
  unsigned int num_lines;				/**< Number of lines. */
  char* filename;					/**< File name. */
  char* mode;						/**< Opening mode ("r", "w"). */
  FILE* fd;						/**< File descriptor. */
} aligner_dataset_file_t;

/**
* @brief Aligner dataset batch
* 
* Batch of aligner dataset lines
*/
typedef struct aligner_dataset_batch {
    int num_lines;					/**< Number of lines. */
    aligner_dataset_line_t** aligner_dataset_lines_p;	/**< Pointers to the lines. */
    struct aligner_dataset_batch* prev_p;		/**< Pointer to the next batch in the list. */
    struct aligner_dataset_batch* next_p;		/**< Pointer to the previous batch in the list. */
} aligner_dataset_batch_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
 *  @brief Open an aligner dataset file in read mode
 *  @param filename filename containing the dataset to read
 *  @return aligner_dataset_file_t
 *  
 *  This function opens a dataset file and creates the structure for
 *  handling the file in read mode
 */
aligner_dataset_file_t* aligner_dataset_fopen(char* filename);

/**
 *  @brief Open an aligner dataset file in the specified mode
 *  @param filename filename containing the dataset to be opened
 *  @param mode open mode ("r", "w", "a", ...) 
 *  @return aligner_dataset_file_t
 *  
 *  This function opens a dataset file and creates the structure for
 *  handling the file in the specified mode
 */
aligner_dataset_file_t* aligner_dataset_fopen(char* filename, char* mode);

/**
 *  @brief Close an aligner dataset file and free the file handler
 *  @param aligner_dataset_file file handler
 *  @return void
 *  
 *  This function closes an aligner dataset file and free its associated
 *  handler structure
 */
void aligner_dataset_fclose(aligner_dataset_file_t* aligner_dataset_file);

/**
 *  @brief Fills an aligner dataset batch
 *  @param dataset_file_p handler of the datase file
 *  @param batch_p pointer to the batch to be filled
 *  @param num_lines number of dataset lines to got (0 if all)
 *  @return void
 *  
 *  This function sorts a dataset file and puts the resulting sorted dataset in the given 
 *  ouput directory
 */
unsigned int aligner_dataset_read_batch(aligner_dataset_file_t* dataset_file_p, aligner_dataset_batch_t* batch_p, int num_lines);

/**
 *  @brief Fills an aligner dataset list 
 *  @param dataset_file_p pointer to the dataset file handler
 *  @param list_p pointer to the aligner dataset list to be filled
 *  @param num_lines number of dataset lines to got (0 if all)
 *  @return void
 *  
 *  Fills an aligner dataset list with the content of the specified dataset file
 */
unsigned int aligner_dataset_read_list(aligner_dataset_file_t* dataset_file_p, aligner_dataset_list_t* list_p, int num_lines);

/**
 *  @brief Creates an aligner dataset batch
 *  @param num_lines number of lines of the batch
 *  @return aligner_dataset_batch_t batch
 *  
 *  Creates and returns an aligner dataset batch
 */
aligner_dataset_batch_t* aligner_dataset_batch_new(int num_lines);

/**
 *  @brief Frees an aligner dataset batch
 *  @param batch_p pointer to the batch to be freed
 *  
 *  Frees memory allocation for aligner dataset batch
 */
void aligner_dataset_batch_free(aligner_dataset_batch_t* batch_p);

/**
 *  @brief Prints the content of the aligner dataset batch
 *  @param batch_p pointer to the batch to be printed
 *  
 *  Prints in the standard output the content of the aligner dataset batch line by line
 */
void aligner_dataset_batch_print(aligner_dataset_batch_t* batch_p);

#endif    /* ALIGNER_DATASET_FILE_H */