
#ifndef ALIGNER_DATASET_H
#define ALIGNER_DATASET_H

#include <stdint.h>
#include <stdio.h>

#include "commons.h"
#include "log.h"

#define INIT_ALIGNER_DATASET_LIST_LINES		20000000
#define INCR_ALIGNER_DATASET_LIST_LINES		10000000
#define ALIGNER_DATASET_SORTED_FILE_SUFFIX	".sorted"

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Aligner dataset line
* 
* Line of an aligner dataset file
*/
typedef struct aligner_dataset_line {
    short int chromosome;		/**< Chromosome. */
    short int strand;			/**< Strand (0 or 1). */
    short int num_errors;		/**< Number of errors. */
    short int num_indels;		/**< Number of insertions&deletions. */
    int start;				/**< Start position of the alignment. */
    int end;				/**< End position of the alignment. */
    char* seq_id;  			/**< Sequence id. */
} aligner_dataset_line_t;

/**
* @brief List of aligner dataset lines
* 
* List containing aligner dataset lines from a dataset file
*/
typedef struct aligner_dataset_list {
    size_t num_lines;					/**< Number of lines in the list. */
    size_t allocated_lines;				/**< Number of allocated lines in the list. */
    aligner_dataset_line_t** aligner_dataset_lines_p;	/**< Pointers to aligner dataset lines. */
    char** seq_id_p;					/**< Pointers to sequence ids. */
    int* indices_p;					/**< Indices to obtain the seq ids . */
} aligner_dataset_list_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Creates a new aligner dataset line
*  @param seq_id_length Length of the sequence id
*  @return pointer to the created aligner dataset line
*  
*  Creates and returns an aligner dataset line
*/
aligner_dataset_line_t* aligner_dataset_line_new(short int seq_id_length);

/**
*  @brief Frees an aligner dataset line
*  @param aligner_dataset_line_p Pointer to the aligner dataset line to be freed
*  @param all Flag to indicate that inner mallocs must or must not be freed
*  @return void
*  
*  Frees an aligner dataset line
*/
void aligner_dataset_line_free(aligner_dataset_line_t* aligner_dataset_line_p, int all);

/**
*  @brief Creates a list of aligner dataset lines
*  @param num_lines Number of lines allocated in the list
*  @return pointer to the aligner dataset list
*  
*  Creates a list of aligner dataset lines
*/
aligner_dataset_list_t* aligner_dataset_list_new(size_t num_lines);

/**
*  @brief Frees a list of aligner dataset lines
*  @param list_p Pointer to the list to free
*  @return void
*  
*  Creates a list of aligner dataset lines
*/
void aligner_dataset_list_free(aligner_dataset_list_t* list_p);

/**
*  @brief Inserts a given aligner dataset line in the list
*  @param list_p pointer to the list in which line will be inserted
*  @param[in,out] line_p pointer to the line to insert
*  @return void
*  
*  Inserts a given aligner dataset line into a list
*/
void aligner_dataset_list_insert_line(aligner_dataset_list_t* list_p, aligner_dataset_line_t* line_p);

/**
 *  @brief Resize a dataset line list
 *  @param list_p pointer to the list to be resized
 *  @param length additional length for the list to be resized
 *  @return void
 *  
 *  This function resize a dataset line list with the additional indicated length
 */
aligner_dataset_list_t* aligner_dataset_list_realloc(aligner_dataset_list_t* list_p, size_t num_lines);

/**
 *  @brief Sorts a dataset list
 *  @param list_p pointer to the list to be sorted
 *  @return void
 *  
 *  This function sorts dataset lines by id in the list
 */
void aligner_dataset_list_sort_by_id(aligner_dataset_list_t* list_p);

/**
*  @brief Writes aligner dataset lines in a list to a file
*  @param list_p pointer to the list to write
*  @param filename name of the file in which lines are written
*  @return void
*  
*  Writes to disk aligner dataset lines in the list
*/
void aligner_dataset_list_write(aligner_dataset_list_t* list_p, char* filename);

#endif  /* ALIGNER_DATASET_H */