
#ifndef SORT_H
#define SORT_H

#include "cuda_commons.h"

#define SORTED_FILE_SUFFIX ".sorted"

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Sorts a BAM file by chromosome/position
*  @param batch_size max. size of the bam batch
*  @param input_filename bam file name
*  @param output_directory output directory where sorted file will be written
*  @return void
*  
*  Sorts a BAM file by chromosome/position
*/
void sort_bam_file(size_t batch_size, char* input_filename, char* output_directory);

/**
*  @brief Sorts a BAM file by sequence id
*  @param batch_size max. size of the bam batch
*  @param input_filename bam file name
*  @param output_directory output directory where sorted file will be written
*  @return void
*  
*  Sorts a BAM file by sequence id
*/
void sort_bam_file_by_id(size_t batch_size, char* input_filename, char* output_directory);

/**
*  @brief Sorts a read dataset file by sequence id
*  @param input_filename dataset file name
*  @param output_directory output directory where  sorted file will be written
*  @return void
*  
*  Sorts a read dataset file by sequence id
*  It is assumed that the file fits in memory 
*/
void sort_dataset_by_id(char* dataset_input, char* output_directory);

#endif  /* SORT_H */
