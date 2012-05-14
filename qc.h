
#ifndef QC_H
#define QC_H

#include <stdlib.h>
#include <stdio.h>

#include "cuda_commons.h"

#define VALID_ALIGNMENT_FILE_SUFFIX	".valid"
#define INVALID_ALIGNMENT_FILE_SUFFIX	".invalid"
#define NO_ALIGNMENT_FILE_SUFFIX

#define MAX_GFF_FILE_LINES		5000

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Performs quality control of a BAM file
*  @param batch_size max. size of the bam batch 
*  @param batch_list_size max. size of the bam batch list
*  @param gpu_num_threads number of gpu threads
*  @param gpu_num_blocks number of gpu blocks
*  @param cpu_num_threads number of cpu threads
*  @param base_quality base quality for quality values normalization
*  @param max_distance_size max distance between paired-end mappings to consider
*  @param input_filename bam file name
*  @param output_directory output directory where output report/files will be written
*  @param gff_filename gff file name
*  @return void
*  
*  Performs quality control of a BAM file (highest level function)
*/
void qc_bam_file(size_t batch_size, int batch_list_size, int gpu_num_threads, int gpu_num_blocks, int cpu_num_threads, int base_quality, int max_distance_size, char* input_filename, char* output_directory, char* gff_filename);

#endif  /* QC_H */
