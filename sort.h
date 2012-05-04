
#ifndef SORT_H_
#define SORT_H_

#include "bam.h"

#include "aligner_dataset.h"
#include "aligner_dataset_file.h"
#include "commons.h"
#include "file_utils.h"
#include "log.h"
#include "system_utils.h"
#include "cuda_commons.h"

#define SORTED_FILE_SUFFIX	".sorted"

/*
      SORT BAM FILE
*/

void sort_bam_file(size_t batch_size, char* input_filename, char* output_directory);
void sort_bam_file_by_id(size_t batch_size, char* input_filename, char* output_directory);

/*
      SORT SAM FILE
*/

void sort_sam_file(char* input_filename, char* output_directory);

/*
      SORT DATASET FILE
*/

void sort_dataset_by_id(char* dataset_input, char* output_directory);

#endif
