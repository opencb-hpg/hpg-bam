
#ifndef QC_H_
#define QC_H_

#include <stdlib.h>
#include <stdio.h>

#include "bam_coverage.h"
#include "bam_data_batch.h"
#include "bam_data_batch_list.h"
#include "bam_qc_batch.h"
#include "bam_qc_report.h"
#include "bam_reader.h"
#include "commons.h"
#include "cuda_commons.h"
#include "file_utils.h"
#include "gff_data.h"
#include "gff_reader.h"
#include "list.h"
#include "log.h"
#include "qc_hash.h"
#include "qc_kernel_omp.h"
#include "sam.h"
#include "system_utils.h"


#define VALID_ALIGNMENT_FILE_SUFFIX	".valid"
#define INVALID_ALIGNMENT_FILE_SUFFIX	".invalid"
#define NO_ALIGNMENT_FILE_SUFFIX

#define MAX_GFF_FILE_LINES		5000


/*
      SORT BAM FILE
*/

void qc_bam_file(size_t batch_size, int batch_list_size, int gpu_num_threads, int gpu_num_blocks, int cpu_num_threads, int base_quality, int max_distance_size, char* input_filename, char* output_directory, char* gff_filename);


#endif