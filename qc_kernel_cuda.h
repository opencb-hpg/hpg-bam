
#ifndef QC_KERNEL_CUDA_H
#define QC_KERNEL_CUDA_H

#include "bam_qc_batch.h"
#include "bam_data_batch.h"
#include "commons.h"
#include "cuda_commons.h"
//#include "log.h"

#define NUMBER_OF_NTS_IN_GENOME 	3200000000
#define WINDOW_SIZE_FOR_COVERAGE_UNIFORMITY	1000

/*
	kernel_bam_qc implementation
*/

__global__ void kernel_bam_qc_basic_stats(bam_data_core_t* d_core_data_p, qc_info_t* d_qc_info_p, int* d_strand_counter_p, int* d_map_quality_p, int* d_alignment_length_p, int num_alignments);
__global__ void kernel_bam_qc_map_errors(bam_data_core_t* d_core_data_p, uint32_t* d_cigar_data_p, qc_alignment_t* d_qc_alignment_p, int num_alignments);

#endif