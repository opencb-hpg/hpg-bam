
#ifndef QC_KERNEL_OMP_H
#define QC_KERNEL_OMP_H

#include "bam_data_batch.h"
#include "bam_qc_batch.h"
#include "commons.h"
#include "cuda_commons.h"
#include "log.h"

#define NUMBER_OF_NTS_IN_GENOME  		3200000000
#define WINDOW_SIZE_FOR_COVERAGE_UNIFORMITY 	1000

/* **************************************
 *    		Functions  		*
 * *************************************/

void cpu_bam_qc_basic_stats(bam_data_core_t* core_data_p, int* strand_counter_p, int* map_quality_p, int* alignment_length_p, int num_alignments, int cpu_num_threads);
void cpu_bam_qc_map_errors(bam_data_core_t* core_data_p, uint32_t* cigar_data_p, qc_alignment_t* qc_alignment_p, int num_alignments);

#endif  /* QC_KERNEL_OMP_H */
