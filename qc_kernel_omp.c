
#include <omp.h>

#include "bam.h"
#include "qc_kernel_omp.h"

/* ******************************************************************************
 *    		Kernel functions implementations (CPU - OpenMP)  		*
 * *****************************************************************************/

void cpu_bam_qc_basic_stats(bam_data_core_t* core_data_p, int* strand_counter_p, int* map_quality_p, int* alignment_length_p, int num_alignments, int cpu_num_threads) {

#pragma omp parallel for num_threads(cpu_num_threads) shared(core_data_p, strand_counter_p, map_quality_p, alignment_length_p)
    for (int k = 0; k < num_alignments; k++) {
#pragma omp critical
        {
            strand_counter_p[0] += core_data_p[k].strand;
            map_quality_p[0] += core_data_p[k].map_quality;
            alignment_length_p[0] += core_data_p[k].alignment_length;
        }
    }

}

void cpu_bam_qc_map_errors(bam_data_core_t* core_data_p, uint32_t* cigar_data_p, qc_alignment_t* qc_alignment_p, int num_alignments) {
    int cigar_start_pos, cigar_end_pos, cigar_operation, cigar_num;
    uint32_t cigar_position;

    for (int k = 0; k < num_alignments; k++) {
        cigar_start_pos = core_data_p[k].cigar_index;
        cigar_end_pos = core_data_p[k+1].cigar_index;

        for (int i = cigar_start_pos; i < cigar_end_pos; i++) {
            cigar_position = cigar_data_p[i];
            cigar_operation = (cigar_position & BAM_CIGAR_MASK);
            cigar_num = cigar_position >> BAM_CIGAR_SHIFT;

            switch (cigar_operation) {
                case BAM_CMATCH:
                    qc_alignment_p[k].counters[M] += cigar_num;
                    break;
                case BAM_CINS:
                    qc_alignment_p[k].counters[I] += cigar_num;
                    break;
                case BAM_CDEL:
                    qc_alignment_p[k].counters[D] += cigar_num;
                    break;
                case BAM_CEQUAL:
                    qc_alignment_p[k].counters[EQUAL] += cigar_num;
                    break;
                case BAM_CDIFF:
                    qc_alignment_p[k].counters[X] += cigar_num;
                    break;
            }
        }

        qc_alignment_p[k].counters[MISMATCHES] = qc_alignment_p[k].counters[D] + qc_alignment_p[k].counters[I] + qc_alignment_p[k].counters[X];
    }
}





