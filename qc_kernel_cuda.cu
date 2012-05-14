
#include <omp.h>

#include "bam.h"
#include "qc_kernel_cuda.h"

extern "C" {
    #include "log.h"
}

//------------------------------------------------------------------------------------
//  kernels call functions
//------------------------------------------------------------------------------------

void call_kernel_basic_stats(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, int* d_strand_counter_p, int* d_map_quality_p,  int* d_alignment_length_p, int num_alignments) {
    kernel_bam_qc_basic_stats <<< dimGrid, dimBlock>>>(d_core_data_p, d_strand_counter_p, d_map_quality_p, d_alignment_length_p, num_alignments);
}

void call_kernel_map_errors(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, uint32_t* d_cigar_data_p, qc_alignment_t* d_qc_alignment_p, int num_alignments) {
    kernel_bam_qc_map_errors <<< dimGrid, dimBlock>>>(d_core_data_p, d_cigar_data_p, d_qc_alignment_p, num_alignments);
}

//------------------------------------------------------------------------------------
//  kernel for bam qc basic statistics ( G P U  implementation)
//------------------------------------------------------------------------------------

__global__ void kernel_bam_qc_basic_stats(bam_data_core_t* d_core_data_p, int* d_strand_counter_p, int* d_map_quality_p, int* d_alignment_length_p, int num_alignments) {

    __shared__ int s_strand_counter_p[OPTIMAL_THREADS_FOR_COMPUTE_CAPABILITY_20];
    __shared__ int s_map_quality_p[OPTIMAL_THREADS_FOR_COMPUTE_CAPABILITY_20];
    __shared__ int s_alignment_length_p[OPTIMAL_THREADS_FOR_COMPUTE_CAPABILITY_20];

    int k;
    unsigned int tid = threadIdx.x;

    k = blockIdx.x * blockDim.x + threadIdx.x;

//   if (threadIdx.x == 0 && blockIdx.x == 0) {
//     int n, s = 0, q = 0;
//     for(n=0 ; n<blockDim.x ; n++) {
//       s += d_bam_data_batch_p[n + k].strand;
//       q += d_bam_data_batch_p[n + k].map_quality;
//     }
//     printf("total: strand (+) = %i, quality = %i\n", s, q);
//   }
//   __syncthreads();

    if (k < num_alignments) {
        s_strand_counter_p[tid] = d_core_data_p[k].strand;
        s_map_quality_p[tid] = d_core_data_p[k].map_quality;
        s_alignment_length_p[tid] = d_core_data_p[k].alignment_length;
        //if (threadIdx.x == 0 && blockIdx.x == 0) { printf("init: bid = %i, tid = %i: counter[0] = %i\n", blockIdx.x, threadIdx.x, s_strand_counter_p[0]);}
    } else {
        s_strand_counter_p[tid] = 0;
        s_map_quality_p[tid] = 0;
        s_alignment_length_p[tid] = 0;
    }
    __syncthreads();

    for (unsigned int s = 1; s < blockDim.x; s *= 2) {
        if (tid % (2*s) == 0) {
            s_strand_counter_p[tid] += s_strand_counter_p[tid + s];
            s_map_quality_p[tid] += s_map_quality_p[tid + s];
            s_alignment_length_p[tid] += s_alignment_length_p[tid + s];
        }
        //if (threadIdx.x == 0 && blockIdx.x == 0) { printf("step %i: bid = %i, tid = %i: counter[0] = %i\n", s, blockIdx.x, threadIdx.x, s_strand_counter_p[0]);}
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) {
        d_strand_counter_p[blockIdx.x] = s_strand_counter_p[0];
        d_map_quality_p[blockIdx.x] = s_map_quality_p[0];
        d_alignment_length_p[blockIdx.x] = s_alignment_length_p[0];
        //if (threadIdx.x == 0 && blockIdx.x == 0) { printf("done: bid = %i, tid = %i: counter[0] = %i\n", blockIdx.x, threadIdx.x, s_strand_counter_p[0]);}
    }
}

//------------------------------------------------------------------------------------
//  kernel for bam qc alignment mismatch count ( G P U implementation)
//------------------------------------------------------------------------------------

__global__ void kernel_bam_qc_map_errors(bam_data_core_t* d_core_data_p, uint32_t* d_cigar_data_p, qc_alignment_t* d_qc_alignment_p, int num_alignments) {

    int k;
    int cigar_start_pos, cigar_end_pos, cigar_operation, cigar_num;
    uint32_t cigar_position;

    k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k < num_alignments) {
        cigar_start_pos = d_core_data_p[k].cigar_index;
        //cigar_end_pos = cigar_start_pos + 1;
        cigar_end_pos = d_core_data_p[k+1].cigar_index;

        for (int i = cigar_start_pos; i < cigar_end_pos; i++) {
            cigar_position = d_cigar_data_p[i];
            cigar_operation = (cigar_position & BAM_CIGAR_MASK);
            cigar_num = cigar_position >> BAM_CIGAR_SHIFT;
            //printf("cigar operation: %i, cigar num nts: %i\n", cigar_operation, cigar_num);

            switch (cigar_operation) {
                case BAM_CMATCH:
                    d_qc_alignment_p[k].counters[M] += cigar_num;
                    break;  //M
                case BAM_CINS:
                    d_qc_alignment_p[k].counters[I] += cigar_num;
                    break;  //I
                case BAM_CDEL:
                    d_qc_alignment_p[k].counters[D] += cigar_num;
                    break;  //D
                case BAM_CEQUAL:
                    d_qc_alignment_p[k].counters[EQUAL] += cigar_num;
                    break; //=
                case BAM_CDIFF:
                    d_qc_alignment_p[k].counters[X] += cigar_num;
                    break;  //X
            }
        }

        d_qc_alignment_p[k].counters[MISMATCHES] = d_qc_alignment_p[k].counters[D] + d_qc_alignment_p[k].counters[I] + d_qc_alignment_p[k].counters[X];
    }
}
