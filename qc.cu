
#ifndef QC_CU
#define QC_CU

extern "C" {
    #include "bam_coverage.h"
    #include "bam_data_batch.h"
    #include "bam_data_batch_list.h"
    #include "bam_qc_batch.h"
    #include "bam_qc_report.h"
    #include "bam_reader.h"
    #include "commons.h"
    #include "file_utils.h"
    #include "gff_data.h"
    #include "gff_reader.h"
    #include "log.h"
    #include "qc.h"
    #include "qc_kernel_omp.h"
    #include "sam.h"
    #include "system_utils.h"
}

#include "qc_kernel_cuda.h"

/* ******************************************************
 *    		Private thread functions 		*
 * *****************************************************/

void* qc_calc_server(void* params_p);
void* cpus_server(void* params_p);
void* results_server(void* params_p);

/* **********************************************
 *    		Global variables 		*
 * *********************************************/

list_t bam_qc_batch_list;

int bam_batch_reader_alive = 1;
int gpus_thread_alive = 1;
int cpus_thread_alive = 1;

pthread_mutex_t gpus_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t cpus_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;

/* **********************************************
 *    		Extern calls	 		*
 * *********************************************/

extern void call_kernel_basic_stats(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, int* d_strand_counter_p, int* d_alignment_length_p, int* d_map_quality_p, int num_alignments);
extern void call_kernel_map_errors(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, uint32_t* d_cigar_data_p, qc_alignment_t* d_qc_alignment_p, int num_alignments);

/* **********************************************************************
 *    		Private thread functions implementations 		*
 * *********************************************************************/

/* ******************************************************
 *    		QC calc server thread	 		*
 * *****************************************************/

void* qc_calc_server(void* params_p) {
    LOG_DEBUG("Thread-GPU: START\n");

    if (time_flag) {
        start_timer(t1_qc_calc_server);
    }

    qc_calc_server_input_t* input_p = (qc_calc_server_input_t*) params_p;

    int cpu_num_threads = input_p->cpu_num_threads;
    list_item_t* bam_data_batch_list_item_p = NULL;
    list_t* gpu_batch_list_p = input_p->gpu_batch_list_p;
    list_t* cpu_batch_list_p = input_p->cpu_batch_list_p;
    bam_data_batch_t* bam_data_batch_p = NULL;
    bam_qc_batch_t* bam_qc_batch_p = NULL;

    // variables for store output results in both CPU and GPU
    qc_alignment_t* qc_alignment_p;
    int* strand_counter_p;
    int* map_quality_p;
    int* alignment_length_p;

    bam_data_core_t* d_core_data_p;
    uint32_t* d_cigar_data_p;
    qc_alignment_t* d_qc_alignment_p;
    int* d_strand_counter_p;
    int* d_map_quality_p;
    int* d_alignment_length_p;

    // selecting GPU device
    CUDA_SAFE_CALL(cudaSetDevice(input_p->gpu_device_id[0]));

    int reads_alive;

    reads_alive = list_get_writers(gpu_batch_list_p);


    while ((bam_data_batch_list_item_p = list_remove_item(gpu_batch_list_p)) != NULL) {
        LOG_DEBUG("Thread-GPU: waiting for batch....\n");

        //if (time_flag) { 
        //    start_timer(t1_gpu); 
        //}

        if ((time_flag) && (gpus_standby_time == 0.0)) {
            stop_timer(t1_active_reader, t1_active_gpus, gpus_standby_time);
        }

        char log_message[50];
        sprintf(log_message, "Thread-GPU: processing for batch %i....\n", bam_data_batch_list_item_p->id);
        LOG_DEBUG(log_message);

        number_of_batchs++;

        // allocation memory for output results
        bam_data_batch_p = (bam_data_batch_t*) bam_data_batch_list_item_p->data_p;
        int num_alignments = bam_data_batch_p->num_alignments;
        int num_blocks;

//cpu_num_threads = 1;
        if (cpu_num_threads == 0) {  // GPU implementation
            num_blocks = (num_alignments / input_p->gpu_num_threads) + 1;

            dim3 dimBlock(input_p->gpu_num_threads, 1, 1);
            dim3 dimGrid(num_blocks, 1, 1);

            strand_counter_p = (int*) calloc(num_blocks, sizeof(int));
            map_quality_p = (int*) calloc(num_blocks, sizeof(int));
            alignment_length_p = (int*) calloc(num_blocks, sizeof(int));
            qc_alignment_p = (qc_alignment_t*) calloc(bam_data_batch_p->num_alignments, sizeof(qc_alignment_t));

            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_core_data_p, (unsigned int)(num_alignments + 1) * sizeof(bam_data_core_t), 0) );
            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_strand_counter_p, (unsigned int) num_blocks * sizeof(int), 0) );
            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_map_quality_p, (unsigned int) num_blocks * sizeof(int), 0) );
            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_alignment_length_p, (unsigned int) num_blocks * sizeof(int), 0) );
            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_cigar_data_p, (unsigned int) bam_data_batch_p->num_cigar_operations * sizeof(uint32_t), 0) );
            CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_qc_alignment_p, (unsigned int) num_alignments * sizeof(qc_alignment_t), 0) );

            CUDA_SAFE_CALL( cudaMemcpy(d_core_data_p, bam_data_batch_p->core_data_p, (num_alignments + 1) * sizeof(bam_data_core_t), cudaMemcpyHostToDevice) );
            CUDA_SAFE_CALL( cudaMemcpy(d_cigar_data_p, bam_data_batch_p->cigar_data_p, bam_data_batch_p->num_cigar_operations * sizeof(uint32_t), cudaMemcpyHostToDevice) );

            CUDA_START_TIMER();
            call_kernel_basic_stats(dimGrid, dimBlock, d_core_data_p, d_strand_counter_p, d_map_quality_p, d_alignment_length_p, num_alignments);
            call_kernel_map_errors(dimGrid, dimBlock, d_core_data_p, d_cigar_data_p, d_qc_alignment_p, num_alignments);
            CUDA_STOP_TIMER();

            // copy result from GPU (GPU -> CPU)
            CUDA_SAFE_CALL( cudaMemcpy(strand_counter_p, d_strand_counter_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
            CUDA_SAFE_CALL( cudaMemcpy(map_quality_p, d_map_quality_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
            CUDA_SAFE_CALL( cudaMemcpy(alignment_length_p, d_alignment_length_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
            CUDA_SAFE_CALL( cudaMemcpy(qc_alignment_p, d_qc_alignment_p, num_alignments * sizeof(qc_alignment_t), cudaMemcpyDeviceToHost) );

            // free device memory
            CUDA_SAFE_CALL( cudaFreeHost(d_core_data_p) );
            CUDA_SAFE_CALL( cudaFreeHost(d_cigar_data_p) );
            CUDA_SAFE_CALL( cudaFreeHost(d_strand_counter_p) );
            CUDA_SAFE_CALL( cudaFreeHost(d_map_quality_p) );
            CUDA_SAFE_CALL( cudaFreeHost(d_qc_alignment_p) );
        } else {
            // accumulation of partial results is only made once
            num_blocks = 1;

            strand_counter_p = (int*) calloc(num_blocks, sizeof(int));
            map_quality_p = (int*) calloc(num_blocks, sizeof(int));
            alignment_length_p = (int*) calloc(num_blocks, sizeof(int));
            qc_alignment_p = (qc_alignment_t*) calloc(bam_data_batch_p->num_alignments, sizeof(qc_alignment_t));

            if (time_flag) {
                start_timer(t1_gpu);
            }
            cpu_bam_qc_basic_stats(bam_data_batch_p->core_data_p, strand_counter_p, map_quality_p, alignment_length_p, num_alignments, cpu_num_threads);
            cpu_bam_qc_map_errors(bam_data_batch_p->core_data_p, bam_data_batch_p->cigar_data_p, qc_alignment_p, num_alignments);
            if (time_flag) {
                stop_timer(t1_gpu, t2_gpu, gpu_time);
            }
        }

        //if (time_flag) { 
        //    stop_timer(t1_gpu, t2_gpu, gpu_time); 
        //}

        // create a new qc_batch object
        bam_qc_batch_p = (bam_qc_batch_t*) malloc(sizeof(bam_qc_batch_t));
        bam_qc_batch_p->id = bam_data_batch_list_item_p->id;
        bam_qc_batch_p->num_alignments = bam_data_batch_p->num_alignments;
        bam_qc_batch_p->num_blocks = num_blocks;
        bam_qc_batch_p->qc_alignment_p = qc_alignment_p;
        bam_qc_batch_p->strand_counter_p = strand_counter_p;
        bam_qc_batch_p->map_quality_p = map_quality_p;
        bam_qc_batch_p->alignment_length_p = alignment_length_p;
        //bam_qc_batch_p->alignments_p = bam_data_batch_list_item_p->alignments_p;

        // and insert it into the bam_qc_batch_list
        list_item_t* item_p = list_item_new(bam_qc_batch_p->id, 0, bam_qc_batch_p);
        list_insert_item(item_p, &bam_qc_batch_list);

        // copy the the current batch item to the cpu batch list in order to perform CPU qc operations
        list_insert_item(bam_data_batch_list_item_p, cpu_batch_list_p);

        sprintf(log_message, "Thread-GPU:...processing for batch %i done !\n", bam_data_batch_list_item_p->id);
        LOG_DEBUG(log_message);

        //if (time_flag) { stop_timer(t1_gpu, t2_gpu, gpu_time); }
    } // end of external while loop

    pthread_mutex_lock(&gpus_thread_alive_lock);
    gpus_thread_alive--;
    pthread_mutex_unlock(&gpus_thread_alive_lock);

    list_decr_writers(cpu_batch_list_p);
    list_decr_writers(&bam_qc_batch_list);

    if (time_flag) {
        stop_timer(t1_qc_calc_server, t2_qc_calc_server, qc_calc_server_time);
    }
    LOG_DEBUG("Thread-GPU: END\n");

    // exiting...
    pthread_exit(0);
}

/* ******************************************************
 *    		QC calc server thread	 		*
 * *****************************************************/

void* cpus_server(void* params_p) {
    double coverage_time = 0.0;
    struct timeval t1_coverage, t2_coverage;

    LOG_DEBUG("Thread-CPU: START\n");

    if (time_flag) {
        start_timer(t1_cpus_server);
    }

    //initialize str_coverage_matrix
    str_coverage_matrix_init();

    cpus_server_input_t* input_p = (cpus_server_input_t*) params_p;

    qc_mapping_counter_t* qc_mapping_counter_p = (qc_mapping_counter_t*) input_p->qc_mapping_counter;
    int max_distance_size = input_p->max_distance_size;
    int cpu_num_threads =  input_p->cpu_num_threads;
    bam_data_batch_t* bam_data_batch_p = NULL;
    list_item_t* bam_data_batch_list_item_p = NULL;
    list_t* cpu_batch_list_p = input_p->cpu_batch_list_p;
    char* gff_filename = input_p->gff_filename;
    char* output_directory = input_p->output_directory;
    char* input_filename = input_p->input_filename;

    // delete previous coverage file to append new data
    bam_coverage_counter_delete_file(output_directory, input_filename);

    // variables for store intermediate and output results in both CPU and GPU
    qc_hash_t* qc_hash_p = (qc_hash_t*) qc_hash_new(QC_HASH_LENGTH);

    // variables for coverage (regions data)
    bam_chromosome_coverage_t bam_chromosome_coverage[num_of_chromosomes];

    for (int j = 0; j < num_of_chromosomes; j++) {
        bam_chromosome_coverage_init(&bam_chromosome_coverage[j]);
    }

    gff_data_t* gff_data_p = gff_data_new(gff_filename);

    int gpus_alive;

    gpus_alive = list_get_writers(cpu_batch_list_p);

    while ((bam_data_batch_list_item_p = list_remove_item(cpu_batch_list_p)) != NULL) {
        if ((time_flag) && (cpus_standby_time == 0.0)) {
            stop_timer(t1_active_reader, t1_active_cpus, cpus_standby_time);
        }

        char log_message[50];
        sprintf(log_message, "Thread-CPU: processing for batch %i....\n", bam_data_batch_list_item_p->id);
        LOG_DEBUG(log_message);

        // allocation memory for output results
        bam_data_batch_p = (bam_data_batch_t*) bam_data_batch_list_item_p->data_p;
        int num_alignments = bam_data_batch_p->num_alignments;
        int cpu_num_threads = input_p->cpu_num_threads;

        if (time_flag) {
            start_timer(t1_cpu);
        }

        char* id_seq;
        int tid, start_coordinate, seq_length;
        short int paired_end;
        bam_data_core_t* core_data_p;

        for (int i = 0; i < bam_data_batch_p->num_alignments; i++) {
            id_seq = &(bam_data_batch_p->id_seq_data_p[bam_data_batch_p->core_data_p[i].id_seq_index]);

            core_data_p = &(bam_data_batch_p->core_data_p[i]);

            tid = core_data_p->chromosome;
            start_coordinate = core_data_p->start_coordinate;
            seq_length = core_data_p->alignment_length;
            paired_end = core_data_p->paired_end;

            qc_hash_insert_alignment(qc_hash_p, id_seq, tid, start_coordinate, seq_length, paired_end);
        }

        //if (time_flag) { 
        //    start_timer(t1_coverage); 
        //}

        if (gff_data_batch_in_region(bam_data_batch_p, gff_data_p) != 0) {
            bam_coverage_compute(bam_data_batch_p, bam_chromosome_coverage, gff_data_p, output_directory, input_filename, cpu_num_threads);
        }

        //if (time_flag) { 
        //    stop_timer(t1_coverage, t2_coverage, coverage_time); 
        //}

        if (time_flag) {
            stop_timer(t1_cpu, t2_cpu, cpu_time);
        }

        sprintf(log_message, "Thread-CPU:...processing for batch %i done !\n", bam_data_batch_list_item_p->id);
        LOG_DEBUG(log_message);
 
        // free the current batch item if all processing with the batch is performed
        bam_data_batch_free((bam_data_batch_t*) bam_data_batch_list_item_p->data_p);
        list_item_free(bam_data_batch_list_item_p);
 
        // ask again for reads server status
        gpus_alive = list_get_writers(cpu_batch_list_p);
    } // end of external while loop

    // print the last counters
    bam_coverage_counter_mark_to_print(bam_chromosome_coverage, true);
    bam_coverage_counter_print(bam_chromosome_coverage, output_directory, input_filename);

    //qc_hash_list_print(qc_hash_p->qc_hash_list_p);

    //calculate over the qc hash table to obtain:
    //    - Mean distance between paired ends
    //    - Histogram of mappings per reads

    unsigned long mean_paired_end_distance = 0;

    if (time_flag) {
        start_timer(t1_cpu);
    }

    qc_hash_perform_calculations(qc_hash_p, qc_mapping_counter_p, &mean_paired_end_distance, max_distance_size, cpu_num_threads); //--

    qc_mapping_counter_p->mean_paired_end_distance = mean_paired_end_distance;
    if (time_flag) {
        stop_timer(t1_cpu, t2_cpu, cpu_time);
    }

    //free qc hash structure, gff data and chromosome coverage
    for (int j = 0; j < num_of_chromosomes; j++) {
        bam_chromosome_coverage_clear(&bam_chromosome_coverage[j]);
    }

    qc_hash_free(qc_hash_p, true);
    gff_data_free(gff_data_p);

    // --------------- D E B U G ----------------

    printf("--------------- D E B U G ----------------\n");

    for (int i = 0; i <= (MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
        printf("qc_mapping_counter_p->num_mappings_histogram[%i]: %i\n", i, qc_mapping_counter_p->num_mappings_histogram[i]);
    }
    printf("mean_paired_end_distance: %ld\n\n", qc_mapping_counter_p->mean_paired_end_distance);

    printf("--------------- D E B U G ----------------\n");

    // --------------- D E B U G ----------------

    pthread_mutex_lock(&cpus_thread_alive_lock);
    cpus_thread_alive--;
    pthread_mutex_unlock(&cpus_thread_alive_lock);

    if (time_flag) {
        stop_timer(t1_cpus_server, t2_cpus_server, cpus_server_time);
    }

    LOG_DEBUG("Thread-CPU: END\n");

    // exiting...
    pthread_exit(0);
}

/* ******************************************************
 *    		Results server thread	 		*
 * *****************************************************/

void* results_server(void* params_p) {
    LOG_DEBUG("Thread-RESULTS: START\n");

    if (time_flag) {
        start_timer(t1_results_server);
    }

    results_server_input_t* input_p = (results_server_input_t*) params_p;

    // variables for storing qc report information
    bam_qc_report_t bam_qc_report;
    memset(&bam_qc_report, 0, sizeof(bam_qc_report_t));

    qc_mapping_counter_t* qc_mapping_counter_p = (qc_mapping_counter_t*) input_p->qc_mapping_counter;
    int nb_total_threads = input_p->gpu_num_blocks * input_p->gpu_num_threads;
    int base_quality = input_p->base_quality;
    int i, alignments;

    // go through the results_batch list, and process it
    // take and remove the first item, and so on...
    list_item_t* item_p = NULL;
    bam_qc_batch_t* bam_qc_batch_p = NULL;
    int gpus_alive, cpus_alive;

    // getting gpus thread status
    pthread_mutex_lock(&gpus_thread_alive_lock);
    gpus_alive = gpus_thread_alive;
    pthread_mutex_unlock(&gpus_thread_alive_lock);

    pthread_mutex_lock(&cpus_thread_alive_lock);
    cpus_alive = cpus_thread_alive;
    pthread_mutex_unlock(&cpus_thread_alive_lock);

    //iteration until not NULL is returned, then process batch
    while ((item_p = list_remove_item(&bam_qc_batch_list)) != NULL) {
        bam_qc_batch_p = (bam_qc_batch_t*) item_p->data_p;
        //printf("while... gpus_alive: %i, cpus_alive: %i, bam_qc_batch_p is NULL: %i\n", gpus_alive, cpus_alive, (bam_qc_batch_p == NULL) ? 1:0);
        if ((time_flag) && (results_standby_time == 0.0)) {
            stop_timer(t1_active_reader, t1_active_results, results_standby_time);
        }
        if (time_flag) {
            start_timer(t1_result);
        }

        char log_message[50];
        sprintf(log_message, "Thread-RESULTS: processing for bam batch %i....\n", bam_qc_batch_p->id);
        LOG_DEBUG(log_message);

        // result processing batch per batch
        alignments = bam_qc_batch_p->num_alignments;
        bam_qc_report.num_alignments += alignments;

        for (int k = 0; k < bam_qc_batch_p->num_blocks; k++) {
            bam_qc_report.strand_counter += bam_qc_batch_p->strand_counter_p[k];
            bam_qc_report.mean_map_quality += bam_qc_batch_p->map_quality_p[k];
            bam_qc_report.mean_alignment_length += bam_qc_batch_p->alignment_length_p[k];
        }

        for (int k = 0; k < bam_qc_batch_p->num_alignments; k++) {
            bam_qc_report.map_error_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[MISMATCHES] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[MISMATCHES] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
            bam_qc_report.map_deletion_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[D] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[D] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
            bam_qc_report.map_insertion_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[I] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[I] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
            bam_qc_report.map_matching_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[EQUAL] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[EQUAL] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
        }

        sprintf(log_message, "Thread-RESULTS: ....processing for batch %i done !\n", bam_qc_batch_p->id);
        LOG_DEBUG(log_message);

        // free ALL memory
        bam_qc_batch_free(bam_qc_batch_p, false);

        if (time_flag) {
            stop_timer(t1_result, t2_result, result_time);
        }


        // getting gpus and cpus thread status
        pthread_mutex_lock(&gpus_thread_alive_lock);
        gpus_alive = gpus_thread_alive;
        pthread_mutex_unlock(&gpus_thread_alive_lock);

        pthread_mutex_lock(&cpus_thread_alive_lock);
        cpus_alive = cpus_thread_alive;
        pthread_mutex_unlock(&cpus_thread_alive_lock);
    } // end of batch loop

    printf("bam_qc_report.num_alignments: %lu, strand (+): %i, strand (-): %i\n", bam_qc_report.num_alignments, bam_qc_report.strand_counter, (bam_qc_report.num_alignments - bam_qc_report.strand_counter));

    if (time_flag) {
        start_timer(t1_result);
    }

    //calculate mean quality and mean length per alignment
    if (bam_qc_report.num_alignments > 0) {
        printf("bam_qc_report.mean_read_quality: %lu, num_alignments: %lu, mean_quality: %i\n", bam_qc_report.mean_map_quality, bam_qc_report.num_alignments, (bam_qc_report.mean_map_quality / bam_qc_report.num_alignments));
        printf("bam_qc_report.mean_alignment_length: %lu, num_alignments: %lu, mean_alignment_length: %i\n", bam_qc_report.mean_alignment_length, bam_qc_report.num_alignments, (bam_qc_report.mean_alignment_length / bam_qc_report.num_alignments));
        bam_qc_report.mean_map_quality /= bam_qc_report.num_alignments;
        bam_qc_report.mean_alignment_length /= bam_qc_report.num_alignments;
    } else {
        printf("bam_qc_report.mean_read_quality: %lu, num_alignments: %lu, mean_quality: 0\n", bam_qc_report.mean_map_quality, bam_qc_report.num_alignments);
        printf("bam_qc_report.mean_alignment_length: %lu, num_alignments: %lu, mean_alignment_length: 0\n", bam_qc_report.mean_alignment_length, bam_qc_report.num_alignments);
        bam_qc_report.mean_map_quality = 0;
        bam_qc_report.mean_alignment_length = 0;
    }

    if (time_flag) {
        stop_timer(t1_result, t2_result, result_time);
    }

    // and finally, print qc report, data files and graphs
    // when cpu data is ready (cpus_alive = 0)
    while (cpus_alive > 0) {
        sched_yield();

        usleep(10000);

        pthread_mutex_lock(&cpus_thread_alive_lock);
        cpus_alive = cpus_thread_alive;
        pthread_mutex_unlock(&cpus_thread_alive_lock);
    }

    if (time_flag) {
        start_timer(t1_result);
    }

    bam_qc_report.num_mappings_histogram = qc_mapping_counter_p->num_mappings_histogram;
    bam_qc_report.mean_paired_end_distance = qc_mapping_counter_p->mean_paired_end_distance;

    if (time_flag) {
        stop_timer(t1_result, t2_result, result_time);
    }

    if (time_flag) {
        start_timer(t1_reporting);
    }
    generate_report(bam_qc_report, input_p->filename, input_p->base_quality, input_p->report_directory, 1);
    if (time_flag) {
        stop_timer(t1_reporting, t2_reporting, reporting_time);
    }

    if (time_flag) {
        stop_timer(t1_results_server, t2_results_server, results_server_time);
    }
    LOG_DEBUG("Thread-RESULTS: END\n");

    // exiting...
    pthread_exit(0);
}

/* **************************************************************
 *    		Public functions implementations 		*
 * *************************************************************/

void qc_bam_file(size_t batch_size, int batch_list_size, int gpu_num_threads, int gpu_num_blocks, int cpu_num_threads, int base_quality, int max_distance_size, char* input_filename, char* output_directory, char* gff_filename) {
    // number of GPUs is obtained, and initializes the number of GPU threads 'alive'
    int num_gpu_devices;
    cudaError_t cudaResultCode = cudaGetDeviceCount(&num_gpu_devices);
    if (cudaResultCode != cudaSuccess) {
        num_gpu_devices = 0;
    }
    gpus_thread_alive = num_gpu_devices;

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (!prop.canMapHostMemory) {
        LOG_FATAL("device does not support MapHostMemory\n");
    }

    //initializing bam_data_batch_list_gpu, bam_data_batch_list_cpu and bam_qc_batch_list
    list_t bam_data_batch_list_gpu;
    list_t bam_data_batch_list_cpu;

    if (num_gpu_devices > 0) {
        list_init("bam_data_batch_list_cpu", num_gpu_devices, batch_list_size, &bam_data_batch_list_cpu);
    } else {
        list_init("bam_data_batch_list_cpu", 1, batch_list_size, &bam_data_batch_list_cpu);
    }
    list_init("bam_data_batch_list_gpu", 1, batch_list_size, &bam_data_batch_list_gpu);
    list_init("bam_qc_batch_list", ((num_gpu_devices > 0) ? num_gpu_devices : 1), batch_list_size, &bam_qc_batch_list);

    //initializing qc_mapping_counter
    qc_mapping_counter_t qc_mapping_counter;
    qc_mapping_counter_init(&qc_mapping_counter);

    //multi-threads
    bam_reader_t* bam_reader_p = bam_reader_by_batch_new(input_filename, batch_size, base_quality, &bam_data_batch_list_gpu, LIST_INSERT_MODE);
    bam_reader_start(bam_reader_p);

    //some local variables
    void* r;

    // multi-threads
    pthread_t* qc_calc_server_thread_p = (pthread_t*) malloc(((num_gpu_devices == 0) ? 1 : num_gpu_devices) * sizeof(pthread_t));
    pthread_t* cpus_server_thread_p = (pthread_t*) malloc(cpu_num_threads * sizeof(pthread_t));
    pthread_t results_server_thread;

    //calling GPU threads to process the bam data,
    //but first, prepare input parameter
    int i;
    qc_calc_server_input_t** qc_calc_server_input_p = (qc_calc_server_input_t**) calloc(num_gpu_devices, sizeof(qc_calc_server_input_t*));

    if (num_gpu_devices > 0) {
        for (i = 0; i < num_gpu_devices; i++) {
            qc_calc_server_input_p[i] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
        }
    } else {
        qc_calc_server_input_p[0] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
    }

    if (num_gpu_devices > 0) { //GPU implementacion
        for (i = 0; i < num_gpu_devices; i++) {
            qc_calc_server_input_p[i]->num_gpu_devices = num_gpu_devices;
            qc_calc_server_input_p[i]->cpu_num_threads = 0;
            qc_calc_server_input_p[i]->gpu_device_id[0] = i;
            qc_calc_server_input_p[i]->gpu_num_blocks = gpu_num_blocks;
            qc_calc_server_input_p[i]->gpu_num_threads = gpu_num_threads;
            qc_calc_server_input_p[i]->gpu_batch_list_p = &bam_data_batch_list_gpu;
            qc_calc_server_input_p[i]->cpu_batch_list_p = &bam_data_batch_list_cpu;
            pthread_create(&qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
        }
    } else {   //CPU implementacion
        qc_calc_server_input_p[0]->num_gpu_devices = 0;
        qc_calc_server_input_p[0]->cpu_num_threads = cpu_num_threads;
        qc_calc_server_input_p[0]->gpu_device_id[0] = 0;
        qc_calc_server_input_p[0]->gpu_num_blocks = 0;
        qc_calc_server_input_p[0]->gpu_num_threads = 0;
        qc_calc_server_input_p[0]->gpu_batch_list_p = &bam_data_batch_list_gpu;
        qc_calc_server_input_p[0]->cpu_batch_list_p = &bam_data_batch_list_cpu;
        pthread_create(&qc_calc_server_thread_p[0], NULL, qc_calc_server, (void*) &qc_calc_server_input_p[0]);
    }

    // calling CPU threads to process the bam data,
    cpus_server_input_t* cpus_server_input_p = (cpus_server_input_t*) calloc(1, sizeof(cpus_server_input_t));

    for (i = 0; i < 1; i++) {
        cpus_server_input_p->cpu_num_threads = cpu_num_threads;
        cpus_server_input_p->max_distance_size = max_distance_size;
        cpus_server_input_p->cpu_batch_list_p = &bam_data_batch_list_cpu;
        cpus_server_input_p->qc_mapping_counter = &qc_mapping_counter;
        cpus_server_input_p->gff_filename = gff_filename;
        cpus_server_input_p->output_directory = output_directory;
        cpus_server_input_p->input_filename = input_filename;

        pthread_create(&cpus_server_thread_p[i], NULL, cpus_server, (void*) cpus_server_input_p);
    }

    //calling thread to process results from GPU,
    results_server_input_t results_server_input;
    results_server_input.gpu_num_blocks = gpu_num_blocks;
    results_server_input.gpu_num_threads = gpu_num_threads;
    results_server_input.base_quality = base_quality;
    results_server_input.qc_mapping_counter = &qc_mapping_counter;
    results_server_input.filename = input_filename;
    results_server_input.report_directory = output_directory;
    pthread_create(&results_server_thread, NULL, results_server, (void*) &results_server_input);
    num_alignments = bam_reader_join(bam_reader_p);

    for (int i = 0; i < num_gpu_devices; i++) {
        pthread_join(qc_calc_server_thread_p[i], &r);
    }
    free(qc_calc_server_thread_p);

    for (int i = 0; i < 1; i++) {
        pthread_join(cpus_server_thread_p[i], &r);
    }

    free(cpus_server_thread_p);
    pthread_join(results_server_thread, &r);

    //free thread stuff and parameters
    if (num_gpu_devices > 0) {
        for (i = 0; i < num_gpu_devices; i++) {
            free(qc_calc_server_input_p[i]);
        }
    } else {
        free(qc_calc_server_input_p[0]);
    }

    free(qc_calc_server_input_p);
    free(cpus_server_input_p);
    bam_reader_free(bam_reader_p);
}

#endif /* QC_CU */
