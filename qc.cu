/*
 *
 *  Created on: Aug 4, 2011
 *      Author: victor
 */

#ifndef QC_CU_
#define QC_CU_

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
    #include "list.h"    
    #include "log.h"
    #include "qc.h"
    #include "qc_hash.h"
    #include "qc_kernel_omp.h"
    #include "sam.h"
    #include "system_utils.h"
}

#include "qc_kernel_cuda.h"

//------------------------------------------------------------------------------------

// global variables for qc process

// bam_qc_batch_list_t bam_qc_batch_list;
list_t bam_qc_batch_list;

int bam_batch_reader_alive = 1;
int gpus_thread_alive = 1;
int cpus_thread_alive = 1;

pthread_mutex_t gpus_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t cpus_thread_alive_lock = PTHREAD_MUTEX_INITIALIZER;

//====================================================================================
// thread functions
//====================================================================================

// structures for gpu server thread
//
typedef struct qc_calc_server_input {
  int num_gpu_devices;
  int cpu_num_threads;
  int gpu_device_id[256];
  int nb_reads_per_batch;
  int gpu_num_blocks;
  int gpu_num_threads;
  list_t* gpu_batch_list_p;
  list_t* cpu_batch_list_p;
} qc_calc_server_input_t;

// structures for cpu server thread
//
typedef struct cpus_server_input {
  int cpu_num_threads;
  int max_distance_size;
  list* cpu_batch_list_p;
  qc_mapping_counter_t* qc_mapping_counter;
  char* gff_filename;
  char* output_directory;
  char* input_filename;
} cpus_server_input_t;

// structures for results server thread
//
typedef struct results_server_input {
  int gpu_num_blocks;
  int gpu_num_threads;
  int base_quality;
  qc_mapping_counter_t* qc_mapping_counter;
  char* filename;
  char* report_directory;
} results_server_input_t;

// threads functions

void* qc_calc_server(void* params_p);
void* cpus_server(void* params_p);
void* results_server(void* params_p);


//-----------------------------------------------------
// qc_calc_server,
//
// this thread gets bam data from the bam data batch list,
// copy them to GPU if exists, execute kernel and insert the 
// results into the qc batch list. If not exist, the
// same calculations are performed in CPU using OMP
//-----------------------------------------------------

extern void call_kernel_basic_stats(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, qc_info_t* d_qc_info_p, int* d_strand_counter_p, int* d_alignment_length_p, int* d_map_quality_p, int num_alignments);
extern void call_kernel_map_errors(dim3 dimGrid, dim3 dimBlock, bam_data_core_t* d_core_data_p, uint32_t* d_cigar_data_p, qc_alignment_t* d_qc_alignment_p, int num_alignments);

void* qc_calc_server(void* params_p) {
	
  LOG_DEBUG("Thread-GPU: START\n");
  
  if (time_flag) { start_timer(t1_qc_calc_server); }
  
  qc_calc_server_input_t* input_p = (qc_calc_server_input_t*) params_p;
  
  int cpu_num_threads = input_p->cpu_num_threads;
  list_item_t* bam_data_batch_list_item_p = NULL;
  list_t* gpu_batch_list_p = input_p->gpu_batch_list_p;
  list_t* cpu_batch_list_p = input_p->cpu_batch_list_p;
  bam_data_batch_t* bam_data_batch_p = NULL;
  bam_qc_batch_t* bam_qc_batch_p = NULL;
  
  // variables for store output results in both CPU and GPU
  //
  qc_info_t* qc_info_p;
  qc_alignment_t* qc_alignment_p;
  int* strand_counter_p;
  int* map_quality_p;
  int* alignment_length_p;

  bam_data_core_t* d_core_data_p;
  uint32_t* d_cigar_data_p;
  qc_info_t* d_qc_info_p;
  qc_alignment_t* d_qc_alignment_p;
  int* d_strand_counter_p;
  int* d_map_quality_p;
  int* d_alignment_length_p;
  
  // selecting GPU device
  //  
  CUDA_SAFE_CALL( cudaSetDevice(input_p->gpu_device_id[0]) );	
  
  int reads_alive;
  
  //reads_alive = bam_data_batch_list_get_producers(gpu_batch_list_p);
  reads_alive = list_get_writers(gpu_batch_list_p);
  
  
  //bam_data_batch_list_print(batch_list_p);
  
  //bam_data_batch_list_item_p = bam_data_batch_list_remove(gpu_batch_list_p);
  bam_data_batch_list_item_p = list_remove_item(gpu_batch_list_p);

  while (reads_alive>0 || bam_data_batch_list_item_p!=NULL) {
    LOG_DEBUG("Thread-GPU: waiting for batch....\n");
    //printf("reads_alive: %i, bam_data_batch_list_item_p is NULL: %i\n", reads_alive, (bam_data_batch_list_item_p == NULL));
    
    if (bam_data_batch_list_item_p==NULL) {
      //printf("Thread-GPU: waiting 0.1 s for reads....\n");
      
      sched_yield();
      
      // Delay for a bit
      //struct timespec ts;
      //ts.tv_sec = 0;
      //ts.tv_nsec = 100000000;
      //nanosleep (&ts, NULL);
      
      usleep(10000);
      
    } else {
printf("1 ---->\n");      
     
      //if (time_flag) { start_timer(t1_gpu); }
      
      if ((time_flag) && (gpus_standby_time == 0.0)) { stop_timer(t1_active_reader, t1_active_gpus, gpus_standby_time); }

      char log_message[50];
      sprintf(log_message, "Thread-GPU: processing for batch %i....\n", bam_data_batch_list_item_p->id);
      LOG_DEBUG(log_message);

      number_of_batchs++;
 
      // allocation memory for output results
      //
      qc_info_p = (qc_info_t*) calloc(1, sizeof(qc_info_t));
      
      //int num_alignments = bam_data_batch_list_item_p->data_p->num_alignments;
      bam_data_batch_p = (bam_data_batch_t*) bam_data_batch_list_item_p->data_p;
      int num_alignments = bam_data_batch_p->num_alignments;      
      int num_blocks;
printf("2 ---->\n");

      if (cpu_num_threads == 0) {		// GPU implementation

	num_blocks = (num_alignments / input_p->gpu_num_threads) + 1;
	
	dim3 dimBlock(input_p->gpu_num_threads, 1, 1);
	dim3 dimGrid(num_blocks, 1, 1);

	strand_counter_p = (int*) calloc(num_blocks, sizeof(int));
	map_quality_p = (int*) calloc(num_blocks, sizeof(int));
	alignment_length_p = (int*) calloc(num_blocks, sizeof(int));
	qc_alignment_p = (qc_alignment_t*) calloc(bam_data_batch_p->num_alignments, sizeof(qc_alignment_t));
printf("3 ---->\n");
printf("3.0 ---->\n");
	//printf("bam data batch id %i, num_alignments: %i\n", bam_data_batch_list_item_p->id, bam_data_batch_list_item_p->batch_p->num_alignments);
	
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_core_data_p, (unsigned int) (num_alignments + 1) * sizeof(bam_data_core_t), 0) );
printf("3.1 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_qc_info_p, (unsigned int) sizeof(qc_info_t), 0) );
printf("3.2 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_strand_counter_p, (unsigned int) num_blocks * sizeof(int), 0) );
printf("3.3 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_map_quality_p, (unsigned int) num_blocks * sizeof(int), 0) );
printf("3.4 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_alignment_length_p, (unsigned int) num_blocks * sizeof(int), 0) );
printf("3.5 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_cigar_data_p, (unsigned int) bam_data_batch_p->num_cigar_operations * sizeof(uint32_t), 0) );
printf("3.6 ---->\n");
	CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_qc_alignment_p, (unsigned int) num_alignments * sizeof(qc_alignment_t), 0) );
printf("4 ---->\n");
	//printf("memory usage: data_size = %.2f MB, data_indices_size = %.2f MB, gpu_result = %.2f MB, gpu_kmers = %.2f MB\n", fastq_batch_list_item_p->batch_p->data_size / 1e6, fastq_batch_list_item_p->batch_p->data_indices_size / 1e6,  ((unsigned int) fastq_batch_list_item_p->batch_p->num_reads * sizeof(qc_read_t)) / 1e6, ((unsigned int) fastq_batch_list_item_p->batch_p->num_reads * KMERS_COMBINATIONS * sizeof(qc_kmers_t)) / 1e6);
	
	CUDA_SAFE_CALL( cudaMemset((void*) d_qc_info_p, 0, (unsigned int) sizeof(qc_info_t)) );

	CUDA_SAFE_CALL( cudaMemcpy(d_core_data_p, bam_data_batch_p->core_data_p, (num_alignments + 1) * sizeof(bam_data_core_t), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL( cudaMemcpy(d_cigar_data_p, bam_data_batch_p->cigar_data_p, bam_data_batch_p->num_cigar_operations * sizeof(uint32_t), cudaMemcpyHostToDevice) );
	
	CUDA_START_TIMER();
	call_kernel_basic_stats(dimGrid, dimBlock, d_core_data_p, d_qc_info_p, d_strand_counter_p, d_map_quality_p, d_alignment_length_p, num_alignments);
	call_kernel_map_errors(dimGrid, dimBlock, d_core_data_p, d_cigar_data_p, d_qc_alignment_p, num_alignments);
	CUDA_STOP_TIMER();
printf("5 ---->\n");	
  //       for (int j=0; j < bam_data_batch_list_item_p->batch_p->num_cigar_operations; j++) {
  // 	printf("cigar operation: %i, num nts: %i\n", (bam_data_batch_list_item_p->batch_p->cigar_data_p[j])&BAM_CIGAR_MASK, (bam_data_batch_list_item_p->batch_p->cigar_data_p[j])>>BAM_CIGAR_SHIFT);
  //       }
 
	// copy result from GPU (GPU -> CPU)
	//
	CUDA_SAFE_CALL( cudaMemcpy(qc_info_p, d_qc_info_p,  sizeof(qc_info_t), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(strand_counter_p, d_strand_counter_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(map_quality_p, d_map_quality_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(alignment_length_p, d_alignment_length_p, num_blocks * sizeof(int), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(qc_alignment_p, d_qc_alignment_p, num_alignments * sizeof(qc_alignment_t), cudaMemcpyDeviceToHost) );
printf("6 ---->\n");
	// free memory
	//
	CUDA_SAFE_CALL( cudaFreeHost(d_core_data_p) );
	CUDA_SAFE_CALL( cudaFreeHost(d_cigar_data_p) );
	CUDA_SAFE_CALL( cudaFreeHost(d_qc_info_p) );
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

	if (time_flag) { start_timer(t1_gpu); }
	cpu_bam_qc_basic_stats(bam_data_batch_p->core_data_p, strand_counter_p, map_quality_p, alignment_length_p, num_alignments, cpu_num_threads);
	cpu_bam_qc_map_errors(bam_data_batch_p->core_data_p, bam_data_batch_p->cigar_data_p, qc_alignment_p, num_alignments);
	if (time_flag) { stop_timer(t1_gpu, t2_gpu, gpu_time); }
      }
printf("7 ---->\n");
	//if (time_flag) { stop_timer(t1_gpu, t2_gpu, gpu_time); }
       
      // create a new qc_batch object
      //
      bam_qc_batch_p = (bam_qc_batch_t*) malloc(sizeof(bam_qc_batch_t));
      bam_qc_batch_p->id = bam_data_batch_list_item_p->id;
      bam_qc_batch_p->num_alignments = bam_data_batch_p->num_alignments;
      bam_qc_batch_p->num_blocks = num_blocks;
      bam_qc_batch_p->qc_info_p = qc_info_p;
      bam_qc_batch_p->qc_alignment_p = qc_alignment_p;
      bam_qc_batch_p->strand_counter_p = strand_counter_p;
      bam_qc_batch_p->map_quality_p = map_quality_p;
      bam_qc_batch_p->alignment_length_p = alignment_length_p;
      //bam_qc_batch_p->alignments_p = bam_data_batch_list_item_p->alignments_p;
printf("8 ---->\n");
      // and insert it into the bam_qc_batch_list
      //
      //bam_qc_batch_list_insert(bam_qc_batch_p, &bam_qc_batch_list);      
      list_item_t* item_p = list_item_new(bam_qc_batch_p->id, 0, bam_qc_batch_p);
      list_insert_item(item_p, &bam_qc_batch_list);
printf("-----> inserting in bam_qc_batch_list id: %i ---->\n", bam_qc_batch_p->id);
      
      //qc_batch_list_print(&qc_batch_list);
    
      // copy the the current batch item to the cpu batch list in order to perform CPU qc operations
      //bam_data_batch_list_insert(bam_data_batch_list_item_p, cpu_batch_list_p);
      list_insert_item(bam_data_batch_list_item_p, cpu_batch_list_p);
      //bam_data_batch_list_item_free(bam_data_batch_list_item_p, true);
printf("-----> inserting in cpu_batch_list_p ---->\n");
      sprintf(log_message, "Thread-GPU:...processing for batch %i done !\n", bam_data_batch_list_item_p->id);
      LOG_DEBUG(log_message);

      //if (time_flag) { stop_timer(t1_gpu, t2_gpu, gpu_time); }
    } // end if-else
    
    // ask again for reads server status
    //		
    //if (time_flag) { start_timer(t1_gpu); }
    //reads_alive = bam_data_batch_list_get_producers(gpu_batch_list_p);
    reads_alive = list_get_writers(gpu_batch_list_p);

    // next batch...: the first in the list
    //    
    //bam_data_batch_list_item_p = bam_data_batch_list_remove(gpu_batch_list_p);
    bam_data_batch_list_item_p = list_remove_item(gpu_batch_list_p);
      
    //if (time_flag) { stop_timer(t1_gpu, t2_gpu, gpu_time); }
  } // end of external while loop

  pthread_mutex_lock(&gpus_thread_alive_lock);
  gpus_thread_alive--;
  pthread_mutex_unlock(&gpus_thread_alive_lock);
  
  //bam_data_batch_list_decr_producers(cpu_batch_list_p);
  list_decr_writers(cpu_batch_list_p);
  list_decr_writers(&bam_qc_batch_list);

  if (time_flag) { stop_timer(t1_qc_calc_server, t2_qc_calc_server, qc_calc_server_time); }
  LOG_DEBUG("Thread-GPU: END\n");
  
  // exiting....
  //
  pthread_exit(0);
}

//---------------------------------------------------------
// cpus_server,
//
// this thread gets bam data from the bam data batch list,
// and calculates number of duplicated alignments,  
// paired end distance and coverage
//---------------------------------------------------------

void* cpus_server(void* params_p) {

  double coverage_time = 0.0;
  struct timeval t1_coverage, t2_coverage;
  
  LOG_DEBUG("Thread-CPU: START\n");
  
  if (time_flag) { start_timer(t1_cpus_server); }
  
  //initialize str_coverage_matrix
  str_coverage_matrix_init();
  
  cpus_server_input_t* input_p = (cpus_server_input_t*) params_p;
  
  qc_mapping_counter_t* qc_mapping_counter_p = (qc_mapping_counter_t*) input_p->qc_mapping_counter;
  int max_distance_size = input_p->max_distance_size;
  int cpu_num_threads =  input_p->cpu_num_threads;
  //bam_data_batch_list_item_t* bam_data_batch_list_item_p = NULL;
  //bam_data_batch_list_t* cpu_batch_list_p = input_p->cpu_batch_list_p;
  bam_data_batch_t* bam_data_batch_p = NULL;
  list_item_t* bam_data_batch_list_item_p = NULL;
  list_t* cpu_batch_list_p = input_p->cpu_batch_list_p;  
  char* gff_filename = input_p->gff_filename;
  char* output_directory = input_p->output_directory;
  char* input_filename = input_p->input_filename;

  // delete previous coverage file to append new data
  bam_coverage_counter_delete_file(output_directory, input_filename);

  // variables for store intermediate and output results in both CPU and GPU
  //
  qc_hash_t* qc_hash_p = (qc_hash_t*) qc_hash_new(QC_HASH_LENGTH);
  
  // variables for coverage (regions data)
  //
  bam_chromosome_coverage_t bam_chromosome_coverage[NUM_OF_CHROMOSOMES];

  int j;
  for (j=0; j < NUM_OF_CHROMOSOMES; j++) {
    bam_chromosome_coverage_init(&bam_chromosome_coverage[j]); 
  }

  gff_data_t* gff_data_p = gff_data_new(gff_filename);
  //printf("num_gff_lines: %i\n", (gff_data_p == NULL) ? 0 : gff_data_p->num_regions);

  //gff_data_print_lines(gff_data_p);
  //gff_data_print_regions(gff_data_p);

  int gpus_alive;

  //gpus_alive = bam_data_batch_list_get_producers(cpu_batch_list_p);
  gpus_alive = list_get_writers(cpu_batch_list_p);
  
  //bam_data_batch_list_print(cpu_batch_list_p);
  
  //bam_data_batch_list_item_p = bam_data_batch_list_remove(cpu_batch_list_p);
  bam_data_batch_list_item_p = list_remove_item(cpu_batch_list_p);
  //bam_data_batch_list_item_p = NULL;
  
  // ---------------- T E S T ------------------
  int test_counter_alignments = 0;

  while (gpus_alive>0 || bam_data_batch_list_item_p!=NULL) {
    //printf("Thread-CPU: waiting for batch....\n");
    //printf("gpus_alive: %i, bam_data_batch_list_item_p is NULL: %i\n", gpus_alive, (bam_data_batch_list_item_p == NULL));

//     if (bam_data_batch_list_item_p==NULL) {
//       //printf("Thread-CPU: waiting 1 s for reads....\n");
// 
//       sched_yield();
//       
//       // Delay for a bit
//       //struct timespec ts;
//       //ts.tv_sec = 0;
//       //ts.tv_nsec = 100000000;
//       //nanosleep (&ts, NULL);
//       
//       usleep(10000);     	    
//     } else {
//        
//       if ((time_flag) && (cpus_standby_time == 0.0)) { stop_timer(t1_active_reader, t1_active_cpus, cpus_standby_time); }
// 		
// 	char log_message[50];
// 	sprintf(log_message, "Thread-CPU: processing for batch %i....\n", bam_data_batch_list_item_p->id);
// 	LOG_DEBUG(log_message);	
// 
// 	// allocation memory for output results
// 	//
// 	bam_data_batch_p = (bam_data_batch_t*) bam_data_batch_list_item_p->data_p;
// 	int num_alignments = bam_data_batch_p->num_alignments;
// 	int cpu_num_threads = input_p->cpu_num_threads;
// 
// 	if (time_flag) { start_timer(t1_cpu); }
// 	
// 	char* id_seq;
// 	int tid, start_coordinate, seq_length;
// 	short int paired_end;
// 	bam_data_core_t* core_data_p;
// 
// 	for (int i=0; i < bam_data_batch_p->num_alignments; i++) {
// 	  id_seq = &(bam_data_batch_p->id_seq_data_p[bam_data_batch_p->core_data_p[i].id_seq_index]);
// 	  
// 	  // ---------------- T E S T ------------------
// 	  //test_counter_alignments++;
// 	  //char* id_seq_test = (char*) calloc(strlen(id_seq) + 2, sizeof(char));
// 	  //strcpy(id_seq_test, id_seq);
// 	  //strcat(id_seq_test, (i%2) == 1 ? "/1" : "/2");
// 	  //printf("id_seq_test: %s\n", id_seq_test);	
// 	  // ---------------- T E S T ------------------
// 	  
// 	  core_data_p = &(bam_data_batch_p->core_data_p[i]);
// 
// 	  tid = core_data_p->chromosome;
// 	  start_coordinate = core_data_p->start_coordinate;
// 	  seq_length = core_data_p->alignment_length;
// 	  paired_end = core_data_p->paired_end;
// 
// 	  //printf("id seq: %s\n", id_seq);
// 	  
// 	  qc_hash_insert_alignment(qc_hash_p, id_seq, tid, start_coordinate, seq_length, paired_end);
// 
// 	  // ---------------- T E S T ------------------	
// 	  //qc_hash_insert_alignment(qc_hash_p, id_seq_test, tid, start_coordinate, seq_length);
// 	  // ---------------- T E S T ------------------	  
// 	}
// 
// 	//if (time_flag) { start_timer(t1_coverage); }
// 	
// 	if (gff_data_batch_in_region(bam_data_batch_p, gff_data_p) != 0) {
// 	    //printf("\nstart_position: %i\n", bam_data_batch_list_item_p->batch_p->start_positions[0]);
// 	    //printf("end_position: %i\n", bam_data_batch_list_item_p->batch_p->last_alignments_position);
// 	    //printf("chromosome: %i, region_start: %i, region_end: %i\n", gff_data_p->gff_regions_p[gff_data_p->actual_region].chromosome, gff_data_p->gff_regions_p[gff_data_p->actual_region].start, gff_data_p->gff_regions_p[gff_data_p->actual_region].end);
// 	    //printf("gff_data_batch_in_region: %i\n", gff_data_batch_in_region(bam_data_batch_list_item_p->batch_p, gff_data_p));
// 	  bam_coverage_compute(bam_data_batch_p, bam_chromosome_coverage, gff_data_p, output_directory, input_filename, cpu_num_threads);
// 	}
// 	
// 	//if (time_flag) { stop_timer(t1_coverage, t2_coverage, coverage_time); }
// 	//printf("total coverage time       (s): \t%10.5f\n", 0.000001 * coverage_time);
// 	//printf("write time                (s): \t%10.5f\n", 0.000001 * write_time);
//       
// 	if (time_flag) { stop_timer(t1_cpu, t2_cpu, cpu_time); }
// 
// 	sprintf(log_message, "Thread-CPU:...processing for batch %i done !\n", bam_data_batch_list_item_p->id);
// 	LOG_DEBUG(log_message);
// 
// 	// F R E E  R E S O U R C E S
// printf("5 ---->\n"); 
// 	// free the current batch item if all processing with the batch is performed
// 	//bam_data_batch_list_item_free(bam_data_batch_list_item_p, true);
// 	bam_data_batch_free((bam_data_batch_t*) bam_data_batch_list_item_p->data_p);
// 	list_item_free(bam_data_batch_list_item_p);
// 	//printf("Thread-CPU:...processing for batch %i done !\n", qc_batch_p->id); 
// 	
// 
//       
//       } // end if-else
//  
      // ask again for reads server status
      //		
      //gpus_alive = bam_data_batch_list_get_producers(cpu_batch_list_p);
      gpus_alive = list_get_writers(cpu_batch_list_p);
printf("6 ---->\n");    
      // next batch...: the first in the list
      //
      //bam_data_batch_list_item_p = bam_data_batch_list_remove(cpu_batch_list_p);
      bam_data_batch_list_item_p = list_remove_item(cpu_batch_list_p);
printf("7 ---->\n");
  } // end of external while loop

  // ---------------- T E S T ------------------	
    
    /*qc_hash_list_t* list_p;
    printf("qc_hash_p->length: %i\n", qc_hash_p->length);
    for (int l=0; l < qc_hash_p->length; l++) {
      list_p = &(qc_hash_p->qc_hash_list_p[l]);
      if (list_p->length == 0) continue;
      printf("list-length: %i\n", list_p->length);
    }*/
    
  // ---------------- T E S T ------------------	
  
  // print the last counters 
  bam_coverage_counter_mark_to_print(bam_chromosome_coverage, true);
  bam_coverage_counter_print(bam_chromosome_coverage, output_directory, input_filename);

  //qc_hash_list_print(qc_hash_p->qc_hash_list_p);

// ---------------- T E S T ------------------	
//   int list_length = 0;
//   int num_lists = 0;
//   int min_list_length = 1;
//   int max_list_length = 1;
//   int mean_list_length = 0;  
//   int count_alignments = 0;
//   
//   qc_hash_list_item_t* item_aux_p;
//   
//   for (int j=0; j < qc_hash_p->length; j++) {
//     //count_alignments += qc_hash_p->qc_hash_list_p[j].length;
//     
//     item_aux_p = qc_hash_p->qc_hash_list_p[j].first_p;
//     list_length = qc_hash_p->qc_hash_list_p[j].length;
//     
//     mean_list_length += list_length;
//     
//     if (list_length > 0) { num_lists++; }
//     
//     if ((list_length < min_list_length) && (list_length != 0)) {
//       min_list_length = list_length;      
//     }
// 
//     if (list_length > max_list_length)  {
//       max_list_length = list_length;
//     }
// 
//     for (int k=0; k < qc_hash_p->qc_hash_list_p[j].length; k++) {
//       count_alignments += item_aux_p->num_pairends1;
//       count_alignments += item_aux_p->num_pairends2;
//       item_aux_p = item_aux_p->next_p;
//     }
//   }
//   printf("\ncount_alignments: %i\n", count_alignments);
//   //printf("test_counter_alignments: %i\n\n", test_counter_alignments);
//   printf("num of lists: %i\n", num_lists);
//   printf("min list length: %i\n", min_list_length);
//   printf("mean list length: %f\n", 1.0 * mean_list_length / num_lists);
//   printf("max list length: %i\n", max_list_length);
// ---------------- T E S T ------------------	  
  
  //calculate over the qc hash table to obtain:
  //    - Mean distance between paired ends
  //    - Histogram of mappings per reads
  
  //unsigned int num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 2];
  unsigned long mean_paired_end_distance = 0;
      
  if (time_flag) { start_timer(t1_cpu); }
  //qc_hash_perform_calculations(qc_hash_p, num_mappings_histogram, &mean_paired_end_distance, max_distance_size, cpu_num_threads);
  qc_hash_perform_calculations(qc_hash_p, qc_mapping_counter_p, &mean_paired_end_distance, max_distance_size, cpu_num_threads); //--

  //memcpy(qc_mapping_counter_p->num_mappings_histogram, num_mappings_histogram, sizeof(num_mappings_histogram));  
  
  //int sum_map = 0;
  //for (int i=1; i<=(MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
  //  sum_map += (qc_mapping_counter_p->num_mappings_histogram[i] * i);
  //}

  //printf("mean_paired_end_distance: %li, sum_map: %i\n", mean_paired_end_distance, sum_map);

  //qc_mapping_counter_p->mean_paired_end_distance = mean_paired_end_distance / sum_map;
  qc_mapping_counter_p->mean_paired_end_distance = mean_paired_end_distance;
  if (time_flag) { stop_timer(t1_cpu, t2_cpu, cpu_time); }
  
  //free qc hash structure, gff data and chromosome coverage
  for (j=0; j < NUM_OF_CHROMOSOMES; j++) {
    bam_chromosome_coverage_clear(&bam_chromosome_coverage[j]); 
  }
  
  //bam_chromosome_coverage_clear(&bam_chromosome_coverage);
//   for (j=0; j < NUM_OF_CHROMOSOMES; j++) {
//     free(&bam_chromosome_coverage[j]); 
//   }

  qc_hash_free(qc_hash_p, true);
  gff_data_free(gff_data_p);  
    
  // --------------- D E B U G ----------------
  
  printf("--------------- D E B U G ----------------\n");
  
  for (int i=0; i<=(MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
    printf("qc_mapping_counter_p->num_mappings_histogram[%i]: %i\n", i, qc_mapping_counter_p->num_mappings_histogram[i]);
  }
  printf("mean_paired_end_distance: %ld\n\n", qc_mapping_counter_p->mean_paired_end_distance);
  //printf("mean_paired_end_distance: %ld\n\n", mean_paired_end_distance / sum_map);
  printf("--------------- D E B U G ----------------\n");
  
  // --------------- D E B U G ----------------
 
  pthread_mutex_lock(&cpus_thread_alive_lock);
  cpus_thread_alive--;
  pthread_mutex_unlock(&cpus_thread_alive_lock);

  if (time_flag) { stop_timer(t1_cpus_server, t2_cpus_server, cpus_server_time); }
  
  LOG_DEBUG("Thread-CPU: END\n");  

  // exiting....
  //
  pthread_exit(0);
}

//---------------------------------------------------------
// results_server,
//
// this thread gets qc results from the qc batch list,
// and makes the post calculations needed for intermediate  
// results. It also launches the reports generation
//---------------------------------------------------------

void* results_server(void* params_p) {

  LOG_DEBUG("Thread-RESULTS: START\n");
  
  if (time_flag) { start_timer(t1_results_server); }

  results_server_input_t* input_p = (results_server_input_t*) params_p;
 
  // variables for storing qc report information
  //
  bam_qc_report_t bam_qc_report;
  memset(&bam_qc_report, 0, sizeof(bam_qc_report_t));

  qc_mapping_counter_t* qc_mapping_counter_p = (qc_mapping_counter_t*) input_p->qc_mapping_counter;  
  int nb_total_threads = input_p->gpu_num_blocks * input_p->gpu_num_threads;
  int base_quality = input_p->base_quality;
  int i, alignments;
    
  // go through the results_batch list, and process it
  // take and remove the first item, and so on...
  //
  list_item_t* item_p = NULL;
  bam_qc_batch_t* bam_qc_batch_p = NULL;
  int gpus_alive, cpus_alive;
  
  // getting gpus thread status
  //
  pthread_mutex_lock(&gpus_thread_alive_lock);
  gpus_alive = gpus_thread_alive;
  pthread_mutex_unlock(&gpus_thread_alive_lock);

  pthread_mutex_lock(&cpus_thread_alive_lock);
  cpus_alive = cpus_thread_alive;
  pthread_mutex_unlock(&cpus_thread_alive_lock);
  
  //bam_qc_batch_list_print(&bam_qc_batch_list);
  
  // get the first element in the list
  //
  //bam_qc_batch_p = bam_qc_batch_list_remove(&bam_qc_batch_list);
  item_p = list_remove_item(&bam_qc_batch_list);
  bam_qc_batch_p = (bam_qc_batch_t*) item_p->data_p;  
  
  //tic("----> processing qc_batch");
  while (gpus_alive>0 || cpus_alive>0) {
    //printf("while... gpus_alive: %i, cpus_alive: %i, bam_qc_batch_p is NULL: %i\n", gpus_alive, cpus_alive, (bam_qc_batch_p == NULL) ? 1:0);
    if (bam_qc_batch_p==NULL) {
      //printf("Thread-RESULTS: waiting 1 s for GPU outputs....\n");
      
      sched_yield();
      
      // Delay for a bit
      //struct timespec ts;
      //ts.tv_sec = 0;
      //ts.tv_nsec = 100000000;
      //nanosleep (&ts, NULL);
      
      usleep(10000);
	    
    } else {
      
      //number_of_batchs++;
      if ((time_flag) && (results_standby_time == 0.0)) { stop_timer(t1_active_reader, t1_active_results, results_standby_time); }
      if (time_flag) { start_timer(t1_result); }
	      
      char log_message[50];
      sprintf(log_message, "Thread-RESULTS: processing for bam batch %i....\n", bam_qc_batch_p->id);
      LOG_DEBUG(log_message);

      // result processing batch per batch
      alignments = bam_qc_batch_p->num_alignments;
      bam_qc_report.num_alignments += alignments;
     
      for (int k=0; k<bam_qc_batch_p->num_blocks; k++) {
	bam_qc_report.strand_counter += bam_qc_batch_p->strand_counter_p[k];
	bam_qc_report.mean_map_quality += bam_qc_batch_p->map_quality_p[k];
	bam_qc_report.mean_alignment_length += bam_qc_batch_p->alignment_length_p[k];
	//printf("processing bam qc batch p id: %i, strand (+): %i, map quality: %i, alignment length: %i\n", bam_qc_batch_p->id, bam_qc_batch_p->strand_counter_p[k], bam_qc_batch_p->map_quality_p[k], bam_qc_batch_p->alignment_length_p[k]);
      }

      for (int k=0; k<bam_qc_batch_p->num_alignments; k++) {
	bam_qc_report.map_error_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[MISMATCHES] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[MISMATCHES] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
	bam_qc_report.map_deletion_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[D] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[D] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
	bam_qc_report.map_insertion_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[I] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[I] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;
	bam_qc_report.map_matching_histogram[(bam_qc_batch_p->qc_alignment_p[k].counters[EQUAL] <= MAX_MAP_ERRORS_IN_HISTOGRAM) ? bam_qc_batch_p->qc_alignment_p[k].counters[EQUAL] : (MAX_MAP_ERRORS_IN_HISTOGRAM + 1)]++;	
      }
      
      sprintf(log_message, "Thread-RESULTS: ....processing for batch %i done !\n",bam_qc_batch_p->id);
      LOG_DEBUG(log_message);
      
      // free ALL memory
      //
      bam_qc_batch_free(bam_qc_batch_p, false);

      if (time_flag) { stop_timer(t1_result, t2_result, result_time); }
    } //end of if-else
  
    // getting gpus and cpus thread status
    //
    pthread_mutex_lock(&gpus_thread_alive_lock);
    gpus_alive = gpus_thread_alive;
    pthread_mutex_unlock(&gpus_thread_alive_lock);
    
    pthread_mutex_lock(&cpus_thread_alive_lock);
    cpus_alive = cpus_thread_alive;
    pthread_mutex_unlock(&cpus_thread_alive_lock);
  
    // next batch...
    //
    //bam_qc_batch_p = bam_qc_batch_list_remove(&bam_qc_batch_list);
    item_p = list_remove_item(&bam_qc_batch_list);
    
    if (item_p != NULL) {
      bam_qc_batch_p = (bam_qc_batch_t*) item_p->data_p;
    } else {
      //bam_qc_batch_p = NULL;
      break;
    }
  } // end of batch loop

  printf("bam_qc_report.num_alignments: %i, strand (+): %i, strand (-): %i\n", bam_qc_report.num_alignments, bam_qc_report.strand_counter, (bam_qc_report.num_alignments - bam_qc_report.strand_counter));

  if (time_flag) { start_timer(t1_result); }

  // calculate mean quality and mean length per alignment
  //
  printf("bam_qc_report.mean_read_quality: %i, num_alignments: %i, mean_quality: %i\n", bam_qc_report.mean_map_quality, bam_qc_report.num_alignments, (bam_qc_report.mean_map_quality / bam_qc_report.num_alignments));
  printf("bam_qc_report.mean_alignment_length: %i, num_alignments: %i, mean_alignment_length: %i\n", bam_qc_report.mean_alignment_length, bam_qc_report.num_alignments, (bam_qc_report.mean_alignment_length / bam_qc_report.num_alignments));
     
  bam_qc_report.mean_map_quality /= bam_qc_report.num_alignments;
  bam_qc_report.mean_alignment_length /= bam_qc_report.num_alignments;
  
  if (time_flag) { stop_timer(t1_result, t2_result, result_time); }		

  // and finally, print qc report, data files and graphs
  // when cpu data is ready (cpus_alive = 0)

  while (cpus_alive>0) {
    //printf("waiting....\n");
    sched_yield();

    // Delay for a bit
    //struct timespec ts;
    //ts.tv_sec = 0;
    //ts.tv_nsec = 100000000;
    //nanosleep (&ts, NULL);	
    
    usleep(10000);
    
    pthread_mutex_lock(&cpus_thread_alive_lock);
    cpus_alive = cpus_thread_alive;
    pthread_mutex_unlock(&cpus_thread_alive_lock);
  }  

  if (time_flag) { start_timer(t1_result); }

  bam_qc_report.num_mappings_histogram = qc_mapping_counter_p->num_mappings_histogram;
  bam_qc_report.mean_paired_end_distance = qc_mapping_counter_p->mean_paired_end_distance;
  
  if (time_flag) { stop_timer(t1_result, t2_result, result_time); }

  // --------------- D E B U G ----------------
  
//   printf("--------------- D E B U G ----------------\n");
//   for (int i=0; i<=(MAX_MAPPING_COUNT_IN_HISTOGRAM + 1); i++) {
//     printf("bam_qc_report.num_mappings_histogram[%i]: %i\n", i, bam_qc_report.num_mappings_histogram[i]);
//   }
//   printf("bam_qc_report.mean_paired_end_distance: %ld\n\n", bam_qc_report.mean_paired_end_distance);
  
  // --------------- D E B U G ----------------

  if (time_flag) { start_timer(t1_reporting); }
  generate_report(bam_qc_report, input_p->filename, input_p->base_quality, input_p->report_directory, 1);
  if (time_flag) { stop_timer(t1_reporting, t2_reporting, reporting_time); }

  if (time_flag) { stop_timer(t1_results_server, t2_results_server, results_server_time); }
  LOG_DEBUG("Thread-RESULTS: END\n");

  // exiting....
  //
  pthread_exit(0);
   
}


/*
      QC BAM FILE
*/

void qc_bam_file(size_t batch_size, int batch_list_size, int gpu_num_threads, int gpu_num_blocks, int cpu_num_threads, int base_quality, int max_distance_size, char* input_filename, char* output_directory, char* gff_filename) {

/*  bam_data_core_t* d_core_data_aux_p;
  CUDA_SAFE_CALL( cudaHostAlloc((void**) &d_core_data_aux_p, (unsigned int) (num_alignments + 1) * sizeof(bam_data_core_t), 0) );*/
  
  // number of GPUs is obtained, and initializes the number of GPU threads 'alive'
  //
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
  
  //initializing bam_data_batch_list, bam_qc_batch_list and qc_mapping_counter
  //
  //bam_data_batch_list_t bam_data_batch_list_gpu;
  //bam_data_batch_list_t bam_data_batch_list_cpu;
 
  //bam_data_batch_list_init(&bam_data_batch_list_gpu, 0);
  //bam_data_batch_list_init(&bam_data_batch_list_cpu, 0);
  //bam_qc_batch_list_init(&bam_qc_batch_list);
  
  list_t bam_data_batch_list_gpu;
  list_t bam_data_batch_list_cpu;

  if (num_gpu_devices > 0) {    
    list_init("bam_data_batch_list_cpu", num_gpu_devices, batch_list_size, &bam_data_batch_list_cpu);
  } else {
    list_init("bam_data_batch_list_cpu", 1, batch_list_size, &bam_data_batch_list_cpu);
  }
  list_init("bam_data_batch_list_gpu", 1, batch_list_size, &bam_data_batch_list_gpu);
  list_init("bam_qc_batch_list", 1, batch_list_size, &bam_qc_batch_list);

  qc_mapping_counter_t qc_mapping_counter;	//PASS TO THE CPU AND RESULTS THREADS AND BUILD THE MUTEX HANDLERS
  qc_mapping_counter_init(&qc_mapping_counter);

  // multi-threads
  //  
  bam_reader_t* bam_reader_p = bam_reader_by_batch_new(input_filename, batch_size, base_quality, &bam_data_batch_list_gpu, LIST_INSERT_MODE);
  bam_reader_start(bam_reader_p);
  //num_alignments = bam_reader_join(bam_reader_p);

  // some local variables
  //
  void* r;
  
  // multi-threads
  //
  pthread_t* qc_calc_server_thread_p = (pthread_t*) malloc(((num_gpu_devices == 0) ? 1 : num_gpu_devices) * sizeof(pthread_t));
  pthread_t* cpus_server_thread_p = (pthread_t*) malloc(cpu_num_threads * sizeof(pthread_t));
  pthread_t results_server_thread;

  // calling GPU threads to process the bam data,
  // but first, prepare input parameter
  //
  int i;
  qc_calc_server_input_t** qc_calc_server_input_p = (qc_calc_server_input_t**) calloc(num_gpu_devices, sizeof(qc_calc_server_input_t*));
  
  if (num_gpu_devices > 0) {
    for (i=0; i < num_gpu_devices; i++) {
      qc_calc_server_input_p[i] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
    }
  } else {
    qc_calc_server_input_p[0] = (qc_calc_server_input_t*) calloc(1, sizeof(qc_calc_server_input_t));
  }

  if (num_gpu_devices > 0) {	// GPU implementacion
    for (i=0; i < num_gpu_devices; i++) {
	qc_calc_server_input_p[i]->num_gpu_devices = num_gpu_devices;
	qc_calc_server_input_p[i]->cpu_num_threads = 0;
	qc_calc_server_input_p[i]->gpu_device_id[0] = i;
	qc_calc_server_input_p[i]->gpu_num_blocks = gpu_num_blocks;
	qc_calc_server_input_p[i]->gpu_num_threads = gpu_num_threads;
	qc_calc_server_input_p[i]->gpu_batch_list_p = &bam_data_batch_list_gpu;
	qc_calc_server_input_p[i]->cpu_batch_list_p = &bam_data_batch_list_cpu;
	//bam_data_batch_list_incr_producers(&bam_data_batch_list_cpu);
	pthread_create(&qc_calc_server_thread_p[i], NULL, qc_calc_server, (void*) qc_calc_server_input_p[i]);
    }
  } else {			// CPU implementacion
    qc_calc_server_input_p[0]->num_gpu_devices = 0;
    qc_calc_server_input_p[0]->cpu_num_threads = cpu_num_threads;
    qc_calc_server_input_p[0]->gpu_device_id[0] = 0;
    qc_calc_server_input_p[0]->gpu_num_blocks = 0;
    qc_calc_server_input_p[0]->gpu_num_threads = 0;
    qc_calc_server_input_p[0]->gpu_batch_list_p = &bam_data_batch_list_gpu;
    qc_calc_server_input_p[0]->cpu_batch_list_p = &bam_data_batch_list_cpu;
    //bam_data_batch_list_incr_producers(&bam_data_batch_list_cpu);
    pthread_create(&qc_calc_server_thread_p[0], NULL, qc_calc_server, (void*) &qc_calc_server_input_p[0]);
  }

  /*for (int i=0; i < num_gpu_devices; i++) {
   pthread_join(qc_calc_server_thread_p[i], &r);
  }*/  
  
  // calling CPU threads to process the bam data,
  //
  //for (int i=0; i < cpu_num_threads; i++) {
  cpus_server_input_t* cpus_server_input_p = (cpus_server_input_t*) calloc(1, sizeof(cpus_server_input_t));
  
  for (i=0; i < 1; i++) {
    cpus_server_input_p->cpu_num_threads = cpu_num_threads;
    cpus_server_input_p->max_distance_size = max_distance_size;
    cpus_server_input_p->cpu_batch_list_p = &bam_data_batch_list_cpu;
    cpus_server_input_p->qc_mapping_counter = &qc_mapping_counter;
    cpus_server_input_p->gff_filename = gff_filename;
    cpus_server_input_p->output_directory = output_directory;
    cpus_server_input_p->input_filename = input_filename;
      
    pthread_create(&cpus_server_thread_p[i], NULL, cpus_server, (void*) cpus_server_input_p);
  }
  
/*  
  // wait for all terminating
  //  
  for (int i=0; i < num_gpu_devices; i++) {
   pthread_join(qc_calc_server_thread_p[i], &r);
  }  
  
  for (int i=0; i < 1; i++) {
   pthread_join(cpus_server_thread_p[i], &r);
  } */

  // calling thread to process results from GPU,
  //
  results_server_input_t results_server_input;
  results_server_input.gpu_num_blocks = gpu_num_blocks;
  results_server_input.gpu_num_threads = gpu_num_threads;
  results_server_input.base_quality = base_quality;
  results_server_input.qc_mapping_counter = &qc_mapping_counter;
  results_server_input.filename = input_filename;
  results_server_input.report_directory = output_directory;
  pthread_create(&results_server_thread, NULL, results_server, (void*) &results_server_input);

  num_alignments = bam_reader_join(bam_reader_p);
  
  for (int i=0; i < num_gpu_devices; i++) {
   pthread_join(qc_calc_server_thread_p[i], &r);
   //pthread_detach(qc_calc_server_thread_p[i]);
  }
  free(qc_calc_server_thread_p);
  
  for (int i=0; i < 1; i++) {
   pthread_join(cpus_server_thread_p[i], &r);
   //pthread_detach(cpus_server_thread_p[i]);
  }
  free(cpus_server_thread_p);
  
  pthread_join(results_server_thread, &r);

  // free thread stuff and parameters
  
  if (num_gpu_devices > 0) {
    for (i=0; i < num_gpu_devices; i++) {
      free(qc_calc_server_input_p[i]);
    }
  } else {
    free(qc_calc_server_input_p[0]);
  }
  
  free(qc_calc_server_input_p);  
  free(cpus_server_input_p);
  bam_reader_free(bam_reader_p);

}

#endif /* QC_CU_ */