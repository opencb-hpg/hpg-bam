
#ifndef QC_H
#define QC_H

#include <stdlib.h>
#include <stdio.h>

#include <cprops/trie.h>

#include "cuda_commons.h"
#include "list.h"
#include "mappings_db.h"
#include "qc_hash.h"

#define VALID_ALIGNMENT_FILE_SUFFIX	".valid"
#define INVALID_ALIGNMENT_FILE_SUFFIX	".invalid"
#define NO_ALIGNMENT_FILE_SUFFIX

#define MAX_GFF_FILE_LINES		5000

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief QC calc server parameters
* 
* Structure containing parameters to pass to the qc calc server
*/
typedef struct qc_calc_server_input {
    int num_gpu_devices;		/**< Number of GPU devices. */
    int cpu_num_threads;		/**< Number of CPU threads. */
    int gpu_device_id[256];		/**< Ids of the GPU devices. */
    int gpu_num_blocks;			/**< Number of GPU blocks. */
    int gpu_num_threads;		/**< Number of GPU threads. */
    list_t* gpu_batch_list_p;		/**< Pointer to the gpu batch list. */
    list_t* cpu_batch_list_p;		/**< Pointer to the cpu batch list. */
} qc_calc_server_input_t;

/**
* @brief cpu server parameters
* 
* Structure containing parameters to pass to the cpu server
*/
typedef struct cpus_server_input {
    int cpu_num_threads;			/**< Number of CPU threads. */
    int max_distance_size;			/**< Max. distance to consider between pairend-end mappings. */
    list_t* cpu_batch_list_p;			/**< Pointer to the cpu batch list. */
    qc_mapping_counter_t* qc_mapping_counter;	/**< Pointer to qc mapping counter. */
    char* gff_filename;				/**< GFF file name. */
    char* output_directory;			/**< Directory to print output files. */
    char* input_filename;			/**< BAM file name. */
} cpus_server_input_t;

/**
* @brief results server parameters
* 
* Structure containing parameters to pass to the results server
*/
typedef struct results_server_input {
    int gpu_num_blocks;				/**< Number of GPU blocks. */
    int gpu_num_threads;			/**< Number of GPU threads. */
    int base_quality;				/**< Base quality for quality normalization. */
    qc_mapping_counter_t* qc_mapping_counter;	/**< Pointer to qc mapping counter. */
    char* filename;				/**< BAM file name. */
    char* report_directory;			/**< Directory to print report. */
} results_server_input_t;

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
