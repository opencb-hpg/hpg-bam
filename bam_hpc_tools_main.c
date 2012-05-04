#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "aligner_dataset.h"
#include "aligner_dataset_file.h"
#include "bam.h"
#include "commons.h"
#include "convert.h"
#include "file_utils.h"
#include "log.h"
#include "qc.h"
#include "sam.h"
#include "sort.h"
#include "string_utils.h"
#include "system_utils.h"

#define DEFAULT_MAX_DISTANCE_SIZE		500

#define DEFAULT_GPU_NUM_BLOCKS 			16
#define DEFAULT_GPU_NUM_THREADS 		512
#define DEFAULT_GPU_NUM_DEVICES 		0
#define DEFAULT_CPU_NUM_THREADS 		2
#define DEFAULT_BATCH_SIZE_MB 			32
#define DEFAULT_BATCH_LIST_SIZE 		4

#define DEFAULT_MAX_NUM_HITS			10
#define DEFAULT_MIN_QUALITY			0
#define DEFAULT_MAX_QUALITY			255
#define DEFAULT_MIN_MISMATCHES	 		0
#define DEFAULT_MAX_MISMATCHES	 		4

#define DEFAULT_BASE_QUALITY			PHRED33

#define BAM_GPU_TOOLS_USAGE_HELP "USAGE: bam-gpu-tools [--qc] [--sort] [--filter] --bam <bam_file> --outdir </path/to/dir> [--sam <sam_file> --to-sam | --to-bam] \
[--max-distance-size <max_distance_in_nts>] [--phred-quality <33|64|sanger|solexa>] [--conf <config_filename>] [--gff <gff_filename>] \
[--gpu-num-blocks <num_blocks>] [--gpu-num-threads <num_threads>] [--gpu-num-devices <num_devices>] [--cpu-num-threads <num_threads>] \
[--batch-size] [--batch-list-size] [--log-level <1-5>] [--log-file <log_filename>] [--verbose] [--t | --time] \n"


/* **********************************************
 *    		Global variables  		*
 * **********************************************/

int time_flag = 0;

double read_time = 0.0;
struct timeval t1_read, t2_read;

double gpu_time = 0.0;
struct timeval t1_gpu, t2_gpu;

double cpu_time = 0.0;
struct timeval t1_cpu, t2_cpu;

double result_time = 0.0;
struct timeval t1_result, t2_result;

double write_time = 0.0;
struct timeval t1_write, t2_write;

double sort_time = 0.0;
struct timeval t1_sort, t2_sort;

double convert_time = 0.0;
struct timeval t1_convert, t2_convert;

double reporting_time = 0.0;
struct timeval t1_reporting, t2_reporting;

double reader_server_time = 0.0;
struct timeval t1_reader_server, t2_reader_server;

double qc_calc_server_time = 0.0;
struct timeval t1_qc_calc_server, t2_qc_calc_server;

double cpus_server_time = 0.0;
struct timeval t1_cpus_server, t2_cpus_server;

double results_server_time = 0.0;
struct timeval t1_results_server, t2_results_server;

double gpus_standby_time = 0.0;
double cpus_standby_time = 0.0;
double results_standby_time = 0.0;
struct timeval t1_active_reader, t1_active_gpus, t1_active_cpus, t1_active_results;

double total_time = 0.0;
struct timeval t1_total, t2_total;

int num_alignments = 0;
int number_of_batchs = 0;
double mean_batch_size = 0;

unsigned int nts_with_coverage = 0;
unsigned long mean_coverage = 0;

int main(int argc, char **argv) {
  // setting global variables for logger
  log_level = LOG_DEFAULT_LEVEL;
  verbose = 1;
  log_filename = NULL;
  
  // variables to store steps to perform
  int qc_step = 0;
  int sort_step = 0; 
  int filter_step = 0; 
  int convert_step = 0;
  int validate_step = 0;
  int sort_dataset_step = 0;

  // variables to store qc options
  int max_distance_size = DEFAULT_MAX_DISTANCE_SIZE;
  int base_quality = DEFAULT_BASE_QUALITY;
  
  // variables to store common parameters
  char* phred_quality = (char*)"";
  
  // variables to store hpc parameters
  int gpu_num_blocks = 	DEFAULT_GPU_NUM_BLOCKS;			// 16
  int gpu_num_threads = DEFAULT_GPU_NUM_THREADS;		// 512
  int gpu_num_devices = DEFAULT_GPU_NUM_DEVICES;		// -1
  int cpu_num_threads = DEFAULT_CPU_NUM_THREADS;		// 2
  size_t batch_size = 	DEFAULT_BATCH_SIZE_MB * 1000000;	// 64MB
  int batch_list_size = DEFAULT_BATCH_LIST_SIZE;		// 4  
  
  // variables to store io parameters
  int to_sam_flag = 0;  
  int to_bam_flag = 0;
  char* sam_input = NULL;
  char* bam_input = NULL;
  char* output_directory = NULL;  
  char* gff_input = NULL;
  char* dataset_input = NULL;  
  
  // variables to store sort parameters
  int sort_by_id_flag;
  
  // variables to store filter parameters
  int max_num_hits = DEFAULT_MAX_NUM_HITS;
  int min_quality = DEFAULT_MIN_QUALITY;
  int max_quality = DEFAULT_MAX_QUALITY;
  int min_mismatches = DEFAULT_MIN_MISMATCHES;
  int max_mismatches = DEFAULT_MAX_MISMATCHES;
  
  //variables for long_options  
  int c;
  int option_index = 0;

  // struct defining options and its associated option internal value
  static struct option long_options[] = {
    /*	QC parameters	*/
    {"qc",  	      		no_argument, 0, 'a'},
    {"quality-control",		no_argument, 0, 'a'},
    {"max-distance-size",   	required_argument, 0, 'b'},
    
    /*	Sort  parameters	*/    
    {"sort",  	      		no_argument, 0, 'c'},
    {"by-id",  	      		no_argument, 0, 'd'},
    
    /*	Commons parameters	*/    
    {"phred-quality",   	required_argument, 0, 'e'},
    
    /*	HPC parameters		*/
    {"gpu-num-blocks",		required_argument, 0, 'f'},	// num-blocks: 16 (dafault 16)
    {"gpu-num-threads",		required_argument, 0, 'g'},	// num-threads: 512
    {"gpu-num-devices",		required_argument, 0, 'h'},	// num-devices: 1 (NOTE: 0 for all of them)
    {"cpu-num-threads",  	required_argument, 0, 'i'},    
    {"batch-size",     		required_argument, 0, 'j'},	 
    {"batch_list_size", 	required_argument, 0, 'k'},
    
    /*	IO parameters	*/
    {"sam",  	      		required_argument, 0, 'l'},
    {"bam",  	      		required_argument, 0, 'm'},
    {"to-sam",  		no_argument, 0, 'n'},
    {"to-bam",			no_argument, 0, 'o'},
    {"o",  			required_argument, 0, 'p'},
    {"outdir",  		required_argument, 0, 'p'},
    {"conf",  			required_argument, 0, 'q'},    
    {"gff",  			required_argument, 0, 'r'},
    
    /*	LOG parameters		*/    
    {"log-level",  		required_argument, 0, 's'},
    {"log-file",		required_argument, 0, 't'},
    {"v",			required_argument, 0, 'u'},	// verbose: if False no 'console' output    
    {"verbose",  		required_argument, 0, 'u'},    
    {"t",  	      		no_argument, 0, 'v'},
    {"time",  	      		no_argument, 0, 'v'},
   
    /*	Filter  parameters	(PENDING OF STUDY AND IMPLEMENTATION)	*/
    {"filter", 	      		no_argument, 0, 'w'},
    {"max-num-hits",  		required_argument, 0, 'x'},
    {"min-quality",  		required_argument, 0, 'y'},
    {"max-quality",  		required_argument, 0, 'z'},
    {"min-mismatches",  	required_argument, 0, 'A'},
    {"max-mismatches",  	required_argument, 0, 'B'}, 
    
    /*	Validate BAM		*/
    {"validate", 	      	no_argument, 0, 'C'},
    {"dataset",  		required_argument, 0, 'D'},

    /*	Sort dataset		*/
    {"sort-dataset", 	      	no_argument, 0, 'E'},
    
    {0, 0, 0, 0}
  }; 
  
  set_log_level(1);

  int argc_with_file_options = 0;
  char** argv_with_file_options = NULL;
  char** argv_from_file_options = NULL;
	
  for(int i=0; i<argc; i++) {
    if(strcmp(argv[i], "--conf") == 0) {
      
      char str[256];
      strcpy(str, "Reading config file: ");
      strcat(str, argv[i+1]);
      LOG_DEBUG(str);
      
      argv_from_file_options = parse_conf_file(argv[i+1]);
      int num_conf_lines = count_lines(argv[i+1]);
      
      argv_with_file_options = (char **)malloc((argc + 2*num_conf_lines) * sizeof(char *));
      array_concat(argv_with_file_options, argc, (const char**)argv, 2*num_conf_lines, (const char**)argv_from_file_options);
      
      char command_line[1024];
      strcpy(command_line,"Command line: ");
      argc_with_file_options = argc + 2 * num_conf_lines;

      for(int i=0; i<argc_with_file_options; i++) {
	strcat(command_line, argv_with_file_options[i]);
	strcat(command_line, " ");
      }
      
      LOG_INFO(command_line);
    }
  }
  
  argc_with_file_options = argc;
  argv_with_file_options = argv;

  // validation of no argument launch
  if (argc < 2) {
    printf(BAM_GPU_TOOLS_USAGE_HELP); 
    exit(0);
  }	

  while((c = getopt_long(argc_with_file_options, argv_with_file_options, "", long_options, &option_index)) != -1) {	
    switch (c) {
      	
	/*	PARSING QC PARAMETERS	*/
	case 'a':
	  //printf("option --qc selected, quality control selected\n");
	  qc_step = true;
	  break;
	  
	case 'b':
	  //printf("option --max-distance-size with value '%s'\n", optarg);
	  if (is_numeric(optarg) == 1) {	
	      sscanf(optarg, "%i", &max_distance_size);
	  } else {
	      LOG_WARN("--max-distance-size is not a valid number, assuming default value 500\n");
	  }
	  break;  
	
	/*	PARSING SORT PARAMETERS	*/
	case 'c':
	  //printf("option --sort selected, sorting enabled\n");
	  sort_step = 1;
	  break;
	  
	case 'd':
	  //printf("option --sort selected, sorting enabled\n");
	  sort_by_id_flag = 1;
	  break; 
      
	/*	PARSING COMMON PARAMETERS	*/
	case 'e':
	  //printf("option --phred-quality with value '%s'\n", optarg);
	  if (base_quality == PHRED33) {
	    if (strcmp(optarg, "33") == 0) {
	      base_quality = PHRED33;
	    } else if (strcmp(optarg, "64") == 0) {
	      base_quality = PHRED64;
	    } else if (strcmp(optarg, "sanger") == 0) {
	      base_quality = PHRED33;
	    } else if (strcmp(optarg, "solexa") == 0) {
	      base_quality = PHRED64;
	    } else {
	      LOG_WARN("Incorrect quality scale (33 or 64). Assuming 33.\n");
	    }
	  }  
	  break;	  
	
	/*	PARSING HPC PARAMETERS	*/  
	case 'f':
	  if (gpu_num_blocks == DEFAULT_GPU_NUM_BLOCKS) {
	    if(is_numeric(optarg) != 0) {
	      sscanf(optarg, "%i", &gpu_num_blocks);
	      
	      if(gpu_num_blocks < 8) {
		gpu_num_blocks = 8;
		LOG_WARN("--gpu-num-blocks is not a valid number, assuming 8\n");
	      }
	      
	    } else {
	      LOG_FATAL("--gpu-num-blocks is not a valid number, aborting execution\n");
	    }
	  }
	  break;

	case 'g':
	  if (gpu_num_threads == DEFAULT_GPU_NUM_THREADS) {
	    if(is_numeric(optarg) != 0) {
	      sscanf(optarg, "%i", &gpu_num_threads);
	      
	      if(gpu_num_threads < 32) {
		gpu_num_threads = 32;
		LOG_WARN("--gpu-num-threads is not a valid number, assuming 32\n");
	      }
	      
	    } else {
	      LOG_FATAL("--gpu-num-threads is not a valid number, aborting execution\n");
	    }
	  }
	  break;

	case 'h':
	  //printf("option --grid-block-size with value '%s'\n", optarg);
	  if (gpu_num_devices == DEFAULT_GPU_NUM_DEVICES) {
	    if(is_numeric(optarg) != 0) {
	      sscanf(optarg, "%i", &gpu_num_devices);
	    } else {
	      LOG_FATAL("--gpu-num-devices is not a valid number, aborting execution\n");
	    }
	  }
	  break;

	case 'i':
	  //printf("option --threads with value '%s'\n", optarg);
	  if (cpu_num_threads == 1) {
	    sscanf(optarg, "%i", &cpu_num_threads);
	  }
	  break;
	  
	case 'j':
	  //printf("option --batch-size with value '%li'\n", batch_size);
	  if (batch_size == DEFAULT_BATCH_SIZE_MB * 1000000) {
	    if(is_numeric(optarg) != 0) {
	      sscanf(optarg, "%i", &batch_size);
	      
	      // batch-size > 16MB 
	      if (batch_size < 16) {
		batch_size = 16000000;
		LOG_WARN("the value --batch-size parameter must be at least 16000000 (16 MB) \n");
	      }
	    } else {
	      batch_size = 64;
	      LOG_FATAL("--batch-size is not a valid number, aborting execution\n");
	    }
	  }      
	  break;
	  
	case 'k':
	  //printf("option --list-length with value '%s'\n", optarg);
	  if (batch_list_size == 10) {
	    sscanf(optarg, "%i", &batch_list_size);
	  }
	  break;
      
	/*	PARSING IO PARAMETERS	*/
	case 'l':
	  //printf("option --sam with value '%s'\n", optarg);
	  if (sam_input == NULL) {
	    sam_input = optarg;
	  }
	  break;
	  
	case 'm':
	  //printf("option --bam with value '%s'\n", optarg);
	  if (bam_input == NULL) {
	    bam_input = optarg;
	  }
	  break;	  
	  
	case 'n':
	  //printf("option --to-sam with value '%s'\n", optarg);
	  if (to_sam_flag == 0) {
	    to_sam_flag = 1;
	    convert_step = 1;
	  }
	  break;  

	case 'o':
	  //printf("option --to-bam with value '%s'\n", optarg);
	  if (to_bam_flag == 0) {
	    to_bam_flag = 1;
	    convert_step = 1;
	  }
	  break;      
           
	case 'p':
	  //printf("option --output-dir with value '%s'\n", optarg);
	  output_directory = optarg;
	  break;
		
	case 'q':
	  //printf("option --conf filled '%s'\n", optarg);
	  break;   
    
	case 'r':
	  //printf("option --gff with value '%s'\n", optarg);
	  if (gff_input == NULL) {
	    gff_input = optarg;
	  }
	  break;
	  
	/*	PARSING LOG PARAMETERS	*/
	case 's':
	  //printf("option --log-level with value '%s'\n", optarg);
	  if(is_numeric(optarg) != 0) {
	    sscanf(optarg, "%i", &log_level);
	    LOG_LEVEL(log_level);
	  } else {
	    LOG_WARN("--log-level is not a valid number, assuming default level ERROR\n");
	  }
	  break;

	case 't':
	  //printf("option --fastq or --fq with value '%s'\n", optarg);
	  if(log_filename == NULL) {
	    log_filename = (char*) calloc(strlen(optarg) + 1, sizeof(char));
	    strcpy(log_filename, optarg);
	  }
	  break;
	  
	case 'u':
	  //printf("option --log-level with value '%s'\n", optarg);
	  if(strcmp(optarg, "true") == 0) {
	    LOG_VERBOSE(1);
	  } else if(strcmp(optarg, "false") == 0) {
	    LOG_VERBOSE(0);
	  } else {
	    LOG_WARN("--verbose parameter must be true or false, assuming default value false");
	  }
	  break;

	case 'v':
	  //printf("option --time selected. timing enabled\n");
	  time_flag = 1;
	  break;
	  
	/*	FILTER PARAMETERS (NOT IMPLEMENTED)	*/    
	case 'w':
	  //printf("option --filter selected, filtering enabled\n");
	  filter_step = 1;
	  break;	      
	  
	case 'x':
	  //printf("option --max-num-hits with value '%s'\n", optarg);
	  if (max_num_hits == 0) {
	    if (is_numeric(optarg) == 1) {	
		sscanf(optarg, "%i", &max_num_hits);
	    } else {
		LOG_WARN("--max-num-hits is not a valid number, assuming default value 10\n");
	    }
	  }            
	  break;      
      
	case 'y':
	  //printf("option --min-quality with value '%s'\n", optarg);
	  if (min_quality == 0) {
	    if (is_numeric(optarg) == 1) {	
		sscanf(optarg, "%i", &min_quality);
	    } else {
		LOG_WARN("--min-quality is not a valid number, assuming default value 0\n"); 
	    }
	  }
	  break;            

	case 'z':
	  //printf("option --max-quality with value '%s'\n", optarg);
	  if (max_quality == 0) {
	    if (is_numeric(optarg) == 1) {	
		sscanf(optarg, "%i", &max_quality);
	    } else {
		LOG_WARN("--max-quality is not a valid number, assuming default value 255\n");  
	    }
	  }
	  break;                 
          
	case 'A':
	  //printf("option --min-mismatches with value '%s'\n", optarg);
	  if (min_mismatches == 0) {
	    if (is_numeric(optarg) == 1) {	
		sscanf(optarg, "%i", &min_mismatches);
	    } else {
		LOG_INFO("--min-mismatches is not a valid number, assuming default value 0\n");
	    }
	  }
	  break;            

	case 'B':
	  //printf("option --max-mismatches with value '%s'\n", optarg);
	  if (max_mismatches == 0) {
	    if (is_numeric(optarg) == 1) {	
		sscanf(optarg, "%i", &max_mismatches);
	    } else {
		LOG_WARN("--max-mismatches is not a valid number, assuming default value 4\n");
	    }
	  }
	  break; 
	  
	/*	VALIDATE BAM PARAMETERS 	*/    
	case 'C':
	  //printf("option --validate selected, performing BAM validation\n");
	  filter_step = 1;
	  break;	      
	  
	case 'D':
	  //printf("option --dataset with value '%s'\n", optarg);
	  if (dataset_input == NULL) {
	    dataset_input = optarg;
	  }         
	  break;
	  
	case 'E':
	  //printf("option --validate selected, performing BAM validation\n");
	  sort_dataset_step = 1;
	  break;	        
      
	case ':':       /* option without mandatory operand */
	  fprintf(stderr, "Option -%c requires an operand\n", optopt);
	  break;
	  
	case '?':
	  printf(BAM_GPU_TOOLS_USAGE_HELP); 
	  break;
	  
	default:
	  printf(BAM_GPU_TOOLS_USAGE_HELP); 
	  LOG_FATAL("Default case of parameters parsing. Aborting program");
	}
  }
 
  // free argv_with_file_options
  if (argc_with_file_options == 0) {
    for(int i=0; i<argc_with_file_options; i++) {
      free(argv_with_file_options[i]);
    }
  }

  // QC step, BAM <-> SAM conversion, bam sorting, dataset sorting, validation and filtering are exclusive options  
  if (convert_step) { 
    qc_step = 0;     
    sort_step = 0;     
    filter_step = 0;
    validate_step = 0;
    sort_dataset_step = 0;
  }
  
  if (qc_step) { 
    convert_step = 0;     
    sort_step = 0;     
    filter_step = 0;
    validate_step = 0;
    sort_dataset_step = 0;    
  }
  
  if (sort_step) { 
    convert_step = 0;     
    qc_step = 0;     
    filter_step = 0;
    validate_step = 0;
    sort_dataset_step = 0;    
  }  
  
  if (sort_dataset_step) { 
    convert_step = 0; 
    qc_step = 0;     
    sort_step = 0;     
    filter_step = 0;
    validate_step = 0;
  }
  
  if (validate_step) { 
    convert_step = 0;
    qc_step = 0;     
    sort_step = 0;     
    filter_step = 0;
    sort_dataset_step = 0;
  }  

  // if no action is specified only quality control is performed  
  if ((qc_step == 0) && (sort_step == 0)  && (filter_step == 0) && (convert_step == 0) && (validate_step == 0) && (sort_dataset_step == 0)) {
      qc_step = 1;
  }

  // validating that minimal input parameters are filled
  // output directory and a bam file are mandatory  
  if ((!convert_step) && (output_directory == NULL)) {
    printf("--outdir option is mandatory\n");    
    printf(BAM_GPU_TOOLS_USAGE_HELP); 
    exit(0);
  }  

  if ((bam_input == NULL) && (!sort_dataset_step)) {
    printf("--bam option is mandatory\n");    
    printf(BAM_GPU_TOOLS_USAGE_HELP); 
    exit(0);
  }
  
  // only one direction can be converted BAM -> SAM or SAM -> BAM, not both  
  if ((to_sam_flag) && (to_bam_flag)) {
    to_sam_flag = 0;
    printf(BAM_GPU_TOOLS_USAGE_HELP);
    LOG_WARN("both conversion directions are selected, assuming --to-bam option\n");
  }
  
  // if BAM <-> SAM conversion is activated both SAM and BAM files must be informed  
  if (convert_step) {
    if ((sam_input == NULL) || (bam_input == NULL)) {
      printf(BAM_GPU_TOOLS_USAGE_HELP); 
      LOG_FATAL("conversion option is activated, --sam and --bam options are both mandatory\n");
    }
  }

  // validating BAM file according to initial dataset
  if (validate_step) {
    if (dataset_input == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP); 
      LOG_FATAL("validate BAM option is activated, --dataset option is mandatory\n");
    } else if (bam_input == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP);
      LOG_FATAL("validate BAM option is activated, --bam option is mandatory\n");      
    } else if (output_directory == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP);
      LOG_FATAL("validate BAM option is activated, --outdir option is mandatory\n");      
    }
  }
  
  // both dataset file and output directory are mandatory for sorting dataset and validate BAM
  if ((sort_dataset_step) || (validate_step)) {
    if (dataset_input == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP); 
      LOG_FATAL("sort dataset option is activated, --dataset option is mandatory\n");
    } else if (output_directory == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP);
      LOG_FATAL("validate BAM option is activated, --outdir option is mandatory\n");      
    }
  }
  
  //bam input file is mandatory for its validation
  if (validate_step) {
    if (bam_input == NULL) {
      printf(BAM_GPU_TOOLS_USAGE_HELP);
      LOG_FATAL("validate BAM option is activated, --bam option is mandatory\n");      
    }   
  }  
   
  // listh length must be at least 4  
  if (batch_list_size < 4) {
    batch_list_size = 4;
    LOG_WARN("--batch_list_size must be at least 4, assuming default value 4\n");
  }   
 
  // print any remaining command line arguments that are not options
  if (optind < argc) {
    printf("no valid options: ");
    while (optind < argc) printf("%s ", argv[optind++]);
    printf("\n");
    printf(BAM_GPU_TOOLS_USAGE_HELP);
  }
  
  // aplying heuristic values if default values have not been modified  
  if (batch_size == 1000000 * DEFAULT_BATCH_SIZE_MB ) {
    batch_size = get_optimal_batch_size(BAM_QC, 0);
  }
   
  if (cpu_num_threads == DEFAULT_CPU_NUM_THREADS) {
    cpu_num_threads = get_optimal_cpu_num_threads();  
  }

  if (gpu_num_threads == DEFAULT_GPU_NUM_THREADS) {
    gpu_num_threads = get_optimal_gpu_num_threads();
  }
 
  // start measuring time after options
  if (time_flag) { start_timer(t1_total); }
  
  // control of the execution flow depending on the filled options
  if (to_sam_flag == 1) {
    convert_bam_to_sam(bam_input, sam_input);
  } else if (to_bam_flag == 1) {
    convert_sam_to_bam(sam_input, bam_input);
  }

  if ((sort_step) && (bam_input != NULL)) {
    #define THRUST-GPU
    if (sort_by_id_flag) {
        sort_bam_file_by_id(batch_size, bam_input, output_directory);
    } else {
        sort_bam_file(batch_size, bam_input, output_directory);
    }
  }

  if (qc_step) {
    qc_bam_file(batch_size, batch_list_size, gpu_num_threads, gpu_num_blocks, cpu_num_threads, base_quality, max_distance_size, bam_input, output_directory, gff_input);
  }

  if (filter_step) {
    //filter_bam_file();
  }
  
  if (validate_step) {
    //validate_bam_file(bam_input, dataset_input, output_directory);
  }
  
  if (sort_dataset_step) {
    sort_dataset_by_id(dataset_input, output_directory);
  }
  
  total_time = 0;
  if (time_flag) { stop_timer(t1_total, t2_total, total_time); }

  if (time_flag) {   
    printf("\n");  
    printf("number of alignments     : \t%10i\n\n", num_alignments);
    printf("number of batches        : \t%10i\n\n", number_of_batchs);
    printf("mean alignments per batch  : \t%10.2f\n", 1.0 * num_alignments / number_of_batchs); 
        
    printf("total time           (s): \t%10.5f\n", 0.000001 * total_time);  
    printf("\n");  
    printf("total read time      (s): \t%10.5f\n", 0.000001 * read_time);
    printf("total gpu time       (s): \t%10.5f\n", 0.000001 * gpu_time);
    printf("total cpu time       (s): \t%10.5f\n", 0.000001 * cpu_time);
    printf("total result time    (s): \t%10.5f\n", 0.000001 * result_time);
    printf("total write time     (s): \t%10.5f\n", 0.000001 * write_time);
    
    if (sort_step) {
      printf("total sort time      (s): \t%10.5f\n", 0.000001 * sort_time);
    }
 	
    if (qc_step) {
	printf("total reporting time (s): \t%10.5f\n", 0.000001 * reporting_time); 
    }
    
    if (convert_step) {
	printf("total convert time (s): \t%10.5f\n", 0.000001 * convert_time);       
    }
    
    printf("\n----------- elapsed server times -----------\n\n");
    printf("reader server time (s): \t%10.5f\n", 0.000001 * reader_server_time);
    printf("qc calc server time(s): \t%10.5f\n", 0.000001 * qc_calc_server_time);
    printf("cpu server time    (s): \t%10.5f\n", 0.000001 * cpus_server_time);
    printf("result server time (s): \t%10.5f\n", 0.000001 * results_server_time);     

    printf("\n------ standby times until first batch -----\n\n");
    printf("gpu standby time      (s): \t%10.5f\n", 0.000001 * gpus_standby_time);
    printf("cpu standby time      (s): \t%10.5f\n", 0.000001 * cpus_standby_time);
    printf("result standby time   (s): \t%10.5f\n", 0.000001 * results_standby_time);
  }
    
  //free memory
  //if (log_filename != NULL) free(log_filename);
  //if (argv_from_file_options != NULL) free(argv_from_file_options);  
  //if (sam_input != NULL) free(sam_input);
  //if (bam_input != NULL) free(bam_input);
  //if (output_directory != NULL) free(output_directory);
  //if (gff_input != NULL) free(gff_input);
  //if (dataset_input != NULL) free(dataset_input);
    
  return 1;
}