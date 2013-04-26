#include "options.h"

//const char DEFAULT_OUTPUT_NAME[30] = "hpg-aligner_output";

//------------------------------------------------------------------------

options_t *options_new(void) {
  options_t *opts = (options_t*) calloc (1, sizeof(options_t));
  
  //======================= COMMON OPTIONS ====================
  opts->in_filename = NULL;
  opts->out_filename = NULL;
  opts->log_level = 0;
  opts->num_threads = 2;
  opts->batch_size = 10000;
  opts->gff_region_filename = NULL;
  opts->regions = NULL;

  /*
  opts->in_filename2 = NULL;
  opts->report_all =  0;
  opts->output_name = strdup(DEFAULT_OUTPUT_NAME);
  opts->num_gpu_threads = DEFAULT_GPU_THREADS;
  //GET Number System Cores
  //----------------------------------------------
  size_t num_cores = 0;
  if (num_cores = get_optimal_cpu_num_threads()) {
    opts->num_cpu_threads = num_cores;
  }else {
    opts->num_cpu_threads = DEFAULT_CPU_THREADS;
  }
  //----------------------------------------------
  opts->max_intron_length = DEFAULT_MAX_INTRON_LENGTH;
  opts->min_num_seeds = DEFAULT_MIN_NUM_SEEDS;
  opts->max_num_seeds = DEFAULT_MAX_NUM_SEEDS;
  opts->cal_seeker_errors = DEFAULT_CAL_SEEKER_ERRORS;
  opts->write_size = DEFAULT_WRITE_BATCH_SIZE;
  opts->bwt_threads = DEFAULT_BWT_THREADS;
  opts->region_threads = DEFAULT_REGION_THREADS;
  opts->num_cal_seekers = DEFAULT_NUM_CAL_SEEKERS;
  opts->num_sw_servers = DEFAULT_NUM_SW_THREADS;
  opts->min_score = DEFAULT_SW_MIN_SCORE;
  opts->match = DEFAULT_SW_MATCH;
  opts->mismatch = DEFAULT_SW_MISMATCH;
  opts->gap_open = DEFAULT_SW_GAP_OPEN;
  opts->gap_extend = DEFAULT_SW_GAP_EXTEND;
  opts->min_intron_length = DEFAULT_MIN_INTRON_LENGTH;
  opts->pair_mode = DEFAULT_PAIR_MODE;
  opts->pair_min_distance = DEFAULT_PAIR_MIN_DISTANCE;
  opts->pair_max_distance = DEFAULT_PAIR_MAX_DISTANCE;
  opts->timming = 0;
  opts->statistics = 0;
  opts->report_best = 0;
  opts->report_n_hits = 0;
  opts->gpu_process = 0;
  opts->bwt_set = 0;
  opts->reg_set = 0;
  opts->cal_set = 0;
  opts->sw_set = 0;

  opts->min_cal_size = 0; 
  opts->seeds_max_distance = 0;
  opts->batch_size = 0;
  opts->min_seed_size = 0;
  opts->seed_size = 0;
  opts->flank_length = 0;
  */
  return opts;
}

//------------------------------------------------------------------------

void validate_options(options_t *opts, char *mode) {

  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n\n");
    usage_cli();
  }

  /*
  int value_dir = exists(opts->output_name);
  int DEFAULT_READ_BATCH_SIZE;
  int DEFAULT_SEED_SIZE;
  int DEFAULT_FLANK_LENGTH;
  int DEFAULT_MIN_SEED_SIZE;
  int DEFAULT_MIN_CAL_SIZE;
  int DEFAULT_SEEDS_MAX_DISTANCE;

  if (strcmp("dna", mode) == 0) {
    DEFAULT_READ_BATCH_SIZE = 20000;
    DEFAULT_SEED_SIZE	= 20;
    DEFAULT_FLANK_LENGTH = 5;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = 30;
    DEFAULT_SEEDS_MAX_DISTANCE = 100;
  }else if (strcmp("rna", mode) == 0) {
    DEFAULT_READ_BATCH_SIZE = 200000;
    DEFAULT_SEED_SIZE = 15;
    DEFAULT_FLANK_LENGTH = 30;
    DEFAULT_MIN_SEED_SIZE = 15;
    DEFAULT_MIN_CAL_SIZE = 20;
    DEFAULT_SEEDS_MAX_DISTANCE = 60;
  }

  if (strcmp("dna", mode) == 0 || strcmp("rna", mode) == 0) {
    if (!value_dir) {
      create_directory(opts->output_name);
    }
    
    if (!opts->in_filename) {
      printf("Not filename input found. Please, insert it with option '-f FILENAME'.\n");
      usage_cli();
    }

    
    if (!opts->bwt_dirname) {
      printf("Not BWT index input found. Please, insert it with option '-i DIRNAME'.\n");
      usage_cli();
    }
  }

  if (!opts->min_cal_size) {
    opts->min_cal_size = DEFAULT_MIN_CAL_SIZE;
  }
 
  if (!opts->seeds_max_distance) {
    opts->seeds_max_distance = DEFAULT_SEEDS_MAX_DISTANCE;
  }
  
  if (!opts->batch_size) {
    opts->batch_size = DEFAULT_READ_BATCH_SIZE;
  }
   
  if (!opts->min_seed_size) {
    opts->min_seed_size = DEFAULT_MIN_SEED_SIZE;
  }
  
  if (!opts->seed_size) {
    opts->seed_size = DEFAULT_SEED_SIZE;
  }

  if (!opts->flank_length) {
    opts->flank_length = DEFAULT_FLANK_LENGTH;
  }
  */
}

//------------------------------------------------------------------------

void options_free(options_t *opts) {
     if(opts == NULL) { return; }

     if (opts->in_filename)	{ free(opts->in_filename); }
     if (opts->out_filename)	{ free(opts->out_filename); }
     if (opts->gff_region_filename)	{ free(opts->gff_region_filename); }
     if (opts->regions) { free(opts->regions); }

     /*
     if (opts->in_filename2  != NULL) { free(opts->in_filename2); }
     if (opts->bwt_dirname  != NULL)	{ free(opts->bwt_dirname); }     
     if (opts->genome_filename  != NULL) { free(opts->genome_filename); }
     if (opts->output_name  != NULL)	{ free(opts->output_name); }
     if (opts->extend_name != NULL) {free(opts->extend_name); }
     */
     free(opts);
}

//------------------------------------------------------------------------

void options_display(options_t *opts) {
     /*
     char* in_filename2 = NULL;
     if (opts->in_filename2 != NULL) {
	  in_filename2 = strdup(opts->in_filename2);
     }
     char* bwt_dirname =  strdup(opts->bwt_dirname);
     char* genome_filename =  NULL;
     if (opts->genome_filename != NULL) {
	  genome_filename =  strdup(opts->genome_filename);
     }
     unsigned int  report_all = (unsigned int)opts->report_all;
     unsigned int  report_best = (unsigned int)opts->report_best;
     unsigned int  report_n_hits = (unsigned int)opts->report_n_hits;
     
     if ((report_best == 0) && (report_n_hits == 0)) {
	  report_all = 1;
     }
     
     unsigned int num_gpu_threads =  (unsigned int)opts->num_gpu_threads;
     unsigned int num_cpu_threads =  (unsigned int)opts->num_cpu_threads;
     unsigned int cal_seeker_errors =  (unsigned int)opts->cal_seeker_errors; 
     unsigned int min_cal_size =  (unsigned int)opts->min_cal_size; 
     unsigned int seeds_max_distance =  (unsigned int)opts->seeds_max_distance; 
     unsigned int bwt_threads =  (unsigned int)opts->bwt_threads; 
     unsigned int batch_size =  (unsigned int)opts->batch_size; 
     unsigned int write_size =  (unsigned int)opts->write_size;  
     unsigned int num_cal_seekers =  (unsigned int)opts->num_cal_seekers;
     unsigned int region_threads =  (unsigned int)opts->region_threads;
     unsigned int num_sw_servers =  (unsigned int)opts->num_sw_servers;
     unsigned int min_seed_size =  (unsigned int)opts->min_seed_size;
     unsigned int seed_size =  (unsigned int)opts->seed_size;
     unsigned int min_num_seeds =  (unsigned int)opts->min_num_seeds;
     unsigned int max_num_seeds =  (unsigned int)opts->max_num_seeds;
     unsigned int max_intron_length =  (unsigned int)opts->max_intron_length;
     unsigned int flank_length =  (unsigned int)opts->flank_length;
     unsigned int pair_mode =  (unsigned int)opts->pair_mode;
     unsigned int pair_min_distance =  (unsigned int)opts->pair_min_distance;
     unsigned int pair_max_distance =  (unsigned int)opts->pair_max_distance;
     unsigned int min_intron_length =  (unsigned int)opts->min_intron_length;
     unsigned int gpu_process = (unsigned int)opts->gpu_process;
     float min_score =  (float)opts->min_score;
     float match =   (float)opts->match;
     float mismatch =   (float)opts->mismatch;
     float gap_open =   (float)opts->gap_open;
     float gap_extend =   (float)opts->gap_extend;
     */  
     printf("PARAMETERS CONFIGURATION\n");
     printf("=================================================\n");
     printf("Main options\n");
     printf("\tBAM input filename  : %s\n", opts->in_filename);
     if (opts->regions) {
       printf("\tRegions             : %s\n", opts->regions);
     } else if (opts->gff_region_filename) {
       printf("\tGFF region filename : %s\n", opts->gff_region_filename);
     }
     //     printf("\tOutput filename: %s\n", out_filename);
     printf("\n");
     printf("Report options\n");
     printf("\tLog level: %d\n",  (int) opts->log_level);
     printf("\n");
     printf("Architecture options\n");
     printf("\tNum. threads: %d\n",  (int) opts->num_threads);
     printf("\tBatch size  : %d alignments\n",  (int) opts->batch_size);
     /*
     printf("Num gpu threads %d\n", num_gpu_threads);
     printf("GPU Process: %s\n",  gpu_process == 0 ? "Disable":"Enable");
     printf("Num cpu threads %d\n",  num_cpu_threads);
     printf("Report all hits: %s\n",  report_all == 0 ? "Disable":"Enable");
     printf("Report best hits: %d\n",  report_best);
     printf("Report n hits: %d\n",  report_n_hits);
     printf("CAL seeker errors: %d\n",  cal_seeker_errors);
     printf("Batch size: %dBytes\n",  batch_size);
     printf("Write size: %dBytes\n",  write_size);
     printf("BWT Threads: %d\n",  bwt_threads);
     printf("Region Threads: %d\n",  region_threads);
     printf("Num CAL seekers: %d\n", num_cal_seekers);
     printf("Num SW servers: %d\n",  num_sw_servers);
     printf("SEEDING and CAL PARAMETERS\n");
     printf("\tMin. number of seeds: %d\n",  min_num_seeds);
     printf("\tMax. number of seeds: %d\n",  max_num_seeds);
     printf("\tSeed size: %d\n",  seed_size);
     printf("\tMin seed size: %d\n",  min_seed_size);
     printf("\tMin CAL size: %d\n",  min_cal_size);
     printf("\tSeeds max distance: %d\n",  seeds_max_distance);
     printf("\tFlank length: %d\n", flank_length);
     printf("RNA PARAMETERS\n");
     printf("\tMax intron length: %d\n", max_intron_length);
     printf("\tMin intron length: %d\n", min_intron_length);
     printf("PAIR-MODE PARAMETERS\n");
     printf("\tPair mode: %d\n", pair_mode);
     printf("\tMin. distance: %d\n", pair_min_distance);
     printf("\tMax. distance: %d\n", pair_max_distance);
     printf("SMITH-WATERMAN PARAMETERS\n");
     printf("\tMin score  : %0.4f\n", min_score);
     printf("\tMatch      : %0.4f\n", match);
     printf("\tMismatch   : %0.4f\n", mismatch);
     printf("\tGap open   : %0.4f\n", gap_open);
     printf("\tGap extend : %0.4f\n", gap_extend);
     */
     printf("=================================================\n");


     /*
     if (in_filename2 != NULL) free(in_filename2);
     free(bwt_dirname);
     free(genome_filename);
     free(output_name);
     */     
}

//--------------------------------------------------------------------

void** argtable_options_new(void) {
     void **argtable = (void**)malloc((NUM_OPTIONS + 1) * sizeof(void*));

     // NOTICE that order cannot be changed as is accessed by index in other functions
     argtable[0] = arg_file0("f", "in-file", NULL, "Input file name (BAM format)");
     argtable[1] = arg_file0("o", "out-file", NULL, "Output file name (BAM format)");
     argtable[2] = arg_int0("l", "log-level", NULL, "Log debug level");
     argtable[3] = arg_lit0("h", "help", "Help option");
     argtable[4] = arg_int0(NULL, "num-threads", NULL, "Number of threads");
     argtable[5] = arg_int0(NULL, "batch-size", NULL, "Batch size (in number of alignments)");
     argtable[6] = arg_file0("r", "gff-refion-file", NULL, "Region file name (GFF format)");
     argtable[7] = arg_file0(NULL, "regions", NULL, "Regions (e.g., 1:3000-3200,4:100-200,...)");
     /*
     argtable[1] = arg_file0("i", "bwt-index", NULL, "BWT directory name");
     argtable[3] = arg_lit0(NULL, "report-all", "Report all alignments");
     argtable[4] = arg_file0("o", "outdir", NULL, "Output directory");
     argtable[5] = arg_int0(NULL, "gpu-threads", NULL, "Number of GPU Threads");
     argtable[6] = arg_int0(NULL, "cpu-threads", NULL, "Number of CPU Threads");
     argtable[7] = arg_int0("r", "index-ratio", NULL, "BWT index compression ratio");
     argtable[8] = arg_int0(NULL, "cal-seeker-errors", NULL, "Number of errors in CAL Seeker");
     argtable[9] = arg_int0(NULL, "min-cal-size", NULL, "Minimum CAL size");
     argtable[10] = arg_int0(NULL, "max-distance-seeds", NULL, "Maximum distance between seeds");
     argtable[11] = arg_int0(NULL, "read-batch-size", NULL, "Batch Size");
     argtable[12] = arg_int0(NULL, "write-batch-size", NULL, "Write Size");
     argtable[13] = arg_int0(NULL, "num-cal-seekers", NULL, "Number of CAL Seekers");
     argtable[14] = arg_int0(NULL, "num-sw-servers", NULL, "Number of Smith-Waterman servers");
     argtable[15] = arg_int0(NULL, "num-bwt-threads", NULL, "Number of BWT threads");
     argtable[16] = arg_int0(NULL, "num-region-threads", NULL, "Number of region threads");
     argtable[17] = arg_int0(NULL, "seed-size", NULL, "Number of nucleotides in a seed");
     argtable[18] = arg_int0(NULL, "min-seed-size", NULL, "Minimum number of nucleotides in a seed");
     argtable[19] = arg_int0(NULL, "cal-flank-size", NULL, "Flank length for CALs");
     argtable[20] = arg_dbl0(NULL, "sw-match", NULL, "Match value for Smith-Waterman algorithm");
     argtable[21] = arg_dbl0(NULL, "sw-mismatch", NULL, "Mismatch value for Smith-Waterman algorithm");
     argtable[22] = arg_dbl0(NULL, "sw-gap-open", NULL, "Gap open penalty for Smith-Waterman algorithm");
     argtable[23] = arg_dbl0(NULL, "sw-gap-extend", NULL, "Gap extend penalty for Smith-Waterman algorithm");
     argtable[24] = arg_dbl0(NULL, "sw-min-score", NULL, "Minimum score for valid mappings");
     argtable[25] = arg_int0(NULL, "max-intron-size", NULL, "Maximum intron size");
     argtable[26] = arg_int0(NULL, "min-intron-size", NULL, "Minimum intron size");
     argtable[27] = arg_lit0("t", "time", "Timming mode active");
     argtable[28] = arg_lit0("s", "stats", "Statistics mode active");
     argtable[30] = arg_file0("e", "ext", NULL, "File extend name");
     argtable[31] = arg_file0("g", "ref-genome", NULL, "Reference genome");
     argtable[32] = arg_file0("j", "fq2,fastq2", NULL, "Reads file input #2 (for paired mode)");
     argtable[33] = arg_int0(NULL, "paired-mode", NULL, "Pair mode: 0 = single-end, 1 = paired-end, 2 = mate-pair [Default 0]");
     argtable[34] = arg_int0(NULL, "paired-min-distance", NULL, "Minimum distance between pairs");
     argtable[35] = arg_int0(NULL, "paired-max-distance", NULL, "Maximum distance between pairs");
     argtable[36] = arg_int0(NULL, "report-n-best", NULL, "Report the <n> best alignments");
     argtable[37] = arg_int0(NULL, "report-n-hits", NULL, "Report <n> hits");
     argtable[38] = arg_int0(NULL, "min-num-seeds", NULL, "Minimum number of seeds per read");
     argtable[39] = arg_int0(NULL, "max-num-seeds", NULL, "Maximum number of seeds per read");
     argtable[40] = arg_lit0(NULL, "gpu-enable", "Enable GPU Process");
     */
     argtable[NUM_OPTIONS] = arg_end(20);
     
     return argtable;
}

//------------------------------------------------------------------------

void argtable_options_free(void **argtable) {
     if(argtable != NULL) {
	  arg_freetable(argtable, NUM_OPTIONS + 1);	// struct end must also be freed
	  free(argtable);
     }
}

//------------------------------------------------------------------------

int read_config_file(const char *filename, options_t *opts) {
	if (filename == NULL || opts == NULL) {
		return -1;
	}

	config_t *config = (config_t*) calloc (1, sizeof(config_t));
	int ret_code = config_read_file(config, filename);
	if (ret_code == CONFIG_FALSE) {
		LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
		return -1;
	}

	const char *tmp_string;
	long tmp_int;

	config_destroy(config);
	free(config);

	return ret_code;
}

//------------------------------------------------------------------------


/**
 * @brief Initializes an options_t structure from argtable parsed CLI with default values. Notice that options are order dependent.
 * @return A new options_t structure initialized with default values.
 *
 * Initializes the only default options from options_t.
 */
options_t *read_CLI_options(void **argtable, options_t *opts) {	
  if (((struct arg_file*)argtable[0])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[0])->filename)); }
  if (((struct arg_file*)argtable[1])->count) { opts->out_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_int*)argtable[2])->count) { opts->log_level = *(((struct arg_int*)argtable[2])->ival); }
  if (((struct arg_int*)argtable[3])->count) { opts->help = ((struct arg_int*)argtable[3])->count; }
  if (((struct arg_int*)argtable[4])->count) { opts->num_threads = *(((struct arg_int*)argtable[4])->ival); }
  if (((struct arg_int*)argtable[5])->count) { opts->batch_size = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_file*)argtable[6])->count) { opts->gff_region_filename = strdup(*(((struct arg_file*)argtable[6])->filename)); }
  if (((struct arg_file*)argtable[7])->count) { opts->regions = strdup(*(((struct arg_file*)argtable[7])->filename)); }

  /*
  if (((struct arg_file*)argtable[1])->count) { opts->bwt_dirname = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[3])->count) { opts->report_all = (((struct arg_int *)argtable[3])->count); }
  if (((struct arg_file*)argtable[4])->count) { free(opts->output_name); opts->output_name = strdup(*(((struct arg_file*)argtable[4])->filename)); }  
  if (((struct arg_int*)argtable[5])->count) { opts->num_gpu_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { opts->num_cpu_threads = *(((struct arg_int*)argtable[6])->ival); }
  if (((struct arg_int*)argtable[7])->count) { opts->index_ratio = *(((struct arg_int*)argtable[7])->ival); }
  if (((struct arg_int*)argtable[8])->count) { opts->cal_seeker_errors = *(((struct arg_int*)argtable[8])->ival); }
  if (((struct arg_int*)argtable[9])->count) { opts->min_cal_size = *(((struct arg_int*)argtable[9])->ival); }
  if (((struct arg_int*)argtable[10])->count) { opts->seeds_max_distance = *(((struct arg_int*)argtable[10])->ival); }
  if (((struct arg_int*)argtable[12])->count) { opts->write_size = *(((struct arg_int*)argtable[12])->ival); }
  if (((struct arg_int*)argtable[13])->count) { opts->cal_set = 1; opts->num_cal_seekers = *(((struct arg_int*)argtable[13])->ival); }
  if (((struct arg_int*)argtable[14])->count) { opts->sw_set = 1; opts->num_sw_servers = *(((struct arg_int*)argtable[14])->ival); }
  if (((struct arg_int*)argtable[15])->count) { opts->bwt_set = 1; opts->bwt_threads = *(((struct arg_int*)argtable[15])->ival); }
  if (((struct arg_int*)argtable[16])->count) { opts->reg_set = 1; opts->region_threads = *(((struct arg_int*)argtable[16])->ival); }
  if (((struct arg_int*)argtable[17])->count) { opts->seed_size = *(((struct arg_int*)argtable[17])->ival); }
  if (((struct arg_int*)argtable[18])->count) { opts->min_seed_size = *(((struct arg_int*)argtable[18])->ival); }
  if (((struct arg_int*)argtable[19])->count) { opts->flank_length = *((struct arg_int*)argtable[19])->ival; }
  if (((struct arg_dbl*)argtable[20])->count) { opts->match = *((struct arg_dbl*)argtable[20])->dval; }
  if (((struct arg_dbl*)argtable[21])->count) { opts->mismatch = *(((struct arg_dbl*)argtable[21])->dval); }
  if (((struct arg_dbl*)argtable[22])->count) { opts->gap_open = *(((struct arg_dbl*)argtable[22])->dval); }
  if (((struct arg_dbl*)argtable[23])->count) { opts->gap_extend = *(((struct arg_dbl*)argtable[23])->dval); }
  if (((struct arg_dbl*)argtable[24])->count) { opts->min_score = *(((struct arg_dbl*)argtable[24])->dval); }
  if (((struct arg_int*)argtable[25])->count) { opts->max_intron_length = *(((struct arg_int*)argtable[25])->ival); }
  if (((struct arg_int*)argtable[26])->count) { opts->min_intron_length = *(((struct arg_int*)argtable[26])->ival); }
  if (((struct arg_int*)argtable[27])->count) { opts->timming = ((struct arg_int*)argtable[27])->count; }
  if (((struct arg_int*)argtable[28])->count) { opts->statistics = ((struct arg_int*)argtable[28])->count; }
  if (((struct arg_file*)argtable[30])->count) { opts->extend_name = strdup(*(((struct arg_file*)argtable[30])->filename)); }
  if (((struct arg_file*)argtable[31])->count) { opts->genome_filename = strdup(*(((struct arg_file*)argtable[31])->filename)); }
  if (((struct arg_file*)argtable[32])->count) { opts->in_filename2 = strdup(*(((struct arg_file*)argtable[32])->filename)); }
  if (((struct arg_int*)argtable[33])->count) { opts->pair_mode = *(((struct arg_int*)argtable[33])->ival); }
  if (((struct arg_int*)argtable[34])->count) { opts->pair_min_distance = *(((struct arg_int*)argtable[34])->ival); }
  if (((struct arg_int*)argtable[35])->count) { opts->pair_max_distance = *(((struct arg_int*)argtable[35])->ival); }
  if (((struct arg_int*)argtable[36])->count) { opts->report_best = *(((struct arg_int*)argtable[36])->ival); }
  if (((struct arg_int*)argtable[37])->count) { opts->report_n_hits = *(((struct arg_int*)argtable[37])->ival); }
  if (((struct arg_int*)argtable[38])->count) { opts->min_num_seeds = *(((struct arg_int*)argtable[38])->ival); }
  if (((struct arg_int*)argtable[39])->count) { opts->max_num_seeds = *(((struct arg_int*)argtable[39])->ival); }
  if (((struct arg_int*)argtable[40])->count) { 
    #ifdef HPG_GPU
       opts->gpu_process = (((struct arg_int *)argtable[40])->count); 
    #else
       opts->gpu_process = 0; 
    #endif
  }
  */
  return opts;
}

//------------------------------------------------------------------------

options_t *parse_options(int argc, char **argv) {
  void **argtable = argtable_options_new();
  
  options_t *opts = options_new();
  if (argc < 2) {
    usage(argtable);
    exit(-1);
  } else {
    
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*)argtable[3])->count) {
      usage(argtable);
	argtable_options_free(argtable);
	options_free(opts);
	exit(0);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_OPTIONS], "hpg-bam");
      usage(argtable);
      exit(-1);
    }else {
      opts = read_CLI_options(argtable, opts);
      if(opts->help) {
	usage(argtable);
	argtable_options_free(argtable);
	options_free(opts);
	exit(0);
      }
      // Check if 'help' option has been provided.
    }
    
  }

  argtable_options_free(argtable);

  return opts;
}

//------------------------------------------------------------------------

void usage(void **argtable) {
  printf("Usage:\n./hpg-bam {stats | compare | sort | filter | realignment}");
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
}

//------------------------------------------------------------------------

void usage_cli() {
  void **argtable = argtable_options_new();
  usage(argtable);
  exit(0);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
