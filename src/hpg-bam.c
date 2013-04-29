/*
 * main.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "containers/khash.h"

#include "bioformats/features/region/region_table.h"
#include "bioformats/features/region/region_table_utils.h"
#include "bioformats/bam/bam_db.h"
#include "bioformats/bam/bam_stats.h"
#include "bioformats/bam/bam_stats_report.h"
#include "bioformats/bam/bam_filter.h"

#include "bioformats/db/db_utils.h"

#include "stats_options.h"
/*
#include "filter_options.h"
#include "index_options.h"
#include "sort_options.h"
*/
//------------------------------------------------------------------------

void usage(char *exec_name);
region_table_t *build_region_table(char *by_string, char *by_gff_file);

extern int bam_index_build2(const char *fn, const char *_fnidx);

//------------------------------------------------------------------------
//                    M A I N     F U N C T I O N
//------------------------------------------------------------------------
//KHASH_MAP_INIT_STR(str, int)

int main (int argc, char *argv[]) {
  /*
  int ret;
  khiter_t k;

  char key[] = "I.1234", key2[] = "I.000";
  khash_t(str) *h;
  h = kh_init(str);


  k = kh_put(str, h, key, &ret);
  if (ret == 0) {
    kh_value(h, k) = (kh_value(h, k) + 1);
  } else if (ret == 1) {
    kh_value(h, k) = 3;
  }
  k = kh_get(str, h, key);
  printf("**** key %s -> value = %i\n", kh_key(h, k), kh_value(h, k));

  k = kh_put(str, h, key, &ret);
  if (ret == 0) {
    kh_value(h, k) = (kh_value(h, k) + 1);
  } else if (ret == 1) {
    kh_value(h, k) = 3;
  }
  k = kh_get(str, h, key);
  printf("**** key %s -> value = %i\n", kh_key(h, k), kh_value(h, k));

  k = kh_put(str, h, key, &ret);
  if (ret == 0) {
    kh_value(h, k) = (kh_value(h, k) + 1);
  } else if (ret == 1) {
    kh_value(h, k) = 3;
  }
  k = kh_get(str, h, key);
  printf("**** key %s -> value = %i\n", kh_key(h, k), kh_value(h, k));




  printf("ret = %i\n", ret);

  k = kh_put(str, h, key, &ret);
  printf("ret = %i\n", ret);

  //  if (!ret) kh_del(str, h, k);

  //  k = kh_put(str, h, key2, &ret);
  //  kh_value(h, k) = 23;

  for (k = kh_begin(h); k != kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      printf("key %s -> value = %i\n", kh_key(h, k), kh_value(h, k));
    }
  }


  k = kh_get(str, h, key2);
  printf("**** key %s -> value = %i\n", kh_key(h, k), kh_value(h, k));


  kh_destroy(str, h); 

  exit(-1);
  */

  // init logs, then after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  char *exec_name = argv[0];  
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    usage(exec_name);
  }

  char *command_name = argv[1];  

  argc--;
  argv++;

  if (strcmp(command_name, "stats") == 0) {

    //--------------------------------------------------------------------
    //                  S T A T S     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display stats options
    stats_options_t *opts = parse_stats_options(exec_name, command_name, 
						argc, argv);
    validate_stats_options(opts);
    display_stats_options(opts);

    // now, we can set logs according to the command-line
    init_log_custom(opts->log_level, 1, "hpg-bam.log", "w");

    // create db table and hash
    khash_t(stats_chunks) *hash = NULL;
    sqlite3 *db = NULL;
    //    sqlite3_stmt *stmt = NULL;
    if (opts->db) {
      hash = kh_init(stats_chunks);
      create_stats_db(opts->out_dbname, BAM_CHUNKSIZE, create_bam_query_fields, &db);
      //      prepare_statement_bam_query_fields(db, &stmt);
    }

    // create the region table
    region_table_t *region_table = build_region_table(opts->region_list,
						      opts->gff_region_filename);
    

    // set parameters
    bam_stats_input_t* input = bam_stats_input_new(opts->in_filename, region_table,
						   opts->num_threads, opts->batch_size,
						   (void *) db, (void *) hash);
    bam_stats_output_t* output = bam_stats_output_new();
    
    // run and display stats
    bam_stats(input, output);
    report_bam_stats_output(opts->in_filename, opts->out_dirname, (void *) db, output);

    // free memory
    bam_stats_input_free(input);
    bam_stats_output_free(output);
    if (region_table) free_table(region_table);
    stats_options_free(opts);

    if (db) {
      // after all records insertions, insert touched chunks from hash,
      // and after all inserts then create index
      insert_chunk_hash(BAM_CHUNKSIZE, hash, db);
      create_stats_index(create_bam_index, db);

      // finally, close db and free hash
      sqlite3_close(db);
      kh_destroy(stats_chunks, hash); 
    }

    /*
  } else if (strcmp(command_name, "filter" ) == 0) {

    //--------------------------------------------------------------------
    //                  F I L T E R     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display filter options
    filter_options_t *opts = parse_filter_options(exec_name, command_name, 
						  argc, argv);
    validate_filter_options(opts);
    display_filter_options(opts);
    
    // create the region table
    region_table_t *region_table = build_region_table(opts->region_list,
						      opts->gff_region_filename);
    
    // set parameters
    int by_num_errors, min_num_errors, max_num_errors;
    if (opts->min_num_errors != -1 ||
	opts->max_num_errors != -1) {
      by_num_errors = 1;
      min_num_errors = opts->min_num_errors; 
      max_num_errors = opts->max_num_errors; 
    } else {
      by_num_errors = 0;
    }
    int by_quality, min_quality, max_quality;
    if (opts->min_quality != -1 ||
	opts->max_quality != -1) {
      by_quality = 1;
      min_quality = opts->min_quality; 
      max_quality = opts->max_quality; 
    } else {
      by_quality = 0;
    }
    int by_length, min_length, max_length;
    if (opts->min_length != -1 ||
	opts->max_length != -1) {
      by_length = 1;
      min_length = opts->min_length; 
      max_length = opts->max_length; 
    } else {
      by_length = 0;
    }

    bam_filter_input_t* input = new_bam_filter_input(opts->in_filename, opts->out_dirname,
						     opts->mapped, opts->unmapped, 
						     opts->proper_pairs, opts->unique, 
						     by_num_errors, min_num_errors, max_num_errors,
						     by_quality, min_quality, max_quality,
						     by_length, min_length, max_length,	
						     region_table,
						     opts->num_threads, opts->batch_size);
    
    // run filter
    filter_bam(input);

    // free memory
    free_bam_filter_input(input);
    if (region_table) free_table(region_table);
    free_filter_options(opts);

  } else if (strcmp(command_name, "index" ) == 0) {

    //--------------------------------------------------------------------
    //                  I N D E X     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display index options
    index_options_t *opts = parse_index_options(exec_name, command_name, 
						argc, argv);
    validate_index_options(opts);
    display_index_options(opts);

    // set parameters
    char idx_filename[strlen(opts->out_dirname) + strlen(opts->in_filename) + 10];
    sprintf(idx_filename, "%s/%s.bai", opts->out_dirname, opts->in_filename);

    // run index
    bamFile bf = bam_open(opts->in_filename, "r");
    bam_index_t *idx = bam_index_core(bf);
    bam_close(bf);
    
    FILE *idxf = fopen(idx_filename, "wb");
    bam_index_save(idx, idxf);
    bam_index_destroy(idx);
    fclose(idxf);

    printf("Index created in %s\n", idx_filename);

    // free memory
    free_index_options(opts);

  } else if (strcmp(command_name, "sort" ) == 0) {

    //--------------------------------------------------------------------
    //                  S O R T     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display sort options
    sort_options_t *opts = parse_sort_options(exec_name, command_name, 
					      argc, argv);
    validate_sort_options(opts);
    display_sort_options(opts);

    // set parameters
    int is_by_qname = 0, is_stdout = 0;
    char *sorted_filename[strlen(opts->in_filename)];
    char *path[strlen(opts->out_dirname) + strlen(opts->in_filename) + 100];

    strcpy(sorted_filename, opts->in_filename);
    char *ext = strstr(sorted_filename, ".bam");
    if (ext) {
      *ext = 0;
    }
    sprintf(path, "%s/%s.sorted", opts->out_dirname, sorted_filename);

    if (strcmp("name", opts->criteria) == 0) {
      is_by_qname = 1;
    }

    // run sort
    bam_sort_core_ext(is_by_qname, opts->in_filename, path, opts->max_memory, is_stdout);

    printf("Sorted BAM file in %s.bam\n", path);

    // free memory
    free_sort_options(opts);
    */
  } else {

    //--------------------------------------------------------------------
    //                  U N K N O W N     C O M M A N D
    //--------------------------------------------------------------------

    usage(exec_name);
  }
  printf("Done !\n");
}

//------------------------------------------------------------------------

void usage(char *exec_name) {
    printf("Program: %s (High-performance tools for handling BAM files)\n", exec_name);
    printf("Version: 1.0.0\n");
    printf("\n");
    printf("Usage: %s <command> [options]\n", exec_name);
    printf("\n");
    printf("Command: stats\t\tstatistics summary\n");
    //    printf("         filter\t\tfilter a BAM file by using advanced criteria\n");
    //    printf("         index\t\tindex a BAM file (using the samtools 0.1.18)\n");
    //    printf("         sort\t\tsort a BAM file (using the samtools 0.1.18)\n");
    //    printf("         compare\tcompare two BAM files\n");
    //    printf("         realignment\trealign locally a BAM file\n");
    //    printf("         recalibrate\trecalibrate a BAM file\n");
    printf("\n");    
    printf("For more information about a certain command, type %s <command> --help\n", exec_name);
    exit(-1);
}

//------------------------------------------------------------------------

region_table_t *build_region_table(char *by_string, char *by_gff_file) {
  region_table_t *region_table = NULL;
  if (by_string) {
    region_table = parse_regions(by_string, 1, 
				 "http://ws.bioinfo.cipf.es/", "cel", "latest");
  } else if (by_gff_file) {
    region_table = parse_regions_from_gff_file(by_gff_file,
					       "http://ws.bioinfo.cipf.es/", "cel", "latest");
  }
  return region_table;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
