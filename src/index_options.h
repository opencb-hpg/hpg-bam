#ifndef INDEX_OPTIONS_H
#define INDEX_OPTIONS_H

/*
 * index_options.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "argtable2.h"
#include "libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//============================ DEFAULT VALUES ============================

//------------------------------------------------------------------------

#define NUM_INDEX_OPTIONS	3

//------------------------------------------------------------------------

typedef struct index_options { 
  //  int log_level;
  //  int verbose;
  int help;
  //  int num_threads;
  //  int batch_size;

  char* in_filename;
  char* out_dirname;
  //  char* gff_region_filename;
  //  char* region_list;

  char *exec_name;
  char *command_name;
} index_options_t;

//------------------------------------------------------------------------

index_options_t *new_index_options(char *exec_name, char *command_nane);

index_options_t *parse_index_options(char *exec_name, char *command_nane,
				     int argc, char **argv);

void free_index_options(index_options_t *opts);

void validate_index_options(index_options_t *opts);

void display_index_options(index_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif