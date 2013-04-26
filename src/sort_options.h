#ifndef SORT_OPTIONS_H
#define SORT_OPTIONS_H

/*
 * sort_options.h
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

#define NUM_SORT_OPTIONS	5

//------------------------------------------------------------------------

typedef struct sort_options { 
  int help;

  size_t max_memory;
  char *criteria;

  char* in_filename;
  char* out_dirname;

  char *exec_name;
  char *command_name;
} sort_options_t;

//------------------------------------------------------------------------

sort_options_t *new_sort_options(char *exec_name, char *command_nane);

sort_options_t *parse_sort_options(char *exec_name, char *command_nane,
				     int argc, char **argv);

void free_sort_options(sort_options_t *opts);

void validate_sort_options(sort_options_t *opts);

void display_sort_options(sort_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
