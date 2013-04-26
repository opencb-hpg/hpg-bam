#ifndef OPTIONS_H
#define OPTIONS_H

/*
 * hpg-bam.h
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

#define NUM_OPTIONS			8

//------------------------------------------------------------------------

typedef struct options {
  int log_level;
  int help;
  int num_threads;
  int batch_size;
  char* in_filename;
  char* out_filename;
  char* gff_region_filename;
  char* regions;
} options_t;

//------------------------------------------------------------------------

options_t *options_new(void);

void options_free(options_t *options);

void options_display(options_t *options);

void usage_cli();

void validate_options(options_t *options, char *mode);

/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 *
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
void** argtable_options_new(void);


/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 *
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void argtable_options_free(void **argtable_options);


/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 *
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
options_t *read_CLI_options(void **argtable_options, options_t *options);



/**
 * @brief Reads the configuration parameters of the application.
 * @param filename file the options data are read from
 * @param options_data options values (vcf filename, output directory...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 *
 * Reads the basic configuration parameters of the application. If the configuration file can't be
 * read, these parameters should be provided via the command-line interface.
 */
int read_config_file(const char *filename, options_t *options);



options_t *parse_options(int argc, char **argv);


void usage(void **argtable);

#endif
