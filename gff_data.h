
#ifndef GFF_DATA_H
#define GFF_DATA_H

#include <stdio.h>

#include "bam_data_batch_list.h"
#include "commons.h"

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for storing a gff file line
*
* Structure containing data of a gff file line (see http://www.sanger.ac.uk/resources/software/gff/spec.html)
*/
typedef struct gff_line {
  char* seqname;   		/**< Name of the sequence (chromosome or scaffold). */
  char* source;   		/**< Program that generated the feature. */
  char* feature;   		/**< Type of feature ("CDS", "start_codon", "stop_codon", and "exon"). */
  int start;   			/**< Start position. */
  int end;   			/**< End position. */
  char* score;   		/**< Score between 0 and 1000. */
  char* strand;   		/**< "+", "-" or "." (unknown). */
  char* frame;   		/**< "0", "1", "2" or ".". */
  char* group;   		/**< Group of lines. */
} gff_line_t;

/**
* @brief Structure for storing a region defined in a gff file
*
* Structure containing data of a region defined in a gff file
*/
typedef struct gff_region {
  short int chromosome;   	/**< Chromosome. */
  int start;   			/**< Start position. */
  int end;   			/**< End position. */
} gff_region_t;

/**
* @brief Structure for storing gff file data
*
* Structure for storing gff file data
*/
typedef struct gff_data {
  int num_regions;   		/**< Number of regions in data. */
  int actual_region;   		/**< Actual region. */
  gff_line_t* gff_lines_p;   	/**< Pointer to gff lines. */
  gff_region_t* gff_regions_p;  /**< Pointer to gff regions. */
  pthread_mutex_t lock;   	/**< Lock for gff data. */
} gff_data_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a new gff_data
*  @param gff_filename gff filename
*  @return pointer to created gff_data
*  
*  Creates a new gff_data
*/
gff_data_t* gff_data_new(char* gff_filename);

/**
*  @brief Frees gff lines
*  @param[in,out] gff_lines_p pointer to gff_lines structure
*  @return void
*  
*  Frees gff lines
*/
void gff_lines_free(gff_line_t* gff_lines_p);

/**
*  @brief Frees gff data
*  @param[in,out] gff_data_p pointer to gff_data structure
*  @return void
*  
*  Frees gff data
*/
void gff_data_free(gff_data_t* gff_data_p);

/**
*  @brief Prints regions in gff data
*  @param gff_data_p pointer to gff_data
*  @return void
*  
*  Prints regions in gff data
*/
void gff_data_print_regions(gff_data_t* gff_data_p);

/**
*  @brief Prints lines in gff data
*  @param gff_data_p pointer to gff_data
*  @return void
*  
*  Prints lines in gff data
*/
void gff_data_print_lines(gff_data_t* gff_data_p);

/**
*  @brief Determines if a given alignment is inside a region in gff_data
*  @param gff_data_p pointer to gff_data
*  @param chromosome chromosome of the alignment
*  @param start_coordinate start coordinate of the alignment
*  @param end_coordinate end coordinate of the alignment
*  @param regions pointer to regions where alignments are located
*  @return 1 if true, 0 if false
*  
*  Determines if a give aligment is inside a region
*/
int gff_data_alignment_in_region(gff_data_t* gff_data_p, int chromosome, int start_coordinate, int end_coordinate, int* regions);

/**
*  @brief Determines if a batch is inside a region
*  @param bam_data_batch_p pointer to the bam data batch 
*  @param gff_data_p pointer to gff_data
*  @return 1 if true, 0 if false
*  
*  Determines if a batch is inside a region
*/
int gff_data_batch_in_region(bam_data_batch_t* bam_data_batch_p, gff_data_t* gff_data_p);

#endif  /* GFF_DATA_H */