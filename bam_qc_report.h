
#ifndef BAM_QC_REPORT_H
#define BAM_QC_REPORT_H

#include <stdio.h>

#include "bam_commons.h"
#include "bam_qc_batch.h"
#include "file_utils.h"

#define QC_SUFFIX   			".qc"
#define MAP_ERRORS_HISTOGRAM_SUFFIX     ".map.error.hist.dat"
#define NUM_MAPPINGS_HISTOGRAM_SUFFIX 	".num.map.hist.dat"
#define HTML_FILE_SUFFIX  		".qc.html"
#define VALID_FILE_SUFFIX  		".valid"
#define INVALID_FILE_SUFFIX  		".invalid"

#define POS_COLUMN         		1
#define MAP_ERRORS_COLUMN 		2
#define DELETION_COLUMN  		3
#define INSERTION_COLUMN 		4
#define MATCHING_COLUMN  		5

#define NUM_ALIGNMENTS_COLUMN 		2

#define GRAPH_MAX_COLUMNS 		10
#define MAX_TITLE_LENGTH 		25
#define MAX_MAP_ERRORS_IN_HISTOGRAM 	10

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for BAM QC report
*
* Structure containing information of the BAM QC report.
*/
typedef struct bam_qc_report {
    unsigned long num_alignments;   						/**< Number of alignments. */
    unsigned long strand_counter;					   	/**< Strand counter (forward), num_aligments - strand_counter (reverse). */
    unsigned long mean_map_quality;   						/**< Mean quality mapping. */
    unsigned long mean_alignment_length;   					/**< Mean length of the alignments. */
    unsigned long mean_paired_end_distance;   					/**< Mean distance between paired ends. */
    unsigned int map_error_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];   	/**< Histogram of errors (D, I or X). */
    unsigned int map_deletion_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];	/**< Histogram of deletions (D). */
    unsigned int map_insertion_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];	/**< Histogram of insertions (I). */
    unsigned int map_matching_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];   	/**< Histogram of matchings (=). */
    unsigned int* num_mappings_histogram;   					/**< Histogram of number of mappings per read. */
} bam_qc_report_t;

/**
* @brief Structure for QC graphs parameters
* 
* Structure containing parameters for graphics representation of QC graphs
*/
typedef struct qc_graph {  
    int x_autoscale;					/**< Autoscale flag for X axis. */
    int x_start;					/**< X axis start coordinate. */
    int x_end;						/**< X axis end coordinate. */
    int y_autoscale;					/**< Autoscale flag for Y axis. */
    int y_start;					/**< Y axis start coordinate. */
    int y_end;						/**< Y axis end coordinate. */
    int lmargin;					/**< Left margin. */
    int rmargin;					/**< Right margin. */	
    int tmargin;					/**< Top margin. */
    int bmargin;					/**< Bottom margin. */
    int num_y_columns;					/**< Number of columns in the Y axis. */
    int x_column;					/**< Number of column in data file with X axis data. */
    int y_columns[GRAPH_MAX_COLUMNS];			/**< Numbers of columns in data file that are represented in Y axis. */
    char* y_titles[GRAPH_MAX_COLUMNS];			/**< Titles of series represented in Y axis. */
    char* title;					/**< Title of the graphic. */
    char* xlabel;					/**< Label of X axis. */	
    char* ylabel;					/**< Label of Y axis. */
    char* type; 					/**< Type of graph. */
} qc_graph_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Generates HTML and text reports and data files
*  @param bam_qc_report pointer to the bam_qc_report structure
*  @param inputfilename filename of the processed bam file 
*  @param base_quality base quality for quality normalization
*  @param report_directory directory where output files will be written
*  @param valid flag of report for valid or invalid alignments (0: invalid, 1: valid)
*  @return void
*  
*  Generates the output of QC process: HTML and text reports and data files
*/
void generate_report(bam_qc_report_t bam_qc_report, char* inputfilename, int base_quality, char* report_directory, int valid);

#endif
