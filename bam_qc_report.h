
#ifndef BAM_QC_REPORT_H
#define BAM_QC_REPORT_H

#include <stdio.h>

#include "bam_commons.h"
#include "bam_qc_batch.h"
#include "file_utils.h"

#define QC_SUFFIX			".qc"
#define MAP_ERRORS_HISTOGRAM_SUFFIX    	".map.error.hist.dat"
#define NUM_MAPPINGS_HISTOGRAM_SUFFIX	".num.map.hist.dat"
#define HTML_FILE_SUFFIX		".qc.html"
#define VALID_FILE_SUFFIX		".valid"
#define INVALID_FILE_SUFFIX		".invalid"



#define POS_COLUMN	        1
#define MAP_ERRORS_COLUMN	2
#define DELETION_COLUMN		3
#define INSERTION_COLUMN	4
#define MATCHING_COLUMN		5

#define NUM_ALIGNMENTS_COLUMN	2

#define GRAPH_MAX_COLUMNS	10
#define MAX_TITLE_LENGTH	25
#define MAX_MAP_ERRORS_IN_HISTOGRAM	10


typedef struct bam_qc_report {
  unsigned long num_alignments;
  unsigned long strand_counter;
  unsigned long mean_map_quality;
  unsigned long mean_alignment_length;
  unsigned long mean_paired_end_distance;
  unsigned int map_error_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];
  unsigned int map_deletion_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];
  unsigned int map_insertion_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];
  unsigned int map_matching_histogram[MAX_MAP_ERRORS_IN_HISTOGRAM + 2];
  unsigned int* num_mappings_histogram;
  //unsigned int num_mappings_histogram[MAX_MAPPING_COUNT_IN_HISTOGRAM + 2];
} bam_qc_report_t;

typedef struct bam_qc_graph {
  int x_autoscale;
  int x_start;
  int x_end;
  int y_autoscale;
  int y_start;
  int y_end;
  int lmargin;
  int rmargin;	
  int tmargin;
  int bmargin;
  int num_y_columns;
  int x_column;
  int y_columns[GRAPH_MAX_COLUMNS];
  char* y_titles[GRAPH_MAX_COLUMNS];
  char* title;
  char* xlabel;
  char* ylabel;
  char* type; 
} qc_graph_t;

// generate report function
void generate_report(bam_qc_report_t bam_qc_report, char* inputfilename, int base_quality, char* report_directory, int valid);



#endif
