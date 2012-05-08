
#ifndef GFF_DATA_H_
#define GFF_DATA_H_

#include <stdio.h>

#include "bam_data_batch_list.h"
#include "commons.h"

//====================================================================================
//  gff_data.h
//
//  structures and methods for storing and managing gff file data
//====================================================================================

typedef struct gff_line {
  char* seqname;
  char* source;
  char* feature;
  int start;
  int end;
  char* score;
  char* strand;
  char* frame;
  char* group;
} gff_line_t;

typedef struct gff_region {
  short int chromosome;
  int start;
  int end;
} gff_region_t;

typedef struct gff_data {
  int num_regions;
  int actual_region;
  gff_line_t* gff_lines_p;
  gff_region_t* gff_regions_p;
  pthread_mutex_t lock;
} gff_data_t;


// ------------------------------------------------
// gff data functions 
// ------------------------------------------------

gff_data_t* gff_data_new(char* gff_filename);
//gff_data_t* gff_data_new();
void gff_lines_free(gff_line_t* gff_lines_p);
void gff_data_free(gff_data_t* gff_data_p);
void gff_data_print_regions(gff_data_t* gff_data_p);
void gff_data_print_lines(gff_data_t* gff_data_p);
int gff_data_alignment_in_region(gff_data_t* gff_data_p, int chromosome, int start_coordinate, int end_coordinate, int* regions);
int gff_data_batch_in_region(bam_data_batch_t* bam_data_batch_p, gff_data_t* gff_data_p);

#endif