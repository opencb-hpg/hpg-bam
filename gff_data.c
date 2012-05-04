
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gff_reader.h"
#include "log.h"


// ------------------------------------------------
//  p r i v a t e  f u n c t i o n s
// ------------------------------------------------

void gff_data_fill_regions_(gff_region_t* gff_regions_p, gff_line_t* gff_lines_p, int num_regions);
void gff_data_overlap_regions_(gff_data_t* gff_data_p);
int gff_data_sort_regions_(const void* gff_region_a, const void* gff_region_b);


//-----------------------------------------------------
// gff_data_new
//-----------------------------------------------------

gff_data_t* gff_data_new(char* gff_filename) {
 
  gff_data_t* gff_data_p = NULL; 
  
  //if (strcmp(gff_filename, "") != 0) {
  if (gff_filename != NULL) {  
    gff_data_p = (gff_data_t*) calloc(1, sizeof(gff_data_t));
    gff_data_p->gff_lines_p = (gff_line_t*) calloc(MAX_GFF_LINES, sizeof(gff_line_t));
    gff_data_p->gff_regions_p = (gff_region_t*) calloc(MAX_GFF_LINES, sizeof(gff_region_t));
    gff_data_p->lock = PTHREAD_MUTEX_INITIALIZER;
    
    gff_data_p->num_regions = gff_file_read(gff_filename, gff_data_p->gff_lines_p);
    gff_data_p->actual_region = 0;
    gff_data_fill_regions_(gff_data_p->gff_regions_p, gff_data_p->gff_lines_p, gff_data_p->num_regions);
    qsort(gff_data_p->gff_regions_p, gff_data_p->num_regions, sizeof(gff_region_t), gff_data_sort_regions_);
    gff_data_overlap_regions_(gff_data_p);
    //gff_data_print_regions(gff_data_p);
  }

  return gff_data_p;
}

//-----------------------------------------------------
// gff_lines_free
//-----------------------------------------------------

void gff_lines_free(gff_line_t* gff_lines_p) {

  int count = 0;
  
  while (gff_lines_p[count].seqname != NULL) {
    
    free(gff_lines_p[count].seqname);
    if (gff_lines_p[count].source != NULL) free(gff_lines_p[count].source);
    if (gff_lines_p[count].feature != NULL) free(gff_lines_p[count].feature);
    if (gff_lines_p[count].score != NULL) free(gff_lines_p[count].score);
    if (gff_lines_p[count].strand != NULL) free(gff_lines_p[count].strand);
    if (gff_lines_p[count].frame != NULL) free(gff_lines_p[count].frame);
    if (gff_lines_p[count].group != NULL) free(gff_lines_p[count].group);
    
    count++;
  }
  
  free(gff_lines_p);
}

//-----------------------------------------------------
// gff_data_free
//-----------------------------------------------------

void gff_data_free(gff_data_t* gff_data_p) { 
  
  if (gff_data_p != NULL) {
    gff_lines_free(gff_data_p->gff_lines_p);
    free(gff_data_p->gff_regions_p);
    free(gff_data_p); 
  }
  
}

//-----------------------------------------------------
// gff_data_print_regions
//-----------------------------------------------------

void gff_data_print_regions(gff_data_t* gff_data_p) {
  
  printf("num_regions: %i\n", gff_data_p->num_regions);
  int count = 0;
  
  while (count < gff_data_p->num_regions) {
    printf("chr: %hi, start: %i, end: %i\n", gff_data_p->gff_regions_p[count].chromosome, gff_data_p->gff_regions_p[count].start, gff_data_p->gff_regions_p[count].end);
    count++;
  }
  
}


//-----------------------------------------------------
// gff_lines_print
//-----------------------------------------------------

void gff_data_print_lines(gff_data_t* gff_data_p) {
  
  printf("num_lines: %i\n", gff_data_p->num_regions);
  int count = 0;
  gff_line_t* gff_line_p;
  
  while (count < gff_data_p->num_regions) {
    gff_line_p = &gff_data_p->gff_lines_p[count];
    printf("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s", gff_line_p->seqname, gff_line_p->source, gff_line_p->feature, gff_line_p->start, gff_line_p->end, gff_line_p->score, gff_line_p->strand, gff_line_p->frame, gff_line_p->group);
    count++;
  }
  
}

//-----------------------------------------------------
// gff_data_alignment_in_region
//-----------------------------------------------------

int gff_data_alignment_in_region(gff_data_t* gff_data_p, int chromosome, int start_coordinate, int end_coordinate, int* regions) {

  int num_regions = 0;
  int first_actual_region = -1;
  int search_region;
  
  pthread_mutex_lock(&(gff_data_p->lock));
  search_region = gff_data_p->actual_region;
  pthread_mutex_unlock(&(gff_data_p->lock));
  
  
  while (search_region < gff_data_p->num_regions) {

    //printf("while loop... region %i [chr: %i, start: %i, end: %i]\n", gff_data_p->actual_region, gff_data_p->gff_regions_p[gff_data_p->actual_region].chromosome, gff_data_p->gff_regions_p[gff_data_p->actual_region].start, gff_data_p->gff_regions_p[gff_data_p->actual_region].end);
    
    gff_region_t* gff_region_p = &(gff_data_p->gff_regions_p[search_region]);
  
    if (chromosome == gff_region_p->chromosome) {
      if ((start_coordinate >= gff_region_p->start) && (end_coordinate <= gff_region_p->end)) {
	//printf("fragment in region...\n");
	if (first_actual_region == -1) { 
	  //pthread_mutex_lock(&(gff_data_p->lock));	  
	  first_actual_region = search_region;
	  //pthread_mutex_unlock(&(gff_data_p->lock));
	}
	
	regions[num_regions++] = search_region;
      } else if (((start_coordinate < gff_region_p->start) && (end_coordinate >= gff_region_p->start)) || 
		 ((start_coordinate <= gff_region_p->end) && (end_coordinate > gff_region_p->end))) {	
	//printf("fragment is overlaped with region (start or end)...\n");
	if (first_actual_region == -1) first_actual_region = search_region;
	regions[num_regions++] = search_region;
      } else if ((start_coordinate < gff_region_p->start) && (end_coordinate > gff_region_p->end)) {
	//printf("fragment greater than region...\n");
	if (first_actual_region == -1) first_actual_region = search_region;
	regions[num_regions++] = search_region;
      } 
    } else if ((chromosome < gff_region_p->chromosome) || (end_coordinate < gff_region_p->start)) {
      break;
    }
      
    search_region++;    
  }
  
  if (first_actual_region != -1) {
    pthread_mutex_lock(&(gff_data_p->lock));
    gff_data_p->actual_region = first_actual_region;
    pthread_mutex_unlock(&(gff_data_p->lock));
  }
  
  return num_regions;
}


//-----------------------------------------------------
// gff_data_batch_in_region
//-----------------------------------------------------

int gff_data_batch_in_region(bam_data_batch_t* batch_p, gff_data_t* gff_data_p) {
    
  if (gff_data_p == NULL) return -1;
  
  int regions[100];
  
  for (int i=0; i < batch_p->num_chromosomes_in_batch; i++) {
    if (gff_data_alignment_in_region(gff_data_p, batch_p->chromosomes[i], batch_p->start_positions[i], batch_p->end_positions[i], regions)) {
      //printf("chromosome: %i, start: %i, end: %i\n", batch_p->chromosomes[i], batch_p->start_positions[i], batch_p->end_positions[i]);
      return 1;
    }
  }
  
  return 0;
}

// ------------------------------------------------
//  p r i v a t e  f u n c t i o n s
// ------------------------------------------------

//-----------------------------------------------------
// gff_data_fill_regions_
//-----------------------------------------------------


void gff_data_fill_regions_(gff_region_t* gff_regions_p, gff_line_t* gff_lines_p, int num_regions) {
  
  int i;
  char chromosome_aux[2];
  
  for (i=0; i < num_regions; i++) {

    //printf("filling with line... %s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s", gff_lines_p[i].seqname, gff_lines_p[i].source, gff_lines_p[i].feature, gff_lines_p[i].start, gff_lines_p[i].end, gff_lines_p[i].score, gff_lines_p[i].strand, gff_lines_p[i].frame, gff_lines_p[i].group);
    
    if (strlen(gff_lines_p[i].seqname) == 4) {
      //chromosome_aux[0] = ' ';
      chromosome_aux[1] = gff_lines_p[i].seqname[3];
    } else if (strlen(gff_lines_p[i].seqname) == 5) {
      chromosome_aux[0] = gff_lines_p[i].seqname[3];
      chromosome_aux[1] = gff_lines_p[i].seqname[4];
    } else {
      LOG_ERROR("gff file sequence does not match 'chrnn' pattern");      
    }
    
    sscanf(chromosome_aux, "%i", &(gff_regions_p[i].chromosome));
    //printf("chromosome: %i\n\n", gff_regions_p[i].chromosome);
    gff_regions_p[i].start = gff_lines_p[i].start - 1;
    gff_regions_p[i].end = gff_lines_p[i].end - 1;
    
    chromosome_aux[0] = ' ';
    chromosome_aux[1] = ' ';
  }
  
}

//-----------------------------------------------------
// gff_data_sort_regions_
//-----------------------------------------------------

int gff_data_sort_regions_(const void* gff_region_a, const void* gff_region_b) {
  
  gff_region_t* region_a = (gff_region_t*) gff_region_a;
  gff_region_t* region_b = (gff_region_t*) gff_region_b;
  
  if (region_a->chromosome < region_b->chromosome) {
    return -1;    
  } else if (region_a->chromosome > region_b->chromosome) {
    return 1;
  } else {
    if (region_a->start < region_b->start) {
      return -1;
    } else if (region_a->start > region_b->start) {
      return 1;
    } else {
      if (region_a->end < region_b->end) {
	return -1;
      } else if (region_a->end > region_b->end) {
	return 1;
      } else {
	return 0;
      }	
    }
  }  
}

//-----------------------------------------------------
// gff_data_overlap_regions_
//-----------------------------------------------------

void gff_data_overlap_regions_(gff_data_t* gff_data_p) {
  
  int num_overlapped_regions = 1;
  gff_region_t* gff_overlapped_regions_p =  (gff_region_t*) calloc(gff_data_p->num_regions, sizeof(gff_region_t));
  gff_region_t* gff_regions_p = gff_data_p->gff_regions_p;
  
  memcpy(&gff_overlapped_regions_p[0], &gff_regions_p[0], sizeof(gff_region_t));
  
  for (int i=1; i<gff_data_p->num_regions; i++) {
    
    if ((gff_overlapped_regions_p[num_overlapped_regions-1].chromosome == gff_regions_p[i].chromosome) && (gff_overlapped_regions_p[num_overlapped_regions-1].start < gff_regions_p[i].end)) {
      gff_overlapped_regions_p[num_overlapped_regions-1].end = gff_regions_p[i].end;  
    } else {      
      memcpy(&gff_overlapped_regions_p[num_overlapped_regions], &gff_regions_p[i], sizeof(gff_region_t));
      num_overlapped_regions++;
    }
    
  }
  
  free(gff_data_p->gff_regions_p);
  gff_data_p->num_regions = num_overlapped_regions;
  gff_data_p->gff_regions_p = gff_overlapped_regions_p;
}



























