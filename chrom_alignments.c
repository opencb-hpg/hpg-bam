
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chrom_alignments.h"


//=====================================================
// chromosome alignment functions
//=====================================================

//--------------------------------------------------------
// chrom_alignments_new
// init and malloc for a chromosome alignment structure
//--------------------------------------------------------

chrom_alignments_t* chrom_alignments_new(int num_alignments) {  
    
  chrom_alignments_t* chrom_alignments_p = (chrom_alignments_t*) calloc(1, sizeof(chrom_alignments_t));
  chrom_alignments_p->bam_alignments_p = (bam1_t**) calloc(num_alignments, sizeof(bam1_t*));

  #ifdef THRUST-GPU
    chrom_alignments_p->bam_alignment_coordinates_p.resize(num_alignments);
    chrom_alignments_p->indices_p.resize(num_alignments);
  #else
    chrom_alignments_p->bam_alignment_coordinates_p = (int*) calloc(num_alignments, sizeof(int));
    chrom_alignments_p->indices_p = (int*) calloc(num_alignments, sizeof(int));
  #endif
  
  chrom_alignments_p->alignment_count_lock = PTHREAD_MUTEX_INITIALIZER;
  chrom_alignments_p->complete_lock = PTHREAD_MUTEX_INITIALIZER;
  chrom_alignments_p->alignment_count = 0;
  chrom_alignments_p->complete = 0;
  chrom_alignments_p->allocated_alignment = num_alignments;
  
  chrom_alignments_p->initial_allocation = num_alignments;

  return chrom_alignments_p;
}

//--------------------------------------------------------

chrom_alignments_t* chrom_alignments_realloc(chrom_alignments_t* chrom_alignments_p, int num_alignments) {  
 
  bam1_t** aux_alignments_p = chrom_alignments_p->bam_alignments_p;

  int* aux_coordinates_p = chrom_alignments_p->bam_alignment_coordinates_p;
  int* aux_indices_p = chrom_alignments_p->indices_p;
  
  chrom_alignments_p->bam_alignments_p = (bam1_t**) calloc(num_alignments, sizeof(bam1_t*));
  
  #ifdef THRUST-GPU
    chrom_alignments_p->bam_alignment_coordinates_p.resize(num_alignments);
    chrom_alignments_p->indices_p.resize(num_alignments);
  #else
    chrom_alignments_p->bam_alignment_coordinates_p = (int*) calloc(num_alignments, sizeof(int));
    chrom_alignments_p->indices_p = (int*) calloc(num_alignments, sizeof(int));  
  #endif
  
  memcpy(chrom_alignments_p->bam_alignments_p, aux_alignments_p, chrom_alignments_p->allocated_alignment * sizeof(bam1_t*));
  memcpy(chrom_alignments_p->bam_alignment_coordinates_p, aux_coordinates_p, chrom_alignments_p->allocated_alignment * sizeof(int));
  memcpy(chrom_alignments_p->indices_p, aux_indices_p, chrom_alignments_p->allocated_alignment * sizeof(int));
  
  //free(aux_alignments_p);
  //free(aux_coordinates_p);
  //free(aux_indices_p);
  
  chrom_alignments_p->allocated_alignment = num_alignments;
  
  //chrom_alignments_p->initial_allocation = num_alignments;

  return chrom_alignments_p;
}

//--------------------------------------------------------
// chrom_alignments_free
// free a chromosome alignment structure
//--------------------------------------------------------

void chrom_alignments_free(chrom_alignments_t* chrom_alignments_p) {  
  if (chrom_alignments_p != NULL) {
    if (chrom_alignments_p->bam_alignments_p != NULL) { free(chrom_alignments_p->bam_alignments_p); }
    
    #ifdef THRUST-GPU
      chrom_alignments_p->bam_alignment_coordinates_p.resize(1);
      chrom_alignments_p->indices_p.resize(1);
    #else
      if (chrom_alignments_p->bam_alignment_coordinates_p != NULL) { free(chrom_alignments_p->bam_alignment_coordinates_p); }
      if (chrom_alignments_p->indices_p != NULL) { free(chrom_alignments_p->indices_p); }
    #endif
    
    free(chrom_alignments_p);
  }
}


//--------------------------------------------------------
// chrom_alignments_get_alignment
// gets bam1_t* from a chromosome alignment structure
//--------------------------------------------------------

bam1_t* chrom_alignments_get_alignment(chrom_alignments_t* chrom_alignments_p, int index) {
  bam1_t* alignment_p = NULL;
  
  pthread_mutex_lock(&chrom_alignments_p->alignment_count_lock);
  
  if (chrom_alignments_p->alignment_count > index) {
    alignment_p = chrom_alignments_p->bam_alignments_p[index];
    //alignment_p = chrom_alignments_p->bam_alignments_p[chrom_alignments_p->alignment_count - index];
  }
  
  pthread_mutex_unlock(&chrom_alignments_p->alignment_count_lock);
  
  return alignment_p;
}

//--------------------------------------------------------
// chrom_alignments_is_complete
// determines if all chromosome alignments have been read
//--------------------------------------------------------

int chrom_alignments_is_complete(chrom_alignments_t* chrom_alignments_p) {
  int complete;
  
  pthread_mutex_lock(&chrom_alignments_p->complete_lock);
  complete = chrom_alignments_p->complete;
  pthread_mutex_unlock(&chrom_alignments_p->complete_lock);
  
  return complete;  
}

//--------------------------------------------------------
// chrom_alignments_set_complete
// set complete if alignments have been read
//--------------------------------------------------------

void chrom_alignments_set_complete(chrom_alignments_t* chrom_alignments_p, int complete) {
  pthread_mutex_lock(&chrom_alignments_p->complete_lock);
  chrom_alignments_p->complete = complete;
  pthread_mutex_unlock(&chrom_alignments_p->complete_lock);
}

//=====================================================
// alignment list functions
//=====================================================

//-----------------------------------------------------
// alignments_list_new
// init and malloc for alignments list
//-----------------------------------------------------

alignments_list_t* alignments_list_new(int num_chromosomes) {
  alignments_list_t* list_p = (alignments_list_t*) calloc(1, sizeof(alignments_list_t));
  list_p->num_chromosomes = num_chromosomes;
  list_p->chromosomes_p = (chrom_alignments_t**) calloc(num_chromosomes, sizeof(chrom_alignments_t*));
  
  return list_p;
}

//-----------------------------------------------------
// alignments_list_free
// free for alignments list
//-----------------------------------------------------

void alignments_list_free(alignments_list_t* list_p) {
  if (list_p != NULL) {
    if (list_p->chromosomes_p != NULL) {
      for(int i=0; i<list_p->num_chromosomes; i++) {
	chrom_alignments_free(list_p->chromosomes_p[i]);
      }
      free(list_p->chromosomes_p);
    }
    free(list_p);
  }
}

//-----------------------------------------------------
// alignments_list_new_chrom_alignment
//-----------------------------------------------------

void alignments_list_new_chrom_alignment(int chromosome, int num_alignments, alignments_list_t* list_p) {
  if (list_p == NULL) {
    LOG_FATAL("ERROR: list is null !!\n");
  }
  
  if (list_p->chromosomes_p[chromosome] != NULL) {
    char log_message[50];
    sprintf(log_message, "ERROR: chromosome %i (alignment list) is already allocated !!\n", chromosome);
    LOG_FATAL(log_message);
  }
  
  list_p->chromosomes_p[chromosome] = chrom_alignments_new(num_alignments);
}


//-----------------------------------------------------
// alignments_list_get_chrom_alignment
//-----------------------------------------------------

chrom_alignments_t* alignments_list_get_chrom_alignment(int chromosome, alignments_list_t* list_p) {
 
  return list_p->chromosomes_p[chromosome];  
}


//-----------------------------------------------------
// alignments_list_insert_batch
//-----------------------------------------------------

void alignments_list_insert_batch(bam_batch_t* batch_p, alignments_list_t* list_p) {
  
  bam1_t* alignment_p;
  chrom_alignments_t* chrom_alignment_p;
    
  if (batch_p->type == SINGLE_CHROM_BATCH) {
    chrom_alignment_p = list_p->chromosomes_p[batch_p->alignments_p[0]->core.tid];

    pthread_mutex_lock(&(chrom_alignment_p->alignment_count_lock));
    for(int i=0; i<batch_p->num_alignments; i++) {
      alignment_p = batch_p->alignments_p[i];
      chrom_alignment_p = list_p->chromosomes_p[alignment_p->core.tid];
    
      if (chrom_alignment_p->alignment_count >= chrom_alignment_p->allocated_alignment) {
	chrom_alignments_realloc(chrom_alignment_p, chrom_alignment_p->allocated_alignment + chrom_alignment_p->initial_allocation);
      }
    
      chrom_alignment_p->bam_alignments_p[chrom_alignment_p->alignment_count] = alignment_p;
      chrom_alignment_p->bam_alignment_coordinates_p[chrom_alignment_p->alignment_count] = (int) alignment_p->core.pos;
      chrom_alignment_p->indices_p[chrom_alignment_p->alignment_count] = chrom_alignment_p->alignment_count;  
      chrom_alignment_p->alignment_count++;
    }  
    pthread_mutex_unlock(&(chrom_alignment_p->alignment_count_lock));
    
  } else {
    
    for(int i=0; i<batch_p->num_alignments; i++) {
      alignment_p = batch_p->alignments_p[i];
      chrom_alignment_p = list_p->chromosomes_p[alignment_p->core.tid];
      
      pthread_mutex_lock(&chrom_alignment_p->alignment_count_lock);
    
      if (chrom_alignment_p->alignment_count >= chrom_alignment_p->allocated_alignment) {
	chrom_alignments_realloc(chrom_alignment_p, chrom_alignment_p->allocated_alignment + chrom_alignment_p->initial_allocation);
      }
    
      chrom_alignment_p->bam_alignments_p[chrom_alignment_p->alignment_count] = alignment_p;
      chrom_alignment_p->bam_alignment_coordinates_p[chrom_alignment_p->alignment_count] = (int) alignment_p->core.pos;
      chrom_alignment_p->indices_p[chrom_alignment_p->alignment_count] = chrom_alignment_p->alignment_count;
      chrom_alignment_p->alignment_count++;
    
      pthread_mutex_unlock(&chrom_alignment_p->alignment_count_lock);
    }  
  }
}

//-----------------------------------------------------
// 
//-----------------------------------------------------

bam_batch_t* get_bam_batch_form_alignments_list(size_t batch_size, alignments_list_t* list_p) {
/*  
  int num_alignments = 0;
  long current_size = 0;
  
  bam_batch_t* batch_p = (bam_batch_t*) calloc(1, sizeof(bam_batch_t));
  chrom_alignments_t* chrom_alignment_p;
  
  while (current_size<batch_size) {
    
    batch_p->alignments_p[num_alignments] = 
    
    num_alignments++;
  }
  
  batch_p->num_alignments = num_alignments;
  
  
  */
}









