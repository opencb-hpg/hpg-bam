
#include <stdlib.h>
#include <string.h>

#include "bam.h"
#include "bam_qc_batch.h"

/* ******************************************************
 *      	Function implementations    		*
 * ******************************************************/

bam_qc_batch_t* bam_qc_batch_new(int id) {
    bam_qc_batch_t* bam_qc_batch_p = (bam_qc_batch_t*) calloc(1, sizeof(bam_qc_batch_t));
    bam_qc_batch_p->id = id;

    return bam_qc_batch_p;
}

//-----------------------------------------------------
// bam_qc_batch_list_init
//
// init to zero the object and initialize the lock
//-----------------------------------------------------

// void bam_qc_batch_list_init(bam_qc_batch_list_t* list_p) {
//    //list_p = (bam_qc_batch_list_t*) calloc(1, sizeof(bam_qc_batch_list_t));
//    memset(list_p, 0, sizeof(bam_qc_batch_list_t));
//    list_p->lock = PTHREAD_MUTEX_INITIALIZER;
// }

//-----------------------------------------------------
// bam_qc_batch_list_insert
//
// insert a qc_batch object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

// void bam_qc_batch_list_insert(bam_qc_batch_t* batch_p, bam_qc_batch_list_t* list_p) {
//
//   if (list_p==NULL) return;
//
//   pthread_mutex_lock(&list_p->lock);
//
//   if (list_p->first_p==NULL) {
//     batch_p->prev_p = NULL;
//     batch_p->next_p = NULL;
//     list_p->first_p = batch_p;
//     list_p->last_p = batch_p;
//   } else {
//     list_p->last_p->next_p = batch_p;
//     batch_p->prev_p = list_p->last_p;
//     batch_p->next_p = NULL;
//     list_p->last_p = batch_p;
//   }
//   list_p->length++;
//
//   pthread_mutex_unlock(&list_p->lock);
// }

//-----------------------------------------------------
// bam_qc_batch_list_remove
//
// remove the first bam_qc_batch object in the begining
// of the list, according to fifo order
//-----------------------------------------------------

// bam_qc_batch_t* bam_qc_batch_list_remove(bam_qc_batch_list_t* list_p) {
//
//   if (list_p==NULL) return NULL;
//
//   pthread_mutex_lock(&list_p->lock);
//
//   // just get the first element, and if is not null update
//   // the first-element pointer
//
//   bam_qc_batch_t* batch_p = list_p->first_p;
//
//   if (batch_p != NULL) {
//     list_p->first_p = batch_p->next_p;
//     list_p->length--;
//   }
//
//   pthread_mutex_unlock(&list_p->lock);
//
//   return batch_p;
// }


//-----------------------------------------------------
// bam_qc_batch_list_print
//
// print the content of qc_batch objects in the list
//-----------------------------------------------------

// void bam_qc_batch_list_print(bam_qc_batch_list_t* list_p) {
//
//   pthread_mutex_lock(&list_p->lock);
//
//   printf("Number of items: %i\n", list_p->length);
//
//   bam_qc_batch_t* batch_p = list_p->first_p;
//
//   while(batch_p != NULL) {
//     int strand_counter = 0;
//     int map_quality = 0;
//     int alignment_length = 0;
//
//     for (int i=0; i<batch_p->num_blocks; i++) {
//  strand_counter += batch_p->strand_counter_p[i];
//  map_quality += batch_p->map_quality_p[i];
//  alignment_length += batch_p->alignment_length_p[i];
//     }
//
//     printf("bam_qc_batch id: %i, strand 0: %i, strand 1: %i, mean alignment quality: %i, mean alignment length: %i\n", batch_p->id, (batch_p->num_alignments - strand_counter), strand_counter, map_quality / batch_p->num_alignments, alignment_length);
//     batch_p = batch_p->next_p;
//   }
//
//   pthread_mutex_unlock(&list_p->lock);
// }

//-----------------------------------------------------
// bam_qc_batch_list_items_free
//
// free list itmes
//-----------------------------------------------------

// void bam_qc_batch_list_items_free(bam_qc_batch_list_t* list_p) {
//   bam_qc_batch_t* item_p;
//
//   while ((item_p = bam_qc_batch_list_remove(list_p)) != NULL) {
//     bam_qc_batch_free(item_p, true);
//   }
// }

//-----------------------------------------------------
// bam_qc_batch_free
//
// free memory for this structure,
// the input parameter 'all' indicates free memory
// associated to the pointers in his structure
//-----------------------------------------------------

void bam_qc_batch_free(bam_qc_batch_t* batch_p, int all) {
    if (batch_p == NULL) {
        return;
    }

    if (batch_p->qc_info_p != NULL) { 
        free(batch_p->qc_info_p);
    }
    
    if (batch_p->strand_counter_p != NULL) {
        free(batch_p->strand_counter_p);
    }
    
    if (batch_p->map_quality_p != NULL) { 
        free(batch_p->map_quality_p);
    }
    
    if (batch_p->alignment_length_p != NULL) { 
        free(batch_p->alignment_length_p);
    }
    
    if (batch_p->qc_alignment_p != NULL) { 
        free(batch_p->qc_alignment_p);
    }

    if (all) {
        free_bam1(batch_p->alignments_p, batch_p->num_alignments);
    }

    free(batch_p);
}