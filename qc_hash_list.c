
#include <stdlib.h>
#include <string.h>

#include "bam.h"
#include "qc_hash.h"
#include "qc_hash_list.h"
#include "GeneralHashFunctions.h"


//=====================================================
// functions to manage qc hash list list item
//=====================================================


qc_hash_list_item_t* qc_hash_list_item_new(char* id_seq, int tid, int start_coordinate, int seq_length, int paired_end) {
    qc_hash_list_item_t* item_p;

    item_p = (qc_hash_list_item_t*) calloc(1, sizeof(qc_hash_list_item_t));
    item_p->pairend1_p = (read_position_t*) calloc(ALLOCATED_MAPPINGS_PER_READ, sizeof(read_position_t));
    item_p->pairend2_p = (read_position_t*) calloc(ALLOCATED_MAPPINGS_PER_READ, sizeof(read_position_t));
    item_p->allocated_pairends1 = ALLOCATED_MAPPINGS_PER_READ;
    item_p->allocated_pairends1 = ALLOCATED_MAPPINGS_PER_READ;

    if (paired_end == PAIRED_END1) {
        item_p->num_pairends1 = 1;
        item_p->pairend1_p->tid = tid;
        item_p->pairend1_p->position = start_coordinate + seq_length;
    } else if (paired_end == PAIRED_END2) {
        item_p->num_pairends2 = 1;
        item_p->pairend2_p->tid = tid;
        item_p->pairend2_p->position = start_coordinate;
    }

    //item_p->id_seq = (char*) calloc(strlen(id_seq) - 2, sizeof(char));
    //strncpy(item_p->id_seq, id_seq, strlen(id_seq) - 2);
    item_p->id_seq = (char*) calloc(strlen(id_seq), sizeof(char));
    strncpy(item_p->id_seq, id_seq, strlen(id_seq));

    //item_p->lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&(item_p->lock), NULL);
    item_p->next_p = NULL;
    item_p->prev_p = NULL;

    //printf("qc_hash_list_item_new : %x\n", item_p);

    return item_p;
}

//-----------------------------------------------------
// qc_hash_list_item_free
//
// free memory for this structure,
// the input parameter 'all' indicates free memory
// associated to the pointers in his structure
//-----------------------------------------------------

void qc_hash_list_item_free(qc_hash_list_item_t* item_p, int all) {

    //printf("qc_hash_list_item_free : %x\n", item_p);
    if (all) {
        //printf("\t\tall\n");
        if (item_p->id_seq != NULL) {
            free(item_p->id_seq);
        }
        if (item_p->pairend1_p != NULL) {
            free(item_p->pairend1_p);
        }
        if (item_p->pairend2_p != NULL) {
            free(item_p->pairend2_p);
        }
    }

    free(item_p);
}

//=====================================================
// functions to manage qc hash list
//=====================================================

//-----------------------------------------------------
// qc_hash_list_init
//
// init to zero the object and initialize the lock
//-----------------------------------------------------

void qc_hash_list_init(qc_hash_list_t* list_p) {
    memset(list_p, 0, sizeof(qc_hash_list_t));
    //list_p->lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&(list_p->lock), NULL);
}

//-----------------------------------------------------
// qc_hash_list_lock
//
// lock the qc hast list
//-----------------------------------------------------

void qc_hash_list_lock(qc_hash_list_t* list_p) {
    pthread_mutex_lock(&list_p->lock);
}

//-----------------------------------------------------
// qc_hash_list_unlock
//
// unlock the qc hast list
//-----------------------------------------------------

void qc_hash_list_unlock(qc_hash_list_t* list_p) {
    pthread_mutex_unlock(&list_p->lock);
}

//-----------------------------------------------------
// qc_hash_list_insert
//
// insert an item in the end of the list,
// according to fifo order,
//-----------------------------------------------------

void qc_hash_list_insert(qc_hash_list_item_t* item_p, qc_hash_list_t* list_p) {

    if (list_p == NULL) return;

    pthread_mutex_lock(&list_p->lock);

    if (list_p->first_p == NULL) {
        item_p->prev_p = NULL;
        item_p->next_p = NULL;
        list_p->first_p = item_p;
        list_p->last_p = item_p;
    } else {
        list_p->last_p->next_p = item_p;
        item_p->prev_p = list_p->last_p;
        item_p->next_p = NULL;
        list_p->last_p = item_p;
    }
    list_p->length++;

    pthread_mutex_unlock(&list_p->lock);
}

//-----------------------------------------------------
// bam_data_batch_list_remove
//
// remove the first batch object in the begining
// of the list,
// according to fifo order,
//-----------------------------------------------------

qc_hash_list_item_t* qc_hash_list_remove(qc_hash_list_t* list_p) {

    if (list_p == NULL) return NULL;

    pthread_mutex_lock(&list_p->lock);

    // just get the first element, and if is not null
    // update the first-element pointer
    //
    qc_hash_list_item_t* item_p = list_p->first_p;
    if (item_p != NULL) {
        list_p->first_p = item_p->next_p;
        list_p->length--;
    }

    pthread_mutex_unlock(&list_p->lock);

    return item_p;
}

//-----------------------------------------------------
// qc_hash_list_length
//
// returns list length
//-----------------------------------------------------

int qc_hash_list_length(qc_hash_list_t* list_p) {

    int length = 0;

    if (list_p == NULL) return length;

    pthread_mutex_lock(&list_p->lock);
    length = list_p->length;
    pthread_mutex_unlock(&list_p->lock);

    return length;
}

//-----------------------------------------------------
// qc hash list print
//
//
//-----------------------------------------------------

void qc_hash_list_print(qc_hash_list_t* list_p) {

    if (list_p == NULL) return;

    pthread_mutex_lock(&list_p->lock);

    printf("Number of items: %i\n", list_p->length);

    qc_hash_list_item_t* item_p = list_p->first_p;
    int id = 0;

    while (item_p != NULL) {
        printf("item num: %i, id seq: %i, num_pairends1: %i, num_pairends2: %i\n", id++, item_p->id_seq, item_p->num_pairends1, item_p->num_pairends2);

        printf("\tpaired end 1 -> ");
        for (int i = 0; i < item_p->num_pairends1; i++) {
            printf("tid: %i pos: %i -", item_p->pairend1_p[i].tid, item_p->pairend1_p[i].position);
        }
        printf("\n");

        printf("\tpaired end 2 -> ");
        for (int i = 0; i < item_p->num_pairends2; i++) {
            printf("tid: %i pos: %i -", item_p->pairend2_p[i].tid, item_p->pairend2_p[i].position);
        }
        printf("\n");

        item_p = item_p->next_p;
    }

    pthread_mutex_unlock(&list_p->lock);
}

//-----------------------------------------------------
// qc_hash_list_items_free
//
// free list itmes
//
//-----------------------------------------------------

void qc_hash_list_items_free(qc_hash_list_t* list_p) {
    qc_hash_list_item_t* item_p;

    //if (list_p->length > 0) printf("-------------> list_p->length: %i\n", list_p->length);
    while ((item_p = qc_hash_list_remove(list_p)) != NULL) {
        qc_hash_list_item_free(item_p, 1);
    }
}