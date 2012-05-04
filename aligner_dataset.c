
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "aligner_dataset.h"
#include "sort_thrust.h"

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

aligner_dataset_line_t* aligner_dataset_line_new(short int seq_id_length) {
    aligner_dataset_line_t* aligner_dataset_line_p = (aligner_dataset_line_t*) calloc(1, sizeof(aligner_dataset_line_t));
    
    if (seq_id_length > 0) {
        aligner_dataset_line_p->seq_id = (char*) calloc(seq_id_length, sizeof(char));
    }
    
    return aligner_dataset_line_p;    
}

void aligner_dataset_line_free(aligner_dataset_line_t* aligner_dataset_line_p, int all) {
    if (aligner_dataset_line_p == NULL) {
        return;
    }
  
    if (all) {
        if (aligner_dataset_line_p->seq_id != NULL) {
            free(aligner_dataset_line_p->seq_id);
	}
    }
  
    free(aligner_dataset_line_p);
}

aligner_dataset_list_t* aligner_dataset_list_new(size_t num_lines) {
    aligner_dataset_list_t* list_p = (aligner_dataset_list_t*) calloc(1, sizeof(aligner_dataset_list_t));
    
    list_p->num_lines = 0;
    list_p->allocated_lines = num_lines;
    list_p->aligner_dataset_lines_p = (aligner_dataset_line_t**) calloc(num_lines, sizeof(aligner_dataset_line_t*));
    list_p->seq_id_p = (char**) calloc(num_lines, sizeof(char*));
    list_p->indices_p = (int*) calloc(num_lines, sizeof(int));
    
    return list_p;
}

void aligner_dataset_list_free(aligner_dataset_list_t* list_p) {
    if (list_p == NULL) {
        return;
    }
  
    for (int i=0; i < list_p->num_lines; i++) {
        if (list_p->aligner_dataset_lines_p[i] != NULL) {
	    free(list_p->aligner_dataset_lines_p[i]);
	}
	
	if (list_p->seq_id_p[i] != NULL) {
	    free(list_p->seq_id_p[i]);
	}
    }
    
    free(list_p->indices_p);  
}

void aligner_dataset_list_insert_line(aligner_dataset_list_t* list_p, aligner_dataset_line_t* line_p) {
    list_p->aligner_dataset_lines_p[list_p->num_lines] = line_p;
    list_p->seq_id_p[list_p->num_lines] = line_p->seq_id;
    list_p->indices_p[list_p->num_lines] = list_p->num_lines++;
}

aligner_dataset_list_t* aligner_dataset_list_realloc(aligner_dataset_list_t* list_p, size_t num_lines) {
    aligner_dataset_line_t** lines_aux_p = list_p->aligner_dataset_lines_p;
    char** seq_id_aux_p = list_p->seq_id_p;
    int* indices_aux_p = list_p->indices_p;
    
    list_p->aligner_dataset_lines_p = (aligner_dataset_line_t**) calloc(num_lines, sizeof(aligner_dataset_line_t*));
    list_p->seq_id_p = (char**) calloc(num_lines, sizeof(char*));
    list_p->indices_p = (int*) calloc(num_lines, sizeof(int));
    
    memcpy(list_p->aligner_dataset_lines_p, lines_aux_p, list_p->num_lines * sizeof(aligner_dataset_line_t*));
    memcpy(list_p->seq_id_p, seq_id_aux_p, list_p->num_lines * sizeof(char));
    memcpy(list_p->indices_p, indices_aux_p, list_p->num_lines * sizeof(int));

    free(lines_aux_p);
    free(seq_id_aux_p);
    free(indices_aux_p);

    list_p->num_lines = num_lines;
}

void aligner_dataset_list_sort_by_id(aligner_dataset_list_t* list_p) {
    sort_dataset_by_id(list_p);
}

void aligner_dataset_list_write(aligner_dataset_list_t* list_p, char* filename) {
    FILE* dataset_fd = fopen(filename, "w");  
    aligner_dataset_line_t* line_p;
    
    for (int i = 0; i < list_p->num_lines; i++) {
        line_p = list_p->aligner_dataset_lines_p[list_p->indices_p[i]];
        fprintf(dataset_fd, "%s\t%hi\t%i\t%i\t%hi\t%hi\t%hi\n", line_p->seq_id, line_p->chromosome, line_p->start, line_p->end, line_p->strand, line_p->num_errors, line_p->num_indels);
    }

    fflush(dataset_fd);
    fclose(dataset_fd);
}










