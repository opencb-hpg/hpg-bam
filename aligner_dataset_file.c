#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "aligner_dataset_file.h"

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

aligner_dataset_file_t* aligner_dataset_fopen(char* filename) {
    return aligner_dataset_fopen(filename, (char*)"r");
}

aligner_dataset_file_t* aligner_dataset_fopen(char* filename, char* mode) {
    FILE* fd = fopen(filename, mode);
    if (fd == NULL) {
        char log_message[200];
        sprintf(log_message, "Error opening file '%.150s' in mode (%s) !!!!!\n", filename, mode);
        LOG_FATAL(log_message);
        return NULL;
    }

    aligner_dataset_file_t* aligner_dataset_file = (aligner_dataset_file_t *) calloc(1, sizeof(aligner_dataset_file_t));

    aligner_dataset_file->filename = filename;
    aligner_dataset_file->mode = mode;
    aligner_dataset_file->fd = fd;
    aligner_dataset_file->num_lines = 0;
   
    return aligner_dataset_file;
}

void aligner_dataset_fclose(aligner_dataset_file_t* aligner_dataset_file) {
    fclose(aligner_dataset_file->fd);
    free(aligner_dataset_file);
}

unsigned int aligner_dataset_read_batch(aligner_dataset_file_t* dataset_file_p, aligner_dataset_batch_t* batch_p, int num_lines) {
    unsigned int count = 0;

    if (num_lines == 0) {
        num_lines = INT_MAX;
    }
  
    char* buffer = (char*) calloc(1, MAX_DATASET_LINE_LENGTH);
    aligner_dataset_line_t** lines_p = batch_p->aligner_dataset_lines_p;
    const char delimiters[] = "\t";
    char* token = NULL;
    int pos;
    
    while((fgets(buffer, MAX_DATASET_LINE_LENGTH, dataset_file_p->fd)) && (count < num_lines)) {
        pos = 0;
        token = strtok(buffer, delimiters);
	
	lines_p[count] = aligner_dataset_line_new(0);

        while (token != NULL) {
            switch(pos) {
	        case 0:
		    lines_p[count]->seq_id = (char*) calloc(strlen(token) + 1, sizeof(char));    
		    strcpy(lines_p[count]->seq_id, token);
		    break;
	        case 1:
		    sscanf(token, "%hi", &lines_p[count]->chromosome);
		    break;
	        case 2:
		    sscanf(token, "%i", &lines_p[count]->start);
		    break;
	        case 3:
		    sscanf(token, "%i", &lines_p[count]->end);
		    break;
	        case 4:
		    sscanf(token, "%hi", &lines_p[count]->strand);
		    break;
	        case 5:
		    sscanf(token, "%hi", &lines_p[count]->num_errors);
		    break;
	        case 6:
		    sscanf(token, "%hi", &lines_p[count]->num_indels);
		    break;
	        default:
		    break;
	    }

            token = strtok(NULL, delimiters);      
            pos++;
	}        

        count++;
    }
    
    batch_p->num_lines = count;
    
    free(buffer);

    return count;
}

unsigned int aligner_dataset_read_list(aligner_dataset_file_t* dataset_file_p, aligner_dataset_list_t* list_p, int num_lines) {
    unsigned int count = 0;
    
    if (num_lines == 0) {
        num_lines = INT_MAX;
    }

    char* buffer = (char*) calloc(1, MAX_DATASET_LINE_LENGTH);
    aligner_dataset_line_t** lines_p = list_p->aligner_dataset_lines_p;
    const char delimiters[] = "\t";
    char* token = NULL;
    int pos;
    
    while((fgets(buffer, MAX_DATASET_LINE_LENGTH, dataset_file_p->fd)) && (count < num_lines)) {
        pos = 0;
        token = strtok(buffer, delimiters);
	
	lines_p[count] = aligner_dataset_line_new(0);

        while (token != NULL) {
            switch(pos) {
	        case 0:
		    lines_p[count]->seq_id = (char*) calloc(strlen(token) + 1, sizeof(char));    
		    strcpy(lines_p[count]->seq_id, token);
		    list_p->seq_id_p[count] = lines_p[count]->seq_id;
		    break;
	        case 1:
		    sscanf(token, "%hi", &lines_p[count]->chromosome);
		    break;
	        case 2:
		    sscanf(token, "%i", &lines_p[count]->start);
		    break;
	        case 3:
		    sscanf(token, "%i", &lines_p[count]->end);
		    break;
	        case 4:
		    sscanf(token, "%hi", &lines_p[count]->strand);
		    break;
	        case 5:
		    sscanf(token, "%hi", &lines_p[count]->num_errors);
		    break;
	        case 6:
		    sscanf(token, "%hi", &lines_p[count]->num_indels);
		    break;
	        default:
		    break;
	    }

            token = strtok(NULL, delimiters);      
            pos++;
	}        

        list_p->indices_p[count] = count;
        count++;
    }
    
    list_p->num_lines = count;
    
    free(buffer);

    return count; 
}

aligner_dataset_batch_t* aligner_dataset_batch_new(int num_lines) {
    aligner_dataset_batch_t* batch_p = (aligner_dataset_batch_t*) calloc(1, sizeof(aligner_dataset_batch_t));    
    batch_p->aligner_dataset_lines_p = (aligner_dataset_line_t**) calloc(num_lines, sizeof(aligner_dataset_line_t*));
    
    return batch_p;  
}

void aligner_dataset_batch_free(aligner_dataset_batch_t* batch_p) {
    if (batch_p == NULL) {
        return;      
    }
  
    for (int i = 0; i < batch_p->num_lines; i++) {
        aligner_dataset_line_free(batch_p->aligner_dataset_lines_p[i], 1);
    }
  
    if (batch_p->aligner_dataset_lines_p != NULL) {
        free(batch_p->aligner_dataset_lines_p);
    }
    
    free(batch_p);
}

void aligner_dataset_batch_print(aligner_dataset_batch_t* batch_p) {
    aligner_dataset_line_t** lines_p = batch_p->aligner_dataset_lines_p;
    
    for (int i = 0; i < batch_p->num_lines; i++) {
        printf("%s\t%hi\t%i\t%i\t%hi\t%hi\t%hi\n", lines_p[i]->seq_id, lines_p[i]->chromosome, lines_p[i]->start, lines_p[i]->end, lines_p[i]->strand, lines_p[i]->num_errors, lines_p[i]->num_indels);
    }
}