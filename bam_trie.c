#include "bam_trie.h"

void generate_logs(const char *file, cp_trie *trie);

cp_trie*  dna_bam_to_trie(char * file) {
    // Create trie
    cp_trie *trie = cp_trie_create(0);

    int read_bytes;
    //  char* bam_string;
    int cont = 0;
    bam1_t* bam_p = bam_init1();

    // Open bam file for read
    bamFile bam_fd = bam_open(file, "r");

    //header for BAM file has been done in the opening
    bam_header_t* bam_header_p = bam_header_read(bam_fd);

    while ((read_bytes = bam_read1(bam_fd, bam_p)) > 0) {
        add_region_to_trie(trie, bam_p, bam_header_p);
        cont++;
    }

    // Free Memory
    bam_close(bam_fd);
    bam_destroy1(bam_p);
    return trie;
}

cp_trie *dna_dataset_to_trie(char * file) {
    char* buffer = (char*) calloc(1, MAX_DATASET_LINE_LENGTH);
    const char delimiters[] = "\t";
    char* token = NULL;
    char* id = "";
    char s;
    int pos;
    short int tam;
    FILE* fd = fopen(file, "r");
    FILE* fb = fopen("id.tmp", "wb");

    // Create trie
    cp_trie* trie = cp_trie_create(0);
    //  cp_trie *trie = cp_trie_create_trie(COLLECTION_MODE_DEEP,NULL,(cp_destructor_fn) region_free);
    dna_map_region_t* region;

    if (fd == NULL || fb == NULL) { // Mejorar esta gestión de errores
        printf("Fallo al abrir el fichero");
        return NULL;
    }

    while (fgets(buffer, MAX_DATASET_LINE_LENGTH, fd)) {
        pos = 0;
        token = strtok(buffer, delimiters);
        region = (dna_map_region_t*) calloc(1, sizeof(dna_map_region_t));

        while (token != NULL) {
            switch (pos) {
                case 0:
                    region->chromosome = (char*) calloc(strlen(token) + 1, sizeof(char));
                    strcpy(region->chromosome, token);
                    break;
                case 1:
                    sscanf(token, "%i", &region->start_position);
                    break;
                case 2:
                    sscanf(token, "%i", &region->end_position);
                    break;
                case 3:
                    sscanf(token, "%c", &s);
                    region->strand = (s == '0' ? 0 : 1);
                    // printf("%c%i\n",s,region->strand);
                    break;
                case 4:
                    id = (char*) calloc(strlen(token) + 1, sizeof(char));
                    sscanf(token, "%s", id);
                    break;
            }
            token = strtok(NULL, delimiters);
            pos++;
        }

        tam = strlen(id) + 1;
        fwrite(&tam, sizeof(short int), 1, fb);
        fwrite(id, tam, 1, fb);
        cp_trie_add(trie, id, region);

        free(id);
    }
    
    //close and free resources
    free(buffer);
    fclose(fd);
    fclose(fb);
    return trie;
}


int add_region_to_trie(cp_trie *trie, bam1_t *bam_line, bam_header_t *header) {
    // Create new region
    dna_map_region_t * region = (dna_map_region_t *) calloc(1, sizeof(dna_map_region_t));

    // Add chromosome info
    int tid = bam_line->core.tid;

    if (tid == -1) {
        region->chromosome = "*";
    } else {
        region->chromosome = (char*) calloc(strlen(header->target_name[tid]) + 1, sizeof(char));
        strcpy(region->chromosome, header->target_name[tid]);
    }

    // Add start position info
    region->start_position = bam_line->core.pos; // TODO mirar el FLAG si es 1-base o 0-base

    // Add end position info
    region->end_position = (int) bam_calend(&bam_line->core, bam1_cigar(bam_line));


    // Add strand info
    uint32_t flag = (uint32_t) bam_line->core.flag;
    region->strand = (flag & BAM_FREVERSE) ? 1 : 0;

    // Get the Seq name Id
    char * id = (char*) calloc(bam_line->core.l_qname, sizeof(char));
    strcpy(id, bam1_qname(bam_line));

    cp_trie_add(trie, id, region);
    //  free(id);
    return 1;
}

void dna_intersection(cp_trie *trie, char * filename, char * log_file, char *unmapped_bam, unsigned char align_bam, unsigned char soft_hard, trie_result_t* result) {
    int read_bytes;
    // int log = 0;
    unsigned char equal_res;

    char *id;
    dna_map_region_t* region;
    dna_map_region_t* region_trie = NULL;

    FILE* ftlog = fopen(log_file, "w");

    bam1_t* bam_line = bam_init1();

    // Open bam file for read
    bamFile bam_fd =  bam_open(filename, "r");
    bamFile f_unmapped = bam_open(unmapped_bam, "w");

    //header for BAM file has been done in the opening
    bam_header_t* bam_header_p = bam_header_read(bam_fd);

    bam_header_write(f_unmapped, bam_header_p);

    // Create new region
    region = (dna_map_region_t *) calloc(1, sizeof(dna_map_region_t));

    while ((read_bytes = bam_read1(bam_fd, bam_line)) > 0) {
        id = (char*) calloc(bam_line->core.l_qname, sizeof(char));
        strcpy(id, bam1_qname(bam_line));

        if (bam_line->core.flag == 4) {
            if (align_bam == 1) {
                if (equal_res) {
                    result->right_mapped++;
                    region_trie->rwmapped = 1;
                } else {
                    printf("entra");
                    result->wrong_mapped++;
                    fprintf(ftlog, "Query: %s\n", id);
                    fprintf(ftlog, "Case 1: ");
                    dna_fprint_region(ftlog, region_trie);
                    fprintf(ftlog, "Case 2: ");
                    dna_fprint_region(ftlog, region);
                }
                result->mapped++;
            } else {
                result->not_mapped++;
            }
        } else if ((region_trie = cp_trie_exact_match(trie, bam1_qname(bam_line))) == NULL) {
            result->not_mapped++;
            bam_write1(f_unmapped, bam_line);
        } else {
            result->mapped++;

            region_trie->mapped = 1;

            int tid = bam_line->core.tid;

            region->chromosome = (char*) calloc(strlen(bam_header_p->target_name[tid]) + 1, sizeof(char));
            strcpy(region->chromosome, bam_header_p->target_name[tid]);

            // Add start position info
            region->start_position = bam_line->core.pos;

            // Add end position info
            region->end_position = (int) bam_calend(&bam_line->core, bam1_cigar(bam_line))   ;
            // region->end_position = region->start_position +  (int) bam_cigar2qlen(&bam_line->core, bam1_cigar(bam_line)) -1;

            //printf("trie: %d - %d\n", region_trie->start_position, region_trie->end_position);
            //printf("bam:  %d - %d\n\n", region->start_position, region->end_position);

            // Add strand info
            uint32_t flag = (uint32_t) bam_line->core.flag;
            region->strand = (flag & BAM_FREVERSE) ? 1 : 0;

            equal_res = (soft_hard == 0) ? dna_map_region_equal_soft(region_trie, region) : dna_map_region_equal_hard(region_trie, region);

            if (equal_res) {
                result->right_mapped++;
                region_trie->rwmapped = 1;
            } else {
                result->wrong_mapped++;
                fprintf(ftlog, "Query: %s\n", id);
                fprintf(ftlog, "Case 1: ");
                dna_fprint_region(ftlog, region_trie);
                fprintf(ftlog, "Case 2: ");
                dna_fprint_region(ftlog, region);
            }
            
            free(region->chromosome);
        }
        
        free(id);
    }
    result->filename = (char*) calloc(strlen(filename) + 1, sizeof(char));

    strcpy(result->filename, filename);

    // Free Memory    
    free(region);
    //  free(region_trie);
    bam_close(bam_fd);
    bam_close(f_unmapped);
    bam_destroy1(bam_line);
    fclose(ftlog);
}


void print_result(trie_result_t *result, int log) {
    double total = result->mapped + (double)result->not_mapped;
    if (log == 0) {
        printf("File: %s\n", result->filename);
        printf("Mapped: %d\n", result->mapped);
        printf("\tRight mapped: %d\n", (result->right_mapped));
        printf("\tWrong mapped: %d\n", (result->wrong_mapped));
        printf("Not mapped: %d \n\n", result->not_mapped);
    } else {
        printf("%s\t%d\t%lf\t%d\t%lf\t%d\t%lf\t%d\t%lf\n", result->filename,
               (result->mapped),
               (result->mapped / total),
               (result->not_mapped),
               (result->not_mapped / total),
               (result->right_mapped),
               (result->right_mapped / total),
               (result->wrong_mapped),
               (result->wrong_mapped / total)
              );
    }
}

