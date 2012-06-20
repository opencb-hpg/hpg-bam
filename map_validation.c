
#include "bam_trie.h"
#include "map_validation.h"

/* ******************************************************
 *      	Functions implementations   		*
 * *****************************************************/

void bam_map_validate(int dna_rna, int align_bam, int soft_hard, char* bam_files, char* ref_file, char* wrong_mapped_filename) {
    LOG_DEBUG("BAM-MAP-VALIDATE: START\n");
    
    int pos = 0;
    int log = 1;
    cp_trie* trie;
    trie_result_t* result;
    char* bam_file_token[MAX_NUM_OF_BAM_FILES_TO_VALIDATE];
    
    //parsing the bam file to compare
    bam_file_token[pos++] = strtok(bam_files, ",");
    
    while (bam_file_token[pos - 1] != NULL) {    
        bam_file_token[pos++] = strtok(NULL, ",");
    }

    //trie structure is filled with alignment data from reference file
    if (align_bam == 0 && dna_rna == 1) {         //dataset reference file
        trie = dna_dataset_to_trie(ref_file);
    } else if (align_bam == 1 && dna_rna == 1) {  //bam reference file
        trie = dna_bam_to_trie(ref_file);
    }

    if (dna_rna == 1) {
        for (int i = 0; i < (pos - 1); i++) {
            result = (trie_result_t*) calloc(1, sizeof(trie_result_t));
            dna_intersection(trie, bam_file_token[i], wrong_mapped_filename, "unmapped", align_bam, soft_hard, result);
            print_result(result, log);             
            free(result);
        }
    }

    free(ref_file);
    cp_trie_destroy(trie);

    LOG_DEBUG("BAM-MAP-VALIDATE: END\n");
}
