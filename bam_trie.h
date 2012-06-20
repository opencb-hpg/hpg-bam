
#ifndef BAM_TRIE_H
#define BAM_TRIE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cprops/trie.h>

#include "alignment.h"
#include "bam.h"
#include "bam_commons.h"
#include "bam_file.h"
#include "commons.h"
#include "dna_map_region.h"

#define MAX_DATASET_LINE_LENGTH		512

/* **************************************
 *      	Structures    		*
 * *************************************/

/**
* @brief Structure for storing results of aligment comparison
*
* Structure for storing results of aligment comparison
*/
typedef struct {
    char* filename;		/**< Chromosome. */
    uint32_t mapped;		/**< Number of mapped alignments. */
    uint32_t right_mapped;	/**< Number of right mapped alignments. */
    uint32_t wrong_mapped;	/**< Number of wrong mapped alignments. */
    uint32_t not_mapped;	/**< Number of not mapped alignments. */
} trie_result_t;

/* **************************************
 *      	Functions    		*
 * *************************************/

/**
 *  @brief Converts bam file into trie
 *  @param filename filename containing the bam
 *  @return cp_trie pointer to the trie structure
 *
 *  Convert bam file into trie
 **/
cp_trie* dna_bam_to_trie(char* filename);

/**
 *  @brief Converts simulated file (dwgsim) into trie
 *  @param filename filename containing the simulated dataset
 *  @return cp_trie pointer to the trie structure
 *
 *  Convert simulated file (dwgsim) into trie
**/
cp_trie* dna_dataset_to_trie(char* filename);

/**
 *  @brief Calculates the intersection between the trie and the bam_file
 *  @param trie the Trie calculated with dataset_to_trie or bam_to_trie
 *  @param filename the BAM filename to calculate the intersection
 *  @return trie_result_t
 *
 *  Calculates the intersection between the trie and the bam_file
 **/
void dna_intersection(cp_trie* trie, char* filename, char* log_file, char* unmapped_bam, unsigned char align_bam, unsigned char soft_hard, trie_result_t* result);

int add_region_to_trie(cp_trie* trie, bam1_t* bam_line, bam_header_t* header);

void print_region(dna_map_region_t* reg);
void print_result(trie_result_t* result, int log);

#endif /* BAM_TRIE_H */
