
#ifndef MAP_VALIDATION_H
#define MAP_VALIDATION_H

#include <stdlib.h>
#include <stdio.h>

#include "commons.h"
#include "log.h"

#define MAX_NUM_OF_BAM_FILES_TO_VALIDATE	100

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Performs validation of mappings in a BAM file 
*  @param dna_rna option to indicate dna or rna validation (1: dna, 2: rna)
*  @param align_bam option to indicate if validation comparison is made against an alignment 
*                   dataset file or a BAM file (0: alignment dataset, 1 => BAM)
*  @param soft_hard option to indicate soft or hard/full comparison between alignments (0: soft, 1: hard)
*  @param bam_files bam files to validate
*  @param ref_file reference file to compare
*  @param wrong_mapped_filename filename of the 
*  @return void
*  
*  Performs validation of mappings in a BAM file by comparison with the reference file 
*/

void bam_map_validate(int dna_rna, int align_bam, int soft_hard, char* bam_files, char* ref_file, char* wrong_mapped_filename);

#endif  /* MAP_VALIDATION_H */
