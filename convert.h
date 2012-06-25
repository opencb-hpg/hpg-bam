
#ifndef CONVERT_H
#define CONVERT_H

#include "bam.h"
#include "bam_file.h"
#include "commons.h"
#include "log.h"

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Converts a BAM file into a SAM file
*  @param bam_input bam input filename
*  @param sam_input sam input filename
*  @return void
*  
*  Converts a BAM file into a SAM file
*/
void convert_bam_to_sam(char* bam_input, char* sam_input);

/**
*  @brief Converts a SAM file into a BAM file
*  @param sam_input sam input filename
*  @param bam_input bam input filename
*  @return void
*  
*  Converts a SAM file into a BAM file
*/
void convert_sam_to_bam(char* sam_input, char* bam_input);

#endif /* CONVERT_H */
