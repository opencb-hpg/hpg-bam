
#ifndef FILTER_H
#define FILTER_H

#include "bam.h"
#include "bam_file.h"
#include "commons.h"
#include "log.h"

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Filter a BAM file by a given chromosome
*  @param bam_input bam input filename
*  @param chromosome chromosome to keep in the filtered BAM
*  @return void
*  
*  Filter a BAM file by a given chromosome, output BAM with given chromosome is written to disk
*/
void filter_bam_by_chromosome(char* bam_input, char* output_directory, short int chromosome);

/**
*  @brief Filter a BAM file using input criteria
*  @param bam_input bam input filename
*  @param chromosome chromosome to keep in the filtered BAM
*  @param min_length  min length between paired ends to keep in the filtered BAM
*  @param min_quality min quality between paired ends to keep in the filtered BAM
*  @param max_distance max distance between paired ends to keep in the filtered BAM
*  @return void
*  
*  Filter a BAM file using input criteria (chromosome, alignment length, quality and/or 
*  distance between paired ends), output BAM is written to disk
*/
void filter_bam_by_criteria(char* bam_input, char* output_directory, int max_mismatches, short int chromosome, int min_length, int min_quality, int max_quality, int max_distance);

#endif /* FILTER_H */