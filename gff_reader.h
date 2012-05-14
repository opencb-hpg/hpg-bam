
#ifndef GFF_READER_H
#define GFF_READER_H

#include <stdio.h>
#include "gff_data.h"

#define MAX_GFF_LINES   		25000
#define MAX_GFF_FILE_LINE_LENGTH 	2048

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Reads a gff file and obtains its lines
*  @param filename gff file name
*  @param[in,out] gff_lines_p pointer to gff_lines 
*  @return number of lines read
*  
*  Reads a gff file and obtains its lines
*/
unsigned int gff_file_read(char* filename, gff_line_t* gff_lines_p);

#endif
