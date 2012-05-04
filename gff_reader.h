
#ifndef GFF_READER_H
#define GFF_READER_H

#include <stdio.h>
#include "gff_data.h"

#define MAX_GFF_LINES			25000
#define MAX_GFF_FILE_LINE_LENGTH	2048

//====================================================================================
//  gff_reader.h
//
//  structures and methods for reading gff files
//====================================================================================

// ------------------------------------------------
// gff reader functions 
// ------------------------------------------------

unsigned int gff_file_read(char* filename, gff_line_t* gff_lines_p);

#endif
