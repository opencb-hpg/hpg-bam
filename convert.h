
#ifndef CONVERT_H_
#define CONVERT_H_

#include "bam.h"
#include "bam_file.h"
#include "commons.h"
#include "log.h"


/*
      CONVERT BAM TO SAM
*/

void convert_bam_to_sam(char* bam_input, char* sam_input);

/*
      CONVERT SAM TO BAM
*/

void convert_sam_to_bam(char* sam_input, char* bam_input);


#endif