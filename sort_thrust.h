
#ifndef SORT_THRUST_H
#define SORT_THRUST_H

#include "aligner_dataset.h"
#include "chrom_alignments.h"
#include "commons.h"
#include "log.h"
#include "cuda_commons.h"

/* **************************************
 *    		Functions  		*
 * *************************************/

void sort_alignments_by_position(alignments_list_t* alignments_list_p, int chromosome);
void sort_dataset_list_by_id(aligner_dataset_list_t* aligner_dataset_list_p);

#endif  /* SORT_THRUST_H */
