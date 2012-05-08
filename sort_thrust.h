
#ifndef SORT_THRUST_H
#define SORT_THRUST_H

#include "aligner_dataset.h"
#include "chrom_alignments.h"
#include "commons.h"
#include "log.h"
#include "cuda_commons.h"


/**
*  public functions to sort structures
*/

void sort_alignments_by_position(alignments_list_t* alignments_list_p, int chromosome);
void sort_dataset_list_by_id(aligner_dataset_list_t* aligner_dataset_list_p);

/*
	function to call thrust sorting for a single one-dimension vector
*/

//void sort_vector(int* vector, int length);
  

/*
	function to call thrust sorting for a key-value pair, element count is needed
*/

//void sort_key_value(int* keys, int* values, int length);
//void sort_key_value(thrust::host_vector<int> h_key, thrust::host_vector<int> h_value, int length);


#endif