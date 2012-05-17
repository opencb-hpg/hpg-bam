
#include <iostream>
#include <algorithm>
#include <vector>

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

extern "C" {
    #include "bam_commons.h"
    #include "sort_thrust.h"
}

/* ******************************************************
 *    		Private functions (CPU/GPU) 		*
 * *****************************************************/

void sort_key_value(int* keys, int* values, int length);

void sort_key_value_char(char** keys, int* values, int length);

/* **********************************************************************
 *    		Public functions implementations (CPU/GPU) 		*
 * *********************************************************************/

void sort_alignments_by_position(alignments_list_t* list_p, int chromosome) {
    if (chromosome == ALL_CHROMOSOMES) {
        for (int i = 0; i < NUM_OF_CHROMOSOMES; i++) {
            sort_key_value(list_p->chromosomes_p[i]->bam_alignment_coordinates_p, list_p->chromosomes_p[i]->indices_p, list_p->chromosomes_p[i]->alignment_count);
        }
    } else {
        sort_key_value(list_p->chromosomes_p[chromosome]->bam_alignment_coordinates_p, list_p->chromosomes_p[chromosome]->indices_p, list_p->chromosomes_p[chromosome]->alignment_count);
    }
}

void sort_dataset_list_by_id(aligner_dataset_list_t* list_p) {
    sort_key_value_char(list_p->seq_id_p, list_p->indices_p, list_p->num_lines);
}

/* **********************************************************************
 *    		Private functions implementations (CPU/GPU) 		*
 * *********************************************************************/

void sort_key_value(int* keys, int* values, int length) {
    thrust::host_vector<int> h_key(length);
    thrust::host_vector<int> h_value(length);

    //fill the host vector
    for (int i = 0 ; i < length ; i++) {
        h_key[i] = keys[i];
        h_value[i] = values[i];
    }

    if (time_flag) {
        start_timer(t1_sort);
    }
    thrust::sort_by_key(h_key.begin(), h_key.end(), h_value.begin());
    if (time_flag) {
        stop_timer(t1_sort, t2_sort, sort_time);
    }

    //fills back the original vector
    for (int i = 0 ; i < length ; i++) {
        keys[i] = h_key[i];
        values[i] = h_value[i];
    }
}

void sort_key_value_char(char** keys, int* values, int length) {
    thrust::host_vector<std::string> h_key(length);
    thrust::host_vector<int> h_value(length);

    for (int i = 0 ; i < length ; i++) {
        std::string key(keys[i]);
        h_key[i] = key;
        h_value[i] = values[i];
    }

    if (time_flag) {
        start_timer(t1_sort);
    }
    thrust::sort_by_key(h_key.begin(), h_key.end(), h_value.begin());
    if (time_flag) {
        stop_timer(t1_sort, t2_sort, sort_time);
    }

    for (int i = 0 ; i < length ; i++) {
        //keys[i] = h_key[i].c_str();
        values[i] = h_value[i];
    }
}