
#include <iostream>
#include <algorithm>
#include <vector>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

#include "bam_commons.h"

extern "C" {
    #include "sort_thrust.h"
}

/* ******************************************************
 *    		Private functions (CPU/GPU) 		*
 * *****************************************************/

#ifdef THRUST-GPU	//Compute in GPU

void sort_vector(int* vector, int length);
void sort_key_value(thrust::host_vector<int>* h_key_p, thrust::host_vector<int>* h_value_p, int length);

#else			//Compute in CPU

void sort_key_value(int* keys, int* values, int length);

#endif

void sort_key_value_char(char** keys, int* values, int length);

/* **********************************************************************
 *    		Public functions implementations (CPU/GPU) 		*
 * *********************************************************************/

#ifdef THRUST-GPU	//Compute in GPU

void sort_alignments_by_position(alignments_list_t* list_p, int chromosome) {
    if (chromosome == ALL_CHROMOSOMES) {
        for (int i = 0; i < num_of_chromosomes; i++) {
            sort_key_value(&(list_p->chromosomes_p[i]->bam_alignment_coordinates_p), &(list_p->chromosomes_p[i]->indices_p), list_p->chromosomes_p[i]->alignment_count);
        }
    } else {
        sort_key_value(&(list_p->chromosomes_p[chromosome]->bam_alignment_coordinates_p), &(list_p->chromosomes_p[chromosome]->indices_p), list_p->chromosomes_p[chromosome]->alignment_count);
    }
}

#else			//Compute in CPU

void sort_alignments_by_position(alignments_list_t* list_p, int chromosome) {
    if (chromosome == ALL_CHROMOSOMES) {
        for (int i = 0; i < num_of_chromosomes; i++) {
            sort_key_value(list_p->chromosomes_p[i]->bam_alignment_coordinates_p, list_p->chromosomes_p[i]->indices_p, list_p->chromosomes_p[i]->alignment_count);
        }
    } else {
        sort_key_value(list_p->chromosomes_p[chromosome]->bam_alignment_coordinates_p, list_p->chromosomes_p[chromosome]->indices_p, list_p->chromosomes_p[chromosome]->alignment_count);
    }
}

#endif

void sort_dataset_list_by_id(aligner_dataset_list_t* list_p) {
    sort_key_value_char(list_p->seq_id_p, list_p->indices_p, list_p->num_lines);
}

/* **********************************************************************
 *    		Private functions implementations (CPU/GPU) 		*
 * *********************************************************************/

#ifdef THRUST-GPU	//GPU implementation

void sort_vector(int* vector, int length) {
    thrust::host_vector<int> h_vec(length);
    for (int i = 0 ; i < length ; i++) h_vec[i] = vector[i];

    //copy to device vector
    thrust::device_vector<int> d_vec = h_vec;

    if (time_flag) {
        start_timer(t1_sort);
    }
    thrust::sort(d_vec.begin(), d_vec.end());
    if (time_flag) {
        stop_timer(t1_sort, t2_sort, sort_time);
    }
    
    //copy back to host vector
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

//overloaded method: implementation with vectors
void sort_key_value(int* keys, int* values, int length) {
    thrust::host_vector<int> h_key(length);
    thrust::host_vector<int> h_value(length);
    for (int i = 0 ; i < length ; i++) {
        h_key[i] = keys[i];
        h_value[i] = values[i];
    }

    thrust::device_vector<int> d_key = h_key;
    thrust::device_vector<int> d_value = h_value;

    if (time_flag) {
        start_timer(t1_sort);
    }
    thrust::sort_by_key(d_key.begin(), d_key.end(), d_value.begin());
    if (time_flag) {
        stop_timer(t1_sort, t2_sort, sort_time);
    }

    thrust::copy(d_key.begin(), d_key.end(), h_key.begin());
    thrust::copy(d_value.begin(), d_value.end(), h_value.begin());
}

//overloaded method: implementation with thrust::host_vector
void sort_key_value(thrust::host_vector<int>* h_key, thrust::host_vector<int>* h_value, int length) {
    h_key->resize(length);
    h_value->resize(length);

    thrust::device_vector<int> d_key = *h_key;
    thrust::device_vector<int> d_value = *h_value;

    printf("********** h_key.size() = %i, length = %i\n", h_key->size(), length);

    int key, value;

    if (time_flag) {
        start_timer(t1_sort);
    }
    thrust::stable_sort_by_key(d_key.begin(), d_key.end(), d_value.begin());
    if (time_flag) {
        stop_timer(t1_sort, t2_sort, sort_time);
    }

    thrust::copy(d_key.begin(), d_key.end(), (*h_key).begin());
    thrust::copy(d_value.begin(), d_value.end(), (*h_value).begin());
}

#else			//CPU implementation

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

#endif

#ifdef THRUST-GPU


#else

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

#endif
