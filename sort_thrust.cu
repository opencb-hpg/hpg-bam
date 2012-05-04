/*
  *
 *  Created on: Aug 4, 2011
 *      Author: victor
 */


#include <iostream>
#include <algorithm>
#include <vector>

//#ifdef THRUST-GPU
  #include <thrust/host_vector.h>
  #include <thrust/device_vector.h>
  #include <thrust/sort.h>
  #include <thrust/copy.h>
//#endif
  
#include "bam_commons.h"
#include "sort_thrust.h"

//------------------------------------------------------------------------------------
// function prototypes
//------------------------------------------------------------------------------------

#ifdef THRUST-GPU
  void sort_vector(int* vector, int length);
  void sort_key_value(thrust::host_vector<int>* h_key_p, thrust::host_vector<int>* h_value_p, int length);
#else
  void sort_key_value(int* keys, int* values, int length);
#endif

//------------------------------------------------------------------------------------
// public function to call thrust sorting for a key-value pair 
//------------------------------------------------------------------------------------

#ifdef THRUST-GPU
  
void sort_alignments_by_position(alignments_list_t* list_p, int chromosome) {  
  if (chromosome == ALL_CHROMOSOMES) {
    for (int i=0; i<NUM_OF_CHROMOSOMES; i++) {
      sort_key_value(&(list_p->chromosomes_p[i]->bam_alignment_coordinates_p), &(list_p->chromosomes_p[i]->indices_p), list_p->chromosomes_p[i]->alignment_count);    
    }    
  } else {
      sort_key_value(&(list_p->chromosomes_p[chromosome]->bam_alignment_coordinates_p), &(list_p->chromosomes_p[chromosome]->indices_p), list_p->chromosomes_p[chromosome]->alignment_count);    
  }    
}

#else

void sort_alignments_by_position(alignments_list_t* list_p, int chromosome) {
  if (chromosome == ALL_CHROMOSOMES) {
    for (int i=0; i<NUM_OF_CHROMOSOMES; i++) {
      sort_key_value(list_p->chromosomes_p[i]->bam_alignment_coordinates_p, list_p->chromosomes_p[i]->indices_p, list_p->chromosomes_p[i]->alignment_count);    
    }    
  } else {
      sort_key_value(list_p->chromosomes_p[chromosome]->bam_alignment_coordinates_p, list_p->chromosomes_p[chromosome]->indices_p, list_p->chromosomes_p[chromosome]->alignment_count);    
  }  
}

#endif

#ifdef THRUST-GPU

//------------------------------------------------------------------------------------
// function to call thrust sorting for a single one-dimension vector
//------------------------------------------------------------------------------------

void sort_vector(int* vector, int length) {
    
    thrust::host_vector<int> h_vec(length);
    for(int i=0 ; i<length ; i++) h_vec[i] = vector[i];
    
    thrust::device_vector<int> d_vec = h_vec;
    
    if (time_flag) { start_timer(t1_sort); }
    thrust::sort(d_vec.begin(), d_vec.end());
    if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }
    
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}

//------------------------------------------------------------------------------------
//  function to call thrust sorting for a key-value pair, element count is needed
//------------------------------------------------------------------------------------

void sort_key_value(int* keys, int* values, int length) {    
  thrust::host_vector<int> h_key(length);
  thrust::host_vector<int> h_value(length);
  for(int i=0 ; i<length ; i++) {
    h_key[i] = keys[i];
    h_value[i] = values[i];
  }
    
  thrust::device_vector<int> d_key = h_key;
  thrust::device_vector<int> d_value = h_value;
  
  if (time_flag) { start_timer(t1_sort); }
  thrust::sort_by_key(d_key.begin(), d_key.end(), d_value.begin());
  if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }
  
  thrust::copy(d_key.begin(), d_key.end(), h_key.begin());
  thrust::copy(d_value.begin(), d_value.end(), h_value.begin());  
}

void sort_key_value(thrust::host_vector<int>* h_key, thrust::host_vector<int>* h_value, int length) {
//void sort_key_value(thrust::host_vector<int> h_key, thrust::host_vector<int> h_value, int length) {
    
   h_key->resize(length);
   h_value->resize(length);
   
//  thrust::device_vector<int> d_key(length);// = h_key;
//  thrust::device_vector<int> d_value(length);// = h_value;
  thrust::device_vector<int> d_key = *h_key;
  thrust::device_vector<int> d_value = *h_value;

  printf("********** h_key.size() = %i, length = %i\n", h_key->size(), length);
  
  int key, value;
/*
  if (length>1) 
  {
    int n=length;
    thrust::device_vector<int> d_key_1(n);
    thrust::device_vector<int> d_value_1(n);
    for(int i=0; i<n; i++) {
      d_key_1[i] = h_key[i];//n-i;
      d_value_1[i] = h_value[i];//i;

      d_key[i] = h_key[i];//n-i;
      d_value[i] = h_value[i];//i;
      
    }
    
    thrust::stable_sort_by_key(d_key_1.begin(), d_key_1.end(), d_value_1.begin());
    
    FILE* fd_aux1 = fopen("/tmp/after_sort_1.txt", "w");
    fprintf(fd_aux1, "----- length=%i\n", n);
    for(int i=0 ; i<n ; i++) {
      key = d_key_1[i];
      value = d_value_1[i];
      fprintf(fd_aux1, "key: %i, value: %i\n", key, value);
    }  
    fclose(fd_aux1);
  }
*/  
  
//   if (length>1) {
//     FILE* fd_aux1 = fopen("/tmp/before_sort.txt", "w");
//     fprintf(fd_aux1, "----- length=%i\n", length);
//     for(int i=0 ; i<length ; i++) {
//       key = h_key[i];
//       value = h_value[i];
//       fprintf(fd_aux1, "key: %i, value: %i\n", key, value);
//     }  
//     fclose(fd_aux1);
//   }

//   if (length>1) {
//     FILE* fd_aux3 = fopen("/tmp/before_sort_device.txt", "w");
//     fprintf(fd_aux3, "----- length=%i\n", length);
//     for(int i=0 ; i<length ; i++) {
//       key = d_key[i];
//       value = d_value[i];
//       fprintf(fd_aux3, "key: %i, value: %i\n", key, value);
//     }  
//     fclose(fd_aux3);
//   }

  if (time_flag) { start_timer(t1_sort); }
  thrust::stable_sort_by_key(d_key.begin(), d_key.end(), d_value.begin());
  if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }
  
//   if (length>1) {
//     FILE* fd_aux4 = fopen("/tmp/after_sort_device.txt", "w");
//     fprintf(fd_aux4, "----- length=%i\n", length);
//     for(int i=0 ; i<length ; i++) {
//       key = d_key[i];
//       value = d_value[i];
//       fprintf(fd_aux4, "key: %i, value: %i\n", key, value);
//     }  
//     fclose(fd_aux4);
//   }
   
  thrust::copy(d_key.begin(), d_key.end(), (*h_key).begin());
  thrust::copy(d_value.begin(), d_value.end(), (*h_value).begin());
  
//   if (length>1) {
//     FILE* fd_aux2 = fopen("/tmp/after_sort.txt", "w");
//     fprintf(fd_aux2, "----- length=%i\n", length);
//     for(int i=0 ; i<length ; i++) {
//       key = (*h_key)[i];
//       value = (*h_value)[i];
//       fprintf(fd_aux2, "key: %i, value: %i\n", key, value);
//     }
//     fclose(fd_aux2);
//   }

}

#else

void sort_key_value(int* keys, int* values, int length) {    
  thrust::host_vector<int> h_key(length);
  thrust::host_vector<int> h_value(length);
  for(int i=0 ; i<length ; i++) {
    h_key[i] = keys[i];
    h_value[i] = values[i];
  }
  
  if (time_flag) { start_timer(t1_sort); }
  thrust::sort_by_key(h_key.begin(), h_key.end(), h_value.begin());
  if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }  
  
  for(int i=0 ; i<length ; i++) {
    keys[i] = h_key[i];
    values[i] = h_value[i];
  }
}

// void sort_key_value(int* keys, int* values, int length) {
//   std::sort(keys, keys + length);  // using default comparison (operator <)
// }

#endif

#ifdef THRUST-GPU


#else

void sort_key_value_char(char** keys, int* values, int length) {    
  thrust::host_vector<std::string> h_key(length);
  thrust::host_vector<int> h_value(length);
  
  for(int i=0 ; i<length ; i++) {
    std::string key(keys[i]);
    h_key[i] = key;
    h_value[i] = values[i];
  }
  
  if (time_flag) { start_timer(t1_sort); }
  thrust::sort_by_key(h_key.begin(), h_key.end(), h_value.begin());
  if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }  
  
  for(int i=0 ; i<length ; i++) {
    //keys[i] = h_key[i].c_str();
    values[i] = h_value[i];
  }
}

// void sort_key_value(int* keys, int* values, int length) {
//   std::sort(keys, keys + length);  // using default comparison (operator <)
// }

#endif

void sort_dataset_by_id(aligner_dataset_list_t* list_p) {
    sort_key_value_char(list_p->seq_id_p, list_p->indices_p, list_p->num_lines);
}