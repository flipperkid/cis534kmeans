// Kmeans Framework Header
#ifndef KMEANS_H
#define KMEANS_H
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <limits>
#include "timer.h"

#define uchar unsigned char

/**
 * Images are arranged such that packed 8-bit RGB values
 * That is { {H, S, V}, {H, S, V}, ...}
 *
 * @return array[particle_num*num_of_dimensions+which_dimension]
 */
template <typename T> inline T &array_store(T *arr, const int main_iter, const int sub_iter, const int sub_count = 3) {
    return arr[main_iter*sub_count+sub_iter];
}

template <typename T> inline T array_load(const T *arr, const int main_iter, const int sub_iter, const int sub_count = 3) {
    return arr[main_iter*sub_count+sub_iter];
}

int kmeans_serial(uchar *data, float *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments);

int select_centers_serial(uchar *data, int particle_count, int dimensions, int cluster_count, float *centers);
int select_centerspp_serial(uchar *data, int particle_count, int dimensions, int cluster_count, float *centers);

/**
 * Computes
 * (H_1 - H_C)^2 + (S_1 - S_C)^2 + (V_1 - V_C)^2
 *
 * Maybe should be square root of above?
 */
template <typename T1, typename T2> inline float compute_distance(T1 *particle_data, const int particle_iter, T2 *centers, const int center_iter, const int dimensions) {
    float dist = 0;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        dist += pow( (float)array_store<T1>(particle_data, particle_iter, dim_iter, dimensions) - (float)array_load<T2>(centers, center_iter, dim_iter, dimensions), 2 );
    }
    return dist;
}

#endif

