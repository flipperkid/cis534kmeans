// Kmeans Framework Header
#ifndef KMEANS_OMP_H
#define KMEANS_OMP_H
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <limits>
#include "timer.h"
#include <omp.h>

#define uchar unsigned char

/**
 * Images are arranged such that packed 8-bit RGB values
 * That is { {H, S, V}, {H, S, V}, ...}
 *
 * @return array[particle_num*num_of_dimensions+which_dimension]
 */
/*
template <typename T> 
inline T &array_store(T *arr, const int main_iter, const int sub_iter, const int sub_count = 3)
{
    return arr[main_iter*sub_count+sub_iter];
}

template <typename T> 
inline T array_load(const T *arr, const int main_iter, const int sub_iter, const int sub_count = 3)
{
    return arr[main_iter*sub_count+sub_iter];
}
*/

int kmeans_omp(uchar *data, double *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize, int thread_count);

// int select_centers_omp(uchar *data, int particle_count, int dimensions, int cluster_count, double *centers);
int select_centerspp_omp(uchar *data, int particle_count, int dimensions, int cluster_count, double *centers);

// double compute_distance(uchar *particle_data, const int particle_iter, double *centers, const int center_iter, const int dimensions);

#endif

