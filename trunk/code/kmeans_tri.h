// Kmeans Framework Header
#ifndef KMEANS_TRI_H
#define KMEANS_TRI_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits>
#include "timer.h"
#include "kmeans.h"

template <typename T> inline T _MAX(T val1, T val2) {
    if(val1 > val2) return val1;
    return val2;
}

int kmeans_tri(uchar *data, float *&centers, int particle_count, int dimensions, int cluster_count, uchar *assignments);
void step1_tri( uchar *particle_data, float *centers, int particle_count, int cluster_count, uchar *assignments, int dimensions, bool &assignment_change, 
                float *lbound, float *ubound, bool *ubound_invalid, float **center_dists, float *short_dist );
void step2_tri( uchar *particle_data, float *&centers, float *&means_new, float *center_shift, int *cluster_sizes, int particle_count, int cluster_count, uchar *assignments, int dimensions, 
                float *lbound, float *ubound, bool *ubound_invalid );

#endif

