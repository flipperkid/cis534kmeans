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

int kmeans_tri(uchar *data, double *&centers, int particle_count, int dimensions, int cluster_count, uchar *assignments);
void step1_tri( uchar *particle_data, double *centers, int particle_count, int cluster_count, uchar *assignments, int dimensions, bool &assignment_change, 
                double *lbound, double *ubound, bool *ubound_invalid, double **center_dists, double *short_dist );
void step2_tri( uchar *particle_data, double *&centers, double *&means_new, double *center_shift, int *cluster_sizes, int particle_count, int cluster_count, uchar *assignments, int dimensions, 
                double *lbound, double *ubound, bool *ubound_invalid );

#endif

