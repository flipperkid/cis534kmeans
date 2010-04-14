// Kmeans Framework Header
#ifndef KMEANS_H
#define KMEANS_H
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <limits.h>

#define uchar unsigned char
int kmeans_serial(uchar *data, int particle_count, int dimensions, int cluster_count, uchar *assignments);

int select_centers_serial(uchar *data, int particle_count, int dimensions, int cluster_count, double *centers);

int compute_distance(uchar *particle_data, int particle_iter, double *centers, int center_iter, int dimensions);

#endif

