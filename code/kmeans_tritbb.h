// Kmeans Framework Header
#ifndef KMEANS_TRITBB_H
#define KMEANS_TRITBB_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "timer.h"
#include "kmeans_tri.h"
#include "kmeans_tbb.h"
#include "kmeans.h"

class KmeansStep1tri {
    uchar *particle_data;
    double *centers;
    double *ubound;
    bool *ubound_invalid;
    double *lbound;
    double **center_dists;
    double *short_dist;
    uchar *assignments;
    bool *assignment_change;
    int cluster_count;
    int dimensions;
public:
    void operator() (const blocked_range<size_t>& r) const {
        for (size_t particle_iter = r.begin(); particle_iter != r.end(); ++particle_iter) {
            int cluster_assignment = (int)(assignments[particle_iter]);            
            if( cluster_assignment >= cluster_count || ubound[particle_iter] > short_dist[cluster_assignment] ) {
                // find closest center
                for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
                    if( center_iter != cluster_assignment && ( ubound[particle_iter] > array_load(lbound, particle_iter, center_iter, cluster_count) ||
                        ubound[particle_iter] > center_dists[center_iter][cluster_assignment] ) ) {
                        double dist;
                        if( ubound_invalid[particle_iter] ) {
                            dist = compute_distance<uchar, double>(particle_data, particle_iter, centers, cluster_assignment, dimensions);
                            ubound[particle_iter] = dist;
                            array_store(lbound, particle_iter, cluster_assignment, cluster_count) = dist;
                            ubound_invalid[particle_iter] = false;
                        }
                        else {
                            dist = ubound[particle_iter];
                        }
                        if( dist > array_load(lbound, particle_iter, center_iter, cluster_count) || dist > center_dists[center_iter][cluster_assignment] ) {
                            double dist_new = compute_distance<uchar, double>(particle_data, particle_iter, centers, center_iter, dimensions);
                            array_store(lbound, particle_iter, center_iter, cluster_count) = dist_new;
                            if( dist_new < dist ) {
                                cluster_assignment = center_iter;
                                ubound[particle_iter] = dist_new;
                            }
                        }
                    }
                }
                // assign to closest center
                if( cluster_assignment != (int)(assignments[particle_iter]) ) {
                    assignments[particle_iter] = (uchar)cluster_assignment;
                    *assignment_change = true;
                }
            }
        }    
        return;
    }
    KmeansStep1tri(uchar *_particle_data, double *_centers, double *_ubound, bool *_ubound_invalid, double *_lbound, double **_center_dists, 
        double *_short_dist, uchar *_assignments, int _cluster_count, int _dimensions, bool *_assignment_change) :
        particle_data(_particle_data), centers(_centers), ubound(_ubound), ubound_invalid(_ubound_invalid), lbound(_lbound), center_dists(_center_dists), 
        short_dist(_short_dist), assignments(_assignments), cluster_count(_cluster_count), dimensions(_dimensions), assignment_change(_assignment_change) {}
};            

class KmeansStep2c {
    double *means_new;
    double *centers;
    double *center_shift;
    int dimensions;
public:
    void operator() (const blocked_range<size_t>& r) const {
        for (size_t center_iter = r.begin(); center_iter != r.end(); ++center_iter) {
            center_shift[center_iter] = compute_distance<double, double>(means_new, center_iter, centers, center_iter, dimensions);
        }
        return;
    }
    KmeansStep2c(double *_means_new, double *_centers, double *_center_shift, int _dimensions) :
        means_new(_means_new), centers(_centers), center_shift(_center_shift), dimensions(_dimensions) {}
};            

class KmeansStep2d {
    double *center_shift;
    double *lbound;
    double *ubound;
    bool *ubound_invalid;
    uchar *assignments;
    int cluster_count;
public:
    void operator() (const blocked_range<size_t>& r) const {
        for (size_t particle_iter = r.begin(); particle_iter != r.end(); ++particle_iter) {
            for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
                array_store(lbound, particle_iter, center_iter, cluster_count) = _MAX<double>(0, array_load(lbound, particle_iter, center_iter, cluster_count) - 
                    center_shift[center_iter] - 2*sqrt(array_load(lbound, particle_iter, center_iter, cluster_count)*center_shift[center_iter]));
            }
            ubound[particle_iter] = ubound[particle_iter] + center_shift[(int)(assignments[particle_iter])] + 2*sqrt(ubound[particle_iter]*center_shift[(int)(assignments[particle_iter])]);
            ubound_invalid[particle_iter] = true;
        }
        return;
    }
    KmeansStep2d(double *_center_shift, double *_lbound, double *_ubound, bool *_ubound_invalid, uchar *_assignments, int _cluster_count) :
        center_shift(_center_shift), lbound(_lbound), ubound(_ubound), ubound_invalid(_ubound_invalid), assignments(_assignments), cluster_count(_cluster_count) {}
};            


int kmeans_tritbb(uchar *data, double *&centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize, int thread_count );
void step1_tritbb( uchar *particle_data, double *centers, int particle_count, int cluster_count, uchar *assignments, int dimensions, bool &assignment_change, 
                double *lbound, double *ubound, bool *ubound_invalid, double **center_dists, double *short_dist, int grainSize );
void step2_tritbb( uchar *particle_data, double *&centers, double *&means_new, double *center_shift, int *cluster_sizes, int particle_count, int cluster_count, uchar *assignments, int dimensions, 
                double *lbound, double *ubound, bool *ubound_invalid, int grainSize );

#endif

