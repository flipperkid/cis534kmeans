// Kmeans Framework Header
#ifndef KMEANS_TBB_H
#define KMEANS_TBB_H
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <limits>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "timer.h"
#include "kmeans.h"

using namespace tbb;

class KmeansStep1 {
    uchar *particle_data;
    uchar *assignments;
    double *centers;
    const int cluster_count;
    const int dimensions;
    bool *assignment_change;
public:
    void operator()( const blocked_range<size_t>& r ) const {
        for( size_t particle_iter=r.begin(); particle_iter!=r.end(); ++particle_iter ) {
            // find closest center
            double min_dist = std::numeric_limits<double>::max();
            uchar cluster_assignment = std::numeric_limits<uchar>::max();
            for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
                double dist = compute_distance(particle_data, particle_iter, centers, center_iter, dimensions);
                if (dist < min_dist) {
                    min_dist = dist;
                    cluster_assignment = (uchar)center_iter;
                }
            }
            
            // assign to closest center
            if (cluster_assignment != assignments[particle_iter]) {
                assignments[particle_iter] = cluster_assignment;
                *assignment_change = true;
            }
        }
    }
    KmeansStep1(uchar *_particle_data, double *_centers, int _cluster_count, uchar *_assignments, int _dimensions, bool *_assignment_change) : 
                particle_data(_particle_data), cluster_count(_cluster_count), centers(_centers), assignments(_assignments), dimensions(_dimensions), assignment_change(_assignment_change) {}
};

class KmeansStep2a {
    uchar *particle_data;
    double *centers;
    uchar *assignments;
    const int cluster_count;
    const int dimensions;
public:
    double *sumCenters;
    double *sumClusterSizes;
    void operator() (const blocked_range<size_t>& r) const {
        for( size_t particle_iter=r.begin(); particle_iter!=r.end(); ++particle_iter ) {
            int cluster_assignment = (int)assignments[particle_iter];
            sumClusterSizes[cluster_assignment]++;
            for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                array_store<double>(sumCenters, cluster_assignment, dim_iter, dimensions) += (double)array_load<uchar>(particle_data, particle_iter, dim_iter, dimensions);
            }
        }
    }

    KmeansStep2a( KmeansStep2a &x, split ) : 
        particle_data(x.particle_data), centers(x.centers), assignments(x.assignments), cluster_count(x.cluster_count), dimensions(x.dimensions) {
        sumCenters = new double[cluster_count*dimensions];
        sumClusterSizes = new double[cluster_count];
        for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {    
            sumClusterSizes[cluster_iter] = 0;
            for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                array_store<double>(sumCenters, cluster_iter, dim_iter, dimensions) = 0;
            }
        }
    }

    void join( const KmeansStep2a &y ) {
        for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {    
            sumClusterSizes[cluster_iter] += y.sumClusterSizes[cluster_iter];
            for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                array_store<double>(sumCenters, cluster_iter, dim_iter, dimensions) += (double)array_load<double>(y.sumCenters, cluster_iter, dim_iter, dimensions);
            }
        }
    }

    KmeansStep2a(uchar *_particle_data, double *_centers, uchar *_assignments, int _cluster_count, int _dimensions) :
        particle_data(_particle_data), centers(_centers), assignments(_assignments), cluster_count(_cluster_count), dimensions(_dimensions) {
        sumCenters = new double[_cluster_count*_dimensions];
        sumClusterSizes = new double[_cluster_count];
        for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {    
            sumClusterSizes[cluster_iter] = 0;
            for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                array_store<double>(sumCenters, cluster_iter, dim_iter, dimensions) = 0;
            }
        }
    }

    ~KmeansStep2a() {
        delete [] sumCenters;
        delete [] sumClusterSizes;
    }
};            

class KmeansStep2b {
    double *centers;
    double *center_sums;
    double *cluster_sizes;
    int dimensions;
public:
    void operator() (const blocked_range<size_t>& r) const {
        for (size_t cluster_iter = r.begin(); cluster_iter != r.end(); ++cluster_iter) {
            for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                array_store<double>(centers, cluster_iter, dim_iter, dimensions) =
                    floor(array_load<double>(center_sums, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5);
            }
        }
        return;
    }
    KmeansStep2b(double *_centers, double *_center_sums, double *_cluster_sizes, int _dimensions) :
        centers(_centers), center_sums(_center_sums), cluster_sizes(_cluster_sizes), dimensions(_dimensions) {}
};            

int kmeans_tbb( uchar *data, double *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize );

#endif
