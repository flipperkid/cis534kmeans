#include "kmeans.h"

int kmeans_serial(uchar *particle_data, int particle_count, int dimensions, int cluster_count, uchar *assignments) {
    srand( time(NULL) );

    // assign centers
    double *centers = new double[cluster_count*dimensions];
    double *cluster_sizes = new double[cluster_count];
    if (!select_centers_serial(particle_data, particle_count, dimensions, cluster_count, centers)) {
        printf("Error selecting centers");
        return 0;
    }

    // iterate until clusters converge
    bool assignment_change = true;
    while(assignment_change) {
        assignment_change = false;

        // step 1
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
            
            // find closest center
            int min_dist = INT_MAX;
            uchar cluster_assignment = UCHAR_MAX;
            for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
                int dist = compute_distance(particle_data, particle_iter, centers, center_iter, dimensions);
                if (dist < min_dist) {
                    min_dist = dist;
                    cluster_assignment = (uchar)center_iter;
                }
            }
            
            // assign to closest center
            if (cluster_assignment != assignments[particle_iter]) {
                assignments[particle_iter] = cluster_assignment;
                assignment_change = true;
            }
        }

        // step 2
        if (assignment_change) {

            // initialize centers to 0
            for (int center_arr_iter = 0; center_arr_iter < cluster_count*dimensions; center_arr_iter++) {
                centers[center_arr_iter] = 0;
            }
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                cluster_sizes[cluster_iter] = 0;
            }

            // sum over assigned particle positions
            for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
                int cluster_assignment = (int)assignments[particle_iter];
                cluster_sizes[cluster_assignment]++;
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    centers[cluster_assignment*dimensions+dim_iter] += particle_data[particle_iter*dimensions+dim_iter];
                }
            }

            // divide by number of assigned particles to get mean
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    centers[cluster_iter*dimensions+dim_iter] = floor(centers[cluster_iter*dimensions+dim_iter]/cluster_sizes[cluster_iter] + 0.5);
                }
            }
        }
    }

    // release memory
    delete [] centers;
    return 1;
}

int select_centers_serial(uchar *particle_data, int particle_count, int dimensions, int cluster_count, double *centers) {
    for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
        int particle_select = rand()%particle_count;
        for (int dimension_iter = 0; dimension_iter < dimensions; dimension_iter++) {
            centers[center_iter*dimensions+dimension_iter] = (double)particle_data[particle_select*dimensions+dimension_iter];
        }
    }
    return 1;
}

int compute_distance(uchar *particle_data, int particle_iter, double *centers, int center_iter, int dimensions) {
    int dist = 0;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        dist += abs( (double)particle_data[particle_iter*dimensions+dim_iter] - centers[center_iter*dimensions+dim_iter] );
    }
    return dist;
}
