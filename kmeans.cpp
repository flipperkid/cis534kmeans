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
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for
            
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
            for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel reduce or gather for each cluster and sum
                int cluster_assignment = (int)assignments[particle_iter];
                cluster_sizes[cluster_assignment]++;
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_access<double>(centers, cluster_assignment, dim_iter, dimensions) += (double)array_access<uchar>(particle_data, particle_iter, dim_iter, dimensions);

                }
            }

            // divide by number of assigned particles to get mean
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_access<double>(centers, cluster_iter, dim_iter, dimensions) = floor(array_access<double>(centers, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5);
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
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_access<double>(centers, center_iter, dim_iter, dimensions) = (double)array_access<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    return 1;
}

int select_centerspp_serial(uchar *particle_data, int particle_count, int dimensions, int cluster_count, double *centers) {
    double *particle_distr = new double[particle_count];
    int distr_length = RAND_MAX;
    int particle_select = rand()%particle_count;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        array_access<double>(centers, 0, dim_iter, dimensions) = (double)array_access<uchar>(particle_data, particle_select, dim_iter, dimensions);
    }
    for (int center_iter = 1; center_iter < cluster_count; center_iter++) {
        // calc probability for each particle as center 
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for
            // find closest center
            double min_dist = std::numeric_limits<double>::max();
            for (int center_dist_iter = 0; center_dist_iter < center_iter; center_dist_iter++) {
                double dist = compute_distance(particle_data, particle_iter, centers, center_dist_iter, dimensions);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
            particle_distr[particle_iter] = min_dist;
        }
        double particle_sum = 0;
        for (int particle_iter = 1; particle_iter < particle_count; particle_iter++) { // TODO parallel prefix sum
            particle_distr[particle_iter] += particle_distr[particle_iter-1];
        }
        double distr_ratio = particle_distr[particle_count-1] / distr_length;
        int select_pos = rand()%distr_length;
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {  // TODO search - make binary, parallelize (note break, tricky)
            if (particle_distr[particle_iter] > select_pos*distr_ratio) {
                particle_select = particle_iter;
                break;
            }
        }
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_access<double>(centers, center_iter, dim_iter, dimensions) = (double)array_access<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    return 1;
}

double compute_distance(uchar *particle_data, int particle_iter, double *centers, int center_iter, int dimensions) {
    double dist = 0;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        dist += pow( (double)array_access<uchar>(particle_data, particle_iter, dim_iter, dimensions) - array_access<double>(centers, center_iter, dim_iter, dimensions), 2 );
    }
    return dist;
}
