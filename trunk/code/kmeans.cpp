#include "kmeans.h"

int kmeans_serial(uchar *particle_data, float *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
//    srand( time(NULL) );

    // assign centers
    int *cluster_sizes = new int[cluster_count];
    if (select_centers_serial(particle_data, particle_count, dimensions, cluster_count, centers) != 0 ) {
        printf("Error selecting centers");
        return -1;
    }

    // iterate until clusters converge
    start_timer( 3 );
    int iterations = 0;
    bool assignment_change = true;
    while(assignment_change) {
        assignment_change = false;

        // step 1
        start_timer( 1 );
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for            
            // find closest center
            int cluster_assignment = (int)assignments[particle_iter];
            float min_dist;
            if(cluster_assignment < cluster_count) {
                min_dist = compute_distance<uchar, float>(particle_data, particle_iter, centers, cluster_assignment, dimensions);
            }
            else {
                min_dist = std::numeric_limits<float>::max();
            }
            for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
                if( center_iter != cluster_assignment ) {
                    float dist = compute_distance<uchar, float>(particle_data, particle_iter, centers, center_iter, dimensions);
                    if (dist < min_dist) {
                        min_dist = dist;
                        cluster_assignment = center_iter;
                    }
                }
            }
            // assign to closest center
            if (cluster_assignment != (int)assignments[particle_iter]) {
                assignments[particle_iter] = (uchar)cluster_assignment;
                assignment_change = true;
            }
        }
        stop_timer( 1 );
    
        // step 2
        if (assignment_change) {
            start_timer( 2 );
        
            // initialize centers to 0
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_store<float>(centers, cluster_iter, dim_iter, dimensions) = 0;
                }
            }
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                cluster_sizes[cluster_iter] = 0;
            }

            // sum over assigned particle positions
            for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel reduce or gather for each cluster and sum
                int cluster_assignment = (int)assignments[particle_iter];
                cluster_sizes[cluster_assignment]++;
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_store<float>(centers, cluster_assignment, dim_iter, dimensions) += (float)array_load<uchar>(particle_data, particle_iter, dim_iter, dimensions);

                }
            }

            // divide by number of assigned particles to get mean
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_store<float>(centers, cluster_iter, dim_iter, dimensions) =
                        floor(array_load<float>(centers, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5f);
                }
            }
            stop_timer( 2 );
        }
        iterations++;
    }
    stop_timer( 3 );
    printf("Step1: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, 1 threads\n", get_time_elapsed( 1 ), particle_count, cluster_count, dimensions, iterations);
    printf("Step2: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, 1 threads\n", get_time_elapsed( 2 ), particle_count, cluster_count, dimensions, iterations);
    printf("Total: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, 1 threads\n", get_time_elapsed( 3 ), particle_count, cluster_count, dimensions, iterations);

    // release memory
    delete [] cluster_sizes;
    return 0;
}


/**
 * Naive cluster select
 */
int select_centers_serial(uchar *particle_data, int particle_count, int dimensions, int cluster_count, float *centers) {
    for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
        int particle_select = rand() % particle_count;
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<float>(centers, center_iter, dim_iter, dimensions) = (float)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    return 0;
}


/**
 * Smarter cluster select method
 */
int select_centerspp_serial(uchar *particle_data, int particle_count, int dimensions, int cluster_count, float *centers) {
    float *particle_distr = new float[particle_count];
    int distr_length = RAND_MAX;

    // Select the first center randomly
    int particle_select = rand() % particle_count;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        array_store<float>(centers, 0, dim_iter, dimensions) = (float)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
    }

    // Select the other n - 1 centers
    for (int center_iter = 1; center_iter < cluster_count; center_iter++) {
        // calc probability for each particle as center 
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for
            // find closest center based on previously calculated centers
            float min_dist = std::numeric_limits<float>::max();
            for (int center_dist_iter = 0; center_dist_iter < center_iter; center_dist_iter++) {
                float dist = compute_distance<uchar, float>(particle_data, particle_iter, centers, center_dist_iter, dimensions);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
            particle_distr[particle_iter] = min_dist;
        }

        // Prefix sum
        for (int particle_iter = 1; particle_iter < particle_count; particle_iter++) { // TODO parallel prefix sum
            particle_distr[particle_iter] += particle_distr[particle_iter-1];
        }

        // Select new center
        float distr_ratio = particle_distr[particle_count-1] / distr_length;
        int select_pos = rand() % distr_length;
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {  // TODO search - make binary, parallelize (note break, tricky)
            if (particle_distr[particle_iter] > select_pos * distr_ratio) {
                particle_select = particle_iter;
                break;
            }
        }

        // Copy new center to array
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<float>(centers, center_iter, dim_iter, dimensions) = (float)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    delete [] particle_distr;
    return 0;
}
