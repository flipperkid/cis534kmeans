#include "kmeans.h"
#include "kmeans_omp.h"

int kmeans_omp(uchar *particle_data, double *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
    reset_timer( 4 );

    // assign centers
    int *cluster_sizes = new int[cluster_count];
    start_timer( 4 );
    if (select_centerspp_omp(particle_data, particle_count, dimensions, cluster_count, centers) != 0 ) {
        printf("Error selecting centers");
        return -1;
    }
    stop_timer( 4 );

    // iterate until clusters converge
    start_timer( 3 );
    int iterations = 0;
    bool assignment_change = true;
    while(assignment_change) {
        assignment_change = false;

        // step 1
        start_timer( 1 );
        #pragma omp parallel for
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for
            // find closest center
            int cluster_assignment = (int)assignments[particle_iter];
            double min_dist;
            if(cluster_assignment < cluster_count) {
                min_dist = compute_distance<uchar, double>(particle_data, particle_iter, centers, cluster_assignment, dimensions);
            }
            else {
                min_dist = std::numeric_limits<double>::max();
            }
            for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
                if( center_iter != cluster_assignment ) {
                    double dist = compute_distance<uchar, double>(particle_data, particle_iter, centers, center_iter, dimensions);
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
                    array_store<double>(centers, cluster_iter, dim_iter, dimensions) = 0;
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
                    array_store<double>(centers, cluster_assignment, dim_iter, dimensions) += (double)array_load<uchar>(particle_data, particle_iter, dim_iter, dimensions);

                }
            }

            // divide by number of assigned particles to get mean
            #pragma omp parallel for
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_store<double>(centers, cluster_iter, dim_iter, dimensions) =
                        floor(array_load<double>(centers, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5f);
                }
            }
            #pragma omp barrier
            stop_timer( 2 );
        }

        iterations++;
    }
    stop_timer( 3 );
    printf("Step1: %f seconds, %d iterations, %d clusters\n", get_time_elapsed( 1 ), iterations, cluster_count);
    printf("Step2: %f seconds, %d iterations, %d clusters\n", get_time_elapsed( 2 ), iterations, cluster_count);
    printf("Total: %f seconds, %d iterations, %d clusters\n", get_time_elapsed( 3 ), iterations, cluster_count);

    // release memory
    delete [] cluster_sizes;
    return 0;
}


/**
 * Smarter cluster select method
 */
int select_centerspp_omp(uchar *particle_data, int particle_count, int dimensions, int cluster_count, double *centers) {
    double *particle_distr = new double[particle_count];
    int distr_length = RAND_MAX;

    reset_timer(4);
    reset_timer(5);
    reset_timer(6);


    // Select the first center randomly
    int particle_select = rand() % particle_count;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        array_store<double>(centers, 0, dim_iter, dimensions) = (double)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
    }

    // Select the other n - 1 centers
    for (int center_iter = 1; center_iter < cluster_count; center_iter++) {
        // calc probability for each particle as center 
        start_timer(4);
        #pragma omp parallel for
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) { // TODO parallel for
            // find closest center based on previously calculated centers
            double min_dist = std::numeric_limits<double>::max();
            for (int center_dist_iter = 0; center_dist_iter < center_iter; center_dist_iter++) {
                double dist = compute_distance<uchar, double>(particle_data, particle_iter, centers, center_dist_iter, dimensions);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
            particle_distr[particle_iter] = min_dist;
        }
        stop_timer(4);

        start_timer(5);
        // Prefix sum
        for (int particle_iter = 1; particle_iter < particle_count; particle_iter++) { // TODO parallel prefix sum
            particle_distr[particle_iter] += particle_distr[particle_iter-1];
        }
        /*
        for (int ii = 0; ii < ceil(log(particle_count) / log(2.0)); ii++) {
            #pragma omp parallel for
            for (int jj = 0; jj < particle_count; jj++) {
                if (jj >= (1 << ii)) {
                    particle_distr[jj] += particle_distr[jj - (unsigned int)(1 << ii)];
                }
            }
        }
        */
        stop_timer(5);

        start_timer(6);
        // Select new center
        double distr_ratio = particle_distr[particle_count-1] / distr_length;
        int select_pos = rand() % distr_length;
        #pragma omp parallel for
        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {  // TODO search - make binary, parallelize (note break, tricky)
            if (particle_distr[particle_iter] > select_pos * distr_ratio) {
                particle_select = particle_iter;
                break;
            }
        }
        stop_timer(6);

        // Copy new center to array
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<double>(centers, center_iter, dim_iter, dimensions) = (double)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    delete [] particle_distr;

    printf("Cluster Step1: %f seconds\n", get_time_elapsed( 4 ));
    printf("Cluster Step2: %f seconds\n", get_time_elapsed( 5 ));
    printf("Cluster Step3: %f seconds\n", get_time_elapsed( 6 ));

    return 0;
}

