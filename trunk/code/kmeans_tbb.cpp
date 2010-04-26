#include "kmeans_tbb.h"

int kmeans_tbb(uchar *particle_data, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize) {
    reset_timer();
//    srand( time(NULL) );
    task_scheduler_init init;

    // assign centers
    double *centers = new double[cluster_count*dimensions];
    double *cluster_sizes = new double[cluster_count];
    if (select_centers_serial(particle_data, particle_count, dimensions, cluster_count, centers) != 0 ) {
        printf("Error selecting centers");
        return -1;
    }

    // iterate until clusters converge
    int iterations = 0;
    bool assignment_change = true;
    while(assignment_change) {
        assignment_change = false;
        
        // step 1
        start_timer();
        parallel_for(blocked_range<size_t>(0, particle_count, grainSize), KmeansStep1(particle_data, centers, cluster_count, assignments, dimensions, &assignment_change)); // TODO parallel for
        stop_timer();

        // step 2
        if (assignment_change) {

            // initialize centers to 0
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_access<double>(centers, cluster_iter, dim_iter, dimensions) = 0;
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
                    array_access<double>(centers, cluster_assignment, dim_iter, dimensions) += (double)array_access<uchar>(particle_data, particle_iter, dim_iter, dimensions);

                }
            }

            // divide by number of assigned particles to get mean
            parallel_for(blocked_range<size_t>(0, cluster_count, grainSize), KmeansStep2b(centers, cluster_sizes, dimensions));
            for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
                for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
                    array_access<double>(centers, cluster_iter, dim_iter, dimensions) = floor(array_access<double>(centers, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5);
                }
            }
        }
        iterations++;
    }
    printf("Distance Calc: %f seconds, %d iterations, %d grainsize\n", get_time_elapsed(), iterations, grainSize);

    // release memory
    delete [] centers;
    return 0;
}

