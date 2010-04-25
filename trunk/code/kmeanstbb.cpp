#include "kmeanstbb.h"

int kmeans_tbb(uchar *particle_data, int particle_count, int dimensions, int cluster_count, uchar *assignments) {
    reset_timer();
    start_timer();
//    srand( time(NULL) );
    task_scheduler_init init;

    // assign centers
    double *centers = new double[cluster_count*dimensions];
    double *cluster_sizes = new double[cluster_count];
    if (select_centerspp_serial(particle_data, particle_count, dimensions, cluster_count, centers) != 0 ) {
        printf("Error selecting centers");
        return -1;
    }
    stop_timer();
    printf("%f seconds elapsed selecting centers\n", get_time_elapsed());
    reset_timer();
    start_timer();

    // iterate until clusters converge
    int iterations = 0;
    bool assignment_change = true;
    while(assignment_change) {
        assignment_change = false;
        
        // step 1
        int grainSizeStep1 = 10000;
        parallel_for(blocked_range<size_t>(0, particle_count, grainSizeStep1), KmeansStep1(particle_data, centers, cluster_count, assignments, dimensions, &assignment_change)); // TODO parallel for

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
            int grainSizeStep2b = 2;
            parallel_for(blocked_range<size_t>(0, cluster_count, grainSizeStep2b), KmeansStep2b(centers, cluster_sizes, dimensions));
        }
        iterations++;
    }
    stop_timer();
    printf("%f seconds elapsed revising centers over %d iterations\n", get_time_elapsed(), iterations);

    // release memory
    delete [] centers;
    return 0;
}

