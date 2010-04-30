#include "kmeans_tbb.h"

int kmeans_tbb(uchar *particle_data, double *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize, int thread_count) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
    reset_timer( 4 );
    reset_timer( 5 );
    start_timer( 5 );
//    srand( time(NULL) );
    task_scheduler_init init(thread_count);

    // assign centers
    start_timer( 4 );
    if (select_centerspp_tbb(particle_data, particle_count, dimensions, cluster_count, centers, grainSize) != 0 ) {
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
        parallel_for(blocked_range<size_t>(0, particle_count, grainSize), KmeansStep1(particle_data, centers, cluster_count, assignments, dimensions, &assignment_change)); // TODO parallel for
        stop_timer( 1 );

        // step 2
        start_timer( 2 );
        if (assignment_change) {

            // sum over assigned particle positions
            KmeansStep2a step2a_obj(particle_data, centers, assignments, cluster_count, dimensions);
            parallel_reduce(blocked_range<size_t>(0, particle_count, grainSize), step2a_obj );

            // divide by number of assigned particles to get mean
            parallel_for(blocked_range<size_t>(0, cluster_count, grainSize), KmeansStep2b(centers, step2a_obj.sumCenters, step2a_obj.sumClusterSizes, dimensions));
        }
        stop_timer( 2 );
        iterations++;
    }

    stop_timer( 3 );
    stop_timer( 5 );
    printf("Initialization: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, %d threads, %d grainsize\n", 
            get_time_elapsed( 4 ), particle_count, cluster_count, dimensions, iterations, thread_count, grainSize);
    printf("Step1: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, %d threads, %d grainsize\n", 
            get_time_elapsed( 1 ), particle_count, cluster_count, dimensions, iterations, thread_count, grainSize);
    printf("Step2: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, %d threads, %d grainsize\n", 
            get_time_elapsed( 2 ), particle_count, cluster_count, dimensions, iterations, thread_count, grainSize);
    printf("Convergence: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, %d threads, %d grainsize\n", 
            get_time_elapsed( 3 ), particle_count, cluster_count, dimensions, iterations, thread_count, grainSize);
    printf("Total: %f seconds, %d particles, %d clusters, %d dimensions, %d iterations, %d threads, %d grainsize\n", 
            get_time_elapsed( 5 ), particle_count, cluster_count, dimensions, iterations, thread_count, grainSize);
    
    // release memory
    return 0;
}

int select_centerspp_tbb(uchar *particle_data, int particle_count, int dimensions, int cluster_count, double *centers, int grainSize) {

    double *particle_distrIn = new double[particle_count];
    double *particle_distrOut = new double[particle_count];
    int distr_length = RAND_MAX;

    // Select the first center randomly
    int particle_select = rand() % particle_count;
    for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
        array_store<double>(centers, 0, dim_iter, dimensions) = (double)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
    }

    // Select the other n - 1 centers
    for (int center_iter = 1; center_iter < cluster_count; center_iter++) {
        // calc probability for each particle as center 
        parallel_for(blocked_range<size_t>(0, particle_count, grainSize), KmeansppA(particle_data, centers, particle_distrIn, center_iter, dimensions)); // TODO parallel for

        // Prefix sum
        KmeansppPrefixSum prefixSum( particle_distrIn, particle_distrOut );
        parallel_scan( blocked_range<size_t>(0, particle_count, grainSize), prefixSum );

        // Select new center
        double distr_ratio = particle_distrOut[particle_count-1] / distr_length;
        int select_pos = rand() % distr_length;

        for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {  // TODO search - make binary, parallelize (note break, tricky)
            if (particle_distrOut[particle_iter] > select_pos * distr_ratio) {
                particle_select = particle_iter;
                break;
            }
        }
        // Copy new center to array
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<double>(centers, center_iter, dim_iter, dimensions) = (double)array_load<uchar>(particle_data, particle_select, dim_iter, dimensions);
        }
    }
    delete [] particle_distrIn;
    delete [] particle_distrOut;
    return 0;
}
