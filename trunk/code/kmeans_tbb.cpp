#include "kmeans_tbb.h"

int kmeans_tbb(uchar *particle_data, float *centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize, int thread_count) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
//    srand( time(NULL) );
    task_scheduler_init init(thread_count);

    // assign centers
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
    printf("Step1: %f seconds, %d particles, %d clusters, %d iterations, %d dimensions, %d threads, %d grainsize\n", 
            get_time_elapsed( 1 ), particle_count, cluster_count, iterations, dimensions, thread_count, grainSize);
    printf("Step2: %f seconds, %d particles, %d clusters, %d iterations, %d dimensions, %d threads, %d grainsize\n", 
            get_time_elapsed( 2 ), particle_count, cluster_count, iterations, dimensions, thread_count, grainSize);
    printf("Total: %f seconds, %d particles, %d clusters, %d iterations, %d dimensions, %d threads, %d grainsize\n", 
            get_time_elapsed( 3 ), particle_count, cluster_count, iterations, dimensions, thread_count, grainSize);
    
    // release memory
    return 0;
}

