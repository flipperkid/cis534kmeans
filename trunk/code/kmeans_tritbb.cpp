#include "kmeans_tritbb.h"

int kmeans_tritbb(uchar *particle_data, double *&centers, int particle_count, int dimensions, int cluster_count, uchar *assignments, int grainSize, int thread_count ) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
    reset_timer( 4 );
    reset_timer( 5 );
    start_timer( 5 );
//    srand( time(NULL) );
    task_scheduler_init init(thread_count);

    double *lbound = new double[particle_count*cluster_count];
    for( int lbound_iter = 0; lbound_iter < particle_count*cluster_count; lbound_iter++ ) {
        lbound[lbound_iter] = 0;
    }
    double *ubound = new double[particle_count];
    bool *ubound_invalid = new bool[particle_count];
    for( int ubound_iter = 0; ubound_iter < particle_count; ubound_iter++ ) {
        ubound[ubound_iter] = std::numeric_limits<double>::max();
        ubound_invalid[ubound_iter] = false;
    }
    double **center_dists = new double*[cluster_count];
    for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
        center_dists[center_iter] = new double[cluster_count];
    }
    double *short_dist = new double[cluster_count];
    double *means_new = new double[cluster_count*dimensions];
    double *center_shift = new double[cluster_count];

    // assign centers
    int *cluster_sizes = new int[cluster_count];
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
        step1_tritbb( particle_data, centers, particle_count, cluster_count, assignments, dimensions, assignment_change, lbound, ubound, ubound_invalid, center_dists, short_dist, grainSize );
        stop_timer( 1 );
    
        // step 2
        if (assignment_change) {
            start_timer( 2 );
            step2_tritbb( particle_data, centers, means_new, center_shift, cluster_sizes, particle_count, cluster_count, assignments, dimensions, lbound, ubound, ubound_invalid, grainSize );
            stop_timer( 2 );
        }
        iterations++;
    }
    stop_timer( 3 );
    // release memory
    delete [] cluster_sizes;
    delete [] means_new;
    delete [] center_shift;
    delete [] lbound;
    delete [] ubound;
    delete [] ubound_invalid;
    delete [] center_dists;
    delete [] short_dist;
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
    return 0;
}

void step1_tritbb( uchar *particle_data, double *centers, int particle_count, int cluster_count, uchar *assignments, int dimensions, bool &assignment_change, 
                double *lbound, double *ubound, bool *ubound_invalid, double **center_dists, double *short_dist, int grainSize ) {
    // Precompute distance between centers
    for( int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        short_dist[cluster_iter] = std::numeric_limits<double>::max();
    }
    for( int cluster_iter1 = 0; cluster_iter1 < cluster_count; cluster_iter1++) {
        for( int cluster_iter2 = cluster_iter1+1; cluster_iter2 < cluster_count; cluster_iter2++) {
            double dist = 0.25f*compute_distance<double, double>(centers, cluster_iter1, centers, cluster_iter2, dimensions);
            center_dists[cluster_iter1][cluster_iter2] = dist;
            center_dists[cluster_iter2][cluster_iter1] = dist;
            if(dist < short_dist[cluster_iter1]) {
                short_dist[cluster_iter1] = dist;
            }
            if(dist < short_dist[cluster_iter2]) {
                short_dist[cluster_iter2] = dist;
            }
        }
    }
    parallel_for(blocked_range<size_t>(0, particle_count, grainSize), 
        KmeansStep1tri(particle_data, centers, ubound, ubound_invalid, lbound, center_dists, short_dist, assignments, cluster_count, dimensions, &assignment_change));
}

void step2_tritbb( uchar *particle_data, double *&centers, double *&means_new, double *center_shift, int *cluster_sizes, int particle_count, int cluster_count, uchar *assignments, int dimensions, 
                double *lbound, double *ubound, bool *ubound_invalid, int grainSize ) {
    // initialize centers to 0
    for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<double>(means_new, cluster_iter, dim_iter, dimensions) = 0;
        }
    }
    for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        cluster_sizes[cluster_iter] = 0;
    }
    // sum over assigned particle positions
    KmeansStep2a step2a_obj(particle_data, means_new, assignments, cluster_count, dimensions);
    parallel_reduce(blocked_range<size_t>(0, particle_count, grainSize), step2a_obj );
    // divide by number of assigned particles to get mean
    parallel_for(blocked_range<size_t>(0, cluster_count, grainSize), KmeansStep2b(means_new, step2a_obj.sumCenters, step2a_obj.sumClusterSizes, dimensions));

    parallel_for(blocked_range<size_t>(0, cluster_count, grainSize), KmeansStep2c(means_new, centers, center_shift, dimensions));
    parallel_for(blocked_range<size_t>(0, particle_count, grainSize), KmeansStep2d(center_shift, lbound, ubound, ubound_invalid, assignments, cluster_count));

    double *temp = centers;
    centers = means_new;
    means_new = temp;
}
