#include "kmeans_tri.h"

int kmeans_tri(uchar *particle_data, float *&centers, int particle_count, int dimensions, int cluster_count, uchar *assignments) {
    reset_timer( 1 );
    reset_timer( 2 );
    reset_timer( 3 );
//    srand( time(NULL) );
    float *lbound = new float[particle_count*cluster_count];
    for( int lbound_iter = 0; lbound_iter < particle_count*cluster_count; lbound_iter++ ) {
        lbound[lbound_iter] = 0;
    }
    float *ubound = new float[particle_count];
    bool *ubound_invalid = new bool[particle_count];
    for( int ubound_iter = 0; ubound_iter < particle_count; ubound_iter++ ) {
        ubound[ubound_iter] = std::numeric_limits<float>::max();
        ubound_invalid[ubound_iter] = false;
    }
    float **center_dists = new float*[cluster_count];
    for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
        center_dists[center_iter] = new float[cluster_count];
    }
    float *short_dist = new float[cluster_count];
    float *means_new = new float[cluster_count*dimensions];
    float *center_shift = new float[cluster_count];

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
        step1_tri( particle_data, centers, particle_count, cluster_count, assignments, dimensions, assignment_change, lbound, ubound, ubound_invalid, center_dists, short_dist );
        stop_timer( 1 );
    
        // step 2
        if (assignment_change) {
            start_timer( 2 );
            step2_tri( particle_data, centers, means_new, center_shift, cluster_sizes, particle_count, cluster_count, assignments, dimensions, lbound, ubound, ubound_invalid );
            stop_timer( 2 );
        }
        iterations++;
    }
    stop_timer( 3 );

    printf("Step1: %f seconds, %d particles, %d clusters, %d iterations, dimensions, 1 threads\n", get_time_elapsed( 1 ), particle_count, cluster_count, iterations, dimensions);
    printf("Step2: %f seconds, %d particles, %d clusters, %d iterations, dimensions, 1 threads\n", get_time_elapsed( 2 ), particle_count, cluster_count, iterations, dimensions);
    printf("Total: %f seconds, %d particles, %d clusters, %d iterations, dimensions, 1 threads\n", get_time_elapsed( 3 ), particle_count, cluster_count, iterations, dimensions);

    // release memory
    delete [] cluster_sizes;
    delete [] means_new;
    delete [] center_shift;
    delete [] lbound;
    delete [] ubound;
    delete [] ubound_invalid;
    delete [] center_dists;
    delete [] short_dist;
    return 0;
}

void step1_tri( uchar *particle_data, float *centers, int particle_count, int cluster_count, uchar *assignments, int dimensions, bool &assignment_change, 
                float *lbound, float *ubound, bool *ubound_invalid, float **center_dists, float *short_dist ) {
    // Precompute distance between centers
    for( int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        short_dist[cluster_iter] = std::numeric_limits<float>::max();
    }
    for( int cluster_iter1 = 0; cluster_iter1 < cluster_count; cluster_iter1++) {
        for( int cluster_iter2 = cluster_iter1+1; cluster_iter2 < cluster_count; cluster_iter2++) {
            float dist = 0.25f*compute_distance(centers, cluster_iter1, centers, cluster_iter2, dimensions);
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
    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) { // TODO parallel for
        int cluster_assignment = (int)(assignments[particle_iter]);            
        if( cluster_assignment >= cluster_count || ubound[particle_iter] > short_dist[cluster_assignment] ) {
            // find closest center
            for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
                if( center_iter != cluster_assignment && ( ubound[particle_iter] > array_load(lbound, particle_iter, center_iter, cluster_count) ||
                    ubound[particle_iter] > center_dists[center_iter][cluster_assignment] ) ) {
                    float dist;
                    if( ubound_invalid[particle_iter] ) {
                        dist = compute_distance(particle_data, particle_iter, centers, cluster_assignment, dimensions);
                        ubound[particle_iter] = dist;
                        array_store(lbound, particle_iter, cluster_assignment, cluster_count) = dist;
                        ubound_invalid[particle_iter] = false;
                    }
                    else {
                        dist = ubound[particle_iter];
                    }
                    if( dist > array_load(lbound, particle_iter, center_iter, cluster_count) || dist > center_dists[center_iter][cluster_assignment] ) {
                        float dist_new = compute_distance(particle_data, particle_iter, centers, center_iter, dimensions);
                        array_store(lbound, particle_iter, center_iter, cluster_count) = dist_new;
                        if( dist_new < dist ) {
                            cluster_assignment = center_iter;
                            ubound[particle_iter] = dist_new;
                        }
                    }
                }
            }
            // assign to closest center
            if( cluster_assignment != (int)(assignments[particle_iter]) ) {
                assignments[particle_iter] = (uchar)cluster_assignment;
                assignment_change = true;
            }
        }
    }
}

void step2_tri( uchar *particle_data, float *&centers, float *&means_new, float *center_shift, int *cluster_sizes, int particle_count, int cluster_count, uchar *assignments, int dimensions, 
                float *lbound, float *ubound, bool *ubound_invalid ) {
    // initialize centers to 0
    for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<float>(means_new, cluster_iter, dim_iter, dimensions) = 0;
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
            array_store<float>(means_new, cluster_assignment, dim_iter, dimensions) += (float)array_load<uchar>(particle_data, particle_iter, dim_iter, dimensions);
        }
    }
    // divide by number of assigned particles to get mean
    for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        for (int dim_iter = 0; dim_iter < dimensions; dim_iter++) {
            array_store<float>(means_new, cluster_iter, dim_iter, dimensions) = floor(array_load<float>(means_new, cluster_iter, dim_iter, dimensions)/cluster_sizes[cluster_iter] + 0.5f);
        }
    }
    for (int center_iter = 0; center_iter < cluster_count; center_iter++) {
        center_shift[center_iter] = compute_distance(means_new, center_iter, centers, center_iter, dimensions);
    }
    for( int particle_iter = 0; particle_iter < particle_count; particle_iter++ ) {
        for( int center_iter = 0; center_iter < cluster_count; center_iter++ ) {
            array_store(lbound, particle_iter, center_iter, cluster_count) = _MAX<float>(0, array_load(lbound, particle_iter, center_iter, cluster_count) - 
                center_shift[center_iter] - 2*sqrt(array_load(lbound, particle_iter, center_iter, cluster_count))*sqrt(center_shift[center_iter]));
        }
        ubound[particle_iter] = ubound[particle_iter] + center_shift[(int)(assignments[particle_iter])] + 2*sqrt(ubound[particle_iter])*sqrt(center_shift[(int)(assignments[particle_iter])]);
        ubound_invalid[particle_iter] = true;
    }
    float *temp = centers;
    centers = means_new;
    means_new = temp;
}
