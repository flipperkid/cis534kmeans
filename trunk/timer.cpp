#include "timer.h"

double secs = 0.0f;
double start;
bool running = false;

double get_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    double tv_sec = double(ts.tv_sec);
    double tv_nsec = double(ts.tv_nsec);
    double seconds = tv_sec + (tv_nsec / double(1000000000.0));
    return seconds;
}

void reset_timer( ) { 
    running = false;
    start = get_seconds();
    secs = 0.0f; 
}

void start_timer( ) { 
    start = time(NULL); 
    running = true;
}

void stop_timer( ) {
    if( running ) {
        running = false;
        secs += get_seconds() - start;
    }
}

double get_time_elapsed( ) {
    if( running ) {
        return secs + get_seconds() - start;
    }
    return secs;
}

