#include "timer.h"

#define MAX_TIMERS 10

double secs[MAX_TIMERS];
double start[MAX_TIMERS];
bool running[MAX_TIMERS];

double get_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    double tv_sec = double(ts.tv_sec);
    double tv_nsec = double(ts.tv_nsec);
    double seconds = tv_sec + (tv_nsec / double(1000000000.0));
    return seconds;
}

void reset_timer( int i ) { 
    running[i] = false;
    start[i] = get_seconds();
    secs[i] = 0.0f; 
}

void start_timer( int i ) { 
    start[i] = get_seconds(); 
    running[i] = true;
}

void stop_timer( int i ) {
    if( running[i] ) {
        running[i] = false;
        secs[i] += get_seconds() - start[i];
    }
}

double get_time_elapsed( int i) {
    if( running[i] ) {
        return secs[i] + get_seconds() - start[i];
    }
    return secs[i];
}

