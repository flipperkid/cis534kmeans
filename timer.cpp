#include "timer.h"

double secs = 0.0f;
time_t start;
bool running = false;

void reset_timer( ) { 
    running = false;
    start= time(NULL);
    secs = 0.0f; 
}

void start_timer( ) { 
    start = time(NULL); 
    running = true;
}

void stop_timer( ) {
    if( running ) {
        running = false;
        secs += difftime(time(NULL),start);
    }
}

double get_time_elapsed( ) {
    if( running ) {
        return secs + difftime(time(NULL), start);
    }
    return secs;
}

