#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <ctime>

double get_seconds();
void reset_timer( int i );
void start_timer( int i );
void stop_timer( int i );
double get_time_elapsed( int i);

#endif

