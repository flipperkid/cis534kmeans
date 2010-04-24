#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double get_seconds();
void reset_timer( );
void start_timer( );
void stop_timer( );
double get_time_elapsed( );

#endif

