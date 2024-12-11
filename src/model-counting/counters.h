/*========================================================================
  Copyright (c) 2023 Randal E. Bryant, Carnegie Mellon University
  
  Permission is hereby granted, free of
  charge, to any person obtaining a copy of this software and
  associated documentation files (the "Software"), to deal in the
  Software without restriction, including without limitation the
  rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom
  the Software is furnished to do so, subject to the following
  conditions:
  
  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.
  
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
========================================================================*/


#pragma once

// Count all the interesting stuff

typedef enum { 
    COUNT_DATA_VARIABLES,
    COUNT_INPUT_CLAUSE,
    COUNT_SMOOTH_VARIABLES,
    COUNT_EDGES,
    COUNT_OPERATIONS,
    COUNT_NUM
} counter_t;

typedef enum { HISTO_SUMS, HISTO_NODE_PRODUCTS, HISTO_EDGE_PRODUCTS, HISTO_EDGE_SMOOTHS, HISTO_NUM } histogram_t;

typedef enum { TIME_SETUP, TIME_EVAL, TIME_NUM } runtimer_t;


/* Allow this headerfile to define C++ constructs if requested */
#ifdef __cplusplus
#define CPLUSPLUS
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

void incr_count(counter_t counter);
void incr_count_by(counter_t counter, int val);
void set_count(counter_t counter, int val);
void max_count(counter_t counter, int val);
int get_count(counter_t counter);
long get_long_count(counter_t counter);


void incr_timer(runtimer_t timer, double secs);
void reset_timer(runtimer_t timer);
double get_timer(runtimer_t timer);

void incr_histo(histogram_t h, int datum);
void reset_histo(histogram_t h);
int get_histo_min(histogram_t h);
int get_histo_max(histogram_t h);    
int get_histo_count(histogram_t h);
long get_histo_total(histogram_t h);
double get_histo_avg(histogram_t h);


#ifdef CPLUSPLUS
}
#endif


/* EOF */
    
