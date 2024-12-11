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


#include <limits.h>
#include <stdbool.h>
#include "counters.h"
#include "report.h"

typedef struct {
    int min;
    int max;
    int count;
    long total;
} histo_info_t;

static long int counters[COUNT_NUM];
static double timers[TIME_NUM];
static histo_info_t histograms[HISTO_NUM];

static bool initialized = false;

static void test_init() {
    if (initialized)
	return;
    for (int c = 0; c < COUNT_NUM; c++)
	counters[c] = 0;
    for (int t = 0; t < TIME_NUM; t++)
	timers[t] = 0.0;
    for (int h = 0; h < HISTO_NUM; h++) {
	reset_histo(h);
    }
    initialized = true;
}

static bool counter_ok(counter_t counter) {
    test_init();
    bool ok = counter >= 0 && counter < COUNT_NUM;
    if (!ok)
	err(false, "Invalid counter number %d\n", (int) counter);
    return ok;
}

void incr_count_by(counter_t counter, int val) {
    if (!counter_ok(counter))
	return;
    counters[counter] += val;
}

void incr_count(counter_t counter) {
    incr_count_by(counter, 1);
}

void set_count(counter_t counter, int val) {
    if (!counter_ok(counter))
	return;
    counters[counter] = val;
}

void max_count(counter_t counter, int val) {
    if (!counter_ok(counter))
	return;
    if (counters[counter] < val)
	counters[counter] = val;
}


long get_long_count(counter_t counter) {
    if (!counter_ok(counter))
	return -1;
    return counters[counter];

}

int get_count(counter_t counter) {
    return (int) get_long_count(counter);
}

static bool timer_ok(runtimer_t timer) {
    test_init();
    bool ok = timer >= 0 && timer < TIME_NUM;
    if (!ok)
	err(false, "Invalid timer number %d\n", (int) timer);
    return ok;
}

void reset_timer(runtimer_t timer) {
    timers[timer] = 0;
}

void incr_timer(runtimer_t timer, double secs) {
    if (!timer_ok(timer))
	return;
    timers[timer] += secs;
}

double get_timer(runtimer_t timer) {
    if (!timer_ok(timer))
	return -1;
    return timers[timer];

}

static bool histo_ok(histogram_t histo) {
    test_init();
    bool ok = histo >= 0 && histo < HISTO_NUM;
    if (!ok)
	err(false, "Invalid histo number %d\n", (int) histo);
    return ok;
}

void reset_histo(histogram_t histo) {
    histograms[histo].min = INT_MAX;
    histograms[histo].max = INT_MIN;
    histograms[histo].count = 0;
    histograms[histo].total = 0;
}

void incr_histo(histogram_t histo, int datum) {
    if (!histo_ok(histo))
	return;
    histograms[histo].count++;
    histograms[histo].total += datum;
    if (datum < histograms[histo].min)
	histograms[histo].min = datum;
    if (datum > histograms[histo].max)
	histograms[histo].max = datum;
}

int get_histo_min(histogram_t histo) {
    if (!histo_ok(histo))
	return INT_MAX;
    return histograms[histo].min;
}

int get_histo_max(histogram_t histo) {
    if (!histo_ok(histo))
	return INT_MAX;
    return histograms[histo].max;
}

int get_histo_count(histogram_t histo) {
    if (!histo_ok(histo))
	return INT_MAX;
    return histograms[histo].count;
}

double get_histo_avg(histogram_t histo) {
    if (!histo_ok(histo))
	return 0.0;
    if (histograms[histo].count == 0)
	return 0;
    return (double) histograms[histo].total / histograms[histo].count;
}

long get_histo_total(histogram_t histo) {
    if (!histo_ok(histo))
	return 0;
    return histograms[histo].total;
}
