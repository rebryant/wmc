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


#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "q25.h"

/*
   Declarations.  The program assumes that computation is performed in decimal.
   with some number of decimal digits stored in each word
 */

/* 
   Number of decimal digits in word.
   Must fit into uint32_t
*/

#define DEBUG 0

#define Q25_DIGITS 9
#define Q25_RADIX (1000*1000*1000)

// Stress test
//#define Q25_DIGITS 3
//#define Q25_RADIX (1000)


/*
  Maintain working area for building digit representations.
  Have fixed number of words.  Use q25_t for meta information
  but separate extensible arrays for digits.
*/

/* Working area parameters */


/* How many numbers are in working area */
#define DCOUNT 3
/* How many digits are allocated in initial arrays */
#define INIT_DIGITS 100
/* Default ID for working area */
#define WID 0

bool initialized = false;
/* Per-number components */
static q25_t working_val[DCOUNT];
static uint32_t *digit_buffer[DCOUNT];
static unsigned digit_allocated[DCOUNT];
#if RATIONAL
static uint32_t *ddigit_buffer[DCOUNT];
static unsigned ddigit_allocated[DCOUNT];
#endif

/* Lookup table for powers */
static uint32_t power2[Q25_DIGITS+1];
static uint32_t power5[Q25_DIGITS+1];
static uint32_t power10[Q25_DIGITS+1];

/* 
   Static function prototypes.
   Use id to index the different numbers
*/

/* Put into canonical form */
static void q25_canonize(int id);
// Set working value to number x < RADIX
static void q25_set(int id, uint32_t x);
// Make sure enough digits in working space
static void q25_check(int id, unsigned dcount);
static void q25_show_internal(int id, FILE *outfile);

/* Count of operations */
static long operation_counter = 0;
/* Count of number of active numbers (allocated - freed) */
static long active_counter = 0;
static long peak_active_counter = 0;
/* Count of number of Q25 bytes allocated */
static double active_bytes_q25 = 0;
static double peak_active_bytes_q25 = 0;
static double max_bytes_q25 = 0;
/* Approx bytes if had MPQ representation */
static double active_bytes_mpq = 0;
static double peak_active_bytes_mpq = 0;
static double max_bytes_mpq = 0;


/* Factors used for computing MP sizes */
static double mpf_bytes = 40;
static double dbl_bytes = sizeof(double);
/* Computed as log_2(10^9) * 32 / 8 */
static double mpq_bytes_per_dcount = 3.737169106748283;
/* 1/8 */
static double mpq_bytes_per_p2 = 0.125;
/* log_2(5) / 8 */
static double mpq_bytes_per_p5 = 0.2902410118609203;


/* Support for stack-based memory management */
#define QSTACK_INIT_SIZE  100
static q25_ptr *qstack = NULL;
static int qstack_size = 0;
static int qstack_alloc_size = 0;

/* Static functions */

/**** Computing allocations ****/
static double allocation_q25(q25_ptr q) {
#if RATIONAL    
    return 4 * (5 + q->dcount + q->ddcount);
#else
    return 4 * (4 + q->dcount);
#endif
}

static double allocation_mpq(q25_ptr q) {
    double val = mpq_bytes_per_dcount * q->dcount;
#if RATIONAL
    val += mpq_bytes_per_dcount * q->ddcount;
#endif
    val += mpq_bytes_per_p2 * (q->pwr2 > 0 ? q->pwr2 : -q->pwr2);
    val += mpq_bytes_per_p5 * (q->pwr5 > 0 ? q->pwr5 : -q->pwr5);
    /* Round up to multiple of 8 */
    val = 8 * ceil(0.125 * val);
    /* Add overhead */
    val += 32;
    return val;
}

/* Update statistics when q newly allocated */
static void q25_register(q25_ptr q) {
#if METRIC
    active_counter += 1;
    if (active_counter > peak_active_counter)
	peak_active_counter = active_counter;
    double bytes = allocation_q25(q);
    if (bytes > max_bytes_q25)
	max_bytes_q25 = bytes;
    active_bytes_q25 += bytes;
    if (active_bytes_q25 > peak_active_bytes_q25)
	peak_active_bytes_q25 = active_bytes_q25;
    bytes = allocation_mpq(q);
    if (bytes > max_bytes_mpq)
	max_bytes_mpq = bytes;
    active_bytes_mpq += bytes;
    if (active_bytes_mpq > peak_active_bytes_mpq)
	peak_active_bytes_mpq = active_bytes_mpq;
#endif
}

/* Update statistics when q newly allocated */
static void q25_deregister(q25_ptr q) {
#if METRIC
    active_counter -= 1;
    active_bytes_q25 -= allocation_q25(q);
    active_bytes_mpq -= allocation_mpq(q);
#endif
}

/* Initialize data structures */
static void q25_init() {
    if (initialized)
	return;
    q25_reset_counters();
    initialized = true;
    int id;
    for (id = 0; id < DCOUNT; id++) {
	digit_allocated[id] = INIT_DIGITS;
	digit_buffer[id] = (uint32_t *) calloc(INIT_DIGITS, sizeof(uint32_t));
	digit_buffer[id][0] = 0;
	working_val[id].valid = true;
	working_val[id].infinite = false;
	working_val[id].pwr2 = 0;
	working_val[id].pwr5 = 0;
	working_val[id].dcount = 1;
#if RATIONAL
	ddigit_allocated[id] = INIT_DIGITS;
	ddigit_buffer[id] = (uint32_t *) calloc(INIT_DIGITS, sizeof(uint32_t));	
	ddigit_buffer[id][0] = 0;
	working_val[id].ddcount = 0;
#endif
    }
    int i;
    uint64_t p2 = 1;
    uint64_t p5 = 1;
    for (i = 0; i <= Q25_DIGITS; i++) {
	power10[i] = p2 * p5;
	power2[i] = p2;
	p2 *= 2;
	power5[i] = p5;
	p5 *= 5;
    }
    qstack = (q25_ptr *) calloc(QSTACK_INIT_SIZE, sizeof(q25_ptr));
    qstack_size = 0;
    qstack_alloc_size = QSTACK_INIT_SIZE;
}

// Setting working value to number x < RADIX
static void q25_set(int id, uint32_t x) {
    q25_init();
    working_val[id].valid = true;
    working_val[id].infinite = false;
    working_val[id].pwr2 = 0;
    working_val[id].pwr5 = 0;
    working_val[id].dcount = 1;
#if RATIONAL
    working_val[id].ddcount = 0;
#endif
    digit_buffer[id][0] = x;
    q25_canonize(id);
}

// Move value into working space
static void q25_work(int id, q25_ptr q) {
    q25_check(id, q->dcount);
    working_val[id].valid = q->valid;
    working_val[id].infinite = q->infinite;
    working_val[id].negative = q->negative;
    working_val[id].dcount = q->dcount;
    working_val[id].pwr2 = q->pwr2;
    working_val[id].pwr5 = q->pwr5;
    memcpy(digit_buffer[id], q->digit, working_val[id].dcount * sizeof(uint32_t));
#if RATIONAL
    working_val[id].ddcount = q->ddcount;
    if (q->ddcount > 0)
	memcpy(ddigit_buffer[id], q->digit+q->dcount, working_val[id].ddcount * sizeof(uint32_t));
#endif
}

// Make sure enough digits in working space
static void q25_check(int id, unsigned dcount) {
    q25_init();
    if (dcount <= digit_allocated[id])

	return;
    digit_allocated[id] *= 2;
    if (dcount > digit_allocated[id])
	digit_allocated[id] = dcount;
    digit_buffer[id] = (uint32_t *) realloc(digit_buffer[id], digit_allocated[id] * sizeof(uint32_t));
}

// Clear specified number of digits in workspace.  And set as length
static void q25_clear_digits(int id, unsigned len) {
    q25_check(id, len);
    memset(digit_buffer[WID], 0, len * sizeof(uint32_t));
    working_val[id].dcount = len;
}

#if RATIONAL
// Make sure enough digits in working space
static void q25_dcheck(int id, unsigned ddcount) {
    q25_init();
    if (ddcount <= ddigit_allocated[id])
	return;
    ddigit_allocated[id] *= 2;
    if (ddcount > ddigit_allocated[id])
	ddigit_allocated[id] = ddcount;
    ddigit_buffer[id] = (uint32_t *) realloc(ddigit_buffer[id], ddigit_allocated[id] * sizeof(uint32_t));
}

// Clear specified number of digits in workspace.  And set as length
static void q25_clear_ddigits(int id, unsigned len) {
    q25_dcheck(id, len);
    memset(ddigit_buffer[WID], 0, len * sizeof(uint32_t));
    working_val[id].ddcount = len;
}
#endif


// Divide by a number < RADIX
// Assume dividend is valid, finite, and nonzero, and divisor is nonzero
// Return remainder
static uint32_t q25_div_word(int id, uint32_t divisor) {
    if (divisor == 1)
	return 0;
    uint64_t upper = 0;
    int d;
    for (d = working_val[id].dcount-1; d >= 0; d--) {
	uint64_t dividend = (upper * Q25_RADIX) + digit_buffer[id][d];
	digit_buffer[id][d] = dividend/divisor;
	upper = dividend % divisor;
    }
    // See if upper digit set to 0
    if (working_val[id].dcount > 1 && digit_buffer[id][working_val[id].dcount-1] == 0)
	working_val[id].dcount--;
    return upper;
}

/* Take out multiples of n, where n = 2^p2 * 5^p5, and n <= RADIX */
static void old_q25_reduce_multiple(int id, uint32_t p2, uint32_t p5, uint32_t n) {
    uint32_t word;
    while ((word = digit_buffer[id][0])  % n == 0) {
	int pwr = 0;
	uint64_t scale = 1;
	uint64_t nscale = scale * n;
	while (nscale <= Q25_RADIX && Q25_RADIX % nscale == 0 && word % nscale == 0) {
	    pwr ++;
	    scale = nscale;
	    nscale*= n;
	}
	q25_div_word(id, scale);
	working_val[id].pwr2 += p2*pwr;
	working_val[id].pwr5 += p5*pwr;
    }
}

/* Take out multiples of n, where n = 2^p2 * 5^p5, and n <= RADIX */
static void q25_reduce_multiple(int id, uint32_t p2, uint32_t p5, uint32_t n) {
    uint32_t word;
    while ((word = digit_buffer[id][0])  % n == 0) {
	int pwr = 0;
	uint64_t scale = 1;
	uint64_t rradix = Q25_RADIX;
	uint64_t rword = word;
	// Try expanding to two words.  Allows extracting more powers of two
	if (working_val[id].dcount > 1) {
	    rradix *= Q25_RADIX;
	    rword += (uint64_t) Q25_RADIX * digit_buffer[id][1];
	}
	uint64_t nscale = scale * n;
	while (nscale <= Q25_RADIX && rradix % nscale == 0 && rword % nscale == 0) {
	    pwr ++;
	    scale = nscale;
	    nscale *= n;
	}
	q25_div_word(id, scale);
	working_val[id].pwr2 += p2*pwr;
	working_val[id].pwr5 += p5*pwr;
    }
}


/* Take out as many multiples of 10 as possible.  Assume nonzero */
static void q25_reduce10(int id) {
    // Get as many words as possible
    uint32_t wcount = 0;
    while (wcount < working_val[id].dcount && digit_buffer[id][wcount] == 0)
	wcount++;
    // Shift words down
    uint32_t idest = 0;
    uint32_t isrc = wcount;
    while (isrc < working_val[id].dcount) {
	digit_buffer[id][idest++] = digit_buffer[id][isrc++];
    }
    working_val[id].dcount -= wcount;
    working_val[id].pwr2 += Q25_DIGITS * wcount;
    working_val[id].pwr5 += Q25_DIGITS * wcount;
    // Do the final digits
    q25_reduce_multiple(id, 1, 1, 10);
}

// Take out powers of two
static void q25_reduce2(int id) {
    q25_reduce_multiple(id, 1, 0, 2);
}

// Take out powers of five
static void q25_reduce5(int id) {
    q25_reduce_multiple(id, 0, 1, 5);
}

/* Canonize working value */
static void q25_canonize(int id) {
    if (!working_val[id].valid) {
	working_val[id].infinite = false;
	working_val[id].negative = false;
	working_val[id].dcount = 1;
	digit_buffer[id][0] = 0;
	working_val[id].pwr2 = 0;
	working_val[id].pwr5 = 0;
    } else if (working_val[id].infinite) {
	working_val[id].dcount = 1;
	digit_buffer[id][0] = 0;
	working_val[id].pwr2 = 0;
	working_val[id].pwr5 = 0;
    } else {
	// Make sure have the right number of digits
	while (working_val[id].dcount > 1 && digit_buffer[id][working_val[id].dcount-1] == 0)
	    working_val[id].dcount--;
	if (working_val[id].dcount == 1 && digit_buffer[id][0] == 0) {
	    /* Canonize zero */
	    working_val[id].negative = false;
	    working_val[id].pwr2 = 0;
	    working_val[id].pwr5 = 0;
	} else {
	    // Diminish by powers of 10, 2, and 5
	    q25_reduce10(id);
	    q25_reduce2(id);
	    q25_reduce5(id);
	}
    }
}

// Convert the working version into a true q25_t
static q25_ptr q25_build(int id) {
    q25_canonize(id);
    size_t len = sizeof(q25_t) + (working_val[id].dcount - 1) * sizeof(uint32_t);
    q25_ptr result = (q25_ptr) malloc(len);
    if (result == NULL)
	return NULL;
    result->valid = working_val[id].valid;
    result->infinite = working_val[id].infinite;
    result->negative = working_val[id].negative;
    result->dcount = working_val[id].dcount;
    result->pwr2 = working_val[id].pwr2;
    result->pwr5 = working_val[id].pwr5;
    memcpy(result->digit, digit_buffer[id], working_val[id].dcount * sizeof(uint32_t));
#if RATIONAL
    result->ddcount = working_val[id].ddcount;
    if (working_val[id].ddcount > 0)
	memcpy(result->digit+working_val[id].dcount, digit_buffer[id], working_val[id].ddcount * sizeof(uint32_t));
#endif 
    q25_register(result);
    return result;
}

// Multiply by a number < RADIX
// Assume multiplier is nonzero
static void q25_mul_word(int id, uint32_t multiplier) {
#if DEBUG
    printf("  Multiplying by %u\n", multiplier);
#endif
    q25_check(id, working_val[id].dcount+1);
    if (multiplier == 1)
	return;
    uint64_t upper = 0;
    int d;
    for (d = 0 ; d < working_val[id].dcount; d++) {
	uint64_t ndigit = upper + (uint64_t) multiplier * digit_buffer[id][d];
	digit_buffer[id][d] = ndigit % Q25_RADIX;
	upper = ndigit / Q25_RADIX;
    }
    // See if upper digit set to 0
    if (upper > 0) {
	digit_buffer[id][d] = upper;
	working_val[id].dcount++;
    }
}

// Scale number by power of 2, 5, or 10
static void q25_scale_digits(int id, bool p2, int pwr) {
    int p;
    if (p2)
	working_val[id].pwr2 -= pwr;
    else
	working_val[id].pwr5 -= pwr;
    uint32_t multiplier = p2 ? power2[Q25_DIGITS] : power5[Q25_DIGITS];
    while (pwr > Q25_DIGITS) {
	q25_mul_word(id, multiplier);
	pwr -= Q25_DIGITS;
    }
    multiplier = p2 ? power2[pwr] : power5[pwr];
    q25_mul_word(id, multiplier);
}

/* 
   Compare two working numbers.
   Must have already been scaled so that both numbers have same values for pwr2 & pwr5
   Return -1 (q1<q2), 0 (q1=q2), or +1 (q1>q2)
*/
static int q25_compare_working_magnitude(int id1, int id2) {
    if (working_val[id1].dcount < working_val[id2].dcount)
	return -1;
    if (working_val[id1].dcount > working_val[id2].dcount)
	return 1;
    int d;
    for (d = working_val[id1].dcount-1; d >= 0; d--) {
	if (digit_buffer[id1][d] < digit_buffer[id2][d])
	    return -1;
	if (digit_buffer[id1][d] > digit_buffer[id2][d])
	    return 1;
    }
    return 0;
}

/* How many decimal digits are in representation? */
static int q25_length10(int id) {
    if (!working_val[id].valid || working_val[id].infinite)
	return -1;
    int n10 = (working_val[id].dcount-1) * Q25_DIGITS;
    uint32_t word = digit_buffer[id][working_val[id].dcount-1];
    while (word > 0) {
	n10++;
	word = word/10;
    }
    return n10;
}

/* Get individual decimal digit */
static unsigned q25_get_digit10(int id, int index) {
    int digit = index / Q25_DIGITS;
    int offset = index % Q25_DIGITS;
    uint32_t power = power10[offset];
    if (digit < 0 || digit >= working_val[id].dcount)
	return 0;
    uint32_t word = digit_buffer[id][digit];
    return (word / power) % 10;
}

/* Truncate number by setting lower-order digits to 0 */
static void q25_clear_digit10s(int id, int count) {
    int digit = (count-1) / Q25_DIGITS;
    if (digit < 0 || digit >= working_val[id].dcount)
	return;
    int offset = (count-1) % Q25_DIGITS;
    if (offset > 0) {
	uint32_t power = power10[offset];
	digit_buffer[id][digit] = (digit_buffer[id][digit] / power) * power;
    }
    int i;
    for (i = digit-1; i >= 0; i--) {
	digit_buffer[id][i] = 0;
    }
}


/* Show internal representation */
static void q25_show_internal(int id, FILE *outfile) {
    if (!working_val[id].valid)
	fprintf(outfile, "INVALID");
    if (working_val[id].infinite)
	fprintf(outfile, "INFINITE");
    fprintf(outfile, "[%c,p2=%d,p5=%d", working_val[id].negative ? '-' : '+', working_val[id].pwr2, working_val[id].pwr5);
    int d;
    for (d = working_val[id].dcount-1; d >= 0; d--) {
	fprintf(outfile, "|");
	fprintf(outfile, "%u", digit_buffer[id][d]);
    }
#if RATIONAL    
    if (working_val[id].ddcount > 0) {
	int d;
	fprintf(outfile, "/");
	for (d = working_val[id].ddcount-1; d >= 0; d--) {
	    fprintf(outfile, "%u", ddigit_buffer[id][d]);
	}
    }
#endif
    fprintf(outfile, "]");
}

/**** Externally visible functions ****/

void q25_free(q25_ptr q) {
    if (q) {
	q25_deregister(q);
	q->valid = false;
    }
    free((void *) q);
}

/* Convert int64_t to q25 form */
#define I64_DIGITS 20
q25_ptr q25_from_64(int64_t x) {
    int wcount = (I64_DIGITS + Q25_DIGITS-1)/Q25_DIGITS;
    q25_check(WID, wcount);
    q25_set(WID, 0);
    if (x == 0)
	return q25_build(WID);
    if (x < 0) {
	working_val[WID].negative = true;
	x = -x;
    }
    working_val[WID].dcount = 0;
    while (x > 0) {
	digit_buffer[WID][working_val[WID].dcount++] = x % Q25_RADIX;
	x = x / Q25_RADIX;
    }
    return q25_build(WID);
}

/* Convert int32_t to q25 form */
#define I32_DIGITS 10
q25_ptr q25_from_32(int32_t x) {
    int wcount = (I32_DIGITS + Q25_DIGITS-1)/Q25_DIGITS;
    q25_check(WID, wcount);
    q25_set(WID, 0);
    if (x == 0)
	return q25_build(WID);
    if (x < 0) {
	working_val[WID].negative = true;
	x = -x;
    }
    working_val[WID].dcount = 0;
    while (x > 0) {
	digit_buffer[WID][working_val[WID].dcount++] = x % Q25_RADIX;
	x = x / Q25_RADIX;
    }
    return q25_build(WID);
}

q25_ptr q25_invalid() {
    q25_set(WID, 0);
    working_val[WID].valid = false;
    return q25_build(WID);
}

q25_ptr q25_infinity(bool negative) {
    q25_set(WID, 0);
    working_val[WID].infinite = true;
    if (negative)
	working_val[WID].negative = 1;
    return q25_build(WID);
}


q25_ptr q25_copy(q25_ptr q) {
    q25_work(WID, q);
    return q25_build(WID);
}


q25_ptr q25_scale(q25_ptr q, int32_t p2, int32_t p5) {
    q25_work(WID, q);
    int64_t np2 = (int64_t) p2 + working_val[WID].pwr2;
    if (np2 != (int64_t) (int32_t) np2) {
	return p2 > 0 ? q25_infinity(q->negative) : q25_invalid();
    }
    int64_t np5 = (int64_t) p5 + working_val[WID].pwr5;
    if (np5 != (int64_t) (int32_t) np5) {
	return p5 > 0 ? q25_infinity(q->negative) : q25_invalid();
    }
    working_val[WID].pwr2 = np2;
    working_val[WID].pwr5 = np5;
    return q25_build(WID);
}

void q25_inplace_scale(q25_ptr q, int32_t p2, int32_t p5) {
    q25_deregister(q);
    q25_work(WID, q);
    int64_t np2 = (int64_t) p2 + q->pwr2;
    if (np2 != (int64_t) (int32_t) np2) {
	/* This will mess up the allocation */
	if (p2 > 0)
	    q->infinite = true;
	else
	    q->valid = false;
	return;
    }
    int64_t np5 = (int64_t) p5 + q->pwr5;
    if (np5 != (int64_t) (int32_t) np5) {
	/* This will mess up the allocation */
	if (p5 > 0)
	    q->infinite = true;
	else
	    q->valid = false;
	return;
    }
    q->pwr2 = np2;
    q->pwr5 = np5;
    q25_register(q);
}

q25_ptr q25_negate(q25_ptr q) {
    q25_work(WID, q);
    working_val[WID].negative = !working_val[WID].negative;
    return q25_build(WID);
}

void q25_inplace_negate(q25_ptr q) {
    q->negative = !q->negative;
}

// Can only compute reciprocal when d == 1
// Otherwise invalid
q25_ptr q25_recip(q25_ptr q) {
    q25_set(WID, 1);
    if (!q->valid || q->dcount > 1 || q->digit[0] != 1) {
	working_val[WID].valid = false;
    } else {
	working_val[WID].pwr2 = -q->pwr2;
	working_val[WID].pwr5 = -q->pwr5;
    }
    return q25_build(WID);
}

bool q25_is_valid(q25_ptr q) {
    return q->valid;
}

bool q25_is_zero(q25_ptr q) {
    return q->valid && !q->infinite && q->dcount == 1 && q->digit[0] == 0;
}

bool q25_is_one(q25_ptr q) {
    return q->valid && !q->infinite && q->dcount == 1 && q->digit[0] == 1 
	&& q->pwr2 == 0 && q->pwr5 == 0;
}

bool q25_is_infinite(q25_ptr q, bool *negativep) {
    if (negativep)
	*negativep = q->negative == 1;
    return q->infinite;
}

bool q25_is_negative(q25_ptr q) {
    return q->negative;
}


/* 
   Compare two numbers.  Return -1 (q1<q2), 0 (q1=q2), or +1 (q1>q2)
   Return -2 if either invalid, or comparing same infinities
*/
int q25_compare(q25_ptr q1, q25_ptr q2) {
    if (q1->valid != q2->valid)
	return -2;
    if (q1->infinite) {
	if (q2->infinite) {
	    if (q1->negative == q2->negative)
		return -2;
	    else
		return q1->negative ? -1 : 1;
	} else
	    return q1->negative ? -1 : 1;
    } else if (q2->infinite) {
	return q2->negative ? 1 : -1;
    }
    if (q1->negative && !q2->negative)
	return -1;
    if (!q1->negative && q2->negative)
	return 1;
    if (q1->negative) {
	// Swap two, so that can compare digits
	q25_ptr qt = q1; q1 = q2; q2 = qt;
    }
    /* Must move arguments into working area so that can scale */
    q25_work(1, q1);
    q25_work(2, q2);
    int diff2 = working_val[1].pwr2 - working_val[2].pwr2;
    if (diff2 > 0) {
	q25_scale_digits(1, true, diff2);
    } else if (diff2 < 0) {
	q25_scale_digits(2, true, -diff2);
    }
    int diff5 = working_val[1].pwr5 - working_val[2].pwr5;
    if (diff5 > 0) {
	q25_scale_digits(1, false, diff5);
    } else if (diff5 < 0) {
	q25_scale_digits(2, false, -diff5);
    }
    return q25_compare_working_magnitude(1, 2);
}


q25_ptr q25_add(q25_ptr q1, q25_ptr q2) {
    if (!q1->valid || !q2->valid)
	return q25_invalid();
    if (q1->infinite) {
	if (q2->infinite)
	    return q1->negative == q2->negative ? q25_copy(q1) : q25_invalid();
	else return q25_copy(q1);
    } else if (q2->infinite)
	return q25_copy(q2);
    if (q25_is_zero(q1))
	return q25_copy(q2);
    if (q25_is_zero(q2))
	return q25_copy(q1);


    /* Must move arguments into working area.  Build result with id 0 */
    q25_work(1, q1);
    q25_work(2, q2);
#if DEBUG
    printf("  Working argument 1:");
    q25_show_internal(1, stdout);
    printf("\n  Working argument 2:");
    q25_show_internal(2, stdout);
    printf("\n");
#endif
    int diff2 = working_val[1].pwr2 - working_val[2].pwr2;
    if (diff2 > 0) {
	q25_scale_digits(1, true, diff2);
    } else if (diff2 < 0) {
	q25_scale_digits(2, true, -diff2);
    }
    int diff5 = working_val[1].pwr5 - working_val[2].pwr5;
    if (diff5 > 0) {
	q25_scale_digits(1, false, diff5);
    } else if (diff5 < 0) {
	q25_scale_digits(2, false, -diff5);
    }
#if DEBUG
    printf("  Scaled working argument 1:");
    q25_show_internal(1, stdout);
    printf("\n  Scaled working argument 2:");
    q25_show_internal(2, stdout);
    printf("\n");
#endif
    if (working_val[1].negative == working_val[2].negative) {
	unsigned ndcount = working_val[1].dcount;
	if (working_val[2].dcount > ndcount)
	    ndcount = working_val[2].dcount;
	ndcount += 1;
	q25_set(WID, 0);
	q25_check(WID, ndcount);
	working_val[WID].negative = working_val[1].negative;
	working_val[WID].pwr2 = working_val[1].pwr2;
	working_val[WID].pwr5 = working_val[1].pwr5;
	working_val[WID].dcount = ndcount;
	q25_clear_digits(WID, ndcount);
	uint32_t carry = 0;
	int d;
	for (d = 0; d < ndcount; d++) {
	    uint64_t digit = carry;
	    if (d < working_val[1].dcount)
		digit += digit_buffer[1][d];
	    if (d < working_val[2].dcount)
		digit += digit_buffer[2][d];
	    digit_buffer[WID][d] = digit % Q25_RADIX;
	    carry = digit / Q25_RADIX;
	}
    } else {
	int diff = q25_compare_working_magnitude(1, 2);
	q25_set(WID, 0);
	if (diff != 0) {
	    int tid = diff < 0 ? 2 : 1;
	    int bid = diff < 0 ? 1 : 2;
	    working_val[WID].negative = working_val[tid].negative;
	    working_val[WID].pwr2 = working_val[1].pwr2;
	    working_val[WID].pwr5 = working_val[1].pwr5;
	    working_val[WID].dcount = working_val[tid].dcount;
	    q25_check(WID, working_val[tid].dcount);
	    q25_clear_digits(WID, working_val[tid].dcount);
	    int32_t borrow = 0;
	    int d;
	    for (d = 0; d < working_val[tid].dcount; d++) {
		int64_t digit = -borrow;
		digit += digit_buffer[tid][d];
		if (d < working_val[bid].dcount)
		    digit -= digit_buffer[bid][d];
		if (digit < 0) {
		    digit += Q25_RADIX;
		    borrow = 1;
		} else 
		    borrow = 0;
		digit_buffer[WID][d] = digit;
	    }
	}
    }
#if DEBUG
    printf("  Working Sum:");
    q25_show_internal(WID, stdout);
    printf("\n");
#endif
    operation_counter++;
    return q25_build(WID);
}

q25_ptr q25_one_minus(q25_ptr q) {
    if (!q->valid)
	return q25_copy(q);
    if (q->infinite)
	return q25_negate(q);
    q25_ptr minus_one = q25_from_32(-1);
    q25_ptr sum = q25_add(q, minus_one);
    q25_inplace_negate(sum);
    q25_free(minus_one);
    return sum;
}

q25_ptr q25_mul(q25_ptr q1, q25_ptr q2) {
    if (!q1->valid || !q2->valid)
	return q25_invalid();
    if (q1->infinite) {
	if (q2->infinite)
	    return q25_infinity(q1->negative != q2->negative);
	else return q25_is_zero(q2) ? q25_invalid() : q25_copy(q1);
    } else if (q2->infinite)
	return q25_is_zero(q1) ? q25_invalid() : q25_copy(q2);
    if (q25_is_zero(q1))
	return q25_copy(q1);
    if (q25_is_zero(q2))
	return q25_copy(q2);
    if (q1->dcount == 1 && q1->digit[0] == 1) {
	q25_ptr result = q25_scale(q2, q1->pwr2, q1->pwr5);
	if (q1->negative)
	    result->negative = !result->negative;
	return result;
    }
    if (q2->dcount == 1 && q2->digit[0] == 1) {
	q25_ptr result = q25_scale(q1, q2->pwr2, q2->pwr5);
	if (q2->negative)
	    result->negative = !result->negative;
	return result;
    }
    q25_set(WID, 0);
    // Figure out sign
    working_val[WID].negative = (q1->negative != q2->negative);
    // Set powers
    working_val[WID].pwr2 = q1->pwr2 + q2->pwr2;
    working_val[WID].pwr5 = q1->pwr5 + q2->pwr5;
    // Clear out space for the product
    unsigned len = q1->dcount + q2->dcount + 1;
    q25_clear_digits(WID, len);
    // Make sure q1 is longer
    if (q1->dcount < q2->dcount) {
	q25_ptr qt = q1; q1 = q2; q2 = qt;
    }
    unsigned d1, d2;
    for (d2 = 0; d2 < q2->dcount; d2++) {
	uint64_t digit2 = q2->digit[d2];
	uint64_t carry = 0;
	for (d1 = 0; d1 < q1->dcount; d1++) {
	    uint64_t ndigit = q1->digit[d1] * digit2 + carry + digit_buffer[WID][d1+d2];
	    digit_buffer[WID][d1+d2] = ndigit % Q25_RADIX;
	    carry = ndigit / Q25_RADIX;
	}
	digit_buffer[WID][d1+d2] = carry;
    }
    operation_counter++;
    return q25_build(WID);
}

q25_ptr q25_read(FILE *infile) {
    /* Fill up digit buffer in reverse order */
    int d = 0;
    q25_check(1, d+1);
    digit_buffer[1][d] = 0;
    bool negative = false;
    int pwr10 = 0;
    bool got_point = false;
    /* Number of base 10 digits read */
    int n10 = 0;
    bool first = true;
    while (true) {
	int c = fgetc(infile);
	if (c == '-') {
	    if (first) {
		negative = true;
		first = false;
		continue;
	    }
	    else {
		ungetc(c, infile);
		break;
	    }
	} else if (c == '.') {
	    if (got_point) {
		ungetc(c, infile);
		break;
	    } else
		got_point = true;
	} else if (isdigit(c)) {
	    n10++;
	    if (got_point)
		pwr10--;
	    if (n10 > Q25_DIGITS && (n10-1) % Q25_DIGITS == 0) {
		// Time to start new word
		d++;
		q25_check(1, d+1);
		digit_buffer[1][d] = 0;
	    }
	    unsigned dig = c - '0';
	    digit_buffer[1][d] = 10 * digit_buffer[1][d] + dig;
	} else {
	    ungetc(c, infile);
	    break;
	}
	first = false;
    }
    bool valid = n10 > 0;
    if (valid) {
	// See if there's an exponent
	int c = fgetc(infile);
	if (c == 'e') {
	    // Deal with exponent
	    bool exp_negative = false;
	    int nexp = 0;
	    int exponent = 0;
	    bool exp_first = true;
	    while (true) {
		c = fgetc(infile);
		if (c == '-') {
		    if (exp_first)
			exp_negative = true;
		    else {
			ungetc(c, infile);
			valid = false;
			break;
		    }
		} else if (isdigit(c)) {
		    nexp++;
		    unsigned dig = c - '0';
		    exponent = 10 * exponent + dig;
		} else {
		    ungetc(c, infile);
		    break;
		}
		exp_first = false;
	    }
	    valid = valid && nexp > 0;
	    if (exp_negative)
		exponent = -exponent;
	    pwr10 += exponent;
	} else
	    ungetc(c, infile);
    }
    if (!valid) {
	q25_set(WID, 0);
	working_val[WID].valid = false;
	return q25_build(WID);
    }
    q25_set(WID, 0);
    working_val[WID].negative = negative;
    // Reverse the digits
    unsigned dcount = (n10 + Q25_DIGITS-1) / Q25_DIGITS;
    q25_check(WID, dcount);
    for (d = 0; d < dcount; d++) {
	digit_buffer[WID][d] = digit_buffer[1][dcount - 1 - d];
    }
    // Now could have a problem with the bottom word
    // Slide up to top and let the canonizer fix things
    unsigned extra_count = n10 % Q25_DIGITS;
    if (extra_count > 0) {
	unsigned scale = Q25_DIGITS-extra_count;
	unsigned multiplier = power10[scale];
	digit_buffer[WID][0] *= multiplier;
	pwr10 -= scale;
    }
    working_val[WID].dcount = dcount;
    working_val[WID].pwr2 = pwr10;
    working_val[WID].pwr5 = pwr10;
#if DEBUG
    printf("  Read value before canonizing: ");
    q25_show_internal(WID, stdout);
    printf("\n");
#endif
    return q25_build(WID);
}

q25_ptr q25_from_string(const char *sq) {
    int pos = 0;
    /* Fill up digit buffer in reverse order */
    int d = 0;
    q25_check(1, d+1);
    digit_buffer[1][d] = 0;
    bool negative = false;
    int pwr10 = 0;
    bool got_point = false;
    /* Number of base 10 digits read */
    int n10 = 0;
    bool first = true;
    while (true) {
	int c = sq[pos++];
	if (c == '-') {
	    if (first) {
		negative = true;
		first = false;
		continue;
	    }
	    else {
		--pos;
		break;
	    }
	} else if (c == '.') {
	    if (got_point) {
		--pos;
		break;
	    } else
		got_point = true;
	} else if (isdigit(c)) {
	    n10++;
	    if (got_point)
		pwr10--;
	    if (n10 > Q25_DIGITS && (n10-1) % Q25_DIGITS == 0) {
		// Time to start new word
		d++;
		q25_check(1, d+1);
		digit_buffer[1][d] = 0;
	    }
	    unsigned dig = c - '0';
	    digit_buffer[1][d] = 10 * digit_buffer[1][d] + dig;
	} else {
	    --pos;
	    break;
	}
	first = false;
    }
    bool valid = n10 > 0;
    if (valid) {
	// See if there's an exponent
	int c = sq[pos++];
	if (c == 'e') {
	    // Deal with exponent
	    bool exp_negative = false;
	    int nexp = 0;
	    int exponent = 0;
	    bool exp_first = true;
	    while (true) {
		c = sq[pos++];
		if (c == '-') {
		    if (exp_first)
			exp_negative = true;
		    else {
			--pos;
			valid = false;
			break;
		    }
		} else if (isdigit(c)) {
		    nexp++;
		    unsigned dig = c - '0';
		    exponent = 10 * exponent + dig;
		} else {
		    --pos;
		    break;
		}
		exp_first = false;
	    }
	    valid = valid && nexp > 0;
	    if (exp_negative)
		exponent = -exponent;
	    pwr10 += exponent;
	} else
	    --pos;
    }
    if (!valid) {
	q25_set(WID, 0);
	working_val[WID].valid = false;
	return q25_build(WID);
    }
    q25_set(WID, 0);
    working_val[WID].negative = negative;
    // Reverse the digits
    unsigned dcount = (n10 + Q25_DIGITS-1) / Q25_DIGITS;
    q25_check(WID, dcount);
    for (d = 0; d < dcount; d++) {
	digit_buffer[WID][d] = digit_buffer[1][dcount - 1 - d];
    }
    // Now could have a problem with the bottom word
    // Slide up to top and let the canonizer fix things
    unsigned extra_count = n10 % Q25_DIGITS;
    if (extra_count > 0) {
	unsigned scale = Q25_DIGITS-extra_count;
	unsigned multiplier = power10[scale];
	digit_buffer[WID][0] *= multiplier;
	pwr10 -= scale;
    }
    working_val[WID].dcount = dcount;
    working_val[WID].pwr2 = pwr10;
    working_val[WID].pwr5 = pwr10;
#if DEBUG
    printf("  Read value before canonizing: ");
    q25_show_internal(WID, stdout);
    printf("\n");
#endif
    return q25_build(WID);
}

/* Get approx log10 */
int q25_magnitude(q25_ptr q) {
    if (!q->valid)
	return INT_MAX;
    if (q->infinite)
	return q->negative ? INT_MIN : INT_MAX;
    if (q25_is_zero(q))
	return 0;
    q25_work(WID, q);
    int pwr10 = working_val[WID].pwr5;
    // Scale so that pwr2 = pwr5
    int diff = working_val[WID].pwr2 - working_val[WID].pwr5;
    if (diff > 0) {
	q25_scale_digits(WID, true, diff);
	pwr10 = working_val[WID].pwr5;
    } else if (diff < 0) {
	q25_scale_digits(WID, false, -diff);
	pwr10 = working_val[WID].pwr2;
    }
    int n10 = q25_length10(WID);
    return pwr10+n10-1;
}

q25_ptr q25_round(q25_ptr q, int digits) {
    if (!q->valid || q->infinite || q25_is_zero(q))
	return q25_copy(q);
    q25_work(WID, q);
    int pwr10 = working_val[WID].pwr5;
    // Scale so that pwr2 = pwr5
    int diff = working_val[WID].pwr2 - working_val[WID].pwr5;
    if (diff > 0) {
	q25_scale_digits(WID, true, diff);
	pwr10 = working_val[WID].pwr5;
    } else if (diff < 0) {
	q25_scale_digits(WID, false, -diff);
	pwr10 = working_val[WID].pwr2;
    }
    int n10 = q25_length10(WID);
    if (n10 <= digits)
	return q25_build(WID);
    bool roundup = false;
    int rounding_digit = q25_get_digit10(WID, n10-digits-1);
    //  printf("    Got pwr10 = %d, n10 = %d, rounding digit[%d] = %d\n", pwr10, n10, n10-digits-1, rounding_digit);
    if (n10-digits == 1) {
	/* Round to even happens only if only have one digit being rounded */
	if (rounding_digit == 5) {
	    int last_digit = q25_get_digit10(WID, 1);
	    //	    printf("      RTE: Checking last digit[1] = %d\n", last_digit);
	    roundup = last_digit % 2 == 1;
	} else
	    roundup = rounding_digit > 5;
    } else {
	roundup = rounding_digit >= 5;
    }
    /* Zero out the lower digits */
    q25_clear_digit10s(WID, n10-digits+1);
    q25_ptr interim = q25_build(WID);
#if 0
    /* DEBUG */
    char *si = q25_string(interim);
    printf("      Before rounding up: %s\n", si);
    free(si);
    /* DEBUG */
#endif

    if (roundup) {
	q25_ptr rval = q25_from_32(1);
	int scale = pwr10+n10-digits;
	q25_inplace_scale(rval, scale, scale);
	if (q25_is_negative(q))
	    q25_inplace_negate(rval);
	//	printf("      Round up by adding 10^%d\n", scale);
	q25_ptr result = q25_add(interim, rval);
	q25_free(interim); q25_free(rval);
	return result;
    } else
	return interim;
}


void q25_write(q25_ptr q, FILE *outfile) {
    if (!q->valid) {
	fprintf(outfile, "INVALID");
	return;
    }
    if (q->dcount == 1 && q->digit[0] == 0) {
	fprintf(outfile, "0");
	return;
    }    
    if (q->negative)
	fputc('-', outfile);

    if (q->infinite) {
	fprintf(outfile, "INF");
	return;
    }

    q25_work(WID, q);

    // Scale so that pwr2 = pwr5
    int diff = working_val[WID].pwr2 - working_val[WID].pwr5;
    if (diff > 0) {
	q25_scale_digits(WID, true, diff);
    } else if (diff < 0) {
	q25_scale_digits(WID, false, -diff);
    }
#if DEBUG
    printf("  Scaled for printing: ");
    q25_show_internal(WID, stdout);
    printf("\n");
#endif
    int n10 = q25_length10(WID);
    int p10 = working_val[WID].pwr2;
    int i;
    if (p10 >= 0) {
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    fputc(d, outfile);
	}
	while (p10-- > 0)
	    fputc('0', outfile);
    } else if (-p10 >= n10) {
	fputc('0', outfile);
	fputc('.', outfile);
	while (-p10 > n10) {
	    fputc('0', outfile);
	    p10++;
	}
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    fputc(d, outfile);
	}
    } else {
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    fputc(d, outfile);
	    if (i == -p10)
		fputc('.', outfile);
	}
    }
}

/* Support for dynamically allocated string representation of q25 value */
static int string_allocated = 0;
static int string_position = 0;
static char *string_buffer = NULL;

#define INIT_STRING 250

static void string_init() {
    string_buffer = malloc(INIT_STRING);
    string_allocated = INIT_STRING;
    string_position = 0;
}

static void string_check_length(int len) {
    if (len >= string_allocated) {
	string_allocated *=2;
	string_buffer = realloc(string_buffer, string_allocated);
    }
}

static void string_check_more(int count) {
    string_check_length(string_position+count);
}

static void string_append_char(char c) {
    string_check_more(1);
    string_buffer[string_position++] = c;
    string_buffer[string_position] = 0;
}

static void string_append_string(char *s) {
    int c;
    while ((c = *s++) != 0)
	string_append_char(c);
}

static void string_append_number(int val) {
    if (val == 0) {
	string_append_char('0');
	return;
    }
    if (val < 0) {
	string_append_char('-');
	val = -val;
    }
    int pos = Q25_DIGITS;
    while (power10[pos] > val)
	pos--;
    while (pos > 0) {
	int digit = val / power10[pos];
	char d = digit + '0';
	string_append_char(d);
	val = val % power10[pos];
	pos--;
    }
    char d = val + '0';
    string_append_char(d);
}

char *q25_string(q25_ptr q) {
    string_init();
    if (!q->valid) {
	string_append_string("INVALID");
	return string_buffer;
    }
    if (q->dcount == 1 && q->digit[0] == 0) {
	string_append_char('0');
	return string_buffer;
    }    

    if (q->negative)
	string_append_char('-');

    if (q->infinite) {
	string_append_string("INF");
	return string_buffer;
    }

    q25_work(WID, q);
    // Scale so that pwr2 = pwr5
    int diff = working_val[WID].pwr2 - working_val[WID].pwr5;
    if (diff > 0) {
	q25_scale_digits(WID, true, diff);
    } else if (diff < 0) {
	q25_scale_digits(WID, false, -diff);
    }
    int n10 = q25_length10(WID);
    int p10 = working_val[WID].pwr2;
    int i;
    if (p10 >= 0) {
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    string_append_char(d);
	}
	while (p10-- > 0)
	    string_append_char('0');
    } else if (-p10 >= n10) {
	string_append_char('0');
	string_append_char('.');
	while (-p10 > n10) {
	    string_append_char('0');
	    p10++;
	}
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    string_append_char(d);
	}
    } else {
	for (i = n10-1; i >= 0; i--) {
	    int d10 = q25_get_digit10(WID, i);
	    char d = '0' + d10;
	    string_append_char(d);
	    if (i == -p10)
		string_append_char('.');
	}
    }
    return string_buffer;
}

char *q25_scientific_string(q25_ptr q) {
    string_init();
    if (!q->valid) {
	string_append_string("INVALID");
	return string_buffer;
    }
    if (q->dcount == 1 && q->digit[0] == 0) {
	string_append_string("0.0");
	return string_buffer;
    }    

    if (q->negative)
	string_append_char('-');

    if (q->infinite) {
	string_append_string("INF");
	return string_buffer;
    }

    q25_work(WID, q);
    // Scale so that pwr2 = pwr5
    int diff = working_val[WID].pwr2 - working_val[WID].pwr5;
    if (diff > 0) {
	q25_scale_digits(WID, true, diff);
    } else if (diff < 0) {
	q25_scale_digits(WID, false, -diff);
    }
    int n10 = q25_length10(WID);
    int p10 = working_val[WID].pwr2 + n10 - 1;
    
    /* Leading digit */
    int d10 = q25_get_digit10(WID, n10-1);
    char d = '0' + d10;
    string_append_char(d);
    string_append_char('.');
    int i;
    if (n10 == 1) {
	string_append_char('0');
	return string_buffer;
    }
    for (i = n10-2; i >= 0; i--) {
	int d10 = q25_get_digit10(WID, i);
	char d = '0' + d10;
	string_append_char(d);
    }
    if (p10 != 0) {
	string_append_char('e');
	string_append_number(p10);
    }
    return string_buffer;
}

char *q25_best_string(q25_ptr q) {
    char *fs = q25_string(q);
    char *ss = q25_scientific_string(q);
    if (strlen(ss)+2 <= strlen(fs)) {
	free(fs);
	return ss;
    } else {
	free(ss);
	return fs;
    }
}



/* Show value in terms of its representation */
void q25_show(q25_ptr q, FILE *outfile) {
    q25_work(WID, q);
    q25_show_internal(WID, outfile);
}

/* Try converting to int64_t.  Indicate success / failure */
bool get_int64(q25_ptr q, int64_t *ip) {
    if (!q->valid || q->pwr2 < 0 || q->pwr5 < 0)
	return false;
    if (q->negative) {
	q25_ptr qmin = q25_from_64(INT64_MIN);
	if (q25_compare(q, qmin) < 0)
	    return false;
    } else {
	q25_ptr qmax = q25_from_64(INT64_MAX);
	if (q25_compare(q, qmax) > 0)
	    return false;
    }
    int64_t val = 0;
    int d;
    if (q->negative) {
	for (d = q->dcount-1; d >= 0; d--) {
	    val = val * Q25_RADIX + -q->digit[d];
	}
    } else {
	for (d = q->dcount-1; d >= 0; d--) {
	    val = val * Q25_RADIX + q->digit[d];
	}
    }
    int i;
    for (i = 0; i < q->pwr2; i++)
	val *= 2;
    for (i = 0; i < q->pwr5; i++)
	val *= 5;
    *ip = val;
    return true;
}

q25_ptr q25_from_double(double x) {
    union {
	double dx;
	uint64_t bx;
    } u;
    /* Get bit representation of x */
    u.dx = x;
    uint64_t bits = u.bx;
    unsigned sign = bits >> 63;
    int  biased_exp = (bits >> 52) & 0x7FF;
    int exp = biased_exp - (int) 0x3FF;
    int64_t frac = bits & 0xFFFFFFFFFFFFFL;
    if (biased_exp == 0) {
	/* Denormalized */
	exp++;
#if DEBUG
	printf("x = %.20f.  sign = %d, biased_exp = %u, frac = %ld.  Denormalized.  Exponent = %d\n",
	       x, sign, biased_exp, (long) frac, exp);
#endif
    } else if (biased_exp == 0x7FF) {
	/* Special */
	if (frac == 0) {
#if DEBUG
	    printf("x = %.20f.  sign = %d, biased_exp = %u, Infinity\n",
		   x, sign, biased_exp);
#endif
	    return q25_infinity(sign == 1);
	} else {
#if DEBUG
	    printf("x = %.20f.  sign = %d, biased_exp = %u, Special\n",
		   x, sign, biased_exp);
#endif
	    return q25_invalid();
	}
    } else {
	/* Normalized */
	frac += 1L << 52;
#if DEBUG
	printf("x = %.20f.  sign = %d, biased_exp = %d, frac = %ld.  Normalized.  Exponent = %d\n",
	       x, sign, biased_exp, (long) frac, exp);
#endif
    }
    /* Shift binary point to the right */
    exp -=  52;
    if (sign)
	frac = -frac;
    q25_ptr ifrac = q25_from_64(frac);
    ifrac->pwr2 += exp;
    {
	char *sfrac = q25_string(ifrac);
#if DEBUG
	printf("Q25 representation: %s\n", sfrac);
#endif
	free(sfrac);
    }
    return ifrac;
}

static double bits2double(uint64_t bits) {
    union {
	double dx;
	uint64_t bx;
    } u;
    u.bx = bits;
    return u.dx;
}    

static uint64_t double2bits(double x) {
    union {
	double dx;
	uint64_t bx;
    } u;
    u.dx = x;
    return u.bx;
}    

static double build_infinity(bool neg) {
    uint64_t bits = 0x7FF0000000000000L;
    if (neg)
	bits += (1L << 63);
    return bits2double(bits);
}

static double build_nan() {
    return bits2double(0x7FFfffffffffffffL);
}

static bool double_is_infinite(double x, bool *negativep) {
    uint64_t bits = double2bits(x);
    unsigned sign = bits >> 63;
    int  biased_exp = (bits >> 52) & 0x7FF;
    int64_t frac = bits & 0xFFFFFFFFFFFFFL;
    if (negativep)
	*negativep = (sign == 1);
    return biased_exp == 0x7FF &&  frac == 0;
}

static bool double_is_nan(double x) {
    uint64_t bits = double2bits(x);
    int  biased_exp = (bits >> 52) & 0x7FF;
    int64_t frac = bits & 0xFFFFFFFFFFFFFL;
    return biased_exp == 0x7FF &&  frac != 0;
}

double q25_to_double(q25_ptr q) {
    double x;
    bool negative;
    if (!q25_is_valid(q))
	return build_nan();
    if (q25_is_infinite(q, &negative))
	return build_infinity(negative);
    char *sq = q25_string(q);
    if (sscanf(sq, "%lf", &x) != 1) 
	return build_nan();
    free(sq);
    return x;
}

void q25_reset_counters() {
    operation_counter = 0;
    active_counter = 0;
    peak_active_counter = 0;
    active_bytes_q25 = 0;
    peak_active_bytes_q25 = 0;
    active_bytes_mpq = 0;
    peak_active_bytes_mpq = 0;
    max_bytes_q25 = 0;
    max_bytes_mpq = 0;
}

long q25_operation_count() {
    return operation_counter;
}

double q25_peak_allocation_fp(bool is_mpf) {
    return peak_active_counter * (is_mpf ? mpf_bytes : dbl_bytes);
}

double q25_peak_allocation_q25() {
    return peak_active_bytes_q25;
}

double q25_peak_allocation_mpq() {
    return peak_active_bytes_mpq;
}

double q25_max_allocation_q25() {
    return max_bytes_q25;
}

double q25_max_allocation_mpq() {
    return max_bytes_mpq;
}


/* Stack management */
int q25_enter() {
    return qstack_size;
}
void q25_leave(int pos) {
    while (qstack_size > pos)
	q25_free(qstack[--qstack_size]);
}

q25_ptr q25_mark(q25_ptr q) {
    q25_init();
    if (qstack_size >= qstack_alloc_size) {
	qstack_alloc_size *= 2;
	qstack = (q25_ptr *) realloc(qstack, qstack_alloc_size * sizeof(q25_ptr));
    }
    qstack[qstack_size++] = q;
    return q;
}

/*********************** GMP Code ********************************/
#if ENABLE_GMP

/* Table of powers of 5.  Entry i contains 5**(2**i) */
static mpz_t p5_table[8 * sizeof(uint32_t)];
/* How many entries are valid */
static uint32_t p5_count = 0;

static void generate_p5_entry(uint32_t i) {
    if (i < p5_count)
	return;
    /* First entry in table is 5 */
    if (p5_count == 0) {
	mpz_init(p5_table[0]);
	mpz_set_ui(p5_table[0], 5);
	p5_count++;
    }
    /* Fill with successive squares */
    while (p5_count <= i) {
	uint32_t fromi = p5_count-1;
	uint32_t toi = p5_count;
	p5_count++;
	mpz_init(p5_table[toi]);
	mpz_mul(p5_table[toi], p5_table[fromi], p5_table[fromi]);
    }
}

static void mpz_pow5(mpz_ptr z, uint32_t a) {
    mpz_init(z);
    mpz_set_ui(z, 1);
    uint32_t i = 0;
    while (a > 0) {
	if (a & 0x1) {
	    generate_p5_entry(i);
	    mpz_mul(z, z, p5_table[i]);
	}
	a >>= 1;
	i++;
    }
}

bool q25_to_mpq(mpq_ptr dest, q25_ptr q) {
    if (!q25_is_valid(q))
	return false;
    if (q25_is_infinite(q, NULL))
	return false;
    if (q25_is_zero(q)) {
	mpq_init(dest);
	return true;
    }
    mpz_t num;
    mpz_t den;
    mpz_init(num);
    mpz_init(den);
    mpz_set_ui(num, q->digit[q->dcount-1]);
    mpz_set_ui(den, 1);
    if (q->dcount > 1) {
	mpz_t radix;
	mpz_init(radix);
	mpz_set_ui(radix, Q25_RADIX);
	int d;
	for (d = q->dcount-2; d >= 0; d--) {
	    mpz_mul(num, num, radix);
	    mpz_add_ui(num, num, q->digit[d]);
	}
	mpz_clear(radix);
    }

    /* Scale by powers of 2 & 5 */
    if (q->pwr2 > 0)
	mpz_mul_2exp(num, num, q->pwr2);
    else if (q->pwr2 < 0)
	mpz_mul_2exp(den, den, -q->pwr2);

    if (q->pwr5 > 0) {
	mpz_t scale;
	mpz_pow5(scale, q->pwr5);
	mpz_mul(num, num, scale);
	mpz_clear(scale);
    } else if (q->pwr5 < 0) {
	mpz_t scale;
	mpz_pow5(scale, -q->pwr5);
	mpz_mul(den, den, scale);
	mpz_clear(scale);
    }

    /* Assemble result */
    mpq_init(dest);
    mpq_set_num(dest, num);
    mpq_set_den(dest, den);
    mpz_clear(num);
    mpz_clear(den);
    if (q25_is_negative(q))
	mpq_neg(dest, dest);

    return true;
}

q25_ptr q25_from_mpq(mpq_ptr z) {
    bool is_negative = false;
    switch (mpq_sgn(z)) {
    case -1:
	is_negative = true;
	break;
    case 0:
	return q25_from_32(0);
	break;
    case +1:
    default:
	break;
    }
    mpz_t num;
    mpz_init(num);
    mpq_get_num(num, z);
    mpz_abs(num, num);

    mpz_t two;
    mpz_init(two);
    mpz_set_ui(two, 2);
    mpz_t five;
    mpz_init(five);
    mpz_set_ui(five, 5);

    int p2 = mpz_remove(num, num, two);
    int p5 = mpz_remove(num, num, five);
    mpz_t den;
    mpz_init(den);
    mpq_get_den(den, z);

    int dp2 = mpz_remove(den, den, two);
    p2 -= dp2;
    int dp5 = mpz_remove(den, den, five);
    p5 -= dp5;
    if (mpz_cmp_ui(den, 1) != 0) {
	printf("Denominator not unit\n");
	/* Can't represent this number as a q25 */
	mpz_clears(num, den, two, five, NULL);
	return q25_invalid();
    }
    char *decimal = mpz_get_str(NULL, 10, num);
    q25_ptr result = q25_from_string(decimal);
    q25_inplace_scale(result, p2, p5);
    if (is_negative) 
	q25_inplace_negate(result);
    mpz_clears(num, den, two, five, NULL);
    free(decimal);
    return result;
}

/* Only fails for infinite and special values */
bool q25_to_mpf(mpf_ptr dest, q25_ptr q) {
    mpq_t mz;
    mpq_init(mz);
    if (!q25_to_mpq(mz, q)) {
	mpq_clear(mz);
	return false;
    }
    mpf_set_q(dest, mz);
    mpq_clear(mz);
    return true;
}

q25_ptr q25_from_mpf(mpf_ptr z) {
    mpq_t mz;
    mpq_init(mz);
    mpq_set_f(mz, z);
    q25_ptr result = q25_from_mpq(mz);
    mpq_clear(mz);
    return result;
}

/* Will return false if not integer when rounding disabled */
bool q25_to_mpz(mpz_ptr dest, q25_ptr q, bool round) {
    if (!round && (q->pwr2 < 0 || q->pwr5 < 0))
	return false;
    mpq_t mz;
    mpq_init(mz);
    if (!q25_to_mpq(mz, q)) {
	mpq_clear(mz);
	return false;
    }
    mpz_set_q(dest, mz);
    mpq_clear(mz);
    return true;
}

q25_ptr q25_from_mpz(mpz_ptr z) {
    mpq_t mz;
    mpq_init(mz);
    mpq_set_z(mz, z);
    q25_ptr result = q25_from_mpq(mz);
    mpq_clear(mz);
    return result;
}


#endif /* ENABLE_GMP */
