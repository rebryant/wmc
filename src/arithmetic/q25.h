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

/* Allow functions that make use of the Gnu Multiprecision (GMP) library */
#define ENABLE_GMP 1

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#if ENABLE_GMP
#include <gmp.h>
#endif

/* Allow this headerfile to define C++ constructs if requested */
#ifdef __cplusplus
#define CPLUSPLUS
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif


/* Representation of a number of form -1^(sign) * d * 2^p2 * 5 ^p5
   where:
       d is arbitrary integer, represented as set of digits
          with (RADIX = 10**k for some k)
       p2 and p5 encode positive or negative powers of 2.

   Values are assumed to be immutable.
   All externally visible values stored in canonical form:
       If invalid, then d=0, p2=0, p5=0
       If zero then not negative and d=0, p2=0, p5=0
       If nonzero then not divisible by power of 2 or 5
*/
typedef struct {
    bool valid : 1;       // Is this a valid number
    bool negative: 1;     // Is it negative
    bool infinite: 1;     // Is it + or - infinity
    unsigned dcount : 29; // How many digits does it have (must be at least 1)
    int32_t pwr2;         // Power of 2
    int32_t pwr5;         // Power of 5
    uint32_t digit[1];    // Sequence of digits, each between 0 and RADIX-1
} q25_t, *q25_ptr;

void q25_free(q25_ptr q);

/* Make a fresh copy of number */
q25_ptr q25_copy(q25_ptr q);

/* Convert numbers to q25 form */
q25_ptr q25_from_64(int64_t x);
q25_ptr q25_from_32(int32_t x);
q25_ptr q25_invalid();
q25_ptr q25_infinity(bool negative);


/* Convert from/to double-precision FP.  Assume IEEE representation */
q25_ptr q25_from_double(double x);
double q25_to_double(q25_ptr q);

/* Scale by powers of 2 & 5 */
q25_ptr q25_scale(q25_ptr q, int32_t p2, int32_t p5);
void q25_inplace_scale(q25_ptr q, int32_t p2, int32_t p5);

/* Negative value */
q25_ptr q25_negate(q25_ptr q);
void q25_inplace_negate(q25_ptr q);

/* Absolute value */
q25_ptr q25_abs(q25_ptr q);
void q25_inplace_abs(q25_ptr q);

/* 
   Compute reciprocal 
   Can only compute reciprocal when d == 1
   Otherwise invalid
*/
q25_ptr q25_recip(q25_ptr q);

/* Is it valid */
bool q25_is_valid(q25_ptr q);

/* Is it zero */
bool q25_is_zero(q25_ptr q);

/* Is it one */
bool q25_is_one(q25_ptr q);

/* Is it infinite */
bool q25_is_infinite(q25_ptr q, bool *negativep);

/* Is it negative */
bool q25_is_negative(q25_ptr q);

/* 
   Compare two numbers.  Return -1 (q1<q2), 0 (q1=q2), or +1 (q1>q2)
   Return -2 if either invalid, or comparing two infinities of the same sign
*/
int q25_compare(q25_ptr q1, q25_ptr q2);

/* Addition */
q25_ptr q25_add(q25_ptr q1, q25_ptr q2);

/* Compute 1-x */
q25_ptr q25_one_minus(q25_ptr q);

/* Multiplication */
q25_ptr q25_mul(q25_ptr q1, q25_ptr q2);

/* Get approx log10 of number */
int q25_magnitude(q25_ptr q);

/* Round to specified number of decimal digits */
q25_ptr q25_round(q25_ptr q, int digits);

/* Read from file */
q25_ptr q25_read(FILE *infile);

/* Write decimal representation to file */
void q25_write(q25_ptr q, FILE *outfile);

/* Read from string */
q25_ptr q25_from_string(const char *sq);

/* Generate dynamically allocated string.  Should free() when done */
char *q25_string(q25_ptr q);

/* Generate string representation of form D.DD....DeNNN */
char *q25_scientific_string(q25_ptr q);

/* Choose shorter of 2 sring representations */
char *q25_best_string(q25_ptr q);


/* Set max_digits to <= 0 to print entire representation */
/* Show value in terms of its representation */
void q25_show(q25_ptr q, FILE *outfile);

/* Try converting to int64_t.  Indicate success / failure */
/* Fails if number out of range, or nonintegral */
bool get_int64(q25_ptr q, int64_t *ip);

/* 
   Reset all instrumentation counters and set future level:
   0     None
   1     Estimate MPQ sizes
   2     Find exact MPQ sizes
 */
void q25_reset_counters(int level);
/* Count of number of non-trivial operations since reset */
long q25_operation_count();
/* Get the peak allocated bytes since reset, according to different size models */
double q25_peak_allocation_fp(bool is_mpf);
double q25_peak_allocation_q25();
double q25_peak_allocation_mpq();
/* Get the largest allocation */
double q25_max_allocation_q25();
double q25_max_allocation_mpq();


/* Stack-based memory management.  Call q25_enter() when enter context, q25_exit() when leave */
int q25_enter();
void q25_leave(int pos);
q25_ptr q25_mark(q25_ptr q);

/* Extensions make use of GMP */
#if ENABLE_GMP

/* 
   Convert from q25 to GMP rational.  
   Sets *ok to true if successful.
   Will fail if q is special value
*/
bool q25_to_mpq(mpq_ptr dest, q25_ptr q);
q25_ptr q25_from_mpq(mpq_ptr z);

/* Only fails for infinite and special values */
bool q25_to_mpf(mpf_ptr dest, q25_ptr q);
q25_ptr q25_from_mpf(mpf_ptr z);

/* Will return false if rounding disabled and not integer */
bool q25_to_mpz(mpz_ptr dest, q25_ptr q, bool round);
q25_ptr q25_from_mpz(mpz_ptr z);

#endif /* INCLUDE_GMP */

#ifdef CPLUSPLUS
}
#endif
