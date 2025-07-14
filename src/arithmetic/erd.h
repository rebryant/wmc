/*========================================================================
  Copyright (c) 2025 Randal E. Bryant, Carnegie Mellon University
  
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

#include <stdbool.h>
#include <gmp.h>

#include "erd-header.h"

/* 
   Support two versions: 
   DEFAULT:     0.0 has exp = INT64_MIN
   ERDZ:        0.0 has exp = 0

*/
#define ERDZ 0

#if ERDZ
#define ZEXP 0
#else
#define ZEXP INT64_MIN
#endif

/* Max number of times fractions can be multiplied without overflowing exponent */
#define MAX_MUL 1000

/********************* Double **********************/

#define DBL_EXP_OFFSET 52
#define DBL_SIGN_OFFSET 63
#define DBL_EXP_MASK ((uint64_t) 0x7ff)
#define DBL_MAX_PREC 54
#define DBL_BIAS ((int64_t) 0x3ff)

static uint64_t dbl_get_bits(double x) {
    union {
	double   d;
	uint64_t b;
    } u;
    u.d = x;
    return u.b;
}

static double dbl_from_bits(uint64_t bx) {
    union {
	double   d;
	uint64_t b;
    } u;
    u.b = bx;
    return u.d;
}

/* Get exponent as unsigned integer */
static uint64_t dbl_get_biased_exponent(double x) {
    uint64_t bx = dbl_get_bits(x);
    return (bx >> DBL_EXP_OFFSET) & DBL_EXP_MASK;
}

/* Get exponent as signed integer */
static int64_t dbl_get_exponent(double x) {
    int64_t bexp = (int64_t) dbl_get_biased_exponent(x);
    return bexp - DBL_BIAS;
}

static uint64_t dbl_get_sign(double x) {
    uint64_t bx = dbl_get_bits(x);
    return (bx >> DBL_SIGN_OFFSET) & 0x1;
}

static uint64_t dbl_get_fraction(double x) {
    uint64_t bx = dbl_get_bits(x);
    uint64_t umask = ((int64_t) 1 << DBL_EXP_OFFSET) - 1;
    return bx & umask;
}

/* Signed exponent too small */
static bool dbl_exponent_below(int64_t exp) {
    return exp <= -(int64_t) DBL_BIAS;
}

/* Signed exponent too large */
static bool dbl_exponent_above(int64_t exp) {
    return exp >= (int64_t) DBL_EXP_MASK - DBL_BIAS;
}

static double dbl_assemble(uint64_t sign, int64_t exp, uint64_t frac) {
    int64_t bexp = exp + DBL_BIAS;
    uint64_t bx = frac;
    bx += bexp << DBL_EXP_OFFSET;
    bx += sign << DBL_SIGN_OFFSET;
    return dbl_from_bits(bx);
}

static double dbl_replace_exponent(double x, int64_t exp) {
    uint64_t bexp = (uint64_t) ((exp + DBL_BIAS) & DBL_EXP_MASK) << DBL_EXP_OFFSET;
    uint64_t bx = dbl_get_bits(x);
    uint64_t mask = ~(DBL_EXP_MASK << DBL_EXP_OFFSET);
    bx &= mask;
    bx += bexp;
    return dbl_from_bits(bx);
}

static double dbl_infinity(int sign) {
    return dbl_assemble(sign, DBL_EXP_MASK - DBL_BIAS, 0);
}


/********************* ERD *************************/

bool erd_is_zero(erd_t a) {
    return a.dbl == 0.0;
}

erd_t erd_zero() {
    erd_t nval;
    nval.dbl = 0.0;
    nval.exp = ZEXP;
    return nval;
}


static erd_t erd_normalize(erd_t a) {
    if (erd_is_zero(a))
	return erd_zero();
    erd_t nval;
    nval.exp = a.exp + dbl_get_exponent(a.dbl);
    nval.dbl = dbl_replace_exponent(a.dbl, 0);
    return nval;
}

erd_t erd_from_double(double dval) {
    erd_t nval;
#if ERDZ
    nval.exp = 0;
#else
    nval.exp = dval == 0 ? ZEXP : 0;
#endif
    nval.dbl = dval;
    return erd_normalize(nval);
}

erd_t erd_from_mpf(mpf_srcptr fval) {
    erd_t nval;
    long int exp;
    nval.dbl = mpf_get_d_2exp(&exp, fval);
    if (nval.dbl == 0)
	return erd_zero();
    nval.exp = (int64_t) exp;
    return erd_normalize(nval);
}

void erd_to_mpf(mpf_ptr dest, erd_t eval) {
    mpf_set_d(dest, eval.dbl);
#if !ERDZ
    if (erd_is_zero(eval))
	return;
#endif
    if (eval.exp < 0)
	mpf_div_2exp(dest, dest, -eval.exp);
    else if (eval.exp > 0)
	mpf_mul_2exp(dest, dest, eval.exp);
}

double erd_to_double(erd_t eval) {
    if (eval.dbl == 0)
	return 0.0;
    if (dbl_exponent_below(eval.exp))
	return 0.0;
    if (dbl_exponent_above(eval.exp)) {
	int sign = dbl_get_sign(eval.dbl);
	return dbl_infinity(sign);
    }
    return dbl_replace_exponent(eval.dbl, eval.exp);
}

erd_t erd_negate(erd_t a) {
    erd_t nval;
    if (erd_is_zero(a))
	return a;
    nval.exp = a.exp;
    nval.dbl = -a.dbl;
    return nval;
}

erd_t erd_add(erd_t a, erd_t b) {
#if ERDZ
    if (erd_is_zero(a))
	return b;
    if (erd_is_zero(b))
	return a;
#endif
    if (a.exp > b.exp + DBL_MAX_PREC)
	return a;
    if (b.exp > a.exp + DBL_MAX_PREC)
	return b;
    erd_t nval;
    int64_t ediff = a.exp - b.exp;
    double ad = dbl_replace_exponent(a.dbl, ediff);
    nval.dbl = ad + b.dbl;
    nval.exp = b.exp;
    return erd_normalize(nval);
}

static erd_t erd_quick_mul(erd_t a, erd_t b) {
    erd_t nval;
    nval.exp = a.exp + b.exp;
    nval.dbl = a.dbl * b.dbl;
    return nval;
}


erd_t erd_mul(erd_t a, erd_t b) {
    return erd_normalize(erd_quick_mul(a, b));
}

static erd_t erd_mul_seq_slow(erd_t *val, int len) {
    erd_t result = erd_from_double(1.0);
    int i;
    for (i = 0; i < len; i++)
	result = erd_mul(result, val[i]);
    return result;
}

static erd_t erd_mul_seq_x1(erd_t *val, int len) {
    erd_t result = len == 0 ? erd_from_double(1.0) : val[0];
    int i;
    int count = 1;
    for (i = 1; i < len; i++) {
	erd_t arg = val[i];
	result = erd_quick_mul(result, arg);
	if (++count > MAX_MUL) {
	    count = 0;
	    result = erd_normalize(result);
	}
    }
    return erd_normalize(result);
}

static erd_t erd_mul_seq_x4(erd_t *val, int len) {
    // Assume len >= 4
    erd_t prod[4];
    int i, j;
    for (j = 0; j < 4; j++) 
	prod[j] = val[j];
    int count = 0;
    for (i = 4; i <= len-4; i+= 4) {
	for (j = 0; j < 4; j++)
	    prod[j] = erd_quick_mul(prod[j], val[i+j]);
	if (++count > MAX_MUL) {
	    count = 0;
	    for (j = 0; j < 4; j++)
		prod[j] = erd_normalize(prod[j]);
	}
    }
    if (count * 4 > MAX_MUL) {
	for (j = 0; j < 4; j++)
	    prod[j] = erd_normalize(prod[j]);
    }

    prod[2] = erd_quick_mul(prod[2], prod[3]);
    prod[1] = erd_quick_mul(prod[1], prod[2]);
    prod[0] = erd_quick_mul(prod[0], prod[1]);
    for (; i < len; i++)
	prod[0] = erd_quick_mul(prod[0], val[i]);
    return erd_normalize(prod[0]);
}

/* Compute product of sequence of values */
erd_t erd_mul_seq(erd_t *val, int len) {
    if (len < 4)
	return erd_mul_seq_x1(val, len);
    return erd_mul_seq_x4(val, len);
}

erd_t erd_div(erd_t a, erd_t b) {
    erd_t nval;
    nval.dbl = a.dbl / b.dbl;
    nval.exp = a.exp - b.exp;
    return erd_normalize(nval);
}

