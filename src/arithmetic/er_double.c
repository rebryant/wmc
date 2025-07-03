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
#include <stdio.h>

#include "er_double.h"


/*
  Properties of double-precision representation
*/
#define DBL_EXP_OFFSET 52
#define DBL_SIGN_OFFSET 63
#define DBL_EXP_MASK 0x7ff
#define DBL_MAX_PREC 54
#define DBL_BIAS 0x3ff

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

/* Get exponent as signed integer */
static int dbl_get_exponent(double x) {
    uint64_t bx = dbl_get_bits(x);
    int bexp = (bx >> DBL_EXP_OFFSET) & DBL_EXP_MASK;
    return bexp - DBL_BIAS;
}

static int dbl_get_sign(double x) {
    uint64_t bx = dbl_get_bits(x);
    return (bx >> DBL_SIGN_OFFSET) & 0x1f;
}

static uint64_t dbl_get_fraction(double x) {
    uint64_t bx = dbl_get_bits(x);
    uint64_t umask = ((int64_t) 1 << DBL_EXP_OFFSET) - 1;
    return bx & umask;
}

/* Signed exponent too small */
static bool dbl_exponent_below(int exp) {
    return exp <= -((int) DBL_BIAS);
}

/* Signed exponent too large */
static bool dbl_exponent_above(int exp) {
    return exp >= (int) (DBL_EXP_MASK - DBL_BIAS);
}

static double dbl_assemble(int sign, int exp, uint64_t frac) {
    int bexp = exp + DBL_BIAS;
    uint64_t bx = frac;
    bx += ((uint64_t) bexp) << DBL_EXP_OFFSET;
    bx += ((uint64_t) sign) << DBL_SIGN_OFFSET;
    return dbl_from_bits(bx);
}

static double dbl_replace_exponent(double x, int exp) {
    int sign = dbl_get_sign(x);
    uint64_t frac = dbl_get_fraction(x);
    return dbl_assemble(sign, exp, frac);
}

static double dbl_infinity(int sign) {
    return dbl_assemble(sign, DBL_EXP_MASK - DBL_BIAS, 0);
}

static bool erd_is_zero(erd_t a) {
    return a.dbl == 0;
}

static erd_t erd_zero() {
    erd_t nval;
    nval.exp = 0;
    nval.dbl = 0;
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
    nval.dbl = dval;
    nval.exp = 0;
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
    if (erd_is_zero(a))
	return b;
    if (erd_is_zero(b))
	return a;
    erd_t nval;
    if (a.exp - b.exp > DBL_MAX_PREC)
	return a;
    if (b.exp - a.exp > DBL_MAX_PREC)
	return b;
    int ediff = (int) (a.exp - b.exp);
    double ad = dbl_replace_exponent(a.dbl, ediff);
    nval.dbl = ad + b.dbl;
    nval.exp = b.exp;
    return erd_normalize(nval);
}

erd_t erd_mul(erd_t a, erd_t b) {
#if 0
    /* Not needed */
    if (erd_is_zero(a) || erd_is_zero(b))
	return erd_zero();
#endif
    erd_t nval;
    nval.exp = a.exp + b.exp;
    nval.dbl = a.dbl * b.dbl;
    return erd_normalize(nval);
}

erd_t erd_recip(erd_t a) {
    if (erd_is_zero(a))
	return a;
    erd_t nval;
    nval.exp = -a.exp;
    nval.dbl = 1.0/a.dbl;
    return erd_normalize(nval);
}

int erd_cmp(erd_t a, erd_t b) {
    int sa = dbl_get_sign(a.dbl);
    int sb = dbl_get_sign(b.dbl);
    int za = erd_is_zero(a);
    int zb = erd_is_zero(b);
    if (za) {
	if (zb)
	    return 0;
	else if (sb)
	    return -1;
	return 1;
    }
    if (zb)
	return sa ? 1 : -1;
    int factor = 1;
    if (sa) {
	// a < 0
	if (sb)
	    // b < 0
	    factor = -1;
	else
	    // b > 0
	    return -1;
    } else {
	/* a > 0 */
	if (sb)
	    // b < 0
	    return 1;
    }
    if (a.exp > b.exp)
	return factor;
    else if (a.exp < b.exp)
	return -factor;
    else {
	if (a.dbl < b.dbl)
	    return -1;
	else if (a.dbl > b.dbl)
	    return 1;
	else
	    return 0;
    }
}


/* Debugging support */

static void show_double(double d) {
    int sign = dbl_get_sign(d);
    int exp = dbl_get_exponent(d);
    uint64_t frac = dbl_get_fraction(d);
    printf("Sign=%d, Exp=%d, Frac=0x%llx, Val=%.8f", sign, exp, (unsigned long long) frac, d);
}

static void show_erd(erd_t a) {
    int sign = dbl_get_sign(a.dbl);
    uint64_t frac = dbl_get_fraction(a.dbl);
    printf("Sign=%d, Exp=%lld, Frac=0x%llx", sign, (long long) a.exp,
	   (unsigned long long) frac);
}

