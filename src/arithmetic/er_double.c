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

#include "er_double.h"

/* Important constants */

/* 
   Range of exponents stored in double.
   Should be power of 2 between 64 and 512
*/
#define ER_MODULUS 512

/*
  Properties of double-precision representation
*/
#define DBL_EXP_OFFSET 52
#define DBL_SIGN_OFFSET 63
#define DBL_EXP_MASK 0x7ff
#define DBL_MAX_PREC 54

#define DBL_BIAS 0x3ff;

int get_sign(int64_t val) {
    return val < 0;
}

/* Integer arithmetic that preserves sign.  Assume den > 0 */
int64_t signed_divide(int64_t num, int64_t den) {
    int64_t div;
    if (get_sign(num))
	div = -((-num) / den);
    else
	div = num / den;
    return div;
}

/* Integer arithmetic that preserves sign.  Assume den > 0 */
int64_t signed_remainder(int64_t num, int64_t den) {
    int64_t rem;
    if (get_sign(num))
	rem = -((-num) % den);
    else
	rem = num % den;
    return rem;
}

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

static double dbl_assemble(int sign, int exp, uint64_t frac) {
    int bexp = exp + DBL_BIAS;
    uint64_t bx = frac;
    bx += ((uint64_t) bexp) << DBL_EXP_OFFSET;
    bx += ((uint64_t) sign) << DBL_SIGN_OFFSET;
    return dbl_from_bits(bx);
}

static double dbl_replace_exponent(double x, double exp) {
    int sign = dbl_get_sign(x);
    uint64_t frac = dbl_get_fraction(x);
    return dbl_assemble(sign, exp, frac);
}

static int64_t erd_get_full_exponent(erd_t a) {
    int64_t ehigh = a.exh;
    int64_t elow = (int64_t) dbl_get_exponent(a.dbl);
    return ehigh * ER_MODULUS + elow;
}


static erd_t erd_normalize(erd_t a) {
    erd_t nval;
    int64_t texp = erd_get_full_exponent(a);
    int64_t nehigh = signed_divide(texp, ER_MODULUS);
    int nelow = (int) signed_remainder(texp, ER_MODULUS);
    nval.exh = nehigh;
    nval.dbl = dbl_replace_exponent(a.dbl, nelow);
    return nval;
}

erd_t erd_from_double(double dval) {
    erd_t nval;
    nval.dbl = dval;
    nval.exh = 0;
    return erd_normalize(nval);
}

erd_t erd_from_mpf(mpf_srcptr fval) {
    erd_t nval;
    mpf_t mfval;
    mpf_init(mfval);
    mpf_set(mfval, fval);
    int64_t mpf_exp = mfval[0]._mp_exp;
    mfval[0]._mp_exp = 0;
    double d = mpf_get_d(mfval);
    mpf_exp += dbl_get_exponent(d);
    nval.exh = signed_divide(mpf_exp, ER_MODULUS);
    int nexp = (int) signed_remainder(mpf_exp, ER_MODULUS);
    nval.dbl  = dbl_replace_exponent(d, nexp);
    return erd_normalize(nval);
}

void erd_to_mpf(mpf_ptr dest, erd_t eval) {
    mpf_set_d(dest, eval.dbl);
    if (eval.exh < 0)
	mpf_div_2exp(dest, dest, -eval.exh * ER_MODULUS);
    else
	mpf_mul_2exp(dest, dest, eval.exh * ER_MODULUS);
    
}

erd_t erd_add(erd_t a, erd_t b) {
    erd_t nval;
    /* Sort by exponent */
    erd_t eh, el;
    int64_t th, tl;
    int64_t ta = erd_get_full_exponent(a);
    int64_t tb = erd_get_full_exponent(b);
    if (ta > tb) {
	eh = a; th = ta;
	el = b; tl = tb;
    } else {
	eh = b; th = tb;
	el = a; tl = ta;
    }

    if (th - tl > DBL_MAX_PREC)
	// Smaller value will not affect sum
	return eh;
    // Must equalize extended exponents
    int nexp = dbl_get_exponent(el.dbl) - (eh.exh - el.exh) * ER_MODULUS;
    double nl = dbl_replace_exponent(el.dbl, nexp);

    nval.dbl = eh.dbl + nl;
    nval.exh = eh.exh;
    return erd_normalize(nval);
}    

erd_t erd_mul(erd_t a, erd_t b) {
    erd_t nval;
    nval.exh = a.exh + b.exh;
    nval.dbl = a.dbl * b.dbl;
    return erd_normalize(nval);
}

erd_t erd_recip(erd_t a) {
    erd_t nval;
    nval.exh = -a.exh;
    nval.dbl = 1.0/a.dbl;
    return erd_normalize(nval);
}

int erd_cmp(erd_t a, erd_t b) {
    int sa = dbl_get_sign(a.dbl);
    int sb = dbl_get_sign(b.dbl);
    int factor = 1;
    if (sa) {
	// a < 0
	if (sb)
	    // b < 0
	    factor = -1;
	else
	    // b >= 0
	    return -1;
    } else {
	/* a >= 0 */
	if (sb)
	    // b < 0
	    return 1;
    }
    if (a.exh > b.exh)
	return factor;
    else if (a.exh < b.exh)
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


