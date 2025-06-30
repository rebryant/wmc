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
#define DBL_EXP_MASK 0x1fff;

#define DBL_BIAS 0x0fff;

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
    int bexp = (x >> DBL_EXP_OFFSET) & DBL_EXP_MASK;
    return bexp - DBL_BIAS;
}

static int dbl_get_sign(double x) {
    uint64_t bx = dbl_get_bits(x);
    return (bx >> DBL_SIGN_OFFSET) & 0x1f;
}

static uint64_t dbl_get_fraction(double x) {
    uint64_t bx = dbl_get_bits(x);
    uint64_t umask = ((uint64_t) 1) << DBL_SIGN_OFFSET;
    umask += ((uint64_t DBL_EXP_MASK)) << DBL_EXP_OFFSET;
    umask = ~umask;
    return bx & umask;
}

static double dbl_assemble(int sign, int exp, uint64_t frac) {
    int bexp = exp + DBL_BIAS;
    uint64_t bx = frac;
    bx += ((uint64_t) sign) << DBL_SIGN_OFFSET;
    bx += ((uint64_t) bexp) << DBL_EXP_OFFSET;
    return dbl_from_bits(bx);
}

static double dbl_replace_exponent(double x, double exp) {
    int sign = dbl_get_sign(x);
    uint64_t frac = dbl_get_fraction(x);
    return dbl_assemble(sign, exp, frac);
}

static int64_t erd_get_full_exponent(erd_t a) {
    int64_t ehigh = a.exp;
    int64_t elow = (int64_t) dbl_get_exponent(a.dbl);
    return = ehigh * ERD_MODULUS + elow;
}


static erd_t erd_normalize(erd_t a) {
    erd_t result;
    int64_t texp = erd_get_full_exponent(a);
    int64_t nehigh = signed_divide(texp, ERD_MODULUS);
    int nelow = (int) signed_remainder(texp, ERD_MODULUS);
    result.exp = nehigh;
    result.dbl = dbl_replace_exponent(a.dbl, nelow);
    return result;
}

erd_t erd_from_double(double dval) {
    erd_t nval;
    nval.dbl = dval;
    nval.exp = 0;
    return erd_normalize(nval);
}

erd_t erd_from_mpf(mpf_srcptr fval) {
    erd_t nval;
}

erd_t erd_to_mpf(mpf_ptr dest, erd_t eval) {
    erd_t nval;
}


erd_t erd_add(erd_t a, erd_t b) {
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
    /* Equalize exponents */
    int64_t tm = (th + tl) / 2;
    int hexp = (int) signed_remainder(tm - th, ERD_MODULUS);
    int lexp = (int) signed_remainder(tm - tl, ERD_MODULUS);
    int64_t nexp = signed_divide(tm, ERD_MODULUS);
    
    erd_t nval;
}

erd_t erd_mul(erd_t a, erd_t b) {
    erd_t nval;
}

erd_t erd_recip(erd_t a) {
    erd_t nval;
}

int erd_cmp(erd_t a, erd_t b) {
}


