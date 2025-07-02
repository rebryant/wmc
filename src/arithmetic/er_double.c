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

#include "er_double.h"

/* Simple debugging */
#ifndef EDEBUG
#define EDEBUG 0
#endif


#if EDEBUG
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Forward Declarations */
static void erd_fail();
static const char *mpf_string(mpf_srcptr val);
static bool erd_compare_mpf(erd_t val, mpf_srcptr mval);
static void show_double(double d);
static void show_erd(erd_t val);
static bool erd_mpf_check(const char *context, erd_t val);
static void erd_to_mpf_with_init(mpf_ptr dest, erd_t eval);
#endif /* EDEBUG */

/* Important constants */

/* 
   Range of exponents stored in double.
   Should be power of 2 between 64 and 512
*/
#ifdef SMALL_MODULUS
#define ER_MODULUS 64
#else
#define ER_MODULUS 256
#endif

/*
  Properties of double-precision representation
*/
#define DBL_EXP_OFFSET 52
#define DBL_SIGN_OFFSET 63
#define DBL_EXP_MASK 0x7ff
#define DBL_MAX_PREC 54

#define DBL_BIAS 0x3ff

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

/* Signed exponent too small */
static bool dbl_exponent_below(int exp) {
    return exp <= -((int) DBL_BIAS);
}

/* Signed exponent too large */
static bool dbl_exponent_above(int exp) {
    return exp >= (int) (DBL_EXP_MASK - DBL_BIAS);
}

static bool dbl_bad_exponent(double x) {
    int exp = dbl_get_exponent(x);
    return dbl_exponent_below(exp) || dbl_exponent_above(exp);
}


static double dbl_assemble(int sign, int exp, uint64_t frac) {
    int bexp = exp + DBL_BIAS;
    uint64_t bx = frac;
    bx += ((uint64_t) bexp) << DBL_EXP_OFFSET;
    bx += ((uint64_t) sign) << DBL_SIGN_OFFSET;
    return dbl_from_bits(bx);
}

static double dbl_infinity(int sign) {
    return dbl_assemble(sign, DBL_EXP_MASK - DBL_BIAS, 0);
}

static double dbl_zero() {
    return dbl_from_bits(0);
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
    long int exp;
    double d = mpf_get_d_2exp(&exp, fval);
    int64_t nexp = (int64_t) exp + dbl_get_exponent(d);
    nval.exh = signed_divide(nexp, ER_MODULUS);
    int dexp = signed_remainder(nexp, ER_MODULUS);
    nval.dbl = dbl_replace_exponent(d, dexp);
    return erd_normalize(nval);
}

void erd_to_mpf(mpf_ptr dest, erd_t eval) {
    mpf_set_d(dest, eval.dbl);
    if (eval.exh < 0)
	mpf_div_2exp(dest, dest, -eval.exh * ER_MODULUS);
    else
	mpf_mul_2exp(dest, dest, eval.exh * ER_MODULUS);
}

double erd_to_double(erd_t eval) {
    int64_t exp = erd_get_full_exponent(eval);
    if (dbl_exponent_below(exp))
	return dbl_zero();
    if (dbl_exponent_above(exp)) {
	int sign = dbl_get_sign(eval.dbl);
	return dbl_infinity(sign);
    }
    return dbl_replace_exponent(eval.dbl, exp);
}

erd_t erd_negate(erd_t a) {
    erd_t nval;
    nval.exh = a.exh;
    nval.dbl = -a.dbl;
    return nval;
}

erd_t erd_add(erd_t a, erd_t b) {
#if EDEBUG
    erd_mpf_check("erd_add a", a);
    erd_mpf_check("erd_add b", b);
    mpf_t ma, mb;
    erd_to_mpf_with_init(ma, a);
    erd_to_mpf_with_init(mb, b);
    double da = erd_to_double(a);
    double db = erd_to_double(b);
#endif
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
    nval = erd_normalize(nval);
#if EDEBUG
    if (!dbl_bad_exponent(da) && !dbl_bad_exponent(db)) {
	double dsum = da + db;
	double dnval = erd_to_double(nval);
	if (dsum != dnval) {
	    fprintf(stderr, "Double sum mismatch:\n");
	    fprintf(stderr, "  Double A = %.10f [", da); show_double(da); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    A =       ["); show_erd(a); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double B = %.10f [", db); show_double(db); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    B =       ["); show_erd(b); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double Sum    = %.10f [", dsum); show_double(dsum); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double Result = %.10f [", dnval); show_double(dnval); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    Sum =       ["); show_erd(nval); fprintf(stderr, "]\n");
	} else
	    fprintf(stderr, "Sum OK\n");
    } else
	fprintf(stderr, "Out of range argument for sum\n");
    erd_mpf_check("erd_add result", nval);
    mpf_t csum;
    mpf_init2(csum, 64);
    mpf_add(csum, ma, mb);
    if (!erd_compare_mpf(nval, csum)) {
	fprintf(stderr, "   Computing sum:\n");
	fprintf(stderr, "   A = %s [", mpf_string(ma)); show_erd(a); fprintf(stderr, "]\n");
	fprintf(stderr, "   B = %s [", mpf_string(ma)); show_erd(b); fprintf(stderr, "]\n");
	fprintf(stderr, "   Expected sum = %s\n", mpf_string(csum));
    }
    mpf_clear(ma); mpf_clear(mb); mpf_clear(csum);
#endif /* EDEBUG */
    return nval;
}    

erd_t erd_mul(erd_t a, erd_t b) {
#if EDEBUG
    erd_mpf_check("erd_mul a", a);
    erd_mpf_check("erd_mul b", b);
    mpf_t ma, mb;
    erd_to_mpf_with_init(ma, a);
    erd_to_mpf_with_init(mb, b);
    double da = erd_to_double(a);
    double db = erd_to_double(b);
#endif
    erd_t nval;
    nval.exh = a.exh + b.exh;
    nval.dbl = a.dbl * b.dbl;
#if EDEBUG
    if (!dbl_bad_exponent(da) && !dbl_bad_exponent(db)) {
	double dprod = da * db;
	double dnval = erd_to_double(nval);
	if (dprod != dnval) {
	    fprintf(stderr, "Double prod mismatch:\n");
	    fprintf(stderr, "  Double A = %.10f [", da); show_double(da); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    A =       ["); show_erd(a); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double B = %.10f [", db); show_double(db); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    B =       ["); show_erd(b); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double Prod    = %.10f [", dprod); show_double(dprod); fprintf(stderr, "]\n");
	    fprintf(stderr, "  Double Result = %.10f [", dnval); show_double(dnval); fprintf(stderr, "]\n");
	    fprintf(stderr, "  ERD    Prod =       ["); show_erd(nval); fprintf(stderr, "]\n");
	} else
	    fprintf(stderr, "Product OK\n");
    } else
	fprintf(stderr, "Out-of-range argument for product\n");
    erd_mpf_check("erd_mul result", erd_normalize(nval));
    mpf_t cprod;
    mpf_init2(cprod, 64);
    mpf_mul(cprod, ma, mb);
    if (!erd_compare_mpf(erd_normalize(nval), cprod)) {
	fprintf(stderr, "   Computing product:\n");
	fprintf(stderr, "   A = %s [", mpf_string(ma)); show_erd(a); fprintf(stderr, "]\n");
	fprintf(stderr, "   B = %s [", mpf_string(ma)); show_erd(b); fprintf(stderr, "]\n");
	fprintf(stderr, "   Expected product = %s\n", mpf_string(cprod));
    }
    mpf_clear(ma); mpf_clear(mb); mpf_clear(cprod);
#endif
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

#if EDEBUG

#define SHADOW_TOLERANCE 0.0001

static char shadow_buf[2048];

static const char *mpf_string(mpf_srcptr val) {
    int digits = 10;
    char boffset = 0;
    mp_exp_t ecount;
    char *sval = mpf_get_str(NULL, &ecount, 10, digits, val);
    if (!sval || strlen(sval) == 0 || sval[0] == '0') {
	strcpy(shadow_buf, "0.0");
    } else {
	int voffset = 0;
	bool neg = sval[0] == '-';
	if (neg) {
	    voffset++;
	    shadow_buf[boffset++] = '-';
	}
	if (ecount == 0) {
	    shadow_buf[boffset++] = '0';
	    shadow_buf[boffset++] = '.';
	} else {
	    shadow_buf[boffset++] = sval[voffset++];
	    shadow_buf[boffset++] = '.';
	    ecount--;
	}
	if (sval[voffset] == 0)
	    shadow_buf[boffset++] = '0';
	else {
	    while(sval[voffset] != 0)
		shadow_buf[boffset++] = sval[voffset++];
	}
	if (ecount != 0) {
	    shadow_buf[boffset++] = 'e';
	    snprintf(&shadow_buf[boffset], 24, "%ld", (long) ecount);
	} else
	    shadow_buf[boffset] = 0;
    }
    free(sval);
    return shadow_buf;
}


static bool erd_compare_mpf(erd_t val, mpf_srcptr mval) {
    mpf_t tval;
    mpf_init2(tval, 64);
    erd_to_mpf(tval, val);
    mpf_sub(tval, tval, mval);
    mpf_div(tval, tval, mval);
    mpf_abs(tval, tval);
    double err = mpf_get_d(tval);
    bool ok = true;
    if (err > SHADOW_TOLERANCE) {
	ok = false;
	erd_to_mpf(tval, val);
	fprintf(stderr, "Shadow mismatch:\n");
	fprintf(stderr, "   ERD Value: %s [", mpf_string(tval)); show_erd(val); fprintf(stderr, "]\n");
	fprintf(stderr, "   MPF Value: %s\n", mpf_string(mval));
    }	    
    mpf_clear(tval);
    return ok;
}



static void show_double(double d) {
    int sign = dbl_get_sign(d);
    int exp = dbl_get_exponent(d);
    uint64_t frac = dbl_get_fraction(d);
    printf("Sign=%d, Exp=%d, Frac=0x%llx, Val=%.8f", sign, exp, (unsigned long long) frac, d);
}

static void show_erd(erd_t a) {
    int sign = dbl_get_sign(a.dbl);
    int la = dbl_get_exponent(a.dbl);
    uint64_t frac = dbl_get_fraction(a.dbl);
    int64_t ha = a.exh;
    int64_t ta = erd_get_full_exponent(a);
    printf("Sign=%d, Exp=(%lld*%d)+%d = %lld, Frac=0x%llx", sign, (long long) ha, ER_MODULUS, la,
	   (long long) ta, (unsigned long long) frac);
}

static int fcmp(double a, double b) {
    if (a < b)
	return -1;
    if (a > b)
	return 1;
    return 0;
}

static void check_mismatch(double a, erd_t xa) {
    int64_t texp = erd_get_full_exponent(xa);
    double na = dbl_replace_exponent(xa.dbl, texp);
    if (a == na)
	printf(" *==* ");
    else
	printf(" *!!* ");
    mpf_t from_d;
    mpf_t from_erd;
    mpf_init2(from_d, 64);
    mpf_init2(from_erd, 64);
    mpf_set_d(from_d, a);
    erd_to_mpf(from_erd, xa);
    erd_t mxa = erd_from_mpf(from_d);
    if (mpf_cmp(from_d, from_erd) == 0)
	printf(" *=DMXM=* ");
    else
	printf(" *!DMXM!* ");
    if (erd_cmp(xa, mxa) == 0)
	printf(" *=DXMD=* ");
    else {
	printf(" *!DXMD!* "); 
	show_erd(mxa);
    }

}

static void erd_to_mpf_with_init(mpf_ptr dest, erd_t eval) {
    mpf_init2(dest, 64);
    erd_to_mpf(dest, eval);
}

static bool erd_mpf_check(const char *context, erd_t val) {
    mpf_t mval;
    erd_to_mpf_with_init(mval, val);
    erd_t nval = erd_from_mpf(mval);
    bool ok  = true;
    if (erd_cmp(val, nval) != 0) {
	ok = false;
	fprintf(stderr, "ERD->MPF->ERD failure (%s)\n", context);
	fprintf(stderr, "   Original  value: \n");
	show_erd(val); fprintf(stderr, "\n");
	fprintf(stderr, "   MPF value: %s\n", mpf_string(mval));
	fprintf(stderr, "   Generated value: \n");
	show_erd(val); fprintf(stderr, "\n");
    }
    mpf_clear(mval);
    return ok;
}

volatile int sink;

static void erd_fail() {
    int *p = NULL;
    sink = *p;
}

void er_check(double a, double b) {
    erd_t xa = erd_from_double(a);
    erd_t xb = erd_from_double(b);
    printf("a =     "); show_double(a); printf("\n");
    printf("    --> "); show_erd(xa); check_mismatch(a, xa); printf("\n");
    printf("b =     "); show_double(b); printf("\n");
    printf("    --> "); show_erd(xb); check_mismatch(b, xb); printf("\n");
    printf("a * b = "); show_double(a*b); printf("\n");
    printf("    --> "); show_erd(erd_mul(xa, xb)); check_mismatch(a*b, erd_mul(xa, xb)); printf("\n");
    printf("a + b = "); show_double(a+b); printf("\n");
    printf("    --> "); show_erd(erd_add(xa, xb)); check_mismatch(a+b, erd_add(xa, xb)); printf("\n");
    printf("1/a   = "); show_double(1.0/a);  printf("\n");
    printf("    --> "); show_erd(erd_recip(xa));  check_mismatch(1.0/a, erd_recip(xa)); printf("\n");
    printf("a:b   = "); printf("%d\n", fcmp(a, b));
    printf("    --> "); printf("%d\n", erd_cmp(xa, xb));

}

#endif
