/* Testing implementations of extended-range double */


#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include <gmp.h>

#include "report.h"

/*
  Representation of floating-point numbers based on double,
  but with additional exponent field to support extended range
 */
typedef struct {
    double dbl; 
    int64_t exp; 
} erd_t;

/* 
   Support two versions: 
   erdz: 0.0 has exp = 0
   erdm: 0.0 has exp = INT64_MIN
*/

#define ERDZ 0

#if ERDZ
#define ZEXP 0
#else
#define ZEXP INT64_MIN
#endif

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


/*** Debugging help ***/

void frees(const char *s) {
    free((void *) s);
}

const char* erd_document(erd_t a) {
    char buf[100];
    int sign = dbl_get_sign(a.dbl);
    int exp = dbl_get_exponent(a.dbl);
    uint64_t frac = dbl_get_fraction(a.dbl);
    snprintf(buf, 100, "Sign=%d, Exp=%d+%lld=%lld, Frac=0x%llx", sign, 
	     exp, (long long) a.exp, (long long) exp + a.exp,
	     (unsigned long long) frac);
    return archive_string(buf);
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

erd_t erd_quick_mul(erd_t a, erd_t b) {
    erd_t nval;
    nval.exp = a.exp + b.exp;
    nval.dbl = a.dbl * b.dbl;
    return nval;
}


erd_t erd_mul(erd_t a, erd_t b) {
    return erd_normalize(erd_quick_mul(a, b));
}

erd_t erd_mul_seq_slow(erd_t *val, int len) {
    erd_t result = erd_from_double(1.0);
    int i;
    for (i = 0; i < len; i++)
	result = erd_mul(result, val[i]);
    return result;
}

/* Max number of times fractions can be multiplied without overflowing exponent */
#define MAX_MUL 1000

erd_t erd_mul_seq_x1(erd_t *val, int len) {
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

erd_t erd_mul_seq_x4(erd_t *val, int len) {
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

erd_t erd_mul_seq(erd_t *val, int len) {
    if (len <= 8)
	return erd_mul_seq_x1(val, len);
    return erd_mul_seq_x4(val, len);
}

erd_t erd_div(erd_t a, erd_t b) {
    erd_t nval;
    nval.dbl = a.dbl / b.dbl;
    nval.exp = a.exp - b.exp;
    return erd_normalize(nval);
}

#if 0
q25_ptr erd_to_q25(erd_t a) {
    q25_ptr q = q25_from_double(a.dbl);
    q25_inplace_scale(q, (int32_t) a.exp, 0);
    return q;
}
#endif





/*********** Useful functions ***********/

const char *mpf_string(mpf_srcptr val, int digits) {
    char buf[2048];
    char boffset = 0;
    mp_exp_t ecount;
    char *sval = mpf_get_str(NULL, &ecount, 10, digits, val);
    if (!sval || strlen(sval) == 0 || sval[0] == '0') {
	strcpy(buf, "0.0");
    } else {
	int voffset = 0;
	bool neg = sval[0] == '-';
	if (neg) {
	    voffset++;
	    buf[boffset++] = '-';
	}
	if (ecount == 0) {
	    buf[boffset++] = '0';
	    buf[boffset++] = '.';
	} else {
	    buf[boffset++] = sval[voffset++];
	    buf[boffset++] = '.';
	    ecount--;
	}
	if (sval[voffset] == 0)
	    buf[boffset++] = '0';
	else {
	    while(sval[voffset] != 0)
		buf[boffset++] = sval[voffset++];
	}
	if (ecount != 0) {
	    buf[boffset++] = 'e';
	    snprintf(&buf[boffset], 24, "%ld", (long) ecount);
	} else
	    buf[boffset] = 0;
    }
    free(sval);
    return archive_string(buf);
}

const char *erd_string(erd_t a) {
    mpf_t mval;
    mpf_init2(mval, 64);
    erd_to_mpf(mval, a);
    return mpf_string(mval, 20);
}

/* Time and run summations  */
double run_sum_dbl(double *result, double *dval, int len, int reps) {
    int i, r;
    double t = tod();
    double s = 0.0;
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    s += dval[i];
    *result = s;
    t = tod() - t;
    return t;
}

/* Time and run summations  */
double run_sum_mpf(mpf_t result, double *dval, int len, int reps) {
    /* Run with MPF */
    mpf_t *mval = calloc(len, sizeof(mpf_t));
    int i;
    for (i = 0; i < len; i++) {
	mpf_init2(mval[i], 64);
	mpf_set_d(mval[i], dval[i]);
    }
    int r;
    double t = tod();
    mpf_set_d(result, 0.0);
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    mpf_add(result, result, mval[i]);
    t = tod() - t;
    for (i = 0; i < len; i++)
	mpf_clear(mval[i]);
    free(mval);
    return t;
}

/* Time and run summations  */
double run_sum_erd(erd_t *result, double *dval, int len, int reps) {
    erd_t *eval = calloc(len, sizeof(erd_t));
    int i;
    for (i = 0; i < len; i++)
	eval[i] = erd_from_double(dval[i]);
    int r;
    double t = tod();
    erd_t s = erd_zero();
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    s = erd_add(s, eval[i]);
    *result = s;
    t = tod() - t;
    free(eval);
    return t;
}

/* Time and run products  */
double run_prod_dbl(double *result, double *dval, int len, int reps) {
    int i, r;
    double t = tod();
    double s = 1.0;
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    s *= dval[i];
    *result = s;
    t = tod() - t;
    return t;
}

/* Time and run products  */
double run_prod_mpf(mpf_t result, double *dval, int len, int reps) {
    /* Run with MPF */
    mpf_t *mval = calloc(len, sizeof(mpf_t));
    int i;
    for (i = 0; i < len; i++) {
	mpf_init2(mval[i], 64);
	mpf_set_d(mval[i], dval[i]);
    }
    int r;
    double t = tod();
    mpf_set_d(result, 1.0);
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    mpf_mul(result, result, mval[i]);
    t = tod() - t;
    for (i = 0; i < len; i++)
	mpf_clear(mval[i]);
    free(mval);
    return t;
}

/* Time and run products  */
double run_prod_erd(erd_t *result, double *dval, int len, int reps) {
    erd_t *eval = calloc(len, sizeof(erd_t));
    int i;
    for (i = 0; i < len; i++) {
	eval[i] = erd_from_double(dval[i]);
    }
    int r;
    double t = tod();
    erd_t s = erd_from_double(1.0);
    for (r = 0; r < reps; r++) {
	erd_t p = erd_mul_seq(eval, len);
	s = erd_mul(s, p);
    }
    *result = s;
    t = tod() - t;
    free(eval);
    return t;
}


double uniform_value(double min, double max, double zpct) {
    double z = (double) random() / (double) ((1L<<31)-1);
    if (z * 100 < zpct)
	return 0.0;
    double u = (double) random() / (double) ((1L<<31)-1);
    return min + u * (max-min);
}

double exponential_value(double base, double minp, double maxp, double zpct) {
    double z = (double) random() / (double) ((1L<<31)-1);
    if (z * 100 < zpct)
	return 0.0;
    double exp = uniform_value(minp, maxp, 0);
    return pow(base, exp);
}

double *uniform_array(int len, double min, double max, double zpct, unsigned seed) {
    srandom(seed);
    double *dval = calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = uniform_value(min, max, zpct);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    
    return dval;
}

double *exponential_array(int len, double base, double minp, double maxp, double zpct, unsigned seed) {
    srandom(seed);
    double *dval = calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = exponential_value(base, minp, maxp, zpct);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    return dval;
}

double digit_precision(mpf_srcptr x_est, mpf_srcptr x) {
    if (mpf_cmp(x_est, x) == 0)
	return 1e6;
    if (mpf_cmp_d(x_est, 0.0) == 0)
	return 0.0;
    if (mpf_cmp_d(x, 0.0) == 0)
	return 0.0;
    mpf_t rel;
    mpf_init2(rel, 128);
    mpf_sub(rel, x_est, x);
    mpf_div(rel, rel, x);
    mpf_abs(rel, rel);
    long int exp;
    double drel = mpf_get_d_2exp(&exp, rel);
    double dp = -(log10(drel) + log10(2.0) * exp);
    return dp < 0 ? 0 : dp;
}

void run_sum(char *prefix, double *data, int len, int reps) {
    double dval;
    mpf_t mval;
    mpf_init2(mval, 64);
    erd_t eval;
    double dt = run_sum_dbl(&dval, data, len, reps);
    double mt = run_sum_mpf(mval, data, len, reps);
    double et = run_sum_erd(&eval, data, len, reps);
    const char *ms = mpf_string(mval, 20);
    const char *es = erd_string(eval);
    mpf_t md;
    mpf_init2(md, 64);
    double dpd;
    if (dbl_exponent_below(dbl_get_exponent(dval)) || dbl_exponent_above(dbl_get_exponent(dval)))
	dpd = 0.0;
    else {
	mpf_set_d(md, dval);
	dpd = digit_precision(md, mval);
    }
    mpf_t me;
    mpf_init2(me, 64);
    erd_to_mpf(me, eval);
    double dpe = digit_precision(me, mval);
    long sums = (long) len * reps;
    report(1, "%s: Len = %d reps = %d sums = %ld\n",
	   prefix, len, reps, sums);
    report(1, "    DBL: Sum = %.20f ns/sum = %.2f precision = %.2f\n",
	   dval, dt * 1e9 / sums, dpd);
    report(1, "    ERD: Sum = %s ns/sum = %.2f precision = %.2f\n",
	   es, et * 1e9 / sums, dpe);
    report(1, "    MPF: Sum = %s ns/sum = %.2f\n",
	   ms, mt * 1e9 / sums);
    mpf_clear(mval); mpf_clear(md); mpf_clear(me);
    frees(ms); frees(es);
}

void run_prod(char *prefix, double *data, int len, int reps) {
    double dval;
    mpf_t mval;
    mpf_init2(mval, 64);
    erd_t eval;
    double dt = run_prod_dbl(&dval, data, len, reps);
    double mt = run_prod_mpf(mval, data, len, reps);
    double et = run_prod_erd(&eval, data, len, reps);
    const char *ms = mpf_string(mval, 20);
    const char *es = erd_string(eval);
    mpf_t md;
    mpf_init2(md, 64);
    double dpd;
    if (dbl_exponent_below(dbl_get_exponent(dval)) || dbl_exponent_above(dbl_get_exponent(dval)))
	dpd = 0.0;
    else {
	mpf_set_d(md, dval);
	dpd = digit_precision(md, mval);
    }
    mpf_t me;
    mpf_init2(me, 64);
    erd_to_mpf(me, eval);
    double dpe = digit_precision(me, mval);
    long prods = (long) len * reps;
    report(1, "%s: Len = %d reps = %d prods = %ld\n",
	   prefix, len, reps, prods);
    report(1, "    DBL: Product = %.20f ns/prod = %.2f precision = %.2f\n",
	   dval, dt * 1e9 / prods, dpd);
    report(1, "    ERD: Product = %s ns/prod = %.2f precision = %.2f\n",
	   es, et * 1e9 / prods, dpe);
    report(1, "    MPF: Product = %s ns/prod = %.2f\n",
	   ms, mt * 1e9 / prods);
    mpf_clear(mval); mpf_clear(md); mpf_clear(me);
    frees(ms); frees(es);
}


void usage(char *name) {
    fprintf(stderr, "Usage: %s [-h] [-v VERB] [-n CNT] [-z ZPCT] [-r REPS] [-s SEED] [-d (u|e)] [-m DMIN] [-M DMAX]\n", name);
    fprintf(stderr, "   -h      Print this message\n");
    fprintf(stderr, "   -n CNT  Data size\n"); 
    fprintf(stderr, "   -r REPS Repetitions\n");
    fprintf(stderr, "   -z ZPCT Set percentage of zeroes\n");
    fprintf(stderr, "   -d DIST Distribution: uniform or exponential\n");
    fprintf(stderr, "   -m MIN  Data minimum (Power of 10 when exponential)\n");
    fprintf(stderr, "   -M MAX  Data maximum (Power of 10 when exponential)\n");
}

int main(int argc, char *argv[]) {
    int len = 1000;
    int reps = 1;
    double zpct = 0;
    bool exponential = false;
    double dmin = 0.0;
    double dmax = 1.0;
    unsigned seed = 12345;
    int c;
    while ((c = getopt(argc, argv, "hv:n:r:d:m:M:z:")) != -1) {
	switch(c) {
	case 'h':
	    usage(argv[0]);
	    return 0;
	    break;
	case 'v':
	    set_verblevel(atoi(optarg));
	    break;
	case 'n':
	    len = atoi(optarg);
	    break;
	case 'r':
	    reps = atoi(optarg);
	    break;
	case 'd':
	    if (optarg[0] == 'e')
		exponential = true;
	    break;
	case 'z':
	    zpct = atof(optarg);
	    break;
	case 'm':
	    dmin = atof(optarg);
	    break;
	case 'M':
	    dmax = atof(optarg);
	    break;
	}
    }
    char buf[100];
    report(1, "Running with %s\n\n", ERDZ ? "ERDZ" : "ERDM");
    snprintf(buf, 100, "%s[%.2f, %.2f, Z=%.1f%%]", exponential ? "Exp" : "Uni",
	     dmin, dmax, zpct);
    double *data = exponential 
	? exponential_array(len, 10, dmin, dmax, zpct, seed)
	: uniform_array(len, dmin, dmax, zpct, seed);
    run_sum(buf, data, len, reps);
    report(1, "\n");
    run_prod(buf, data, len, reps);
    free(data);
    return 0;
}
