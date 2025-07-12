/* Testing implementations of extended-range double */


#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include <gmp.h>

#include "report.h"
#include "q25.h"
#include "analysis.h"


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

static double dbl_clear_exponent(double x) {
    uint64_t bx = dbl_get_bits(x);
    uint64_t mask = ~(DBL_EXP_MASK << DBL_EXP_OFFSET);
    bx &= mask;
    return dbl_from_bits(bx);
}


static double dbl_infinity(int sign) {
    return dbl_assemble(sign, DBL_EXP_MASK - DBL_BIAS, 0);
}


/********************* ERDZ *************************/

bool erdz_is_zero(erd_t a) {
    return a.dbl == 0.0;
}

erd_t erdz_zero() {
    erd_t nval;
    nval.dbl = 0.0;
    nval.exp = 0.0;
    return nval;
}


static erd_t erdz_normalize(erd_t a) {
    if (erdz_is_zero(a))
	return erdz_zero();
    erd_t nval;
    nval.exp = a.exp + dbl_get_exponent(a.dbl);
    nval.dbl = dbl_clear_exponent(a.dbl);
    return nval;
}

erd_t erdz_from_double(double dval) {
    erd_t nval;
    nval.exp = 0;
    nval.dbl = dval;
    return erdz_normalize(nval);
}

erd_t erdz_from_mpf(mpf_srcptr fval) {
    erd_t nval;
    long int exp;
    nval.dbl = mpf_get_d_2exp(&exp, fval);
    if (nval.dbl == 0)
	return erdz_zero();
    nval.exp = (int64_t) exp;
    return erdz_normalize(nval);
}

void erdz_to_mpf(mpf_ptr dest, erd_t eval) {
    mpf_set_d(dest, eval.dbl);
    if (eval.exp < 0)
	mpf_div_2exp(dest, dest, -eval.exp);
    else if (eval.exp > 0)
	mpf_mul_2exp(dest, dest, eval.exp);
}

double erdz_to_double(erd_t eval) {
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

erd_t erdz_negate(erd_t a) {
    erd_t nval;
    if (erdz_is_zero(a))
	return a;
    nval.exp = a.exp;
    nval.dbl = -a.dbl;
    return nval;
}

erd_t erdz_add(erd_t a, erd_t b) {
    if (erdz_is_zero(a))
	return b;
    if (erdz_is_zero(b))
	return a;
    if (a.exp > b.exp + DBL_MAX_PREC)
	return a;
    if (b.exp > a.exp + DBL_MAX_PREC)
	return b;
    erd_t nval;
    int64_t ediff = a.exp - b.exp;
    double ad = dbl_replace_exponent(a.dbl, ediff);
    nval.dbl = ad + b.dbl;
    nval.exp = b.exp;
    return erdz_normalize(nval);
}

erd_t erdz_mul(erd_t a, erd_t b) {
    erd_t nval;
    nval.exp = a.exp + b.exp;
    nval.dbl = a.dbl * b.dbl;
    return erdz_normalize(nval);
}

erd_t erdz_muladd(erd_t a, erd_t b, erd_t c) {
    erd_t prod;
    prod.dbl = a.dbl * b.dbl;
    if (prod.dbl == 0)
	return c;
    prod.exp = a.exp + b.exp;
    prod.exp += dbl_get_exponent(prod.dbl);
    prod.dbl = dbl_clear_exponent(prod.dbl);
    if (c.dbl == 0)
	return prod;
    if (prod.exp > c.exp + DBL_MAX_PREC)
	return prod;
    if (c.exp > prod.exp + DBL_MAX_PREC)
	return c;
    erd_t nval;
    int64_t ediff = prod.exp - c.exp;
    double ad = dbl_replace_exponent(prod.dbl, ediff);
    nval.dbl = ad + c.dbl;
    nval.exp = c.exp;
    return erdz_normalize(nval);
}

erd_t erdz_mul_seq_slow(erd_t *val, int len) {
    erd_t result = erdz_from_double(1.0);
    int i;
    for (i = 0; i < len; i++)
	result = erdz_mul(result, val[i]);
    return result;
}

/* Max number of times fractions can be multiplied without overflowing exponent */
#define MAX_MUL 1000

erd_t erdz_mul_seq_x1(erd_t *val, int len) {
    erd_t result = len == 0 ? erdz_from_double(1.0) : val[0];
    int i;
    int count = 1;
    for (i = 1; i < len; i++) {
	erd_t arg = val[i];
	result.dbl *= arg.dbl;
	result.exp += arg.exp;
	if (++count > MAX_MUL) {
	    count = 0;
	    result = erdz_normalize(result);
	}
    }
    return erdz_normalize(result);
}

erd_t erdz_mul_seq_x4(erd_t *val, int len) {
    erd_t prod[4];
    int i, j;
    for (j = 0; j < 4; j++) 
	prod[j] = erdz_from_double(1.0);
    int count = 0;
    for (i = 0; i <= len-4; i+= 4) {
	for (j = 0; j < 4; j++) {
	    erd_t arg = val[i+j];
	    prod[j].dbl *= arg.dbl;
	    prod[j].exp += arg.exp;
	}
	if (++count > MAX_MUL) {
	    count = 0;
	    for (j = 0; j < 4; j++)
		prod[j] = erdz_normalize(prod[j]);
	}
    }
    erd_t result = erdz_normalize(prod[0]);
    for (; i < len; i++)
	result = erdz_mul(result, val[i]);
    for (j = 1; j < 4; j++)
	result = erdz_mul(result, erdz_normalize(prod[j]));
    return result;
}

erd_t erdz_mul_seq(erd_t *val, int len) {
    if (len <= 8)
	return erdz_mul_seq_x1(val, len);
    return erdz_mul_seq_x4(val, len);
}

erd_t erdz_div(erd_t a, erd_t b) {
    erd_t nval;
    nval.dbl = a.dbl / b.dbl;
    nval.exp = a.exp - b.exp;
    return erdz_normalize(nval);
}

q25_ptr erdz_to_q25(erd_t a) {
    q25_ptr q = q25_from_double(a.dbl);
    q25_inplace_scale(q, (int32_t) a.exp, 0);
    return q;
}

double digit_precision_erdz(erd_t x_est, mpf_srcptr x) {
    int pos = q25_enter();
    q25_ptr q_est = erdz_to_q25(x_est);
    q25_ptr q = q25_from_mpf(x);
    double p = digit_precision_q25(q_est, q);
    q25_leave(pos);
    return p;
}


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

const char *erdz_string(erd_t a) {
    mpf_t mval;
    mpf_init2(mval, 64);
    erdz_to_mpf(mval, a);
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
double run_sum_erdz(erd_t *result, double *dval, int len, int reps) {
    erd_t *eval = calloc(len, sizeof(erd_t));
    int i;
    for (i = 0; i < len; i++)
	eval[i] = erdz_from_double(dval[i]);
    int r;
    double t = tod();
    erd_t s = erdz_zero();
    for (r = 0; r < reps; r++)
	for (i = 0; i < len; i++)
	    s = erdz_add(s, eval[i]);
    *result = s;
    t = tod() - t;
    free(eval);
    return t;
}

double uniform_value(double min, double max) {
    double u = (double) random() / (double) ((1L<<31)-1);
    return min + u * (max-min);
}

double exponential_value(double base, double minp, double maxp) {
    double exp = uniform_value(minp, maxp);
    return pow(base, exp);
}

double *uniform_array(int len, double min, double max, unsigned seed) {
    srandom(seed);
    double *dval = calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = uniform_value(min, max);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    
    return dval;
}

double *exponential_array(int len, double base, double minp, double maxp, unsigned seed) {
    srandom(seed);
    double *dval = calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = exponential_value(base, minp, maxp);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    return dval;
}

void run_sum(char *prefix, double *data, int len, int reps) {
    double dval;
    mpf_t mval;
    mpf_init2(mval, 64);
    erd_t eval;
    double dt = run_sum_dbl(&dval, data, len, reps);
    report(2, "Ran double\n");
    double mt = run_sum_mpf(mval, data, len, reps);
    report(2, "Ran MPF\n");
    double et = run_sum_erdz(&eval, data, len, reps);
    report(2, "Ran ERD\n");
    const char *ms = mpf_string(mval, 20);
    const char *es = erdz_string(eval);

    int pos = q25_enter();
    report(2, "Converting to Q25\n");
    q25_ptr dq = q25_from_double(dval);
    q25_ptr mq = q25_from_mpf(mval);
    q25_ptr eq = erdz_to_q25(eval);
    report(2, "Computing digit precision\n");
    double dpd = digit_precision_q25(dq, mq);
    double dpe = digit_precision_q25(eq, mq);
    long sums = (long) len * reps;
    report(1, "%s: Len = %d reps = %d sums = %ld\n",
	   prefix, len, reps, sums);
    report(1, "    Double: Sum = %.20f ns/sum = %.2f precision = %.2f\n",
	   dval, dt * 1e9 / sums, dpd);
    report(1, "    ERD: Sum = %s ns/sum = %.2f precision = %.2f\n",
	   es, et * 1e9 / sums, dpe);
    report(1, "    MPF: Sum = %s ns/sum = %.2f\n",
	   ms, mt * 1e9 / sums);
    q25_leave(pos);
    free((void *) ms); free((void *) es);
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [-h] [-v VERB] [-n CNT] [-r REPS] [-s SEED] [-d (u|e)] [-m DMIN] [-M DMAX]\n", name);
    fprintf(stderr, "   -h      Print this message\n");
    fprintf(stderr, "   -n CNT  Data size\n"); 
    fprintf(stderr, "   -r REPS Repetitions\n");
    fprintf(stderr, "   -d DIST Distribution: uniform or exponential\n");
    fprintf(stderr, "   -m MIN  Data minimum (Power of 10 when exponential)\n");
    fprintf(stderr, "   -M MAX  Data maximum (Power of 10 when exponential)\n");
}

int main(int argc, char *argv[]) {
    int len = 1000;
    int reps = 1;
    bool exponential = false;
    double dmin = 0.0;
    double dmax = 1.0;
    unsigned seed = 12345;
    int c;
    while ((c = getopt(argc, argv, "hv:n:r:d:m:M:")) != -1) {
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
	case 'm':
	    dmin = atof(optarg);
	    break;
	case 'M':
	    dmax = atof(optarg);
	    break;
	}
    }
    char buf[100];
    snprintf(buf, 100, "%s[%.2f, %.2f]", exponential ? "Exp" : "Uni", dmin, dmax);
    printf("Running %s\n", buf);
    double *data = exponential 
	? exponential_array(len, 10, dmin, dmax, seed)
	: uniform_array(len, dmin, dmax, seed);
    printf("Generated data\n");
    run_sum(buf, data, len, reps);
    free(data);
    return 0;
}
