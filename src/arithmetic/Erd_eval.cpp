/* Testing implementations of extended-range double */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>

#include "report.h"
#include "Erd.hh"

/*********** Useful functions ***********/
const char *mpf_string(mpf_class &val, int digits) {
    char buf[2048];
    char boffset = 0;
    mp_exp_t ecount;
    char *sval = mpf_get_str(NULL, &ecount, 10, digits, val.get_mpf_t());
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

const char *erd_mpf_string(Erd a) {
    mpf_class mval(0.0, 64);
    mval = a.get_mpf();
    return mpf_string(mval, 20);
}

static double dbl_sum_seq_x4(double *val, int len) {
    // Assume that len >= 4
    int i, j;
    double sum[4];
    for (j = 0; j < 4; j++)
	sum[j] = val[j];
    for (i = 4; i <= len-4; i+=4)
	for (j = 0; j < 4; j++)
	    sum[j] += val[i+j];
    double result = (sum[0]+sum[1])+(sum[2]+sum[3]);
    for (; i < len; i++)
	result += val[i];
    return result;
}

static double dbl_sum_seq(double *val, int len) {
    if (len >= 4)
	return dbl_sum_seq_x4(val, len);
    double result = 0.0;
    int i;
    for (i = 0; i < len; i++)
	result += val[i];
    return result;
}

/* Time and run summations  */
double run_sum_dbl(double *result, double *dval, int len, int reps) {
    int i, r;
    double t = tod();
    double s = 0.0;
    for (r = 0; r < reps; r++)
	s += dbl_sum_seq(dval, len);
    *result = s;
    t = tod() - t;
    return t;
}

/* Time and run summations  */
double run_sum_mpf(mpf_class &result, double *dval, int len, int reps) {
    std::vector<mpf_class> mval;
    mval.resize(len);
    for (int i = 0; i < len; i++) 
	mval[i] = dval[i];

    double t = tod();
    result = 0.0;
    for (int r = 0; r < reps; r++)
	for (int i = 0; i < len; i++)
	    result += mval[i];
    t = tod() - t;
    return t;
}

/* Time and run summations  */
double run_sum_erd(Erd &result, double *dval, int len, int reps) {
    std::vector<Erd> eval;
    eval.resize(len);
    for (int i = 0; i < len; i++)
	eval[i] = dval[i];
    double t = tod();
    result = 0.0;
    for (int r = 0; r < reps; r++)
	for (int i = 0; i < len; i++)
	    result += eval[i];
    t = tod() - t;
    return t;
}


static double dbl_prod_seq_x4(double *val, int len) {
    // Assume that len >= 4
    int i, j;
    double prod[4];
    for (j = 0; j < 4; j++)
	prod[j] = val[j];
    for (i = 4; i <= len-4; i+=4)
	for (j = 0; j < 4; j++)
	    prod[j] *= val[i+j];
    double result = (prod[0]*prod[1])*(prod[2]*prod[3]);
    for (; i < len; i++)
	result *= val[i];
    return result;
}

static double dbl_prod_seq(double *val, int len) {
    if (len >= 4)
	return dbl_prod_seq_x4(val, len);
    double result = 1.0;
    int i;
    for (i = 0; i < len; i++)
	result *= val[i];
    return result;
}

/* Time and run products  */
double run_prod_dbl(double *result, double *dval, int len, int reps) {
    int i, r;
    double t = tod();
    double s = 1.0;
    for (r = 0; r < reps; r++)
	s *= dbl_prod_seq(dval, len);
    *result = s;
    t = tod() - t;
    return t;
}

/* Time and run products  */
double run_prod_mpf(mpf_class &result, double *dval, int len, int reps) {
    /* Run with MPF */
    std::vector<mpf_class> mval;
    mval.resize(len);
    for (int i = 0; i < len; i++)
	mval[i] = dval[i];
    double t = tod();

    result = 1.0;
    for (int r = 0; r < reps; r++)
	for (int i = 0; i < len; i++)
	    result *= mval[i];
    t = tod() - t;
    return t;
}

/* Time and run products  */
double run_prod_erd(Erd &result, double *dval, int len, int reps) {
    std::vector<Erd> eval;
    eval.resize(len);
    for (int i = 0; i < len; i++) 
	eval[i] = dval[i];

    double t = tod();
    result = 1.0;
    for (int r = 0; r < reps; r++) 
	result *= product_reduce(eval);
    t = tod() - t;
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
    double *dval = (double *) calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = uniform_value(min, max, zpct);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    
    return dval;
}

double *exponential_array(int len, double base, double minp, double maxp, double zpct, unsigned seed) {
    srandom(seed);
    double *dval = (double *) calloc(len, sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
	dval[i] = exponential_value(base, minp, maxp, zpct);
	report(4, "d[%d] = %.5f\n", i, dval[i]);
    }
    return dval;
}

double digit_precision(mpf_class &x_est, mpf_class &x) {
    if (x_est == x)
	return 1e6;
    if (x_est == 0)
	return 0.0;
    if (x == 0)
	return 0.0;
    mpf_t rel;
    mpf_init2(rel, 128);
    mpf_sub(rel, x_est.get_mpf_t(), x.get_mpf_t());
    mpf_div(rel, rel, x.get_mpf_t());
    mpf_abs(rel, rel);
    long int exp;
    double drel = mpf_get_d_2exp(&exp, rel);
    double dp = -(log10(drel) + log10(2.0) * exp);
    return dp < 0 ? 0 : dp;
}

void run_sum(char *prefix, double *data, int len, int reps) {
    double dval;
    mpf_class mval;
    Erd eval;
    double dt = run_sum_dbl(&dval, data, len, reps);
    double mt = run_sum_mpf(mval, data, len, reps);
    double et = run_sum_erd(eval, data, len, reps);
    const char *ms = mpf_string(mval, 20);
    const char *es = erd_mpf_string(eval);
    Erd log10 = eval.log10();
    mpf_class md(0, 64);
    double dpd;
    if (dbl_exponent_below(dbl_get_exponent(dval)) || dbl_exponent_above(dbl_get_exponent(dval)))
	dpd = 0.0;
    else {
	md = dval;
	dpd = digit_precision(md, mval);
    }
    mpf_class me = eval.get_mpf();
    double dpe = digit_precision(me, mval);
    long sums = (long) len * reps;
    report(1, "%s: Len = %d reps = %d sums = %ld\n",
	   prefix, len, reps, sums);
    report(1, "    DBL: Sum = %.20f ps/sum = %.2f precision = %.2f\n",
	   dval, dt * 1e12 / sums, dpd);
    report(1, "    ERD: Sum = %s ps/sum = %.2f precision = %.2f\n",
	   es, et * 1e12 / sums, dpe);
    std::cout << "c     Cout Sum = " << eval << " log10 = " << log10 << std::endl;
    report(1, "    MPF: Sum = %s ps/sum = %.2f\n",
	   ms, mt * 1e12 / sums);
    report(1, "    MPF:DBL = %f  MPF:ERD = %f ERD:DBL = %f\n",
	   mt/dt, mt/et, et/dt);
}

void run_prod(char *prefix, double *data, int len, int reps) {
    double dval;
    mpf_class mval;
    Erd eval;
    double dt = run_prod_dbl(&dval, data, len, reps);
    double mt = run_prod_mpf(mval, data, len, reps);
    double et = run_prod_erd(eval, data, len, reps);
    report(1, "Times: DBL %f MPF %f ERD %f\n", dt, mt, et);
    const char *ms = mpf_string(mval, 20);
    const char *es = erd_mpf_string(eval);
    Erd log10 = eval.log10();
    mpf_class md(0, 64);
    double dpd;
    if (dbl_exponent_below(dbl_get_exponent(dval)) || dbl_exponent_above(dbl_get_exponent(dval)))
	dpd = 0.0;
    else {
	md = dval;
	dpd = digit_precision(md, mval);
    }
    mpf_class me = eval.get_mpf();
    double dpe = digit_precision(me, mval);
    long prods = (long) len * reps;
    report(1, "%s: Len = %d reps = %d prods = %ld\n",
	   prefix, len, reps, prods);
    report(1, "    DBL: Product = %.20f ps/prod = %.2f precision = %.2f\n",
	   dval, dt * 1e12 / prods, dpd);
    report(1, "    ERD: Product = %s ps/prod = %.2f precision = %.2f\n",
	   es, et * 1e12 / prods, dpe);
    std::cout << "c     Cout Product = " << eval << " log10 = " << log10 << std::endl;
    report(1, "    MPF: Product = %s ps/prod = %.2f\n",
	   ms, mt * 1e12 / prods);
    report(1, "    MPF:DBL = %f  MPF:ERD = %f ERD:DBL = %f\n",
	   mt/dt, mt/et, et/dt);
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
    mpf_set_default_prec(64);
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
