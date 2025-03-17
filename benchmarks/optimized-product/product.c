/* Compute products yielding lowest digit precision */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include <gmp.h>

#include <mpfr.h>
#include <mpfi.h>

/* Grab code from model counting */

/* Upper threshold for digit precision metric */
#define MAX_DIGIT_PRECISION (1000*1000)
/* Number of powers of two */
#define P2COUNT 20
/* Number of powers of ten */
#define P10COUNT 6
/* MPF Precision */
#define MPF_PREC 128
/* What is the denominator of the offset? */
#define OFFSET_DENOM (1000*1000*1000)
/* How many digits to display */
#define DIGITS 40

char *mpf_string(mpf_srcptr val, int digits) {
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
    char *rstring = (char *) malloc(strlen(buf)+1);
    strcpy(rstring, buf);
    return (char *) rstring;
}

double digit_precision_mpfr(mpfr_srcptr x_est, mpq_srcptr x) {
    if (mpfr_cmp_q(x_est, x) == 0)
	return (double) MAX_DIGIT_PRECISION;

    mpfr_prec_t save_prec = mpfr_get_default_prec();
    mpfr_prec_t prec = 3*mpfr_get_prec(x_est);
    mpfr_set_default_prec(prec);

    mpfr_t num;
    mpfr_t den;
    if (mpq_sgn(x) == 0) {
	mpfr_init_set_d(den, 1.0, MPFR_RNDN);
	mpfr_init_set(num, x_est, MPFR_RNDN);
	mpfr_abs(num, num, MPFR_RNDN);
	if (mpfr_cmp_d(num, 1.0) > 0)
	    mpfr_set_d(num, 1.0, MPFR_RNDN);
    } else {
	mpfr_init_set_q(den, x, MPFR_RNDN);
	mpfr_abs(den, den, MPFR_RNDN);
	mpfr_init_set_q(num, x, MPFR_RNDN);
	mpfr_sub(num, num, x_est, MPFR_RNDN);
	mpfr_abs(num, num, MPFR_RNDN);
    }
    mpfr_div(num, num, den, MPFR_RNDN);
    mpfr_log10(num, num, MPFR_RNDN);
    double result = -mpfr_get_d(num, MPFR_RNDN);
    if (result < 0)
	result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
	result = MAX_DIGIT_PRECISION;
    mpfr_clears(num, den, NULL);
    mpfr_set_default_prec(save_prec);
    return result;
}

double digit_precision_mpf(mpf_srcptr x_est, mpq_srcptr x) {
    mpfr_t rx_est;
    int prec = mpf_get_prec(x_est);
    mpfr_init2(rx_est, prec);
    mpfr_set_f(rx_est, x_est, MPFR_RNDN);
    double result = digit_precision_mpfr(rx_est, x);
    mpfr_clear(rx_est);
    return result;
}

/* Tables of powers of two */


mpq_t p2_table_mpq[P2COUNT];
mpf_t p2_table_mpf[P2COUNT];

/* Largest value of k that gave improved entry */
int kmax = 0;

int k2_best[P2COUNT];
double dp2_best[P2COUNT];
char *v2_best[P2COUNT];

int k10_best[P10COUNT];
double dp10_best[P10COUNT];
char *v10_best[P10COUNT];

void init() {
    mpf_set_default_prec(MPF_PREC);
    mpfr_set_default_prec(MPF_PREC);
    int i;
    for (i = 0; i < P2COUNT; i++) {
	k2_best[i] = 0;
	dp2_best[i] = MAX_DIGIT_PRECISION;
	mpq_init(p2_table_mpq[i]);
	mpf_init(p2_table_mpf[i]);
	v2_best[i] = NULL;
    }
    for (i = 0; i < P10COUNT; i++) {
	k10_best[i] = 0;
	dp10_best[i] = MAX_DIGIT_PRECISION;
	v10_best[i] = NULL;
    }
}


void show_entry(int power, int k, double dp, char *sval) {
    printf("%d,%d,%.4f,%s\n",
	   power, k, dp, sval);
}


void fill_table(int k) {
    mpq_set_si(p2_table_mpq[0], OFFSET_DENOM + k, OFFSET_DENOM);

    int i = 0;
    double dp;
    mpf_set_q(p2_table_mpf[i], p2_table_mpq[i]);
    dp = digit_precision_mpf(p2_table_mpf[i], p2_table_mpq[i]);
    if (dp < dp2_best[i]) {
	dp2_best[i] = dp;
	if (v2_best[i])
	    free(v2_best[i]);
	v2_best[i] = mpf_string(p2_table_mpf[i], DIGITS);
	k2_best[i] = k;
	kmax = k;
    }

    for (i = 1; i < P2COUNT; i++) {
	mpq_mul(p2_table_mpq[i], p2_table_mpq[i-1], p2_table_mpq[i-1]);
	mpf_mul(p2_table_mpf[i], p2_table_mpf[i-1], p2_table_mpf[i-1]);
	dp = digit_precision_mpf(p2_table_mpf[i], p2_table_mpq[i]);
	if (dp < dp2_best[i]) {
	    dp2_best[i] = dp;
	    if (v2_best[i])
		free(v2_best[i]);
	    v2_best[i] = mpf_string(p2_table_mpf[i], DIGITS);
	    k2_best[i] = k;
	    kmax = k;
	}
    }
}

double dp_power(int k, int pwr, mpf_ptr fval) {
    mpq_t qval;
    mpq_init(qval);
    mpq_set_si(qval, 1, 1);
    mpf_set_d(fval, 1.0);
    int i;
    int cnt = 0;
    for (i = 0; pwr > 0; i++) {
	if (pwr & (1<<i)) {
	    cnt++;
	    pwr ^= (1<<i);
	    mpq_mul(qval, qval, p2_table_mpq[i]);
	    mpf_mul(fval, fval, p2_table_mpf[i]);
	}
    }
    double dp =  digit_precision_mpf(fval, qval);
    mpq_clear(qval);
    return dp;
}

void sweep(int klimit) {
    int k;
    for (k = 0; k <= klimit; k++) {
	fill_table(k);
	int i;
	int p10 = 1;
	mpf_t fval;
	mpf_init(fval);
	for (i = 0; i < P10COUNT; i++) {
	    p10 *= 10;
	    double dp = dp_power(k, p10, fval);
	    if (dp < dp10_best[i]) {
		k10_best[i] = k;
		kmax = k;
		if (v10_best[i])
		    free(v10_best[i]);
		v10_best[i] = mpf_string(fval, DIGITS);
		dp10_best[i] = dp;
	    }
	}
    }
    int i;
    for (i = 0; i < P2COUNT; i++)
	show_entry((1<<i), k2_best[i], dp2_best[i], v2_best[i]);
    int p10 = 1;
    for (i = 0; i < P10COUNT; i++) {
	p10 *= 10;
	show_entry(p10, k10_best[i], dp10_best[i], v10_best[i]);
    }
    fprintf(stderr, "kmax = %d\n", kmax);
}


int main(int argc, char *argv[]) {
    int klim = 100;
    init();
    if (argc > 1) {
	if (strcmp(argv[1], "-h") == 0) {
	    printf("Usage: %s KLIM\n", argv[0]);
	    return 0;
	}
	klim = atoi(argv[1]);
    }
    sweep(klim);
    return 0;
}
