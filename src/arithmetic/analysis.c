#include <math.h>
#include <stdlib.h>
#include "analysis.h"

#define DEBUG 0

static uint64_t double2bits(double x) {
    union {
	double dx;
	uint64_t bx;
    } u;
    u.dx = x;
    return u.bx;
}    


static bool double_is_infinite(double x, bool *negativep) {
    uint64_t bits = double2bits(x);
    unsigned sign = bits >> 63;
    int  biased_exp = (bits >> 52) & 0x7FF;
    int64_t frac = bits & 0xFFFFFFFFFFFFFL;
    if (negativep)
	*negativep = (sign == 1);
    return biased_exp == 0x7FF &&  frac == 0;
}

static bool double_is_nan(double x) {
    uint64_t bits = double2bits(x);
    int  biased_exp = (bits >> 52) & 0x7FF;
    int64_t frac = bits & 0xFFFFFFFFFFFFFL;
    return biased_exp == 0x7FF &&  frac != 0;
}

double digit_error_q25(q25_ptr qx, q25_ptr qy) {
#if DEBUG
    {
	char *sqx = q25_string(qx);
	char *sqy = q25_string(qy);

	printf("Finding digit error for q25 numbers %s : %s.", sqx, sqy);
	free(sqx); free(sqy);
    }
#endif
    if (!q25_is_valid(qx) || q25_is_infinite(qx, NULL)
	|| !q25_is_valid(qy) || q25_is_infinite(qy, NULL)) {
#if DEBUG
	printf("  Invalid arguments.  Return -1.0\n");
#endif
	return -1.0;
    }
    if (q25_compare(qx, qy) == 0) {
#if DEBUG
	printf("  Strict equality.  Return %d\n", MAX_DIGIT_ERROR);
#endif
	return MAX_DIGIT_ERROR;
    }
    int pos = q25_enter();
    q25_ptr aqx = q25_is_negative(qx) ? q25_mark(q25_negate(qx)) : qx;
    q25_ptr aqy = q25_is_negative(qy) ? q25_mark(q25_negate(qy)) : qy;
    q25_ptr denom = q25_compare(aqx, aqy) < 0 ? aqy : aqx;
    q25_ptr nqx = q25_mark(q25_negate(qx));
    q25_ptr diff = q25_mark(q25_add(nqx, qy));
    q25_ptr num = q25_is_negative(diff) ? q25_mark(q25_negate(diff)) : diff;
    /* See whether over limit */
    q25_ptr dscale = q25_scale(denom, -MAX_DIGIT_ERROR, -MAX_DIGIT_ERROR);
    if (q25_compare(num, dscale) < 0) {
#if DEBUG
	printf("  Very small error.  Return %d\n", MAX_DIGIT_ERROR);
#endif
	q25_leave(pos);
	return (double) MAX_DIGIT_ERROR;
    }
    /* Must convert here.  Risks over/underflow */
    double rnum = q25_to_double(num);
    double rdenom = q25_to_double(denom);
    q25_leave(pos);
    if (rdenom == 0 || double_is_infinite(rnum, NULL) || double_is_infinite(rdenom, NULL)) {
#if DEBUG
	printf("  Encountered problem converting to double.  Got rnum = %f, rdenom = %f\n", rnum, rdenom);
#endif
	return -1.0;
    }
    if (rnum == 0) {
	if (rdenom == 0) {
#if DEBUG
	    printf("  Zero denominator and numerator.  Return -1.0\n");
#endif
	    return -1.0;
	} else {
#if DEBUG
	    printf("  Nonzero denominator and Zero numerator.  Return %d\n", MAX_DIGIT_ERROR);
#endif
	    return MAX_DIGIT_ERROR;
	}
    }
    double div = rnum/rdenom;
    double ldiv = -log10(div);
    if (ldiv < 0.0)
	ldiv = 0.0;
    if (ldiv > MAX_DIGIT_ERROR)
	ldiv = (double) MAX_DIGIT_ERROR;
#if DEBUG
    printf("  Normal case. Num = %g, Denom = %g, Ratio = %g.  Return %.4f\n", 
	   rnum, rdenom, div, ldiv);
#endif
    return ldiv;
}

double digit_error_mix(q25_ptr qx, double y) {
    int pos = q25_enter();
    q25_ptr qy = q25_mark(q25_from_double(y));
    double result = digit_error_q25(qx, qy);
    q25_leave(pos);
    return result;
}

double digit_error(double x, double y) {
    int pos = q25_enter();
    q25_ptr qx = q25_mark(q25_from_double(x));
    q25_ptr qy = q25_mark(q25_from_double(y));
    double result = digit_error_q25(qx, qy);
    q25_leave(pos);
    return result;
}
