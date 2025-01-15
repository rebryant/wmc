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

double digit_precision_q25(q25_ptr x_est, q25_ptr x) {
    if (!q25_is_valid(x_est) || q25_is_infinite(x_est, NULL)
	|| !q25_is_valid(x) || q25_is_infinite(x, NULL))
	return 0.0;
    if (q25_compare(x_est, x) == 0)
	return MAX_DIGIT_PRECISION;
    int pos = q25_enter();
    q25_ptr num, denom;
    if (q25_is_zero(x)) {
	denom = q25_mark(q25_from_32(1));
	num = q25_mark(q25_abs(x_est));
	if (q25_compare(num, denom) >= 0)
	    num = denom;
    } else {
	denom = q25_mark(q25_abs(x));
	q25_ptr x_neg = q25_mark(q25_negate(x));
	num = q25_mark(q25_add(x_est, x_neg));
	q25_inplace_abs(num);
    }
    /* See whether over limit */
    q25_ptr dscale = q25_scale(denom, -MAX_DIGIT_PRECISION, -MAX_DIGIT_PRECISION);
    if (q25_compare(num, dscale) < 0) {
	q25_leave(pos);
	return (double) MAX_DIGIT_PRECISION;
    }
    /* Shift within reasonable range */
    int mag_num = q25_magnitude(num);
    int mag_denom = q25_magnitude(denom);
    int scale10 = 0;
    if (mag_num > 0 && mag_denom > 0)
	scale10 = -(mag_num > mag_denom ? mag_denom : mag_num);
    else if (mag_num < 0 && mag_denom < 0)
	scale10 = -(mag_num < mag_denom ? mag_denom : mag_num);
    q25_inplace_scale(num, scale10, scale10);
    q25_inplace_scale(denom, scale10, scale10);
    /* Must convert here.  Risks over/underflow */
    double rnum = q25_to_double(num);
    double rdenom = q25_to_double(denom);

    q25_leave(pos);
    if (rdenom == 0 || double_is_infinite(rnum, NULL))
	return 0.0;
    if (double_is_infinite(rdenom, NULL))
	return MAX_DIGIT_PRECISION;
    if (rnum == 0) {
	if (rdenom == 0)
	    return 0.0;
	else
	    return MAX_DIGIT_PRECISION;
    }
    double div = rnum/rdenom;
    double ldiv = -log10(div);
    if (ldiv < 0.0)
	ldiv = 0.0;
    if (ldiv > MAX_DIGIT_PRECISION)
	ldiv = (double) MAX_DIGIT_PRECISION;
    return ldiv;
}

double digit_precision_mix(q25_ptr x_est, double x) {
    int pos = q25_enter();
    q25_ptr qx = q25_mark(q25_from_double(x));
    double result = digit_precision_q25(x_est, qx);
    q25_leave(pos);
    return result;
}

double digit_precision(double x_est, double x) {
    int pos = q25_enter();
    q25_ptr qx = q25_mark(q25_from_double(x));
    q25_ptr qx_est = q25_mark(q25_from_double(x_est));
    double result = digit_precision_q25(qx, qx_est);
    q25_leave(pos);
    return result;
}
