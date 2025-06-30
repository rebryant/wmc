#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#include "safe_double.h"

static bool bad_exponent(double val) {
    union {
	double   d;
	uint64_t b;
    } u;
    u.d = val;
    unsigned exponent = (u.b >> 52) & 0x1FFF;
    return exponent == 0 || exponent == 0x1fff;
}

static void force_mpf(safe_double_ptr sd) {
    switch(sd->type) {
    case SD_MPF:
	break;
    case SD_DOUBLE:
    case SD_ZERO:
	mpf_init(sd->fval);
	mpf_set_d(sd->fval, sd->type == SD_ZERO ? 0.0 : sd->dval);
	break;
    default:
	break;
    }
    sd->type = SD_MPF;

}

void safe_double_clear(safe_double_ptr sd) {
    if (sd->type == SD_MPF)
	mpf_clear(sd->fval);
}

void safe_double_zero(safe_double_ptr sd) {
    sd->type = SD_ZERO;
}


void safe_double_from_double(safe_double_ptr sd, double val) {
    sd->dval = val;
    sd->type = SD_DOUBLE;
}

void safe_double_from_mpf(safe_double_ptr sd, mpf_srcptr val) {
    mpf_init(sd->fval);
    mpf_set(sd->fval, val);
    sd->type = SD_MPF;
}

void safe_double_from_mpq(safe_double_ptr sd, mpq_srcptr val) {
    static bool initialized = false;
    static mpq_t mpq_zero;
    if (!initialized) {
	mpq_init(mpq_zero);
	mpq_set_d(mpq_zero, 0.0);
	initialized = true;
    }
    if (mpq_equal(mpq_zero, val)) {
	safe_double_zero(sd);
	return;
    }
    double d = mpq_get_d(val);
    if (bad_exponent(d)) {
	mpf_init(sd->fval);
	mpf_set_q(sd->fval, val);
	sd->type = SD_MPF;
    } else {
	sd->dval = d;
	sd->type = SD_DOUBLE;
    }
}

void safe_double_from_sd(safe_double_ptr sd, safe_double_ptr val) {
    switch (val->type) {
    case SD_ZERO:
	safe_double_zero(sd);
	break;
    case SD_DOUBLE:
	safe_double_from_double(sd, val->dval);
	break;
    case SD_MPF:
	safe_double_from_mpf(sd, val->fval);
	break;
    default:
	break;
    }
}


void safe_double_to_mpf(mpf_ptr dest, safe_double_ptr sd) {
    switch(sd->type) {
    case SD_ZERO:
    case SD_DOUBLE:
	mpf_set_d(dest, sd->type == SD_ZERO ? 0.0 : sd->dval);
	break;
    case SD_MPF:
	mpf_set(dest, sd->fval);
	break;
    default:
	break;
    }
}

void safe_double_negate(safe_double_ptr sd) {
    switch (sd->type) {
    case SD_ZERO:
	break;
    case SD_DOUBLE:
	sd->dval = -sd->dval;
	break;
    case SD_MPF:
	mpf_neg(sd->fval, sd->fval);
	break;
    }
}

void safe_double_add_sd(safe_double_ptr sd, safe_double_ptr val) {
    double sum;
    mpf_t fnew;
    switch (sd->type) {
    case SD_ZERO:
	safe_double_from_sd(sd, val);
	break;
    case SD_DOUBLE:
	switch(val->type) {
	case SD_ZERO:
	    break;
	case SD_DOUBLE:
	    sum = sd->dval + val->dval;
	    if (bad_exponent(sum)) {
		force_mpf(sd);
		mpf_init(fnew);
		mpf_set_d(fnew, val->dval);
		mpf_add(sd->fval, sd->fval, fnew);
		mpf_clear(fnew);
	    } else
		sd->dval = sum;
	    break;
	case SD_MPF:
	    force_mpf(sd);
	    mpf_add(sd->fval, sd->fval, val->fval);
	    break;
	default:
	    break;
	}
	break;
    case SD_MPF:
	switch(val->type) {
	case SD_ZERO:
	    break;
	case SD_DOUBLE:
	    mpf_init(fnew);
	    mpf_set_d(fnew, val->dval);
	    mpf_add(sd->fval, sd->fval, fnew);
	    mpf_clear(fnew);
	    break;
	case SD_MPF:
	    mpf_add(sd->fval, sd->fval, val->fval);
	    break;
	default:
	    break;
	}
	break;
    default:
	break;
    }
}

void safe_double_mul_sd(safe_double_ptr sd, safe_double_ptr val) {
    double prod;
    mpf_t fnew;
    switch (sd->type) {
    case SD_ZERO:
	break;
    case SD_DOUBLE:
	switch(val->type) {
	case SD_ZERO:
	    safe_double_zero(sd);
	    break;
	case SD_DOUBLE:
	    prod = sd->dval * val->dval;
	    if (bad_exponent(prod)) {
		force_mpf(sd);
		mpf_init(fnew);
		mpf_set_d(fnew, val->dval);
		mpf_mul(sd->fval, sd->fval, fnew);
		mpf_clear(fnew);
	    } else
		sd->dval = prod;
	    break;
	case SD_MPF:
	    force_mpf(sd);
	    mpf_mul(sd->fval, sd->fval, val->fval);
	    break;
	default:
	    break;
	}
	break;
    case SD_MPF:
	switch(val->type) {
	case SD_ZERO:
	    mpf_clear(sd->fval);
	    safe_double_zero(sd);
	    break;
	case SD_DOUBLE:
	    mpf_init(fnew);
	    mpf_set_d(fnew, val->dval);
	    mpf_mul(sd->fval, sd->fval, fnew);
	    mpf_clear(fnew);
	    break;
	case SD_MPF:
	    mpf_mul(sd->fval, sd->fval, val->fval);
	    break;
	default:
	    break;
	}
	break;
    default:
	break;
    }
}

void safe_double_recip_sd(safe_double_ptr sd) {
    double recip;
    mpf_t fnew;
    switch (sd->type) {
    case SD_ZERO:
	/* Generate invalid value */
	sd->dval = 0.0;
	sd->type = SD_DOUBLE;
	break;
    case SD_DOUBLE:
	recip = 1.0/sd->dval;
	if (bad_exponent(recip)) {
	    mpf_init(fnew);
	    mpf_init(sd->fval);
	    mpf_set_d(sd->fval, 1.0);
	    mpf_set_d(fnew, sd->dval);
	    mpf_div(sd->fval, sd->fval, fnew);
	    mpf_clear(fnew);
	    sd->type = SD_MPF;
	} else 
	    sd->dval = recip;
	break;
    case SD_MPF:
	mpf_init(fnew);
	mpf_set_d(fnew, 1.0);
	mpf_div(fnew, fnew, sd->fval);
	mpf_swap(fnew, sd->fval);
	mpf_clear(fnew);
	break;
    default:
	break;
    }
}

void safe_double_add_double(safe_double_ptr sd, double val) {
    double sum;
    mpf_t fnew;
    switch (sd->type) {
    case SD_ZERO:
	sd->dval = val;
	sd->type = SD_DOUBLE;
	break;
    case SD_DOUBLE:
	sum = sd->dval + val;
	if (bad_exponent(sum)) {
	    force_mpf(sd);
	    mpf_init(fnew);
	    mpf_set_d(fnew, val);
	    mpf_add(sd->fval, sd->fval, fnew);
	    mpf_clear(fnew);
	} else 
	    sd->dval = sum;
	break;
    case SD_MPF:
	mpf_init(fnew);
	mpf_set_d(fnew, val);
	mpf_add(sd->fval, sd->fval, fnew);
	mpf_clear(fnew);
	break;
    default:
	break;
    }
}

void safe_double_mul_double(safe_double_ptr sd, double val) {
    double prod;
    mpf_t fnew;
    switch (sd->type) {
    case SD_ZERO:
	break;
    case SD_DOUBLE:
	prod = sd->dval * val;
	if (bad_exponent(prod)) {
	    force_mpf(sd);
	    mpf_init(fnew);
	    mpf_set_d(fnew, val);
	    mpf_mul(sd->fval, sd->fval, fnew);
	    mpf_clear(fnew);
	} else 
	    sd->dval = prod;
	break;
    case SD_MPF:
	mpf_init(fnew);
	mpf_set_d(fnew, val);
	mpf_mul(sd->fval, sd->fval, fnew);
	mpf_clear(fnew);
	break;
    default:
	break;
    }
}


void safe_double_add_mpf(safe_double_ptr sd, mpf_srcptr val) {
    switch (sd->type) {
    case SD_ZERO:
	safe_double_from_mpf(sd, val);
	break;
    case SD_DOUBLE:
	force_mpf(sd);
	mpf_add(sd->fval, sd->fval, val);
	break;
    case SD_MPF:
	mpf_add(sd->fval, sd->fval, val);
	break;
    default:
	break;
    }
}


void safe_double_mul_mpf(safe_double_ptr sd, mpf_srcptr val) {
    switch (sd->type) {
    case SD_ZERO:
	break;
    case SD_DOUBLE:
	force_mpf(sd);
	mpf_mul(sd->fval, sd->fval, val);
	break;
    case SD_MPF:
	mpf_mul(sd->fval, sd->fval, val);
	break;
    default:
	break;
    }
}


