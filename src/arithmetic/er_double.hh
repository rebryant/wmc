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

// #include <stdbool.h>
#include <stdint.h>

/*
  Properties of double-precision representation
*/
#define DBL_EXP_OFFSET 52
#define DBL_SIGN_OFFSET 63
#define DBL_EXP_MASK 0x7ff
#define DBL_MAX_PREC 54
#define DBL_BIAS 0x3ff

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

static double dbl_replace_exponent(double x, int exp) {
    int sign = dbl_get_sign(x);
    uint64_t frac = dbl_get_fraction(x);
    return dbl_assemble(sign, exp, frac);
}

class Erd {
private:
    double dbl;
    int64_t exp;
    Erd(double d, int64_t e) { 
	if (d == 0) { dbl = d; exp = 0; } 
	else { exp = e + dbl_get_exponent(d); dbl = dbl_replace_exponent(d, 0); } 
    }

public:

    Erd() { dbl = 0.0; exp = 0; }

    Erd(double d) { 
	if (d == 0) {
	    dbl = d; exp = 0; 
	} else {
	    exp = dbl_get_exponent(d);
	    dbl = dbl_replace_exponent(d, 0);
	}
    }

    Erd(mpf_srcptr mval) { long int e; dbl = mpf_get_d_2exp(&e, mval); exp = (int64_t) e; }

    bool is_zero() { return dbl == 0.0; }

    void get_mpf(mpf_ptr dest) {
	mpf_set_d(dest, dbl); 
	if (exp < 0) mpf_div_2exp(dest, dest, -exp);
	else if (exp > 0)  mpf_mul_2exp(dest, dest, exp); 
    }

    Erd add(Erd other) { 
	if (other.is_zero()) return *this;
	if (is_zero()) return other; 
	if (exp - other.exp > DBL_MAX_PREC) return *this;
	if (other.exp - exp > DBL_MAX_PREC) return other;
	int ediff = exp - other.exp;
	double da = dbl_replace_exponent(dbl, ediff);
	return Erd(da + other.dbl, other.exp); 
    }

    Erd mul(Erd other) {
	double prod = dbl * other.dbl;
	int64_t e = exp + other.exp;
	return Erd(prod, e);
    }

#define PROD_RENORM_COUNT 100

    friend Erd product_reduce_x1(std::vector<Erd> &arguments) {
	double dp = 1.0;
	int64_t ep = 0.0;
	int rcount = 0;
	for (int i = 0; i < arguments.size(); i++) {
	    dp *= arguments[i].dbl;
	    ep += arguments[i].exp;
	    if (++rcount >= PROD_RENORM_COUNT) {
		Erd mid = Erd(dp, ep);
		dp = mid.dbl;
		ep = mid.exp;
		rcount = 0;
	    }

	}
	return Erd(dp, ep);
    }

    friend Erd product_reduce_x2(std::vector<Erd> &arguments) {
	double dp[2] = { 1.0, 1.0 };
	int64_t ep[2] = { 0, 0 };
	int rcount = 0;
	int i = 0;
	int len = (int) arguments.size();
	for (i = 0; i <= len - 2; i+= 2) {
	    for (int j = 0; j < 2; j++) {
		dp[j] *= arguments[i+j].dbl;
		ep[j] += arguments[i+j].exp;
	    }
	    if (++rcount >= PROD_RENORM_COUNT) {
		for (int j = 0; j < 2; j++) {
		    Erd mid = Erd(dp[j], ep[j]);
		    dp[j] = mid.dbl;
		    ep[j] = mid.exp;
		}
		rcount = 0;
	    }
	}
	for (; i < len; i++) {
	    dp[0] *= arguments[i].dbl;
	    ep[0] += arguments[i].exp;
	}
	for (int j = 1; j < 2; j++) {
	    dp[0] *= dp[j];
	    ep[0] += ep[j];
	}
	return Erd(dp[0], ep[0]);
    }

    friend Erd product_reduce_x4(std::vector<Erd> &arguments) {
	double dp[4] = { 1.0, 1.0, 1.0, 1.0 };
	int64_t ep[4] = { 0, 0, 0, 0 };
	int rcount = 0;
	int i = 0;
	int len = (int) arguments.size();
	for (i = 0; i <= len - 4; i+= 4) {
	    for (int j = 0; j < 4; j++) {
		dp[j] *= arguments[i+j].dbl;
		ep[j] += arguments[i+j].exp;
	    }
	    if (++rcount >= PROD_RENORM_COUNT) {
		for (int j = 0; j < 4; j++) {
		    Erd mid = Erd(dp[j], ep[j]);
		    dp[j] = mid.dbl;
		    ep[j] = mid.exp;
		}
		rcount = 0;
	    }
	}
	for (; i < len; i++) {
	    dp[0] *= arguments[i].dbl;
	    ep[0] += arguments[i].exp;
	}
	for (int j = 1; j < 4; j++) {
	    dp[0] *= dp[j];
	    ep[0] += ep[j];
	}
	return Erd(dp[0], ep[0]);
    }



    friend Erd product_reduce_slow(std::vector<Erd> &arguments) {
	Erd result(1.0);
	for (Erd arg : arguments)
	    result = result.mul(arg);
	return result;
    }

    friend Erd product_reduce(std::vector<Erd> &arguments) {
	if (arguments.size() > 8)
	    return product_reduce_x2(arguments);
	else
	    return product_reduce_x1(arguments);
    }

};

