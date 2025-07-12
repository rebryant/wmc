/* Testing implementations of extended-range double */


#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>


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
   one where 0.0 has exp = 0 (erdz)
   and one where it has exp = INT64_MIN (erdm)
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

/* Get exponent as signed integer */
static int64_t dbl_get_exponent(double x) {
    uint64_t bx = dbl_get_bits(x);
    int64_t bexp = (bx >> DBL_EXP_OFFSET) & DBL_EXP_MASK;
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
    return exp <= -DBL_BIAS;
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
    Uint64_t bexp = (uint64_t) (exp + DBL_BIAS) << DBL_EXP_OFFSET;
    uint64_t bx = dbl_get_bits(x);
    uint64_t mask = ~(DBL_EXP_MASK << DBL_EXP_OFFSET);
    bx &= mask;
    bx += bexp;
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
    return {0.0, 0};
}


static erd_t erdz_normalize(erd_t a) {
    if (erdz_is_zero(a))
	return erdz_zero();
    erd_t nval;
    nval.exp = a.exp + dbl_get_exponent(a.dbl);
    nval.dbl = dbl_replace_exponent(a.dbl, 0);
    return nval;
}



erd_t erdz_from_double(double dval) {
}

/********************** ERDM *************************/

erd_t erdm_zero() {
    return {0.0, INT64_MIN};
}




erd_t erd_from_mpf(mpf_srcptr fval);
void erd_to_mpf(mpf_ptr dest, erd_t eval);
double erd_to_double(erd_t eval);

erd_t erd_negate(erd_t a);
erd_t erd_add(erd_t a, erd_t b);
erd_t erd_mul(erd_t a, erd_t b);
erd_t erd_recip(erd_t a);
int erd_cmp(erd_t a, erd_t b);

