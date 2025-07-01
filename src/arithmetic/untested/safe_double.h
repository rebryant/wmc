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

#pragma once


#include <gmp.h>

/* Allow this headerfile to define C++ constructs if requested */
#ifdef __cplusplus
#define CPLUSPLUS
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif


typedef enum { SD_MPF, SD_DOUBLE, SD_ZERO } sd_type_t;

/*
  Representation of floating-point numbers where double provides sufficient precision, 
  but perhaps not enough range.  
  Use double when possible and switch to mpf (64-bit precision) when need more range
 */
typedef struct {
    double dval;  // Error if zero, denormalized, infinite, or NaN
    mpf_t  fval;  // Allocated only when used
    sd_type_t type;
} safe_double_t, *safe_double_ptr;

void safe_double_clear(safe_double_ptr sd);

/* Convert from/to double-precision FP.  Assume IEEE representation */
void safe_double_zero(safe_double_ptr sd);
void safe_double_from_double(safe_double_ptr sd, double val);
void safe_double_from_mpf(safe_double_ptr sd, mpf_srcptr val);
void safe_double_from_mpq(safe_double_ptr sd, mpq_srcptr val);
void safe_double_from_sd(safe_double_ptr sd, safe_double_ptr val);
void safe_double_to_mpf(mpf_ptr dest, safe_double_ptr sd);
void safe_double_negate(safe_double_ptr sd);
void safe_double_add_sd(safe_double_ptr sd, safe_double_ptr val);
void safe_double_mul_sd(safe_double_ptr sd, safe_double_ptr val);
void safe_double_recip_sd(safe_double_ptr sd);
void safe_double_add_double(safe_double_ptr sd, double val);
void safe_double_mul_double(safe_double_ptr sd, double val);
void safe_double_add_mpf(safe_double_ptr sd, mpf_srcptr val);
void safe_double_mul_mpf(safe_double_ptr sd, mpf_srcptr val);

#ifdef CPLUSPLUS
}
#endif
