/*========================================================================
  Copyright (c) 2024 Randal E. Bryant, Carnegie Mellon University
  
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

// Construct evaluation graph (egraph) from DNNF file
// Support (weighted) model counting

#pragma once

#include <cstdio>
#include <cstdlib>

#include <stdint.h>
#include <stdarg.h>

#include <vector>
#include <unordered_set>

#include <gmp.h>
#include <gmpxx.h>

#include <mpfr.h>
#include <mpfi.h>

#include "er_double.hh"
#include "q25.h"

/* Some useful utility functions */

/*
  How many digits of precision can we guarantee when all weights are nonnegative?
  Compute floats with specified bit precision over formula with specified number of variables
  constant = 4 for smoothed evaluation and 7 for unsmoothed
*/
double digit_precision_bound(int bit_precision, int nvar, double constant);

/*
  How many bits of floating-point precision are required to achieve
  target digit precision when all weights are nonnegative?
  constant = 4 for smoothed evaluation and 7 for unsmoothed
 */
int required_bit_precision(double target_precision, int nvar, double constant, bool nonnegative);


const char *mpf_string(mpf_srcptr val, int digits);
const char *mpfr_string(mpfr_srcptr val, int digits);

double digit_precision_mpfr(mpfr_srcptr x_est, mpq_srcptr x);
double digit_precision_mpf(mpf_srcptr x_est, mpq_srcptr x);
double digit_precision_d(double x_est, mpq_srcptr x);
double digit_precision_mpfi(mpfi_srcptr v);

/*******************************************************************************************************************
 Graph representing NNF formula
*******************************************************************************************************************/

typedef enum { NNF_NONE, NNF_TRUE, NNF_FALSE, NNF_AND, NNF_OR, NNF_NUM } nnf_type_t;

struct Egraph_operation {
    int indegree;
    nnf_type_t type;
};

struct Egraph_edge {
    int from_id;
    int to_id;
    // Disabled due to zero weight or smoothing value
    bool has_zero;
    std::vector<int> literals;
    std::vector<int> smoothing_variables;
};

struct Egraph_weights {
    std::unordered_map<int,mpq_class> evaluation_weights;
    std::unordered_map<int,mpq_class> smoothing_weights;
    std::vector<mpq_class> rescale_weights;
    bool all_nonnegative;
};

class Egraph {
public:
    std::vector<Egraph_operation> operations;
    std::vector<Egraph_edge> edges;
    int root_id;
    std::unordered_set<int> *data_variables;
    bool is_smoothed;
    int smooth_variable_count;
    int disabled_edge_count;
    // Count of variables in original formula, including those eliminated by projection
    int nvar;

    Egraph(std::unordered_set<int> *data_variables, int nvar);
    void read_nnf(FILE *infile);
    void write_nnf(FILE *outfile);

    Egraph_weights *prepare_weights(std::unordered_map<int,const char *> *literal_string_weights);
    void smooth();

    // Ability to do partial smoothing
    // Remove all smoothing variables
    void reset_smooth();
    // Put in single smoothing variable.  If is_zero, then disable edge
    void smooth_single(int var, bool is_zero);

    bool is_data_variable(int var) { return data_variables->find(var) != data_variables->end(); }
    bool is_literal(int lit) { return lit < 0 ? is_data_variable(-lit) : is_data_variable(lit); }
    bool is_operation(int id) { return id > 0 && id <= operations.size(); }

    void add_operation(int id, nnf_type_t type);
    int add_edge(int from_id, int to_id);
    void add_edge_literal(int eid, int lit);
    void add_smoothing_variable(int eid, int var);
};

/*******************************************************************************************************************
Evaluation via Q25
*******************************************************************************************************************/

class Evaluator_q25 {
private:
    Egraph *egraph;
    // For evaluation
    std::unordered_map<int,q25_ptr> evaluation_weights;
    std::unordered_map<int,q25_ptr> smoothing_weights;
    q25_ptr rescale;

public:

    Evaluator_q25(Egraph *egraph);
    // literal_weights == NULL for unweighted
    q25_ptr evaluate(std::unordered_map<int,const char *> *literal_string_weights);
    void clear_evaluation();
    int max_size;

    
private:
    void prepare_weights(std::unordered_map<int,const char *> *literal_string_weights);
    q25_ptr evaluate_edge(Egraph_edge &e);

};

/*******************************************************************************************************************
Evaluation via double-precision floating point
*******************************************************************************************************************/

class Evaluator_double {
private:
    Egraph *egraph;
    // For evaluation
    std::unordered_map<int,double> evaluation_weights;
    std::unordered_map<int,double> smoothing_weights;
    double rescale;

public:

    Evaluator_double(Egraph *egraph);
    // literal_weights == NULL for unweighted
    double evaluate(std::unordered_map<int,const char *> *literal_string_weights);
    void clear_evaluation();
    
private:
    void prepare_weights(std::unordered_map<int,const char *> *literal_string_weights);
    double evaluate_edge(Egraph_edge &e);
};

/*******************************************************************************************************************
Evaluation via extended-range double-precision.  Use MPF as way to get weights in and out
*******************************************************************************************************************/

class Evaluator_erd {
private:
    Egraph *egraph;
    // For evaluation
    Egraph_weights *weights;
    Erd rescale;
    // Used for product computations
    std::vector<Erd> arguments;
    


public:

    Evaluator_erd(Egraph *egraph, Egraph_weights *weights);
    // literal_weights == NULL for unweighted
    void evaluate(mpf_class &count);
    void clear_evaluation();

private:
    Erd evaluate_edge(Egraph_edge &e);
};


/*******************************************************************************************************************
Evaluation via Gnu multi-precision floating point
*******************************************************************************************************************/

class Evaluator_mpf {
private:
    Egraph *egraph;
    // For evaluation
    Egraph_weights *weights;
    mpf_class rescale;

public:

    Evaluator_mpf(Egraph *egraph, Egraph_weights *weights);
    // literal_weights == NULL for unweighted
    void evaluate(mpf_class &count);
    void clear_evaluation();

private:
    void evaluate_edge(mpf_class &value, Egraph_edge &e);
};

/*******************************************************************************************************************
Evaluation via Gnu multi-precision rational arithmetic
*******************************************************************************************************************/

class Evaluator_mpq {
private:
    Egraph *egraph;
    // For evaluation
    Egraph_weights *weights;
    mpq_class rescale;

public:

    Evaluator_mpq(Egraph *egraph, Egraph_weights *weights);
    // literal_weights == NULL for unweighted
    void evaluate(mpq_class &count);
    void clear_evaluation();
    // Maximum number of bytes in MPQ representation of any generated value
    size_t max_bytes;
    
private:
    void evaluate_edge(mpq_class &value, Egraph_edge &e);
};

/*******************************************************************************************************************
Evaluation via MPFI interval floating point
*******************************************************************************************************************/

class Evaluator_mpfi {
private:
    Egraph *egraph;
    // For evaluation.  Use MPQ to store weights precisely
    Egraph_weights *weights;
    mpfi_t rescale;
    // Measure precision of intermdiate results
    bool instrument;

public:

    Evaluator_mpfi(Egraph *egraph, Egraph_weights *weights, bool instrument);
    void evaluate(mpfi_ptr count);
    void clear_evaluation();
    // Least digit precision estimate encountered.  Only computed when instrument.
    double min_digit_precision;

private:
    void evaluate_edge(mpfi_ptr value, Egraph_edge &e);
};

/*******************************************************************************************************************
Evaluation.  When no negative weights, use MPI.  Otherwise, start with MFPI and switch to MPQ if needed
*******************************************************************************************************************/

typedef enum { COMPUTE_ERD, COMPUTE_MPF, COMPUTE_MPFI, COMPUTE_MPQ, 
	       COMPUTE_ERD_NOMPQ, COMPUTE_MPF_NOMPQ, COMPUTE_MPFI_NOMPQ, COMPUTE_MPQ_NOMPQ } computed_t;

class Evaluator_combo {
private:

    Egraph *egraph;
    Egraph_weights *weights;
    double target_precision;
    int bit_precision;
    int instrument;

public:

    Evaluator_combo(Egraph *egraph, Egraph_weights *weights, double target_precision, int bit_precision, int instrument);
    // literal_weights == NULL for unweighted
    void evaluate(mpf_class &count, bool no_mpq);

    computed_t computed_method;
    const char *method();
    double guaranteed_precision;
    size_t max_bytes;
    int used_bit_precision() { return bit_precision; }
    // Leftover stuff that can be reused
    // Times for different evaluations.  Set to 0.0 if not used
    double erd_seconds;
    double mpf_seconds;
    double mpfi_seconds;
    double mpq_seconds;
    mpq_class mpq_count;
    mpf_class mpf_count;
    mpf_class erd_count;
    mpfi_t mpfi_count;
    double min_digit_precision;
};

