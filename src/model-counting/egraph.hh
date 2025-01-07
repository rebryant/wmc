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


#include "q25.h"

/* Some useful utility functions */

const char *mpf_string(mpf_srcptr val);
const char *mpfr_string(mpfr_srcptr val);

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
    std::vector<int> literals;
    std::vector<int> smoothing_variables;
};

class Egraph {
public:
    std::vector<Egraph_operation> operations;
    std::vector<Egraph_edge> edges;
    int root_id;
    std::unordered_set<int> *data_variables;
    bool is_smoothed;

    Egraph(std::unordered_set<int> *data_variables);
    void read_nnf(FILE *infile);
    void write_nnf(FILE *outfile);

    void smooth();

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
    q25_ptr evaluate(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void clear_evaluation();
    int max_size;

    
private:
    void prepare_weights(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    q25_ptr evaluate_edge(Egraph_edge &e, bool smoothed);

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
    double evaluate(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void clear_evaluation();
    
private:
    void prepare_weights(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    double evaluate_edge(Egraph_edge &e, bool smoothed);
};

/*******************************************************************************************************************
Evaluation via Gnu multi-precision floating point
*******************************************************************************************************************/

class Evaluator_mpf {
private:
    Egraph *egraph;
    // For evaluation
    std::unordered_map<int,mpf_class> evaluation_weights;
    std::unordered_map<int,mpf_class> smoothing_weights;
    mpf_class rescale;

public:

    Evaluator_mpf(Egraph *egraph);
    // literal_weights == NULL for unweighted
    bool evaluate(mpf_class &count, std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void clear_evaluation();

private:
    bool prepare_weights(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void evaluate_edge(mpf_class &value, Egraph_edge &e, bool smoothed);
};

/*******************************************************************************************************************
Evaluation via Gnu multi-precision rational arithmetic
*******************************************************************************************************************/

class Evaluator_mpq {
private:
    Egraph *egraph;
    // For evaluation
    std::unordered_map<int,mpq_class> evaluation_weights;
    std::unordered_map<int,mpq_class> smoothing_weights;
    mpq_class rescale;

public:

    Evaluator_mpq(Egraph *egraph);
    // literal_weights == NULL for unweighted
    bool evaluate(mpq_class &count, std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void clear_evaluation();
    // Maximum number of bytes in MPQ representation of any generated value
    size_t max_bytes;
    
private:
    bool prepare_weights(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void evaluate_edge(mpq_class &value, Egraph_edge &e, bool smoothed);
};

/*******************************************************************************************************************
Evaluation via MPFI interval floating point
*******************************************************************************************************************/

double mpfi_digit_precision(mpfi_srcptr v);

class Evaluator_mpfi {
private:
    Egraph *egraph;
    // For evaluation.  Use MPQ to store weights precisely
    std::unordered_map<int,mpq_class> evaluation_weights;
    std::unordered_map<int,mpq_class> smoothing_weights;
    mpfi_t rescale;
    int target_digit_precision;
    /* Use recomputation to achieve desired precision */
    bool refine;

public:

    Evaluator_mpfi(Egraph *egraph, int target_digit_precision, bool refine);
    // literal_weights == NULL for unweighted
    bool evaluate(mpfi_ptr count, std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void clear_evaluation();
    
    // Instrumentation
    double min_digit_precision;
    long precision_failure_count;

private:
    bool prepare_weights(std::unordered_map<int,const char *> *literal_string_weights, bool smoothed);
    void evaluate_edge(mpfi_ptr value, Egraph_edge &e, bool smoothed);
};


/*******************************************************************************************************************
Evaluation via Gnu multi-precision integer arithmetic (unweighted counting only)
*******************************************************************************************************************/

/** NOT IMPLEMENTED **/

class Evaluator_mpz {
private:
    Egraph *egraph;

public:

    Evaluator_mpz(Egraph *egraph);
    // Can indicate superset of those variables appearing in the NNF
    bool evaluate(mpz_class &count, std::unordered_set<int> *data_variables);
    void clear_evaluation();
    
private:
    void evaluate_edge(mpz_class &value, Egraph_edge &e);
};

