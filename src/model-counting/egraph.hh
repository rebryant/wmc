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
#include <vector>
#include <unordered_set>


#include "q25.h"

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
    q25_ptr evaluate(std::unordered_map<int,q25_ptr> *literal_weights, bool smoothed);
    void clear_evaluation();
    
private:
    void prepare_weights(std::unordered_map<int,q25_ptr> *literal_weights, bool smoothed);
    q25_ptr evaluate_edge(Egraph_edge &e, bool smoothed);

};

// Support for evaluation
bool read_cnf(FILE *infile, std::unordered_set<int> **data_variables, std::unordered_map<int,q25_ptr> **weights);


