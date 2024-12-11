/*========================================================================
  Copyright (c) 2023 Randal E. Bryant, Carnegie Mellon University
  
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

// This code was copied from an earlier knowledge compiler.
// Not all of the features implemented get used in this context.

#pragma once

#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "q25.h"

#define TAUTOLOGY INT_MAX
#define CONFLICT (-TAUTOLOGY)
// Used to convert literal to variable
#define IABS(x) ((x)<0?-(x):(x))
#define IMIN(x,y) ((x)<(y)?(x):(y))

class Cnf {
private:
    
    int nvar;

    // For each clause, its starting index into the literal_sequence array
    // A final index beyond the last clause to make it easy to compute clause lengths
    std::vector<int> clause_offset;
    // Clause literals combined into single sequence
    std::vector<int> literal_sequence;

public:
    Cnf();
    Cnf(int icount);

    bool import_file(FILE *infile, bool process_comments, bool skip_clauses);

    ~Cnf();

    void initialize(int icount);
    // Delete sets for data variables and input weights
    void deallocate();

    int variable_count() { return nvar; }
    size_t clause_count() { return clause_offset.size() - 1; }

    int maximum_clause_id();
    int clause_length(int id);
    // Index literals from 0
    int get_literal(int cid, int lid);

    // Both of these return true if successful
    bool show(FILE *outfile);
    bool write(FILE *outfile, bool show_data_variables, bool show_forget_variables, bool show_weights);

    // Determine whether the file calls for weighted counting
    bool is_weighted() { return !input_weights || input_weights->size() != 0; }

    // Add new clauses one literal at a time
    // Returns clause ID
    int new_clause();
    void add_literal(int lit);
    void finish();

    // Public access to extra information
    std::unordered_set<int> *data_variables;
    // Variables that were detected to have Forget property during preprocessing
    std::unordered_set<int> *forget_variables;
    // Optional weights of data variables
    std::unordered_map<int,q25_ptr> *input_weights;

    bool is_data_variable(int var) { return data_variables->find(var) != data_variables->end(); }
    bool is_forget_variable(int var) { return forget_variables->find(var) != forget_variables->end(); }

};

