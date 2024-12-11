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

// Read .nnf file describing d-DNNF formula.
// Compute unweighted and (optionally) weighted  model counts.
// Weights come from CNF file

#include <unistd.h>
#include <cstring>
#include <ctype.h>

#include "report.h"
#include "counters.h"
#include "egraph.hh"
#include "cnf_info.hh"

static const char *nnf_type_name[NNF_NUM] = { "NONE", "TRUE", "FALSE", "AND", "OR" };
static const char nnf_type_char[NNF_NUM] = { '\0', 't', 'f', 'a', 'o' };

// Current line number
static int line_number = 0;

Egraph::Egraph(std::unordered_set<int> *dvars) {
    data_variables = dvars;
    is_smoothed = false;
}

void Egraph::add_operation(int id, nnf_type_t type) {
    if (id > operations.size())
	operations.resize(id);
    operations[id-1].indegree = 0;
    if (operations[id-1].type != NNF_NONE)
	err(true, "Line %d.  Operation %d already defined\n", line_number, id);
    operations[id-1].type = type;
    switch(type) {
    case NNF_NONE:
	err(true, "Line %d.  Operation %d declared to have no type\n", line_number, id);
	break;
    case NNF_TRUE:
    case NNF_AND:
    case NNF_FALSE:
    case NNF_OR:
	break;
    default:
	err(true, "Line %d.  Operation %d declared with unknown type %d\n", line_number, id, (int) type);
	break;
    }
    incr_count(COUNT_OPERATIONS);
}

int Egraph::add_edge(int from_id, int to_id) {
    size_t eid = edges.size()+1;
    edges.resize(eid);
    edges[eid-1].from_id = from_id;
    edges[eid-1].to_id = to_id;
    root_id = to_id;
    incr_count(COUNT_EDGES);
    operations[to_id-1].indegree++;
    return eid;
}

void Egraph::add_edge_literal(int eid, int lit) {
    if (!is_literal(lit))
	err(true, "Line %d.  Attempt to add invalid literal %d to edge %d\n", line_number, lit, eid);
    edges[eid-1].literals.push_back(lit);
}

void Egraph::add_smoothing_variable(int eid, int var) {
    if (!is_data_variable(var))
	err(true, "Line %d.  Attempt to add invalid smoothing variable %d to edge %d\n", line_number, var, eid);
    edges[eid-1].smoothing_variables.push_back(var);
    incr_count(COUNT_SMOOTH_VARIABLES);
}


// Try to read single alphabetic character from line
// If not found, then push back unread character and return 0
// If hit EOF, then return this
static int get_token(FILE *infile) {
    int c;
    while (true) {
	c = fgetc(infile);
	if (isalpha(c) || c == EOF)
	    return c;
	else if (isspace(c))
	    continue;
	else {
	    ungetc(c, infile);
	    return 0;
	}
    }
}

// Read sequence of numbers from line of input
// Consume end of line character
// Return false if non-numeric value encountered
static bool read_numbers(FILE *infile, std::vector<int> &vec, int *rc) {
    vec.clear();
    while (true) {
	int c = fgetc(infile);
	int val;
	if (c == '\n' || c == EOF) {
	    *rc = c;
	    return true;
	} else if (isspace(c))
	    continue;
	else {
	    ungetc(c, infile);
	    if (fscanf(infile, "%d", &val) == 1) {
		vec.push_back(val);
	    } else
		return false;
	}
    }
    // Won't hit this
    return false;
}

void Egraph::read_nnf(FILE *infile) {
    operations.clear();
    line_number = 0;
    // Capture arguments for each line
    std::vector<int> largs;

    while (true) {
	nnf_type_t type = NNF_NONE;
	line_number++;
	int c = get_token(infile);
	int rc = 0;
	if (c == EOF)
	    break;

	if (c != 0) {
	    // Operation
	    for (int t = NNF_TRUE; t < NNF_NUM; t++)
		if (c == nnf_type_char[t]) {
		    type = (nnf_type_t) t;
		    break;
		}
	    if (type == NNF_NONE)
		err(true, "Line %d.  Unknown NNF command '%c'\n", line_number, c);
	    bool ok = read_numbers(infile, largs, &rc);
	    if (!ok)
		err(true, "Line %d.  Couldn't parse numbers\n", line_number);
	    else if (largs.size() != 2)
		err(true, "Line %d.  Expected 2 numbers.  Found %d\n", line_number, largs.size());
	    else if (largs.back() != 0)
		err(true, "Line %d.  Line not zero-terminated\n", line_number);
	    else {
		int id = largs[0];
		add_operation(id, type);
		report(4, "Line %d.  Created NNF operation %s.  Id %d\n",
		       line_number, nnf_type_name[type], id);
	    }
	} else {
	    // Edge
	    bool ok = read_numbers(infile, largs, &rc);
	    if (!ok)
		err(true, "Line %d.  Couldn't parse numbers\n", line_number);
	    else if (largs.size() == 0 && rc == EOF)
		break;
	    else if (largs.size() < 3)
		err(true, "Line %d.  Expected at least 3 numbers.  Found %d\n", line_number, largs.size());
	    else if (largs.back() != 0)
		err(true, "Line %d.  Line not zero-terminated\n", line_number);
	    int to_id = largs[0];
	    int from_id = largs[1];
	    int eid = add_edge(from_id, to_id);
	    // Add literals
	    int pos;
	    int lcount = 0;
	    int scount = 0;
	    for (pos = 2; largs[pos] != 0; pos++) {
		add_edge_literal(eid, largs[pos]);
		lcount++;
	    }
	    // Add smoothing variables
	    for (++pos; pos < largs.size()-1; pos++) {
		add_smoothing_variable(eid, largs[pos]);
		scount++;
	    }
	    incr_histo(HISTO_EDGE_PRODUCTS, lcount);
	    if (scount > 0) {
		report(4, "Added edge #%d %d <-- %d.  %d literals, %d smoothing variables\n", eid, to_id, from_id, lcount, scount);
		incr_histo(HISTO_EDGE_SMOOTHS, scount);
		is_smoothed = true;
	    }  else
		report(4, "Added edge #%d %d <-- %d.  %d literals\n", eid, to_id, from_id, lcount);
	}
    }
    // Check over the operations
    for (int id = 1; id <= operations.size(); id++) {
	if (operations[id-1].indegree > 1) {
	    if (operations[id-1].type == NNF_AND)
		incr_histo(HISTO_NODE_PRODUCTS, operations[id-1].indegree-1);
	    else 
		incr_histo(HISTO_SUMS, operations[id-1].indegree-1);
	}
    }
}

void Egraph::write_nnf(FILE *outfile) {
    for (int id = 1; id <= operations.size(); id++) {
	if (operations[id-1].type != NNF_NONE) {
	    fprintf(outfile, "%c %d 0\n", nnf_type_char[operations[id-1].type], id);
	}
    }
    for (int id = 1; id <= edges.size(); id++) {
	fprintf(outfile, "%d %d", edges[id-1].to_id, edges[id-1].from_id);
	for (int lit : edges[id-1].literals)
	    fprintf(outfile, " %d", lit);
	if (edges[id-1].smoothing_variables.size() > 0) {
	    fprintf(outfile, " 0");
	    for (int var: edges[id-1].smoothing_variables)
		fprintf(outfile, " %d", var);
	}
	fprintf(outfile, " 0\n");
    }
    fclose(outfile);
}

void Egraph::smooth() {   
    if (is_smoothed)
	return;
    // For each operation, union of all variables on which it depends
    std::vector<std::unordered_set<int>> operation_dependencies;
    operation_dependencies.resize(operations.size());
    // For each edge, union of all variables for edge
    std::vector<std::unordered_set<int>> edge_variables;
    edge_variables.resize(edges.size());
 
    for (int id = 1; id <= edges.size(); id++) {
	int from_id = edges[id-1].from_id;
	int to_id = edges[id-1].to_id;
	if (operations[from_id-1].type == NNF_FALSE)
	    continue;
	for (int v : operation_dependencies[from_id-1])
	    operation_dependencies[to_id-1].insert(v);
	for (int lit : edges[id-1].literals) {
	    int v = IABS(lit);
	    operation_dependencies[to_id-1].insert(v);
	    edge_variables[id-1].insert(v);
	}
	for (int v : edges[id-1].smoothing_variables)
	    edge_variables[id-1].insert(v);
    }
    for (int id = 1; id <= edges.size(); id++) {
	int from_id = edges[id-1].from_id;
	int to_id = edges[id-1].to_id;
	if (operations[from_id-1].type == NNF_FALSE)
	    continue;
	if (operations[to_id-1].type == NNF_AND)
	    continue;
	int scount = 0;
	for (int v : operation_dependencies[to_id-1]) {
	    if (edge_variables[id-1].find(v) == edge_variables[id-1].end() &&
		operation_dependencies[from_id-1].find(v) == operation_dependencies[from_id-1].end()) {
		report(4, "Adding smoothing variable %d on edge #%d (%d <-- %d)\n", v, id, to_id, from_id);
		add_smoothing_variable(id, v);
		scount++;
	    }
	}
	if (scount > 0)
	    incr_histo(HISTO_EDGE_SMOOTHS, scount);
    }
    // Add final variables to root
    int id = edges.size();
    int child_id = edges[id-1].to_id;
    int scount = 0;
    for (int v : *data_variables) {
	if (operation_dependencies[child_id-1].find(v) == operation_dependencies[child_id-1].end()) {
	    report(4, "Adding smoothing variable %d on root edge #%d (%d --> %d)\n", v, id, root_id, child_id);
	    add_smoothing_variable(id, v);
	}
    }
    if (scount > 0)
	incr_histo(HISTO_EDGE_SMOOTHS, scount);
    is_smoothed = true;
}


bool read_cnf(FILE *infile, std::unordered_set<int> **data_variables, std::unordered_map<int,q25_ptr> **weights) {
    Cnf *cnf = new Cnf();
    if (!cnf->import_file(infile, true, false))
	return false;
    *data_variables = cnf->data_variables;
    *weights = cnf->input_weights;
    //    delete cnf;
    return true;
}

// literal_weights == NULL for unweighted
void Evaluator_q25::prepare_weights(std::unordered_map<int,q25_ptr> *literal_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	q25_ptr pwt = NULL;
	q25_ptr nwt = NULL;
	if (!literal_weights) {
	    // Unweighted counting
	    pwt = q25_from_32(1);
	    nwt = q25_from_32(1);
	} else {
	    if (literal_weights->find(v) != literal_weights->end())
		pwt = q25_copy((*literal_weights)[v]);
	    if (literal_weights->find(-v) != literal_weights->end())
		nwt = q25_copy((*literal_weights)[-v]);
	    if (pwt) {
		if (!nwt)
		    nwt = q25_one_minus(pwt);
	    } else {
		if (nwt)
		    pwt = q25_one_minus(nwt);
		else {
		    nwt = q25_from_32(1);
		    pwt = q25_from_32(1);
		}
	    }
	}
	q25_ptr sum = q25_add(pwt, nwt);
	if (smoothed)
	    smoothing_weights[v] = sum;
	else {
	    int mark = q25_enter();
	    q25_ptr recip = q25_mark(q25_recip(sum));
	    if (!q25_is_valid(recip)) {
		char *srecip = q25_string(sum);
		err(true, "Could not get reciprocal of summed weights for variable %d.  Sum = s\n", v, srecip);
		free(srecip);
	    }
	    rescale = q25_mul(q25_mark(rescale), q25_mark(sum));
	    pwt = q25_mul(q25_mark(pwt), recip);
	    nwt = q25_mul(q25_mark(nwt), recip);
	    q25_leave(mark);
	}

	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
}

q25_ptr Evaluator_q25::evaluate_edge(Egraph_edge &e, bool smoothed) {
    q25_ptr result = q25_from_32(1);
    int mark = q25_enter();
    for (int lit : e.literals) {
	q25_ptr wt = evaluation_weights[lit];
	result = q25_mul(q25_mark(result), wt);
    }
    if (smoothed) {
	for (int v : e.smoothing_variables) {
	    q25_ptr wt = smoothing_weights[v];
	    result = q25_mul(q25_mark(result), wt);
	}
    }
    q25_leave(mark);
    if (verblevel >= 4) {
	char *sresult = q25_string(result);
	report(4, "Evaluating edge (%d <-- %d).  Value = %s\n", e.to_id, e.from_id, sresult);
	free(sresult);
    }
    return result;
}

q25_ptr Evaluator_q25::evaluate(bool smoothed) {
    std::vector<q25_ptr> operation_values;
    operation_values.resize(egraph->operations.size());
    for (int id = 1; id <= egraph->operations.size(); id++) {
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    operation_values[id-1] = q25_from_32(1);
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    operation_values[id-1] = q25_from_32(0);
	}
    }
    for (Egraph_edge e : egraph->edges) {
	int mark = q25_enter();
	q25_ptr edge_val = evaluate_edge(e, smoothed);
	q25_ptr product = q25_mark(q25_mul(q25_mark(edge_val), operation_values[e.from_id-1]));
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	q25_ptr new_val = multiply ? 
	    q25_mul(q25_mark(operation_values[e.to_id-1]), product) :
	    q25_add(q25_mark(operation_values[e.to_id-1]), product);
	if (verblevel >= 4) {
	    char *sfrom = q25_string(operation_values[e.from_id-1]);
	    char *sold = q25_string(operation_values[e.to_id-1]);
	    char *sedge = q25_string(edge_val);
	    char *sproduct = q25_string(product);
	    char *snew_val = q25_string(new_val);
	    report(4, "Density: Updating %d from %d.  %s * %s %c %s --> %s\n",
		   e.to_id, e.from_id, sfrom, sedge, multiply ? '*' : '+', sold, snew_val);
	    free(sfrom); free(sold); free(sedge); free(sproduct); free(snew_val);
	}
	operation_values[e.to_id-1] = new_val;
	q25_leave(mark);
    }
    q25_ptr result = operation_values[egraph->root_id-1];
    for (int id = 1; id <= egraph->operations.size(); id++) {
	if (id != egraph->root_id)
	    q25_free(operation_values[id-1]);
    }
    operation_values.clear();
    if (!smoothed) {
	q25_ptr oresult = result;
	result = q25_mul(rescale, oresult);
	if (verblevel >= 4) {
	    char *srescale = q25_string(rescale);
	    char *soresult = q25_string(oresult);
	    char *sresult = q25_string(result);
	    report(4, "Final result: Rescale = %s.  Density = %s.  Result = %s\n",
		   srescale, soresult, sresult);
	    free(srescale); free(soresult), free(sresult);
	} else if (verblevel >= 4) {
	    char *sresult = q25_string(result);
	    report(4, "Smoothed result = %s\n", sresult);
	    free(sresult);
	}
	q25_free(oresult);
    }
    return result;
}

Evaluator_q25::Evaluator_q25(Egraph *eg) { 
    egraph = eg;
    rescale = NULL;
    clear_evaluation();
}
    
void Evaluator_q25::clear_evaluation() {
    for (auto iter : evaluation_weights)
	q25_free(iter.second);
    evaluation_weights.clear();
    for (auto iter : smoothing_weights)
	q25_free(iter.second);
    smoothing_weights.clear();
    q25_free(rescale);
    rescale = q25_from_32(1);
}

