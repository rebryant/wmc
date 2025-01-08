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
#include "analysis.h"
#include "egraph.hh"
#include "cnf_info.hh"

/*
  Useful functions

 */

const char *mpf_string(mpf_srcptr val) {
    char buf[2048];
    char boffset = 0;
    mp_exp_t ecount;
    char *sval = mpf_get_str(NULL, &ecount, 10, 40, val);
    if (!sval || strlen(sval) == 0 || sval[0] == '0') {
	strcpy(buf, "0.0");
    } else {
	int voffset = 0;
	bool neg = sval[0] == '-';
	if (neg) {
	    voffset++;
	    buf[boffset++] = '-';
	}
	if (ecount == 0) {
	    buf[boffset++] = '0';
	    buf[boffset++] = '.';
	} else {
	    buf[boffset++] = sval[voffset++];
	    buf[boffset++] = '.';
	    ecount--;
	}
	if (sval[voffset] == 0)
	    buf[boffset++] = '0';
	else {
	    while(sval[voffset] != 0)
		buf[boffset++] = sval[voffset++];
	}
	if (ecount != 0) {
	    buf[boffset++] = 'e';
	    snprintf(&buf[boffset], 24, "%ld", (long) ecount);
	} else
	    buf[boffset] = 0;
    }
    free(sval);
    return archive_string(buf);
}


const char *mpfr_string(mpfr_srcptr val) {
    mpf_t fval;
    mpf_init(fval);
    mpfr_get_f(fval, val, MPFR_RNDN);
    const char* result = mpf_string(fval);
    mpf_clear(fval);
    return result;
}

double digit_precision_mpfi(mpfi_srcptr v) {
    mpfr_t diam;
    mpfr_init(diam);
    mpfi_diam_abs(diam, v);
    if (mpfr_sgn(diam) == 0) {
	mpfr_clear(diam);
	return MAX_DIGIT_PRECISION;
    }

    mpfr_t avg;
    mpfr_init(avg);
    mpfi_mid(avg, v);
    mpfr_abs(avg, avg, MPFR_RNDN);

    mpfr_t div;
    mpfr_init(div);
    mpfr_div(div, diam, avg, MPFR_RNDN);
    double ddiv = mpfr_get_d(div, MPFR_RNDN);
    double result = -log10(ddiv);
    if (result < 0.0)
	result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
	result = MAX_DIGIT_PRECISION;
    mpfr_clears(diam, avg, div, NULL);
    return result;
}

double digit_precision_mpfr(mpfr_srcptr val, mpq_srcptr ref) {
    mpfr_prec_t save_prec = mpfr_get_default_prec();
    mpfr_t fref;
    mpfr_prec_t vprec = 3*mpfr_get_prec(val);
    mpfr_set_default_prec(vprec);
    mpfr_init_set_q(fref, ref, MPFR_RNDN);
    if (mpfr_cmp(val, fref) == 0) {
	mpfr_clear(fref);
	mpfr_set_default_prec(save_prec);
	return (double) MAX_DIGIT_PRECISION;
    }

    mpfr_t num;
    mpfr_init_set(num, val, MPFR_RNDN);
    mpfr_sub(num, num, fref, MPFR_RNDN);
    mpfr_abs(num, num, MPFR_RNDN);

    mpfr_t den;
    mpfr_init_set(den, val, MPFR_RNDN);
    mpfr_abs(den, den, MPFR_RNDN);
    mpfr_abs(fref, fref, MPFR_RNDN);
    mpfr_add(den, den, fref, MPFR_RNDN);

    mpfr_mul_d(den, den, 0.5, MPFR_RNDN);
    mpfr_div(num, num, den, MPFR_RNDN);
    mpfr_log10(num, num, MPFR_RNDN);
    double result = -mpfr_get_d(num, MPFR_RNDN);
    if (result < 0)
	result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
	result = MAX_DIGIT_PRECISION;
    mpfr_clears(fref, num, den, NULL);
    mpfr_set_default_prec(save_prec);
    return result;
    
}

double digit_precision_mpf(mpf_srcptr val, mpq_srcptr ref) {
    mpfr_t rval;
    int prec = mpf_get_prec(val);
    mpfr_init2(rval, prec);
    mpfr_set_f(rval, val, MPFR_RNDN);
    double result = digit_precision_mpfr(rval, ref);
    mpfr_clear(rval);
    return result;
}

double digit_precision_d(double val, mpq_srcptr ref) {
    mpfr_t rval;
    int prec = 64;
    mpfr_init2(rval, prec);
    mpfr_set_d(rval, val, MPFR_RNDN);
    double result = digit_precision_mpfr(rval, ref);
    mpfr_clear(rval);
    return result;
}




/*******************************************************************************************************************
 Graph representing NNF formula
*******************************************************************************************************************/

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

/*******************************************************************************************************************
Evaluation via Q25
*******************************************************************************************************************/

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

// literal_string_weights == NULL for unweighted
void Evaluator_q25::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	q25_ptr pwt = NULL;
	q25_ptr nwt = NULL;
	if (!literal_string_weights) {
	    // Unweighted counting
	    pwt = q25_from_32(1);
	    nwt = q25_from_32(1);
	} else {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		pwt = q25_from_string((*literal_string_weights)[v]);
		if (!q25_is_valid(pwt))
		    err(true, "Q25: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		nwt = q25_from_string((*literal_string_weights)[-v]);
		if (!q25_is_valid(nwt))
		    err(true, "Q25: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
	    }
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
		err(true, "Q25: Could not get reciprocal of summed weights for variable %d.  Sum = s\n", v, sum);
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
	report(4, "Q25: Evaluating edge (%d <-- %d).  Value = %s\n", e.to_id, e.from_id, sresult);
	free(sresult);
    }
    return result;
}

q25_ptr Evaluator_q25::evaluate(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    prepare_weights(literal_string_weights, smoothed);
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
	    report(4, "Q25: Density: Updating %d from %d.  %s * %s %c %s --> %s\n",
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
	    report(4, "Q25: Final result: Rescale = %s.  Density = %s.  Result = %s\n",
		   srescale, soresult, sresult);
	    free(srescale); free(soresult), free(sresult);
	} else if (verblevel >= 4) {
	    char *sresult = q25_string(result);
	    report(4, "Q25: Smoothed result = %s\n", sresult);
	    free(sresult);
	}
	q25_free(oresult);
    }
    return result;
}

/*******************************************************************************************************************
Evaluation via DOUBLE
*******************************************************************************************************************/

Evaluator_double::Evaluator_double(Egraph *eg) { 
    egraph = eg;
    rescale = 1.0;
    clear_evaluation();
}
    
void Evaluator_double::clear_evaluation() {
    evaluation_weights.clear();
    smoothing_weights.clear();
    rescale = 1.0;
}

// literal_string_weights == NULL for unweighted
void Evaluator_double::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	double pwt = 0.0;
	bool have_pos = false;
	double nwt = 0.0;
	bool have_neg = false;
	if (!literal_string_weights) {
	    // Unweighted counting
	    pwt = 1.0;
	    nwt = 1.0;
	} else {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		if (sscanf((*literal_string_weights)[v], "%lf", &pwt) != 1)
		    err(true, "DBL: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
		have_pos = true;
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		if (sscanf((*literal_string_weights)[-v], "%lf", &nwt) != 1)
		    err(true, "DBL: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
		have_neg = true;
	    }
	    if (have_pos) {
		if (!have_neg)
		    nwt = 1.0 - pwt;
	    } else {
		if (have_neg)
		    pwt = 1.0 - nwt;
		else {
		    nwt = 1.0;
		    pwt = 1.0;
		}
	    }
	}
	double sum = pwt + nwt;
	if (smoothed)
	    smoothing_weights[v] = sum;
	else {
	    if (sum == 0) {
		err(true, "DBL: Could not get reciprocal of summed weights for variable %d.  Sum = %f\n", v, sum);
	    }
	    rescale *= sum;
	    pwt = pwt/sum;
	    nwt = nwt/sum;
	}
	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
}

double Evaluator_double::evaluate_edge(Egraph_edge &e, bool smoothed) {
    double result = 1.0;
    for (int lit : e.literals) {
	double wt = evaluation_weights[lit];
	result *= wt;
    }
    if (smoothed) {
	for (int v : e.smoothing_variables) {
	    double wt = smoothing_weights[v];
	    result *= wt;
	}
    }
    if (verblevel >= 4) {
	report(4, "DBL: Evaluating edge (%d <-- %d).  Value = %f\n", e.to_id, e.from_id, result);
    }
    return result;
}

double Evaluator_double::evaluate(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    prepare_weights(literal_string_weights, smoothed);
    std::vector<double> operation_values;
    operation_values.resize(egraph->operations.size());
    for (int id = 1; id <= egraph->operations.size(); id++) {
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    operation_values[id-1] = 1.0;
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    operation_values[id-1] = 0.0;
	}
    }
    for (Egraph_edge e : egraph->edges) {
	double edge_val = evaluate_edge(e, smoothed);
	double product = edge_val * operation_values[e.from_id-1];
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	double new_val = multiply ? 
	    operation_values[e.to_id-1] * product:
	    operation_values[e.to_id-1] + product;
	if (verblevel >= 4) {
	    double dfrom = operation_values[e.from_id-1];
	    double dold = operation_values[e.to_id-1];
	    report(4, "DBL: Density: Updating %d from %d.  %f * %f %c %f --> %f\n",
		   e.to_id, e.from_id, dfrom, edge_val, multiply ? '*' : '+', dold, new_val);
	}
	operation_values[e.to_id-1] = new_val;
    }

    double result = operation_values[egraph->root_id-1];
    operation_values.clear();
    if (!smoothed) {
	double oresult = result;
	result *= rescale;
	if (verblevel >= 4) {
	    report(4, "DBL: Final result: Rescale = %f.  Density = %f.  Result = %f\n",
		   rescale, oresult, result);
	}
    } else
	report(4, "DBL: Smoothed result = %f\n", result);

    return result;
}

/*******************************************************************************************************************
Evaluation via Gnu multi-precision floating-point arithmetic
*******************************************************************************************************************/


Evaluator_mpf::Evaluator_mpf(Egraph *eg) { 
    egraph = eg;
    clear_evaluation();
}
    
void Evaluator_mpf::clear_evaluation() {
    //    for (auto iter : evaluation_weights)
    //	mpf_clear(iter.second);
    evaluation_weights.clear();
    //    for (auto iter : smoothing_weights)
    //	mpf_clear(iter.second);
    smoothing_weights.clear();
    rescale = 1;
}

static void mpf_one_minus(mpf_ptr dest, mpf_srcptr val) {
    mpf_t one;
    mpf_init(one);
    mpf_set_ui(one, 1);
    mpf_neg(dest, val);
    mpf_add(dest, dest, one);
    mpf_clear(one);
}

// literal_string_weights == NULL for unweighted
bool Evaluator_mpf::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	mpf_class pwt = 1;
	bool gotp = false;
	mpf_class nwt = 1;
	bool gotn = false;

	if (literal_string_weights) {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		q25_ptr qpwt = q25_from_string((*literal_string_weights)[v]);
		if (!q25_is_valid(qpwt)) {
		    err(false, "MPF: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		if (!q25_to_mpf(pwt.get_mpf_t(), qpwt)) {
		    err(false, "MPF: Couldn't convert from q25 to mpf for literal %d with string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		q25_free(qpwt);
		gotp = true;
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		q25_ptr qnwt = q25_from_string((*literal_string_weights)[-v]);
		if (!q25_is_valid(qnwt)) {
		    err(false, "MPF: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		if (!q25_to_mpf(nwt.get_mpf_t(), qnwt)) {
		    err(false, "MPF: Couldn't convert from q25 to mpf for literal %d with string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		q25_free(qnwt);
		gotn = true;
	    }
	    if (gotp) {
		if (!gotn)
		    mpf_one_minus(nwt.get_mpf_t(), pwt.get_mpf_t());
	    } else {
		if (gotn)
		    mpf_one_minus(pwt.get_mpf_t(), nwt.get_mpf_t());
	    }
	}
	mpf_class sum = nwt+pwt;
	if (smoothed)
	    smoothing_weights[v] = sum;
	else {
	    if (cmp(sum, mpf_class(0)) == 0) {
		err(false, "MPF: Weights for variable %d sum to 0\n", v);
		return false;
	    }
	    rescale *= sum;
	    pwt /= sum;
	    nwt /= sum;
	}
	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
    return true;
}

void Evaluator_mpf::evaluate_edge(mpf_class &value, Egraph_edge &e, bool smoothed) {
    value = 1;
    for (int lit : e.literals)
	value *= evaluation_weights[lit];
    if (smoothed) {
	for (int v : e.smoothing_variables)
	    value *= smoothing_weights[v];
    }
    if (verblevel >= 4) {
	mp_exp_t exp;
	char *svalue = mpf_get_str(NULL, &exp, 10, 40, value.get_mpf_t());
	report(4, "MPF: Evaluating edge (%d <-- %d).  Value = 0.%se%ld\n", e.to_id, e.from_id, svalue, exp);
	free(svalue);
    }
}

bool Evaluator_mpf::evaluate(mpf_class &count, std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    if (!prepare_weights(literal_string_weights, smoothed))
	return false;
    std::vector<mpf_class> operation_values;
    operation_values.resize(egraph->operations.size());
    for (int id = 1; id <= egraph->operations.size(); id++) {
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    operation_values[id-1] = 1;
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    operation_values[id-1] = 0;
	}
    }
    for (Egraph_edge e : egraph->edges) {
	char *sold = NULL;
	char *sedge = NULL;
	mp_exp_t eold, eedge, efrom, eproduct, enew_val;	
	mpf_class product;

	evaluate_edge(product, e, smoothed);

	if (verblevel >= 4) {
	    sedge = mpf_get_str(NULL, &eedge, 10, 40, product.get_mpf_t());
	    sold = mpf_get_str(NULL, &eold, 10, 40, operation_values[e.to_id-1].get_mpf_t());
	}

	product *= operation_values[e.from_id-1];
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	if (multiply)
	    operation_values[e.to_id-1] *= product;
	else
	    operation_values[e.to_id-1] += product;
	if (verblevel >= 4) {
	    char *sfrom = mpf_get_str(NULL, &efrom, 10, 40, operation_values[e.from_id-1].get_mpf_t());
	    char *sproduct = mpf_get_str(NULL, &eproduct, 10, 40, product.get_mpf_t());
	    char *snew_val = mpf_get_str(NULL, &enew_val, 10, 40, operation_values[e.to_id-1].get_mpf_t());
	    report(4, "MPF: Density: Updating %d from %d.  0.%se%ld * 0.%se%ld %c 0.%se%ld --> 0.%se%ld\n",
		   e.to_id, e.from_id, sfrom, efrom, sedge, eedge, multiply ? '*' : '+', sold, eold, snew_val, enew_val);
	    free(sfrom); free(sold); free(sedge); free(sproduct); free(snew_val);
	}
    }
    count = operation_values[egraph->root_id-1];
    //    for (int id = 1; id <= egraph->operations.size(); id++)
    //	mpf_clear(operation_values[id-1]);
    operation_values.clear();

    if (smoothed) {
	if (verblevel >= 4) {
	    mp_exp_t ecount;
	    char *scount = mpf_get_str(NULL, &ecount, 10, 40, count.get_mpf_t());
	    report(4, "MPF: Smoothed count = 0.%se%ld\n", scount, ecount);
	    free(scount);
	}
    } else {
	char *socount = NULL;
	mp_exp_t eocount, erescale, ecount;
	if (verblevel >= 4)
	    socount = mpf_get_str(NULL, &eocount, 10, 40, count.get_mpf_t());
	count *= rescale;
	if (verblevel >= 4) {
	    char *srescale = mpf_get_str(NULL, &erescale, 10, 40, rescale.get_mpf_t());
	    char *scount = mpf_get_str(NULL, &ecount, 10, 40, count.get_mpf_t());
	    report(4, "MPF: Final count: Rescale = 0.%se%ld.  Density = 0.%se%ld.  Count = 0.%se%ld\n",
		   srescale, socount, scount);
	    free(srescale); free(socount), free(scount);
	}
    }
    return true;
}


/*******************************************************************************************************************
Evaluation via Gnu multi-precision rational arithmetic
*******************************************************************************************************************/


static size_t mpq_bytes(mpq_srcptr val) {
    /* Overhead */
    size_t size = 32;
    mpz_t mp_num, mp_den;
    mpz_init(mp_num); mpz_init(mp_den);
    mpq_get_num(mp_num, val);
    size += mpz_size(mp_num) * sizeof(mp_limb_t);
    mpq_get_den(mp_den, val);
    size += mpz_size(mp_den) * sizeof(mp_limb_t);
    size += 8 * ceil(0.125 * size);
    mpz_clear(mp_num); mpz_clear(mp_den);
    return size;
}


Evaluator_mpq::Evaluator_mpq(Egraph *eg) { 
    egraph = eg;
    clear_evaluation();
}
    
void Evaluator_mpq::clear_evaluation() {
    //    for (auto iter : evaluation_weights)
    //	mpq_clear(iter.second);
    evaluation_weights.clear();
    //    for (auto iter : smoothing_weights)
    //	mpq_clear(iter.second);
    smoothing_weights.clear();
    rescale = 1;
    max_bytes = 0;
}

static void mpq_one_minus(mpq_ptr dest, mpq_srcptr val) {
    mpq_t one;
    mpq_init(one);
    mpq_set_ui(one, 1, 1);
    mpq_neg(dest, val);
    mpq_add(dest, dest, one);
    mpq_clear(one);
}

// literal_string_weights == NULL for unweighted
bool Evaluator_mpq::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	mpq_class pwt = 1;
	bool gotp = false;
	mpq_class nwt = 1;
	bool gotn = false;

	if (literal_string_weights) {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		q25_ptr qpwt = q25_from_string((*literal_string_weights)[v]);
		if (!q25_is_valid(qpwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		if (!q25_to_mpq(pwt.get_mpq_t(), qpwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		q25_free(qpwt);
		gotp = true;
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		q25_ptr qnwt = q25_from_string((*literal_string_weights)[-v]);
		if (!q25_is_valid(qnwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		if (!q25_to_mpq(nwt.get_mpq_t(), qnwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		q25_free(qnwt);
		gotn = true;
	    }
	    if (gotp) {
		if (!gotn)
		    mpq_one_minus(nwt.get_mpq_t(), pwt.get_mpq_t());
	    } else {
		if (gotn)
		    mpq_one_minus(pwt.get_mpq_t(), nwt.get_mpq_t());
	    }
	}
	mpq_class sum = nwt+pwt;
	if (smoothed)
	    smoothing_weights[v] = sum;
	else {
	    if (cmp(sum, mpq_class(0)) == 0) {
		err(false, "MPQ: Weights for variable %d sum to 0\n", v);
		return false;
	    }
	    rescale *= sum;
	    pwt /= sum;
	    nwt /= sum;
	}
	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
    return true;
}

void Evaluator_mpq::evaluate_edge(mpq_class &value, Egraph_edge &e, bool smoothed) {
    // Use breadth-first evaluation to form balanced binary tree
    std::vector<mpq_class> eval_queue;
    for (int lit : e.literals)
	eval_queue.push_back(evaluation_weights[lit]);
    if (smoothed)
	for (int v : e.smoothing_variables)
	    eval_queue.push_back(smoothing_weights[v]);
    if (eval_queue.size() == 0)
	value = 1;
    else {
	size_t index = 0;
	while (index < eval_queue.size()-1) {
	    eval_queue.push_back(eval_queue[index] * eval_queue[index+1]);
	    index += 2;
	}
	value = eval_queue[index];
    }
    if (verblevel >= 4) {
	char *svalue = mpq_get_str(NULL, 10, value.get_mpq_t());
	report(4, "MPQ: Evaluating edge (%d <-- %d).  Value = %s\n", e.to_id, e.from_id, svalue);
	free(svalue);
    }
    size_t bytes = mpq_bytes(value.get_mpq_t());
    if (bytes > max_bytes)
	max_bytes = bytes;

}

bool Evaluator_mpq::evaluate(mpq_class &count, std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    if (!prepare_weights(literal_string_weights, smoothed))
	return false;
    std::vector<mpq_class> operation_values;
    operation_values.resize(egraph->operations.size());
    for (int id = 1; id <= egraph->operations.size(); id++) {
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    operation_values[id-1] = 1;
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    operation_values[id-1] = 0;
	}
    }
    for (Egraph_edge e : egraph->edges) {
	char *sold = NULL;
	char *sedge = NULL;

	mpq_class product;
	evaluate_edge(product, e, smoothed);

	if (verblevel >= 4) {
	    sedge = mpq_get_str(NULL, 10, product.get_mpq_t());
	    sold = mpq_get_str(NULL, 10, operation_values[e.to_id-1].get_mpq_t());
	}

	product *= operation_values[e.from_id-1];
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	if (multiply)
	    operation_values[e.to_id-1] *= product;
	else
	    operation_values[e.to_id-1] += product;
	size_t bytes = mpq_bytes(operation_values[e.to_id-1].get_mpq_t());
	if (bytes > max_bytes)
	    max_bytes = bytes;
				 
	if (verblevel >= 4) {
	    char *sfrom = mpq_get_str(NULL, 10, operation_values[e.from_id-1].get_mpq_t());
	    char *sproduct = mpq_get_str(NULL, 10, product.get_mpq_t());
	    char *snew_val = mpq_get_str(NULL, 10, operation_values[e.to_id-1].get_mpq_t());
	    report(4, "MPQ: Density: Updating %d from %d.  %s * %s %c %s --> %s\n",
		   e.to_id, e.from_id, sfrom, sedge, multiply ? '*' : '+', sold, snew_val);
	    free(sfrom); free(sold); free(sedge); free(sproduct); free(snew_val);
	}
    }
    count = operation_values[egraph->root_id-1];
    //    for (int id = 1; id <= egraph->operations.size(); id++)
    //	mpq_clear(operation_values[id-1]);
    operation_values.clear();

    if (smoothed) {
	if (verblevel >= 4) {
	    char *scount = mpq_get_str(NULL, 10, count.get_mpq_t());
	    report(4, "MPQ: Smoothed count = %s\n", scount);
	    free(scount);
	}
    } else {
	char *socount = NULL;
	if (verblevel >= 4)
	    socount = mpq_get_str(NULL, 10, count.get_mpq_t());
	count *= rescale;
	if (verblevel >= 4) {
	    char *srescale = mpq_get_str(NULL, 10, rescale.get_mpq_t());
	    char *scount = mpq_get_str(NULL, 10, count.get_mpq_t());
	    report(4, "MPQ: Final count: Rescale = %s.  Density = %s.  Count = %s\n",
		   srescale, socount, scount);
	    free(srescale); free(socount), free(scount);
	}
    }
    return true;
}

/*******************************************************************************************************************
Evaluation via MPFI
*******************************************************************************************************************/

Evaluator_mpfi::Evaluator_mpfi(Egraph *eg, int tdp, bool ref) { 
    egraph = eg;
    target_digit_precision = tdp;
    refine = ref;
    mpfi_init(rescale);
    clear_evaluation();
}
    
void Evaluator_mpfi::clear_evaluation() {
    evaluation_weights.clear();
    smoothing_weights.clear();
    mpfi_set_d(rescale, 1.0);
    min_digit_precision = (double) MAX_DIGIT_PRECISION;
    precision_failure_count = 0;
}


// literal_string_weights == NULL for unweighted
bool Evaluator_mpfi::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    clear_evaluation();
    for (int v : *egraph->data_variables) {
	mpq_class pwt = 1;
	bool gotp = false;
	mpq_class nwt = 1;
	bool gotn = false;

	if (literal_string_weights) {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		q25_ptr qpwt = q25_from_string((*literal_string_weights)[v]);
		if (!q25_is_valid(qpwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		if (!q25_to_mpq(pwt.get_mpq_t(), qpwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", v, (*literal_string_weights)[v]);
		    return false;
		}
		q25_free(qpwt);
		gotp = true;
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		q25_ptr qnwt = q25_from_string((*literal_string_weights)[-v]);
		if (!q25_is_valid(qnwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		if (!q25_to_mpq(nwt.get_mpq_t(), qnwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", -v, (*literal_string_weights)[-v]);
		    return false;
		}
		q25_free(qnwt);
		gotn = true;
	    }
	    if (gotp) {
		if (!gotn)
		    mpq_one_minus(nwt.get_mpq_t(), pwt.get_mpq_t());
	    } else {
		if (gotn)
		    mpq_one_minus(pwt.get_mpq_t(), nwt.get_mpq_t());
	    }
	}
	mpq_class sum = nwt+pwt;
	if (smoothed)
	    smoothing_weights[v] = sum;
	else {
	    if (cmp(sum, mpq_class(0)) == 0) {
		err(false, "MPQ: Weights for variable %d sum to 0\n", v);
		return false;
	    }
	    mpfi_mul_q(rescale, rescale, sum.get_mpq_t());
	    pwt /= sum;
	    nwt /= sum;
	}
	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
    return true;
}

void Evaluator_mpfi::evaluate_edge(mpfi_ptr value, Egraph_edge &e, bool smoothed) {
    mpfi_set_d(value, 1.0);
    for (int lit : e.literals)
	mpfi_mul_q(value, value, evaluation_weights[lit].get_mpq_t());
    if (smoothed) {
	for (int v : e.smoothing_variables)
	    mpfi_mul_q(value, value, smoothing_weights[v].get_mpq_t());
    }
}

bool Evaluator_mpfi::evaluate(mpfi_ptr count, std::unordered_map<int,const char*> *literal_string_weights, bool smoothed) {
    if (!prepare_weights(literal_string_weights, smoothed))
	return false;
    mpfi_t *operation_values = new mpfi_t[egraph->operations.size()];
    for (int id = 1; id <= egraph->operations.size(); id++) {
	mpfi_init(operation_values[id-1]);
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    mpfi_set_d(operation_values[id-1], 1.0);
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    mpfi_set_d(operation_values[id-1], 0.0);
	}
    }
    for (Egraph_edge e : egraph->edges) {
	mpfi_t product;
	mpfi_init(product);
	evaluate_edge(product, e, smoothed);
	mpfi_mul(product, product, operation_values[e.from_id-1]);
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	if (multiply)
	    mpfi_mul(operation_values[e.to_id-1], operation_values[e.to_id-1], product);
	else {
	    mpfi_add(operation_values[e.to_id-1], operation_values[e.to_id-1], product);
	}
	double dp = digit_precision_mpfi(operation_values[e.to_id-1]);
	if (dp < min_digit_precision)
	    min_digit_precision = dp;
	if (dp < target_digit_precision)
	    precision_failure_count++;
    }
    mpfi_swap(count, operation_values[egraph->root_id-1]);
    for (int id = 1; id <= egraph->operations.size(); id++)
	mpfi_clear(operation_values[id-1]);
    delete[] operation_values;

    if (!smoothed)
	mpfi_mul(count, count, rescale);

    return true;
}
