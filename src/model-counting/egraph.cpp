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
#include <math.h>

#include "report.h"
#include "counters.h"
#include "analysis.h"
#include "egraph.hh"
#include "cnf_info.hh"

/*
  Useful functions
 */

/*
  How many digits of precision can we guarantee when all weights are nonnegative?
  Compute floats with specified bit precision over formula with specified number of variables
  constant = 3 for smoothed evaluation and 5 for unsmoothed
*/
double digit_precision_bound(int bit_precision, int nvar, double constant) {
    return (double) bit_precision * log10(2) - log10(nvar * constant);
}

/*
  How many bits of floating-point precision are required to achieve
  target digit precision when all weights are nonnegative?
  constant = 3 for smoothed evaluation and 5 for unsmoothed
 */
int required_bit_precision(double target_precision, int nvar, double constant, bool nonnegative) {
    double minp = target_precision * log2(10.0) + log2(nvar * constant);
    if (nonnegative && minp <= 52)
	return 52;
    /* Must be multiple of 64 */
    return 64 * ceil(minp/64);
}


const char *mpf_string(mpf_srcptr val, int digits) {
    char buf[2048];
    char boffset = 0;
    mp_exp_t ecount;
    char *sval = mpf_get_str(NULL, &ecount, 10, digits, val);
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


const char *mpfr_string(mpfr_srcptr val, int digits) {
    mpf_t fval;
    mpf_init(fval);
    mpfr_get_f(fval, val, MPFR_RNDN);
    const char* result = mpf_string(fval, digits);
    mpf_clear(fval);
    return result;
}

double digit_precision_mpfr(mpfr_srcptr x_est, mpq_srcptr x) {
    if (mpfr_cmp_q(x_est, x) == 0)
	return (double) MAX_DIGIT_PRECISION;

    mpfr_prec_t save_prec = mpfr_get_default_prec();
    mpfr_prec_t prec = 3*mpfr_get_prec(x_est);
    mpfr_set_default_prec(prec);

    mpfr_t num;
    mpfr_t den;
    if (mpq_sgn(x) == 0) {
	mpfr_init_set_d(den, 1.0, MPFR_RNDN);
	mpfr_init_set(num, x_est, MPFR_RNDN);
	mpfr_abs(num, num, MPFR_RNDN);
	if (mpfr_cmp_d(num, 1.0) > 0)
	    mpfr_set_d(num, 1.0, MPFR_RNDN);
    } else {
	mpfr_init_set_q(den, x, MPFR_RNDN);
	mpfr_abs(den, den, MPFR_RNDN);
	mpfr_init_set_q(num, x, MPFR_RNDN);
	mpfr_sub(num, num, x_est, MPFR_RNDN);
	mpfr_abs(num, num, MPFR_RNDN);
    }
    mpfr_div(num, num, den, MPFR_RNDN);
    mpfr_log10(num, num, MPFR_RNDN);
    double result = -mpfr_get_d(num, MPFR_RNDN);
    if (result < 0)
	result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
	result = MAX_DIGIT_PRECISION;
    mpfr_clears(num, den, NULL);
    mpfr_set_default_prec(save_prec);
    return result;
}

double digit_precision_mpfi(mpfi_srcptr v) {
    mpfr_t left;
    mpfr_init(left);
    mpfr_t right;
    mpfr_init(right);
    mpfi_get_left(left, v);
    mpfi_get_right(right, v);
    if (mpfr_sgn(left) != mpfr_sgn(right))
	return 0.0;
    mpfr_t diam;
    mpfr_init(diam);
    mpfi_diam_rel(diam, v);
    if (mpfr_sgn(diam) == 0) {
	mpfr_clear(diam);
	return MAX_DIGIT_PRECISION;
    }
    mpfr_log10(diam, diam, MPFR_RNDN);
    double result = -mpfr_get_d(diam, MPFR_RNDN);
    if (result < 0)
	result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
	result = MAX_DIGIT_PRECISION;
    mpfr_clear(diam);
    return result;
}


double digit_precision_mpf(mpf_srcptr x_est, mpq_srcptr x) {
    mpfr_t rx_est;
    int prec = mpf_get_prec(x_est);
    mpfr_init2(rx_est, prec);
    mpfr_set_f(rx_est, x_est, MPFR_RNDN);
    double result = digit_precision_mpfr(rx_est, x);
    mpfr_clear(rx_est);
    return result;
}

static uint64_t double2bits(double x) {
    union {
	double dx;
	uint64_t bx;
    } u;
    u.dx = x;
    return u.bx;
}    

static bool double_is_special(double x) {
    uint64_t bits = double2bits(x);
    int  biased_exp = (bits >> 52) & 0x7FF;
    return biased_exp == 0x7FF;
}


double digit_precision_d(double x_est, mpq_srcptr x) {
    if (double_is_special(x_est))
	return 0.0;
    mpfr_t rx_est;
    int prec = 64;
    mpfr_init2(rx_est, prec);
    mpfr_set_d(rx_est, x_est, MPFR_RNDN);
    double result = digit_precision_mpfr(rx_est, x);
    mpfr_clear(rx_est);
    return result;
}

static void mpq_one_minus(mpq_ptr dest, mpq_srcptr val) {
    mpq_t one;
    mpq_init(one);
    mpq_set_ui(one, 1, 1);
    mpq_neg(dest, val);
    mpq_add(dest, dest, one);
    mpq_clear(one);
}

// Form product of values in queue.
// Use breadth-first evaluation to form balanced binary tree
static void reduce_product(mpq_class &product, std::vector<mpq_class> &eval_queue) {
    size_t start_size = eval_queue.size();
    if (eval_queue.size() == 0)
	product = 1.0;
    else if (eval_queue.size() == 1) 
	product = eval_queue[0];
    else {
	size_t index = 0;
	while (index < eval_queue.size()-1) {
	    eval_queue.push_back(eval_queue[index] * eval_queue[index+1]);
	    index += 2;
	}
	product = eval_queue[index];
    }
    eval_queue.resize(start_size);
}



/*******************************************************************************************************************
 Graph representing NNF formula
*******************************************************************************************************************/

static const char *nnf_type_name[NNF_NUM] = { "NONE", "TRUE", "FALSE", "AND", "OR" };
static const char nnf_type_char[NNF_NUM] = { '\0', 't', 'f', 'a', 'o' };

// Current line number
static int line_number = 0;

Egraph::Egraph(std::unordered_set<int> *dvars, int nv) {
    data_variables = dvars;
    is_smoothed = false;
    smooth_variable_count = 0;
    disabled_edge_count = 0;
    nvar = nv;
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
    std::unordered_set<int> smoothing_variables;
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
		int var = largs[pos];
		smoothing_variables.insert(var);
		add_smoothing_variable(eid, var);
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
    smooth_variable_count = smoothing_variables.size();
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
    std::unordered_set<int> smoothed_variables; 
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
		smoothed_variables.insert(v);
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
		smoothed_variables.insert(v);
	}
    }
    if (scount > 0)
	incr_histo(HISTO_EDGE_SMOOTHS, scount);
    is_smoothed = true;
    smooth_variable_count = smoothed_variables.size();
}

void Egraph::reset_smooth() {
    if (is_smoothed || (smooth_variable_count == 0 && disabled_edge_count == 0))
	return;
    for (int id = 1; id <= edges.size(); id++) {
	incr_count_by(COUNT_SMOOTH_VARIABLES, - (int) edges[id-1].smoothing_variables.size());
	if (edges[id-1].smoothing_variables.size() != 0)
	    report(4, "Removing %d variables from edge #%d (%d <-- %d)\n",
		   (int) edges[id-1].smoothing_variables.size(), id, edges[id-1].to_id, edges[id-1].from_id);
	edges[id-1].has_zero = false;
	edges[id-1].smoothing_variables.clear();
    }
    reset_histo(HISTO_EDGE_SMOOTHS);
    smooth_variable_count = 0;
    disabled_edge_count = 0;
}

void Egraph::smooth_single(int var, bool is_zero) {
    int disable_count = 0;
    std::vector<bool> var_found;
    var_found.resize(operations.size(), false);
    std::vector<bool> edge_contains;
    edge_contains.resize(edges.size(), false);
    int ecount = 0;
    for (int id = 1; id <= edges.size(); id++) {
	int from_id = edges[id-1].from_id;
	int to_id = edges[id-1].to_id;
	if (var_found[from_id-1]) {
	    var_found[to_id-1] = true;
	    continue;
	}
	for (int lit : edges[id-1].literals) {
	    int evar = IABS(lit);
	    if (evar == var) {
		var_found[to_id-1] = true;
		edge_contains[id-1] = true;
		break;
	    }
	}
    }
    for (int id = 1; id <= edges.size(); id++) {
	int from_id = edges[id-1].from_id;
	int to_id = edges[id-1].to_id;
	if (operations[from_id-1].type == NNF_FALSE)
	    continue;
	if (operations[to_id-1].type == NNF_AND)
	    continue;
	int scount = 0;
	if (var_found[to_id-1] && !var_found[from_id-1] && !edge_contains[id-1]) {
	    ecount++;
	    if (is_zero) {
		edges[id-1].has_zero = true;
		disable_count++;
		disabled_edge_count++;
		report(4, "Disabling edge due to variable %d.  #%d (%d <-- %d)\n", var, id, to_id, from_id);
	    } else {
		add_smoothing_variable(id, var);
		report(4, "Adding smoothing variable %d on edge #%d (%d <-- %d)\n", var, id, to_id, from_id);
	    }
	}
    }
    // Check at root
    int id = edges.size();
    int child_id = edges[id-1].to_id;
    if (!var_found[child_id-1]) {
	if (is_zero) {
	    edges[id-1].has_zero = true;
	    disable_count++;
	    disabled_edge_count++;
	    report(3, "Disabling root due to smoothing of variable %d\n", var);
	} else {
	    ecount++;
	    add_smoothing_variable(id, var);
	}
    }
    if (disable_count > 0)
	report(3, "Disabled %d edges\n", disable_count);
    if (ecount > 0) {
	report(3, "Added smoothing variable %d to %d edges\n", var, ecount);
	smooth_variable_count++;
	incr_histo(HISTO_EDGE_SMOOTHS, 1);
    } else 
	report(3, "No copies of smoothing variable %d needed\n", var, ecount);

}

// literal_string_weights == NULL for unweighted
Egraph_weights * Egraph::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights) {
    Egraph_weights *weights = new(Egraph_weights);
    reset_smooth();
    weights->all_nonnegative = true;
    for (int v : *data_variables) {
	mpq_class pwt = 1;
	bool gotp = false;
	mpq_class nwt = 1;
	bool gotn = false;

	if (literal_string_weights) {
	    if (literal_string_weights->find(v) != literal_string_weights->end()) {
		q25_ptr qpwt = q25_from_string((*literal_string_weights)[v]);
		if (!q25_is_valid(qpwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", v, (*literal_string_weights)[v]);
		    delete weights;
		    return NULL;
		}
		if (!q25_to_mpq(pwt.get_mpq_t(), qpwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", v, (*literal_string_weights)[v]);
		    delete weights;
		    return NULL;
		}
		q25_free(qpwt);
		gotp = true;
	    }
	    if (literal_string_weights->find(-v) != literal_string_weights->end()) {
		q25_ptr qnwt = q25_from_string((*literal_string_weights)[-v]);
		if (!q25_is_valid(qnwt)) {
		    err(false, "MPQ: Couldn't parse input weight for literal %d from string '%s'\n", -v, (*literal_string_weights)[-v]);
		    delete weights;
		    return NULL;
		}
		if (!q25_to_mpq(nwt.get_mpq_t(), qnwt)) {
		    err(false, "MPQ: Couldn't convert from q25 to mpq for literal %d with string '%s'\n", -v, (*literal_string_weights)[-v]);
		    delete weights;
		    return NULL;
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
	if (is_smoothed)
	    weights->smoothing_weights[v] = sum;
	else if (cmp(sum, mpq_class(0)) == 0) {
	    weights->smoothing_weights[v] = sum;
	    smooth_single(v, true);
	} else if (cmp(sum, mpq_class(1)) != 0) {
	    weights->rescale_weights.push_back(sum);
	    pwt /= sum;
	    nwt /= sum;
	}
	weights->evaluation_weights[v] = pwt;
	weights->evaluation_weights[-v] = nwt;
	if (mpq_sgn(pwt.get_mpq_t()) < 0 || mpq_sgn(nwt.get_mpq_t()) < 0)
	    weights->all_nonnegative = false;
    }
    return weights;
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
    egraph->reset_smooth();
    q25_free(rescale);
    rescale = q25_from_32(1);

}

// literal_string_weights == NULL for unweighted
void Evaluator_q25::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights) {
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
	if (egraph->is_smoothed)
	    smoothing_weights[v] = sum;
	else if (q25_is_zero(sum)) {
	    smoothing_weights[v] = sum;
	    egraph->smooth_single(v, true);

	} else {
	    int mark = q25_enter();
	    q25_ptr recip = q25_mark(q25_recip(sum));
	    if (!q25_is_valid(recip)) {
		char *srecip = q25_string(sum);
		err(true, "Q25: Could not get reciprocal of summed weights for variable %d.  Sum = %s\n", v, srecip);
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

q25_ptr Evaluator_q25::evaluate_edge(Egraph_edge &e) {
    if (e.has_zero)
	return q25_from_32(0);
    q25_ptr result = q25_from_32(1);
    int mark = q25_enter();
    for (int lit : e.literals) {
	q25_ptr wt = evaluation_weights[lit];
	result = q25_mul(q25_mark(result), wt);
    }
    for (int v : e.smoothing_variables) {
	q25_ptr wt = smoothing_weights[v];
	result = q25_mul(q25_mark(result), wt);
    }
    q25_leave(mark);
    if (verblevel >= 4) {
	char *sresult = q25_string(result);
	report(4, "Q25: Evaluating edge (%d <-- %d).  Value = %s\n", e.to_id, e.from_id, sresult);
	free(sresult);
    }
    return result;
}

q25_ptr Evaluator_q25::evaluate(std::unordered_map<int,const char*> *literal_string_weights) {
    prepare_weights(literal_string_weights);
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
	q25_ptr edge_val = evaluate_edge(e);
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
    q25_ptr oresult = result;
    result = q25_mul(oresult, rescale);
    q25_free(oresult);

    if (verblevel >= 4) {
	char *sresult = q25_string(result);
	report(4, "Q25: Result = %s\n", sresult);
	free(sresult);
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
    egraph->reset_smooth();
}


// literal_string_weights == NULL for unweighted
void Evaluator_double::prepare_weights(std::unordered_map<int,const char*> *literal_string_weights) {
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
	if (egraph->is_smoothed)
	    smoothing_weights[v] = sum;
	else if (sum == 0.0) {
	    smoothing_weights[v] = sum;
	    egraph->smooth_single(v, true);
	} else {
	    rescale *= sum;
	    pwt = pwt/sum;
	    nwt = nwt/sum;
	}
	evaluation_weights[v] = pwt;
	evaluation_weights[-v] = nwt;
    }
}

double Evaluator_double::evaluate_edge(Egraph_edge &e) {
    if (e.has_zero)
	return 0.0;
    double result = 1.0;
    for (int lit : e.literals) {
	double wt = evaluation_weights[lit];
	result *= wt;
    }
    for (int v : e.smoothing_variables) {
	double wt = smoothing_weights[v];
	result *= wt;
    }
    if (verblevel >= 4) {
	report(4, "DBL: Evaluating edge (%d <-- %d).  Value = %f\n", e.to_id, e.from_id, result);
    }
    return result;
}

double Evaluator_double::evaluate(std::unordered_map<int,const char*> *literal_string_weights) {
    prepare_weights(literal_string_weights);
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
	double edge_val = evaluate_edge(e);
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
    result *= rescale;
    report(4, "DBL: Result = %f\n", result);

    return result;
}

/*******************************************************************************************************************
Evaluation via extended-range double
*******************************************************************************************************************/


Evaluator_erd::Evaluator_erd(Egraph *eg, Egraph_weights *wts) { 
    
    egraph = eg;

    mpf_t mval;
    mpf_init2(mval, 64);

    /* Convert weight values from mpq to Erd */
    evaluation_weights.clear();
    for (auto iter : wts->evaluation_weights) {
	int lit = iter.first;
	mpf_set_q(mval, iter.second.get_mpq_t());
	evaluation_weights[lit] = Erd(mval);
    }

    smoothing_weights.clear();
    for (auto iter : wts->smoothing_weights) {
	int var = iter.first;
	mpf_set_q(mval, iter.second.get_mpq_t());
	smoothing_weights[var] = Erd(mval);
    }

    std::vector<Erd> rescale_weights;
    for (mpq_class qval : wts->rescale_weights) {
	mpf_set_q(mval, qval.get_mpq_t());
	rescale_weights.push_back(Erd(mval));
    }
    rescale = product_reduce(rescale_weights);
}

Erd Evaluator_erd::evaluate_edge(Egraph_edge &e) {
    if (e.has_zero)
	return Erd();
    arguments.clear();
    for (int lit : e.literals) 
	arguments.push_back(evaluation_weights[lit]);

    for (int v : e.smoothing_variables) 
	arguments.push_back(smoothing_weights[v]);

    Erd eval = product_reduce(arguments);
    if (verblevel >= 4) {
	mpf_t mval;
	mpf_init2(mval, 64);
	eval.get_mpf(mval);
	mp_exp_t exp;
	char *svalue = mpf_get_str(NULL, &exp, 10, 40, mval);
	report(4, "MPF: Evaluating edge (%d <-- %d).  Value = 0.%se%ld\n", e.to_id, e.from_id, svalue, exp);
	free(svalue);
	mpf_clear(mval);
    }
    return eval;
}

void Evaluator_erd::evaluate(mpf_class &count) {
    std::vector<Erd> operation_values;
    operation_values.resize(egraph->operations.size());
    for (int id = 1; id <= egraph->operations.size(); id++) {
	switch (egraph->operations[id-1].type) {
	case NNF_TRUE:
	case NNF_AND:
	    operation_values[id-1] = Erd(1.0);
	    break;
	case NNF_FALSE:
	case NNF_OR:
	default:
	    operation_values[id-1] = Erd(0.0);
	}
    }
    for (Egraph_edge e : egraph->edges) {
	Erd product = evaluate_edge(e);
	product = product.mul(operation_values[e.from_id-1]);
	bool multiply = egraph->operations[e.to_id-1].type == NNF_AND;
	if (multiply)
	    operation_values[e.to_id-1] = operation_values[e.to_id-1].mul(product);
	else
	    operation_values[e.to_id-1] = operation_values[e.to_id-1].add(product);
    }
    Erd ecount = operation_values[egraph->root_id-1];

    operation_values.clear();

    ecount = ecount.mul(rescale);
    ecount.get_mpf(count.get_mpf_t());

}

/*******************************************************************************************************************
Evaluation via Gnu multi-precision floating-point arithmetic
*******************************************************************************************************************/


Evaluator_mpf::Evaluator_mpf(Egraph *eg, Egraph_weights *wts) { 
    egraph = eg;

    /* Convert weight values from mpq to mpf */
    evaluation_weights.clear();
    for (auto iter : wts->evaluation_weights) {
	int lit = iter.first;
	evaluation_weights[lit] = iter.second;
    }

    smoothing_weights.clear();
    for (auto iter : wts->smoothing_weights) {
	int var = iter.first;
	smoothing_weights[var] = iter.second;
    }

    std::vector<Erd> rescale_weights;
    rescale = 1.0;
    for (mpq_class qval : wts->rescale_weights)
	rescale *= qval;

}

void Evaluator_mpf::evaluate_edge(mpf_class &value, Egraph_edge &e) {
    if (e.has_zero) {
	value = 0.0;
	return;
    }
    value = 1;
    for (int lit : e.literals)
	value *= evaluation_weights[lit];
    for (int v : e.smoothing_variables)
	value *= smoothing_weights[v];
    if (verblevel >= 4) {
	mp_exp_t exp;
	char *svalue = mpf_get_str(NULL, &exp, 10, 40, value.get_mpf_t());
	report(4, "MPF: Evaluating edge (%d <-- %d).  Value = 0.%se%ld\n", e.to_id, e.from_id, svalue, exp);
	free(svalue);
    }
}

void Evaluator_mpf::evaluate(mpf_class &count) {

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

	evaluate_edge(product, e);

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

    count *= rescale;

    if (verblevel >= 4) {
	mp_exp_t ecount;
	char *scount = mpf_get_str(NULL, &ecount, 10, 40, count.get_mpf_t());
	report(4, "MPF: Count = 0.%se%ld\n", scount, ecount);
	free(scount);
    }
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

Evaluator_mpq::Evaluator_mpq(Egraph *eg, Egraph_weights *wts) { 
    egraph = eg;
    weights = wts;
}
    
void Evaluator_mpq::clear_evaluation() {
    rescale = 1;
    max_bytes = 0;
}

void Evaluator_mpq::evaluate_edge(mpq_class &value, Egraph_edge &e) {
    if (e.has_zero) {
	value = 0.0;
	return;
    }
    std::vector<mpq_class> eval_queue;
    for (int lit : e.literals)
	eval_queue.push_back(weights->evaluation_weights[lit]);
    for (int v : e.smoothing_variables)
	eval_queue.push_back(weights->smoothing_weights[v]);
    reduce_product(value, eval_queue);
    if (verblevel >= 4) {
	char *svalue = mpq_get_str(NULL, 10, value.get_mpq_t());
	report(4, "MPQ: Evaluating edge (%d <-- %d).  Value = %s\n", e.to_id, e.from_id, svalue);
	free(svalue);
    }
    size_t bytes = mpq_bytes(value.get_mpq_t());
    if (bytes > max_bytes)
	max_bytes = bytes;
}

void Evaluator_mpq::evaluate(mpq_class &count) {
    clear_evaluation();
    reduce_product(rescale, weights->rescale_weights);
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
	evaluate_edge(product, e);

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

    count *= rescale;

    if (verblevel >= 4) {
	char *scount = mpq_get_str(NULL, 10, count.get_mpq_t());
	report(4, "MPQ: count = %s\n", scount);
	free(scount);
    }
}

/*******************************************************************************************************************
Evaluation via MPFI
*******************************************************************************************************************/

Evaluator_mpfi::Evaluator_mpfi(Egraph *eg, Egraph_weights *wts, bool instr) { 
    egraph = eg;

    /* Convert weight values from mpq to mpfi */
    int next_idx = 0;
    weight_count = wts->evaluation_weights.size() + wts->smoothing_weights.size();
    weights = new mpfi_t[weight_count];

    evaluation_index.clear();
    for (auto iter : wts->evaluation_weights) {
	int lit = iter.first;
	int idx = next_idx++;
	evaluation_index[lit] = idx;
	mpfi_init(weights[idx]);
	mpfi_set_q(weights[idx], iter.second.get_mpq_t());
    }

    smoothing_index.clear();
    for (auto iter : wts->smoothing_weights) {
	int var = iter.first;
	int idx = next_idx++;
	smoothing_index[var] = idx;
	mpfi_init(weights[idx]);
	mpfi_set_q(weights[idx], iter.second.get_mpq_t());
    }

    mpfi_init(rescale);
    mpfi_set_d(rescale, 1.0);
    for (mpq_class wt : wts->rescale_weights)
	mpfi_mul_q(rescale, rescale, wt.get_mpq_t());

    instrument = instr;
}
    
void Evaluator_mpfi::clear_evaluation() {
    min_digit_precision = MAX_DIGIT_PRECISION;
}

void Evaluator_mpfi::evaluate_edge(mpfi_ptr value, Egraph_edge &e) {
    if (e.has_zero) {
	mpfi_set_d(value, 0.0);
	return;
    }
    mpfi_set_d(value, 1.0);
    for (int lit : e.literals)
	mpfi_mul(value, value, weights[evaluation_index[lit]]);
    for (int v : e.smoothing_variables)
	mpfi_mul(value, value, weights[smoothing_index[v]]);
}

void Evaluator_mpfi::evaluate(mpfi_ptr count) {
    clear_evaluation();

    mpfi_t *operation_values = new mpfi_t[egraph->operations.size()];
    bool *operation_updated = new bool[egraph->operations.size()];
    for (int id = 1; id <= egraph->operations.size(); id++) {
	operation_updated[id-1] = false;
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
    int id = 0;
    for (Egraph_edge e : egraph->edges) {
	id++;
	mpfi_t product;
	mpfi_init(product);
	evaluate_edge(product, e);
	report(4, "Evaluated edge #%d (%d <-- %d)\n", id, e.to_id, e.from_id);
	mpfi_mul(product, product, operation_values[e.from_id-1]);
	if (operation_updated[e.to_id-1]) {
	    bool add = egraph->operations[e.to_id-1].type == NNF_OR;
	    if (add) {
		mpfi_add(operation_values[e.to_id-1], operation_values[e.to_id-1], product);
		if (instrument && add) {
		    double dp = digit_precision_mpfi(operation_values[e.to_id-1]);
		    if (dp < min_digit_precision)
			min_digit_precision = dp;
		}
	    } else 
		mpfi_mul(operation_values[e.to_id-1], operation_values[e.to_id-1], product);
	} else {
	    operation_updated[e.to_id-1] = true;
	    mpfi_swap(operation_values[e.to_id-1], product);
	}
	mpfi_clear(product);
    }
    mpfi_swap(count, operation_values[egraph->root_id-1]);
    double dp = digit_precision_mpfi(count);
    if (dp < min_digit_precision)
	min_digit_precision = dp;
    for (int id = 1; id <= egraph->operations.size(); id++)
	mpfi_clear(operation_values[id-1]);

    delete[] operation_values;
    delete[] operation_updated;

    for (int i = 0; i < weight_count; i++)
	mpfi_clear(weights[i]);

    delete[] weights;

    mpfi_mul(count, count, rescale);



}

/*******************************************************************************************************************
Evaluation.  When no negative weights, use MPI.  Otherwise, start with MFPI and switch to MPQ if needed
*******************************************************************************************************************/

/* Parameters */

// Don't attempt floating-point if it requires too many bits 
#define MPQ_THRESHOLD 1024

static const char* method_name[8] = 
    {"ERD", "MPF", "MPFI", "MPQ", "ERD_ONLY", "MPF_ONLY", "MPFI_ONLY", "MPQ_ABORT"};

Evaluator_combo::Evaluator_combo(Egraph *eg, Egraph_weights *wts, double tprecision, int bprecision, int instr) {
    egraph = eg;
    weights = wts;
    target_precision = tprecision;
    bit_precision = bprecision;
    instrument = instr;
    max_bytes = 24;
    erd_seconds = 0.0;
    mpf_seconds = 0.0;
    mpfi_seconds = 0.0;
    mpq_seconds = 0.0;
    mpq_count = 0.0;
    mpf_count = 0.0;
    erd_count = 0.0;
    mpfi_init(mpfi_count);
    mpfi_set_d(mpfi_count, 0.0);
    min_digit_precision = 0.0;
}

const char *Evaluator_combo::method() {
    return method_name[computed_method];
}

void Evaluator_combo::evaluate(mpf_class &count, bool no_mpq) {
    int constant = egraph->is_smoothed ? 4 : 7;
    if (bit_precision == 0)
	bit_precision = required_bit_precision(target_precision, egraph->nvar, constant,
					       weights->all_nonnegative);
    if (no_mpq) {
	computed_method = weights->all_nonnegative ? 
	    (bit_precision < 54 ? COMPUTE_ERD_NOMPQ : COMPUTE_MPF_NOMPQ)
	    : COMPUTE_MPFI_NOMPQ;
    } else
	computed_method = weights->all_nonnegative ? 
	    (bit_precision < 54 ? COMPUTE_ERD : COMPUTE_MPF)
	    : COMPUTE_MPFI;
    int save_precision = mpf_get_default_prec();
    max_bytes = 8 + bit_precision/8;
    if (bit_precision > MPQ_THRESHOLD)
	computed_method = COMPUTE_MPQ;
    report(3, "Achieving target precision %.1f with %d variables would require %d bit FP.  Starting with %s\n",
	   target_precision, egraph->nvar, bit_precision, method());

    double start_time = tod();
    switch (computed_method) {
    case COMPUTE_ERD:
    case COMPUTE_ERD_NOMPQ:
	{
	    max_bytes = 8;
	    Evaluator_erd ev = Evaluator_erd(egraph, weights);
	    ev.evaluate(count);
	    guaranteed_precision = digit_precision_bound(bit_precision, egraph->nvar, constant);
	    erd_seconds = tod() - start_time;
	    erd_count = count;
	}
	break;
    case COMPUTE_MPF:
    case COMPUTE_MPF_NOMPQ:
	{
	    mpf_set_default_prec(bit_precision);
	    Evaluator_mpf ev = Evaluator_mpf(egraph, weights);
	    ev.evaluate(count);
	    guaranteed_precision = digit_precision_bound(bit_precision, egraph->nvar, constant);
	    mpf_set_default_prec(save_precision);
	    mpf_seconds = tod() - start_time;
	    mpf_count = count;
	}
	break;
    case COMPUTE_MPFI:
    case COMPUTE_MPFI_NOMPQ:
	{
	    save_precision = mpfr_get_default_prec();
	    mpfr_set_default_prec(bit_precision);
	    mpfi_set_prec(mpfi_count, bit_precision);
	    max_bytes *= 2;
	    Evaluator_mpfi ev = Evaluator_mpfi(egraph, weights, instrument);
	    ev.evaluate(mpfi_count);
	    mpfi_seconds = tod() - start_time;
	    min_digit_precision = ev.min_digit_precision;
	    guaranteed_precision = digit_precision_mpfi(mpfi_count);
	    if (guaranteed_precision >= target_precision) {
		mpfr_t mpfr_count;
		mpfr_init(mpfr_count);
		mpfi_mid(mpfr_count, mpfi_count);
		mpf_t mpf_count;
		mpf_init2(mpf_count, bit_precision);
		mpfr_get_f(mpf_count, mpfr_count, MPFR_RNDN);
		count = (mpf_class) mpf_count;
		mpfr_set_default_prec(save_precision);		
	    } else if (no_mpq) {
		report(1, "After %.2f seconds, MPFI gave only guaranteed precision of %.1f.  Aborting\n",
		       tod() - start_time, guaranteed_precision);
		count = 0.0;
		computed_method = COMPUTE_MPQ_NOMPQ;
	    } else {
		mpfr_set_default_prec(save_precision);
		// Try again
		report(1, "After %.2f seconds, MPFI gave only guaranteed precision of %.1f.  Computing with MPQ\n",
		       tod() - start_time, guaranteed_precision);
		double start_time_mpq = tod();
		Evaluator_mpq ev = Evaluator_mpq(egraph, weights);
		ev.evaluate(mpq_count);
		computed_method = COMPUTE_MPQ;
		guaranteed_precision = MAX_DIGIT_PRECISION;
		mpq_seconds = tod() - start_time_mpq;
		mpf_t mpf_count;
		mpf_init2(mpf_count, bit_precision);
		mpf_set_q(mpf_count, mpq_count.get_mpq_t());
		count = (mpf_class) mpf_count;
		max_bytes = ev.max_bytes;
	    }
	}
	break;
    case COMPUTE_MPQ_NOMPQ:
	guaranteed_precision = 0.0;
	count = 0.0;
	break;
    case COMPUTE_MPQ:
	{
	    Evaluator_mpq ev = Evaluator_mpq(egraph, weights);
	    ev.evaluate(mpq_count);
	    guaranteed_precision = MAX_DIGIT_PRECISION;
	    mpf_t mpf_count;
	    mpf_init2(mpf_count, bit_precision);
	    mpf_set_q(mpf_count, mpq_count.get_mpq_t());
	    count = (mpf_class) mpf_count;
	    max_bytes = ev.max_bytes;
	    mpq_seconds = tod() - start_time;
	}
    }
    report(3, "Total time for evaluation %.2f seconds.  Method %s, Guaranteed precision %.1f\n",
	   tod() - start_time, method(), guaranteed_precision);
	   
}
