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

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <ctype.h>

#include "cnf.hh"
#include "q25.h"
#include "counters.h"
#include "report.h"

const char *prefix = "c: CNT:";

void usage(const char *name) {
    lprintf("Usage: %s [-h] [-v VERB] [-L LOG] FORMULA.cnf FORMULA.nnf\n", name);
    lprintf("  -h          Print this information\n");
    lprintf("  -v VERB     Set verbosity level\n");
    lprintf("  -L LOG      Record all results to file LOG\n");
}

// Evaluation data structures

typedef enum { NNF_NONE, NNF_TRUE, NNF_FALSE, NNF_AND, NNF_OR, NNF_NUM } nnf_type_t;
const char *nnf_type_name[NNF_NUM] = { "NONE", "TRUE", "FALSE", "AND", "OR" };
const char nnf_type_char[NNF_NUM] = { '\0', 't', 'f', 'a', 'o' };

struct nnf_node {
    int indegree;
    nnf_type_t type;
    // Density computation for unweighted count
    q25_ptr density;
    // Normalized weighted count
    q25_ptr weighted;
};

// Information about operation.
// Indexed by operation id-1
std::vector<nnf_node> operations;

Cnf *cnf;

std::vector<q25_ptr> norm_pos_weights;
std::vector<q25_ptr> norm_neg_weights;

// Current line number
int line_number = 0;
// ID of final edge destination
int root_id = 1;

// q25 Operation counts
long weighted_operations = 0;
long unweighted_operations = 0;
long start_count = 0;

// Routine to aid the management of q25_ptr's
std::vector<q25_ptr> q25_buffer;
q25_ptr qmark(q25_ptr q) {
    q25_buffer.push_back(q);
    return q;
}

void qflush() {
    for (q25_ptr q : q25_buffer)
	q25_free(q);
    q25_buffer.clear();
}


bool is_variable(int var) {
    return cnf->is_data_variable(var);
}

bool is_literal(int lit) {
    return lit < 0 ? is_variable(-lit) : is_variable(lit);
}

bool is_operation(int id) {
    return id > 0 && id <= operations.size();
}

void add_operation(int id, nnf_type_t type) {
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
	operations[id-1].density = q25_from_32(1);
	operations[id-1].weighted = q25_from_32(1);
	break;
    case NNF_FALSE:
    case NNF_OR:
	operations[id-1].density = q25_from_32(0);
	operations[id-1].weighted = q25_from_32(0);
	break;
    default:
	err(true, "Line %d.  Operation %d declared with unknown type %d\n", line_number, id, (int) type);
	break;
    }
}

q25_ptr reduce_edge(int from_id, int to_id, std::vector<int> lits, bool weighted) {
    if (!weighted && lits.size() > 0) {
	incr_histo(HISTO_EDGE_PRODUCTS, lits.size());
    }
    q25_ptr result = q25_from_32(1);
    if (lits.size() == 0) {
	// Done
    } else if (weighted) {
	for (int lit: lits) {
	    if (!is_literal(lit))
		err(true, "Line %d.  Invalid literal %d\n", line_number, lit);
	    q25_ptr lwt = lit < 0 ? norm_neg_weights[-lit-1] : norm_pos_weights[lit-1];
	    result = q25_mul(qmark(result), lwt);
	}
    } else
	result = q25_scale(qmark(result), -lits.size(), 0);
    qflush();
    if (verblevel >= 4) {
	char *sresult = q25_string(result);
	report(4, "Reducing edge %d <-- %s to %s\n", to_id, from_id);
	free(sresult);
    }
    return result;
}

void evaluate_edge(int from_id, int to_id, q25_ptr edge_value, bool weighted) {
    root_id = to_id;
    if (!is_operation(from_id))
	err(true, "Line %d.  Invalid source ID %d\n", line_number, from_id);
    if (!is_operation(to_id))
	err(true, "Line %d.  Invalid destination ID %d\n", line_number, to_id);
    if (operations[from_id-1].type == NNF_NONE)
	err(true, "Line %d.  Undefined operation for source operation ID %d\n", line_number, from_id);
    bool multiply = false;
    switch(operations[to_id-1].type) {
    case NNF_AND:
	multiply = true;
	break;
    case NNF_OR:
	break;
    default:
	err(true, "Line %d.  Operation %d is not a sum or product\n", line_number, to_id);
    }
    if (!weighted) {
	operations[to_id-1].indegree++;
	q25_ptr product = qmark(q25_mul(edge_value, operations[from_id-1].density));
	char *sfrom, *sold, *sedge, *sproduct;
	if (verblevel >= 4) {
	    sfrom = q25_string(operations[from_id-1].density);
	    sold = q25_string(operations[to_id-1].density);
	    sedge = q25_string(edge_value);
	    sproduct = q25_string(product);
	}
	start_count = q25_operation_count();
	operations[to_id-1].density =
	    multiply ?
	    q25_mul(product, qmark(operations[to_id-1].density)) :
	    q25_add(product, qmark(operations[to_id-1].density));
	unweighted_operations += q25_operation_count() - start_count;
	if (verblevel >= 4) {
	    char *snew = q25_string(operations[to_id-1].density);
	    report(4, "Density: Updating %d from %d.  %s * %s %c %s --> %s\n", sfrom, sedge, multiply ? '*' : '+', sold, snew);
	    free(sfrom); free(sold); free(sedge); free(sproduct); free(snew);
	}
    } else {
	q25_ptr product = qmark(q25_mul(edge_value, operations[from_id-1].weighted));
	char *sfrom, *sold, *sedge, *sproduct;
	if (verblevel >= 4) {
	    sfrom = q25_string(operations[from_id-1].weighted);
	    sold = q25_string(operations[to_id-1].weighted);
	    sedge = q25_string(edge_value);
	    sproduct = q25_string(product);
	}
	start_count = q25_operation_count();
	operations[to_id-1].weighted =
	    multiply ?
	    q25_mul(product, qmark(operations[to_id-1].weighted)) :
	    q25_add(product, qmark(operations[to_id-1].weighted));
	weighted_operations += q25_operation_count() - start_count;
	if (verblevel >= 4) {
	    char *snew = q25_string(operations[to_id-1].weighted);
	    report(4, "Weighted: Updating %d from %d.  %s * %s %c %s --> %s\n", sfrom, sedge, multiply ? '*' : '+', sold, snew);
	    free(sfrom); free(sold); free(sedge); free(sproduct); free(snew);
	}
    }
    qflush();
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

void read_nnf(FILE *infile) {
    operations.clear();
    line_number = 0;
    bool do_weighted = norm_pos_weights.size() > 0;
    // Capture arguments for each line
    std::vector<int> largs;
    // List of IDs for edge values
    std::vector<int> edge_ids;

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
	    edge_ids.clear();
	    for (int i = 2; i < largs.size()-1; i++)
		edge_ids.push_back(largs[i]);
	    // Unweighted
	    q25_ptr edge_value = reduce_edge(from_id, to_id, edge_ids, false);
	    evaluate_edge(from_id, to_id, edge_value, false);
	    q25_free(edge_value);
	    // Weighted
	    if (do_weighted) {
		edge_value = reduce_edge(from_id, to_id, edge_ids, true);
		evaluate_edge(from_id, to_id, edge_value, true);
		q25_free(edge_value);
	    }
	    report(4, "Evaluated edge from %d to %d\n", from_id, to_id);
	}
    }
    // Check over the operations
    for (int id = 1; id <= operations.size(); id++) {
	if (operations[id-1].indegree > 1) {
	    if (operations[id-1].type == NNF_AND) {
		incr_histo(HISTO_NODE_PRODUCTS, operations[id-1].indegree-1);
	    } else {
		incr_histo(HISTO_SUMS, operations[id-1].indegree-1);
	    }
	}

    }
}

q25_ptr prepare_weights() {
    if (cnf->input_weights->size() == 0)
	return NULL;
    q25_ptr rescale = q25_from_32(1);
    int nvar = cnf->variable_count();
    norm_pos_weights.resize(nvar, NULL);
    norm_neg_weights.resize(nvar, NULL);
    for (auto iter : *cnf->input_weights) {
	int lit = iter.first;
	q25_ptr weight = iter.second;
	if (lit > 0)
	    norm_pos_weights[lit-1] = weight;
	else
	    norm_neg_weights[-lit-1] = weight;
    }
    for (int v = 1; v <= nvar; v++) {
	if (!is_variable(v))
	    continue;
	q25_ptr pwt = norm_pos_weights[v-1];
	q25_ptr nwt = norm_neg_weights[v-1];
	q25_ptr sum = NULL;
	if (!pwt && !nwt) {
	    err(false, "No weights provided from variable %d.  Set both phases to weight 1.0\n", v);
	    pwt = qmark(q25_from_32(1));
	    nwt = qmark(q25_from_32(1));
	    sum = qmark(q25_from_32(2));
	} else if (!pwt) {
	    pwt = q25_one_minus(nwt);
	    sum = qmark(q25_from_32(1));
	} else if (!nwt) {
	    nwt = q25_one_minus(pwt);
	    sum = qmark(q25_from_32(1));
	} else
	    sum = qmark(q25_add(pwt, nwt));
	if (q25_is_one(sum)) {
	    norm_pos_weights[v-1] = pwt;
	    norm_neg_weights[v-1] = nwt;
	} else {
	    q25_ptr recip = qmark(q25_recip(sum));
	    if (!q25_is_valid(recip)) {
		char *srecip = q25_string(recip);
		err(true, "Could not get reciprocal of summed weights for variable %d.  Sum = s\n", v, srecip);
		free(srecip);
	    }
	    rescale = q25_mul(qmark(rescale), sum);
	    norm_pos_weights[v-1] = q25_mul(qmark(pwt), recip);
	    norm_neg_weights[v-1] = q25_mul(qmark(nwt), recip);
	}
    }
    qflush();
    if (verblevel >= 3) {
	char *srescale = q25_string(rescale);
	report(3, "Weighted rescale = %s\n", srescale);
	free(srescale);
    }
    return rescale;
}

q25_ptr prepare_count() {
    int vcount = cnf->data_variables->size();
    q25_ptr rescale = q25_scale(qmark(q25_from_32(1)), vcount, 0);
    qflush();
    if (verblevel >= 3) {
	char *srescale = q25_string(rescale);
	report(3, "Count rescale = %s\n", srescale);
	free(srescale);
    }
    return rescale;
}

void run(FILE *cnf_file, FILE *nnf_file) {
    cnf = new Cnf();
    cnf->import_file(cnf_file, true, false);
    start_count = q25_operation_count();
    q25_ptr wrescale = prepare_weights();
    weighted_operations += q25_operation_count() - start_count;
    read_nnf(nnf_file);
    start_count = q25_operation_count();
    q25_ptr crescale = prepare_count();
    q25_ptr count = q25_mul(qmark(crescale), operations[root_id-1].density);
    unweighted_operations += q25_operation_count() - start_count;
    char *scount = q25_string(count);
    lprintf("%s   UNWEIGHTED COUNT = %s\n", prefix, scount);
    lprintf("%s   Unweighted count required %ld q25 operations\n",
	    prefix, unweighted_operations);
    free(scount);
    if (wrescale) {
	q25_ptr wcount = q25_mul(qmark(wrescale), operations[root_id-1].weighted);
	char *swcount = q25_string(wcount);
	lprintf("%s   WEIGHTED COUNT   = %s\n", prefix, swcount);
	lprintf("%s   Weighted count required %ld q25 operations\n",
		prefix, weighted_operations);
	free(swcount);
    }
    qflush();
}

void report_stats() {
    int ndvar = cnf->data_variables->size();
    int sum_count = get_histo_count(HISTO_SUMS);
    long sum_ops = get_histo_total(HISTO_SUMS);
    int edge_product_count = get_histo_count(HISTO_EDGE_PRODUCTS);
    long edge_product_ops = get_histo_total(HISTO_EDGE_PRODUCTS);
    int node_product_count = get_histo_count(HISTO_NODE_PRODUCTS);
    long node_product_ops = get_histo_total(HISTO_NODE_PRODUCTS);

    lprintf("%s   Data variables    : %d\n", prefix, ndvar);
    lprintf("%s   Operations \n", prefix);
    lprintf("%s     Sums            : %d\n", prefix, sum_count);
    lprintf("%s     Edge products   : %d\n", prefix, edge_product_count);
    lprintf("%s     Node Products   : %d\n", prefix, node_product_count);
    lprintf("%s     Operations TOTAL: %d\n", prefix, sum_count + edge_product_count + node_product_count);
    lprintf("%s   Binary Operations \n", prefix);
    lprintf("%s     Sums            : %ld\n", prefix, sum_ops);
    lprintf("%s     Edge products   : %ld\n", prefix, edge_product_ops);
    lprintf("%s     Node Products   : %ld\n", prefix, node_product_ops);
    lprintf("%s     Binops  TOTAL   : %ld\n", prefix, sum_ops + edge_product_ops + node_product_ops);
}

int main(int argc, char *argv[]) {
    int c;
    while ((c = getopt(argc, argv, "hv:L:")) != -1) {
	switch(c) {
	case 'h':
	    usage(argv[0]);
	    exit(0);
	case 'v':
	    set_verblevel(atoi(optarg));
	    break;
	case 'L':
	    set_logname(optarg);
	    break;
	default:
	    lprintf("Unknown commandline option '%c'\n", c);
	    usage(argv[0]);
	    break;
	}
    }
    int argi = optind;
    if (argi >= argc) {
	lprintf("Name of input CNF file required\n");
	usage(argv[0]);
	return 1;
    }
    const char *cnf_name = argv[argi++];
    if (argi >= argc) {
	lprintf("Name of output CNF file required\n");
	usage(argv[0]);
	return 1;
    }
    FILE *cnf_file = fopen(cnf_name, "r");
    if (!cnf_file)
	err(true, "Couldn't open CNF file %s", cnf_name);

    const char *nnf_name = argv[argi++];
    if (argi < argc) {
	lprintf("Unknown argument '%s'\n", argv[argi]);
	usage(argv[0]);
	return 1;
    }

    FILE *nnf_file = fopen(nnf_name, "r");
    if (!nnf_file)
	err(true, "Couldn't open NNF file %s", nnf_name);

    double start = tod();
    run(cnf_file, nnf_file);
    double elapsed = tod() - start;
    report_stats();
    lprintf("%s   Elapsed seconds   : %.3f\n", prefix, elapsed);
    return 0;
}
