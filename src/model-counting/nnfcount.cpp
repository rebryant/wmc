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

#include "cnf_info.hh"
#include "egraph.hh"
#include "counters.h"
#include "report.h"
#include "analysis.h"

void usage(const char *name) {
    lprintf("Usage: %s [-h] [-s] [-v VERB] [-L LOG] [-o OUT.nnf] FORMULA.cnf FORMULA.nnf\n", name);
    lprintf("  -h          Print this information\n");
    lprintf("  -s          Use smoothing, rather than ring evaluation\n");
    lprintf("  -o OUT.nnf  Save copy of formula (including possible smoothing)\n");
    lprintf("  -v VERB     Set verbosity level\n");
    lprintf("  -L LOG      Record all results to file LOG\n");
}

const char *prefix = "c: CNT:";
bool smooth = false;
Cnf *cnf;

void run(FILE *cnf_file, FILE *nnf_file, FILE *out_file) {
    double start_time = tod();
    double end_time;
    cnf = new Cnf();
    cnf->import_file(cnf_file, true, false);
    Egraph eg = Egraph(cnf->data_variables);
    eg.read_nnf(nnf_file);
    if (smooth)
	eg.smooth();
    lprintf("%s     Reading files and constructing graph required %.3f seconds\n",
	    prefix, tod() - start_time);
    if (out_file)
	eg.write_nnf(out_file);
    if (!cnf->is_weighted()) {
	// Regular count
	q25_reset_counters();
	start_time = tod();
	Evaluator_q25 qev = Evaluator_q25(&eg);
	q25_ptr count = qev.evaluate(NULL, smooth);
	end_time = tod();
	long unweighted_operations = q25_operation_count();
	double peak_dbl_bytes = q25_peak_allocation_fp(false);
	double peak_q25_bytes = q25_peak_allocation_q25();
	double peak_mpq_bytes = q25_peak_allocation_mpq();
	char *scount = q25_string(count);
	lprintf("%s   UNWEIGHTED Q25 COUNT    = %s\n", prefix, scount);
	lprintf("%s     Unweighted Q25 count required %ld q25 operations, %.3f seconds, %.0f peak bytes\n",
		prefix, unweighted_operations, end_time - start_time, peak_q25_bytes);
	free(scount);
	qev.clear_evaluation();

	start_time = tod();
	Evaluator_mpq mpqev = Evaluator_mpq(&eg);
	mpq_class qcount = 0;
	if (mpqev.evaluate(qcount, NULL, smooth)) {
	    end_time = tod();
	    q25_ptr ccount = q25_from_mpq(qcount.get_mpq_t());
	    if (q25_compare(count, ccount) == 0) 
		lprintf("%s Q25 Count != MPQ Count\n", prefix);
	    else
		err(false, "Q25 Count != MPQ Count\n");
	    q25_free(ccount);
	    lprintf("%s     Unweighted MPQ count required %.3f seconds, %.0f peak bytes\n",
		    prefix, end_time - start_time, peak_mpq_bytes);
	    mpqev.clear_evaluation();
	} else {
	    lprintf("%s Calculation of unweighted count using mpq failed\n", prefix);
	}

	start_time = tod();
	Evaluator_double dev = Evaluator_double(&eg);
	double dcount = dev.evaluate(NULL, smooth);
	end_time = tod();
	double precision = digit_error_mix(count, dcount);
	lprintf("%s   UNWEIGHTED DBL COUNT    = %.1f (precision = %.3f)\n", prefix, dcount, precision);
	lprintf("%s     Unweighted DBL count required %.3f seconds, %.0f peak bytes\n",
		prefix, end_time - start_time, peak_dbl_bytes);
	q25_free(count);
	dev.clear_evaluation();

    } else {
	q25_reset_counters();
	start_time = tod();
	Evaluator_q25 qev = Evaluator_q25(&eg);
	q25_ptr wcount = qev.evaluate(cnf->input_weights, smooth);
	end_time = tod();
	long weighted_operations = q25_operation_count();
	double peak_dbl_bytes = q25_peak_allocation_fp(false);
	double peak_q25_bytes = q25_peak_allocation_q25();
	double peak_mpq_bytes = q25_peak_allocation_mpq();
	char *swcount = q25_string(wcount);
	lprintf("%s   WEIGHTED Q25 COUNT    = %s\n", prefix, swcount);
	lprintf("%s     Weighted Q25 count required %ld q25 operations, %.3f seconds, %.0f peak bytes\n",
		prefix, weighted_operations, end_time - start_time, peak_q25_bytes);
	free(swcount);
	qev.clear_evaluation();

	start_time = tod();
	Evaluator_mpq mpqev = Evaluator_mpq(&eg);
	mpq_class qwcount = 0;
	if (mpqev.evaluate(qwcount, cnf->input_weights, smooth)) {
	    end_time = tod();
	    q25_ptr cwcount = q25_from_mpq(qwcount.get_mpq_t());
	    if (q25_compare(wcount, cwcount) == 0) 
		lprintf("%s Q25 weighted count == MPQ weighted count\n", prefix);
	    else
		err(false, "Q25 weighted count != MPQ weighted count\n");
	    lprintf("%s     Weighted MPQ count required %.3f seconds, %.0f peak bytes\n",
		    prefix, end_time - start_time, peak_mpq_bytes);
	    q25_free(cwcount);
	    mpqev.clear_evaluation();
	} else {
	    lprintf("%s Calculation of weighted count using mpq failed\n", prefix);
	}
	start_time = tod();
	Evaluator_double dev = Evaluator_double(&eg);
	double dwcount = dev.evaluate(cnf->input_weights, smooth);
	double wprecision = digit_error_mix(wcount, dwcount);
	lprintf("%s   WEIGHTED DBL COUNT    = %.20g (precision = %.3f)\n", prefix, dwcount, wprecision);
	lprintf("%s     Weighted DBL count required %.3f seconds, %.0f peak bytes\n",
		prefix, tod() - start_time, peak_dbl_bytes);
	q25_free(wcount);
	dev.clear_evaluation();

    }
}

void report_stats() {
    int ndvar = cnf->data_variables->size();
    int sum_count = get_histo_count(HISTO_SUMS);
    long sum_ops = get_histo_total(HISTO_SUMS);
    int edge_product_count = get_histo_count(HISTO_EDGE_PRODUCTS);
    long edge_product_ops = get_histo_total(HISTO_EDGE_PRODUCTS);
    int node_product_count = get_histo_count(HISTO_NODE_PRODUCTS);
    long node_product_ops = get_histo_total(HISTO_NODE_PRODUCTS);
    int smoothing_count = get_histo_count(HISTO_EDGE_SMOOTHS);
    long smoothing_ops = get_histo_total(HISTO_EDGE_SMOOTHS);

    lprintf("%s   Data variables    : %d\n", prefix, ndvar);
    lprintf("%s   Operations \n", prefix);
    lprintf("%s     Sums            : %d\n", prefix, sum_count);
    lprintf("%s     Edge products   : %d\n", prefix, edge_product_count);
    lprintf("%s     Node Products   : %d\n", prefix, node_product_count);
    lprintf("%s     Smooth prods    : %d\n", prefix, smoothing_count);
    lprintf("%s     Operations TOTAL: %d\n", prefix, sum_count + edge_product_count + node_product_count + smoothing_count);
    lprintf("%s   Binary Operations \n", prefix);
    lprintf("%s     Sum ops         : %ld\n", prefix, sum_ops);
    lprintf("%s     Edge product ops: %ld\n", prefix, edge_product_ops);
    lprintf("%s     Node product ops: %ld\n", prefix, node_product_ops);
    lprintf("%s     Smooth prod ops : %ld\n", prefix, smoothing_ops);
    lprintf("%s     Binops  TOTAL   : %ld\n", prefix, sum_ops + edge_product_ops + node_product_ops + smoothing_ops);
    lprintf("%s   Graph bytes       : %ld\n", prefix,
	    sum_count + node_product_count + 8 * edge_product_count + 4 * edge_product_ops + 4 * smoothing_ops);
}

int main(int argc, char *argv[]) {
    int c;
    FILE *out_file = NULL;
    while ((c = getopt(argc, argv, "hsv:L:o:")) != -1) {
	switch(c) {
	case 'h':
	    usage(argv[0]);
	    exit(0);
	case 's':
	    smooth = true;
	    break;
	case 'v':
	    set_verblevel(atoi(optarg));
	    break;
	case 'L':
	    set_logname(optarg);
	    break;
	case 'o':
	    out_file = fopen(optarg, "w");
	    if (!out_file) {
		lprintf("Couldn't open output file '%s'\n", optarg);
		return 1;
	    }
	    break;
	default:
	    lprintf("Unknown commandline option '%c'\n", c);
	    usage(argv[0]);
	    return 1;
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
    run(cnf_file, nnf_file, out_file);
    double elapsed = tod() - start;
    report_stats();
    lprintf("%s   Elapsed seconds   : %.3f\n", prefix, elapsed);
    return 0;
}
