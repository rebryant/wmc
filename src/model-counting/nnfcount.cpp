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
    lprintf("Usage: %s [-h] [-s] [-Q] [-v VERB] [-o OUT.nnf] FORMULA.nnf FORMULA_1.cnf ... FORMULA_k.cnf\n", name);
    lprintf("  -h          Print this information\n");
    lprintf("  -Q          Include evaluation using q25\n");
    lprintf("  -s          Use smoothing, rather than ring evaluation\n");
    lprintf("  -v VERB     Set verbosity level\n");
    lprintf("  -o OUT.nnf  Save copy of formula (including possible smoothing)\n");

}

#define BUFLEN 1024
char namebuf[BUFLEN];

char *change_extension(const char *fname, const char *ext) {
    strncpy(namebuf, fname, BUFLEN-1);
    namebuf[BUFLEN-1] = 0;
    /* Find last slash */
    int lpos = strnlen(namebuf, BUFLEN);
    for (; lpos > 0 && namebuf[lpos-1] != '/'; lpos--)
	;
    /* Find last dot */
    int rpos = strnlen(namebuf, BUFLEN)-1;
    for (; rpos > lpos && namebuf[rpos] != '.'; rpos--)
	;
    if (rpos > lpos)
	namebuf[rpos] = 0;
    strncpy(&namebuf[rpos], ext, 10);
    return &namebuf[lpos];
}

const char *prefix = "c: CNT:";
bool smooth = false;
bool use_q25 = false;
Egraph *eg;
Cnf *core_cnf = NULL;
double setup_time = 0;

void setup(FILE *cnf_file, FILE *nnf_file, FILE *out_file) {
    double start_time = tod();
    core_cnf = new Cnf();
    core_cnf->import_file(cnf_file, true, false);
    eg = new Egraph(core_cnf->data_variables);
    eg->read_nnf(nnf_file);
    if (smooth)
	eg->smooth();
    setup_time = tod() - start_time;
    if (out_file)
	eg->write_nnf(out_file);
}

void run(const char *cnf_name) {
    FILE *cnf_file = fopen(cnf_name, "r");
    if (!cnf_file) {
	err(false, "Couldn't open file '%s'.  Skipping\n", cnf_name);
	return;
    }
    lprintf("%s     Reading files and constructing graph required %.3f seconds\n",
	    prefix, setup_time);
    lprintf("%s     Using weights from file '%s'\n", prefix, cnf_name);
    double start_time = tod();
    double end_time;
    Cnf *local_cnf = new Cnf();
    local_cnf->import_file(cnf_file, true, false);
    fclose(cnf_file);
    
    std::unordered_map<int,const char*> *input_weights = NULL;
    const char *wlabel = "UNWEIGHTED";
    if (local_cnf->is_weighted()) {
	input_weights = local_cnf->input_weights;
	wlabel = "WEIGHTED";
    }

    q25_ptr wcount = q25_from_32(0);
    mpq_class qwcount = 0;
    if (use_q25) {
	Evaluator_q25 qev = Evaluator_q25(eg);
	long weighted_operations = 0;
	double peak_q25_bytes = 0;
	double max_q25_bytes = 0;
	start_time = tod();
	wcount = qev.evaluate(input_weights, smooth);
	end_time = tod();
	q25_ptr rwcount = q25_round(wcount, 50);
	char *swcount = q25_best_string(rwcount);
	char cmp = q25_compare(wcount, rwcount) == 0 ? ' ' : '~';
	q25_free(rwcount);
	lprintf("%s   %s Q25 COUNT   %c= %s\n", prefix, wlabel, cmp, swcount);
	lprintf("%s     Q25 required %ld q25 operations, %.3f seconds, %.0f peak (%.0f max) bytes\n",
		prefix, weighted_operations, end_time - start_time, peak_q25_bytes, max_q25_bytes);
	free(swcount);
	qev.clear_evaluation();
    }
    start_time = tod();
    Evaluator_mpq mpqev = Evaluator_mpq(eg);
    if (mpqev.evaluate(qwcount, input_weights, smooth)) {
	end_time = tod();
	if (use_q25) {
	    q25_ptr cwcount = q25_from_mpq(qwcount.get_mpq_t());
	    if (q25_compare(wcount, cwcount) == 0) 
		lprintf("%s   MPQ weighted count == Q25 weighted count\n", prefix);
	    else
		err(false, "Q25 weighted count != MPQ weighted count\n");
	    q25_free(cwcount);
	} else {
	    mpf_t fw;
	    mpf_init(fw);
	    mpf_set_q(fw, qwcount.get_mpq_t());
	    const char *swcount = mpf_string(fw);
	    lprintf("%s   %s MPQ COUNT    = %s\n", prefix, wlabel, swcount);
	    mpf_clear(fw);
	}
	lprintf("%s     MPQ required %.3f seconds, %d max bytes\n",
		prefix, end_time - start_time, mpqev.max_bytes);
	mpqev.clear_evaluation();
    } else {
	lprintf("%s Calculation of weighted count using mpq failed\n", prefix);
    }
    q25_free(wcount);

    start_time = tod();
    Evaluator_mpf mpfev = Evaluator_mpf(eg);
    mpf_class fcount = 0;
    if (mpfev.evaluate(fcount, input_weights, smooth)) {
	end_time = tod();
	double precision = digit_precision_mpf(fcount.get_mpf_t(), qwcount.get_mpq_t());
	const char *sfcount = mpf_string(fcount.get_mpf_t());
	lprintf("%s   %s MPF COUNT    = %s   precision = %.3f\n", prefix, wlabel, sfcount, precision);
	lprintf("%s     MPF required %.3f seconds\n",
		prefix, end_time - start_time);
	mpfev.clear_evaluation();
    } else {
	lprintf("%s Calculation of weighted count using mpf failed\n", prefix);
    }
    start_time = tod();
    Evaluator_double dev = Evaluator_double(eg);
    double dwcount = dev.evaluate(input_weights, smooth);
    end_time = tod();
    double wprecision = digit_precision_d(dwcount, qwcount.get_mpq_t());
    lprintf("%s   %s DBL COUNT    = %.20g   precision = %.3f\n", prefix, wlabel, dwcount, wprecision);
    lprintf("%s     DBL required %.3f seconds\n",
	    prefix, end_time - start_time);
    dev.clear_evaluation();

    start_time = tod();
    Evaluator_mpfi mpfiev = Evaluator_mpfi(eg, 20, false);
    mpfi_t icount;
    mpfi_init(icount);
    if (mpfiev.evaluate(icount, input_weights, smooth)) {
	end_time = tod();
	double est_precision = digit_precision_mpfi(icount);
	mpfr_t mid;
	mpfr_init(mid);
	mpfi_mid(mid, icount);
	double actual_precision = digit_precision_mpfr(mid, qwcount.get_mpq_t());
	const char *sicount = mpfr_string(mid);
	lprintf("%s   %s MPFI COUNT   = %s   precision est = %.3f actual = %.3f\n", prefix, wlabel, sicount,
		est_precision, actual_precision);
	lprintf("%s     MPFI required %.3f seconds\n",
		prefix, end_time - start_time);
	lprintf("%s     MPFI had %ld precision failures and a minimum precision of %.3f\n",
		prefix, mpfiev.precision_failure_count, mpfiev.min_digit_precision);
	mpfiev.clear_evaluation();
    } else {
	lprintf("%s Calculation of weighted count using mpfi failed\n", prefix);
    }
    delete local_cnf;
}

void report_stats() {
    int ndvar = core_cnf->data_variables->size();
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
    int mpf_precision = 128;
    FILE *out_file = NULL;
    while ((c = getopt(argc, argv, "hQsv:o:")) != -1) {
	switch(c) {
	case 'h':
	    usage(argv[0]);
	    exit(0);
	case 'Q':
	    use_q25 = true;
	    break;
	case 's':
	    smooth = true;
	    break;
	case 'v':
	    set_verblevel(atoi(optarg));
	    break;
	case 'o':
	    out_file = fopen(optarg, "w");
	    if (!out_file) {
		printf("Couldn't open output file '%s'\n", optarg);
		return 1;
	    }
	    break;
	default:
	    printf("Unknown commandline option '%c'\n", c);
	    usage(argv[0]);
	    return 1;
	    break;
	}
    }
    int argi = optind;
    const char *nnf_name = argv[argi++];
    if (argi >= argc) {
	printf("Name of input NNF file required\n");
	usage(argv[0]);
	return 1;
    }
    FILE *nnf_file = fopen(nnf_name, "r");
    if (!nnf_file)
	err(true, "Couldn't open NNF file '%s'\n", nnf_name);

    const char *cnf_name = argv[argi];
    if (argi >= argc) {
	printf("Name of input CNF file required\n");
	usage(argv[0]);
	return 1;
    }
    FILE *cnf_file = fopen(cnf_name, "r");
    if (!cnf_file)
	err(true, "Couldn't open CNF file '%s'\n", cnf_name);

    mpf_set_default_prec(mpf_precision);
    mpfr_set_default_prec(mpf_precision);

    double start = tod();
    setup(cnf_file, nnf_file, out_file);
    fclose(cnf_file);

    while (argi < argc) {
	const char *cnf_name = argv[argi++];
	printf("\n");
	char *lname = change_extension(cnf_name, smooth ? ".scount" : ".count");
	printf("Saving results in '%s'\n", lname);
	set_logname(lname);
	run(cnf_name);
	report_stats();
	set_logname(NULL);
    }
    double elapsed = tod() - start;
    printf("%s   Elapsed seconds   : %.3f\n", prefix, elapsed);
    return 0;
}
