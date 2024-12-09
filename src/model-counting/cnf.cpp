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

#include <cstdlib>
#include <cstdio>
#include <ctype.h>
#include <cstring>
#include <queue>
#include <algorithm>

#include "report.h"
#include "counters.h"
#include "cnf.hh"
// From glucose
#include "Solver.h"

const char *flag_char[2] = {
    "axbsp",
    "AXBSP"
};

int parse_classify_flag(char *str) {
    int flag = 0;
    int c;
    while ((c = *str++) != 0) {
	bool found = false;
	for (int pos = 0; pos < classify_pos_count; pos++) {
	    if (c == flag_char[0][pos])
		found = true;
	    else if (c == flag_char[1][pos]) {
		found = true;
		flag |= 1 << pos;
	    }
	}
	if (!found)
	    return -1;
    }
    return flag;
}

void gen_flag_string(int flag, char *dest) {
    for (int pos = 0; pos < classify_pos_count; pos++) {
	char c = flag & (1 << pos) ? flag_char[1][pos] : flag_char[0][pos];
	*dest++ = c;
    }
    *dest = '\0';
}

// Experimental version of BVE used randomization to break ties with variable selection
// Interesting and elegant, but not an improvement.
// Disabled here

#define RANDOM_BVE 0
#define STACK_BVE 0

// Implementation of FIFO queues that don't store duplicates
template <typename T> class unique_queue {
private:
    std::queue<T>q;
    std::unordered_set<T>elements;

    void quick_push(T val) { q.push(val); elements.insert(val); }

public:

    unique_queue() {}

    unique_queue(std::set<T> &vals) {
	for (T val : vals)
	    quick_push(val);
    }

    unique_queue(std::unordered_set<T> &vals) {
	for (T val : vals)
	    quick_push(val);
    }

    unique_queue(std::set<T> *vals) {
	for (T val : *vals)
	    quick_push(val);
    }

    unique_queue(std::unordered_set<T> *vals) {
	for (T val : *vals)
	    quick_push(val);
    }

    
    bool push(T val) {
	bool new_val = !is_member(val);
	if (new_val)
	    quick_push(val);
	return new_val;
    }

    bool empty() {
	return q.empty();
    }

    bool is_member(T val) {
	return elements.find(val) != elements.end();
    }

    T get_and_pop() {
	T val = q.front();
	elements.erase(val);
	q.pop();
	return val;
    }

};

#if RANDOM_BVE
///////////////////////////////////////////////////////////////////
// Managing pairs of ints packed into 64-bit word
// Use these both as edge identifier
// and as cost+unique value
//
// When combining cost functions with unique identifiers
// Have upper 32 bits represent cost
// and lower 32 bits represent unique value, 
//   assigned either sequentially or via a pseudo-random sequence
///////////////////////////////////////////////////////////////////

/*
  32-bit words packed into pairs
 */
static int64_t pack(int upper, int lower) {
    return ((int64_t) upper << 32) | lower;
}


static int64_t ordered_pack(int x1, int x2) {
    return x1 < x2 ? pack(x1, x2) : pack(x2, x1);
}

static int upper(int64_t pair) {
    return (int) (pair>>32);
}

static int lower(int64_t pair) {
    return (int) (pair & 0xFFFFFFFF);
}
#endif

// A handy class for generating random numbers with own seed
// Pseudo-random number generator based on the Lehmer MINSTD RNG
// Use to both randomize selection and to provide a unique value
// for each cost.
class Sequencer {

public:
    Sequencer(int s) { seed = s; }

    Sequencer(void) { seed = default_seed; }

    Sequencer(Sequencer &s) { seed = s.seed; }

    void set_seed(uint64_t s) { seed = s == 0 ? 1 : s; next() ; next(); }

    uint32_t next() { seed = (seed * mval) % groupsize; return seed; }

    // Return next pseudo random value in interval [0.0, 1.0)
    double pseudo_double() { uint32_t val = next(); return (double) val / (double) groupsize; }

    // Return next pseudo random value in range [0, m)
    int pseudo_int(int m) { return (int) (m * pseudo_double()); }

private:
    uint64_t seed;
    const uint64_t mval = 48271;
    const uint64_t groupsize = 2147483647LL;
    const uint64_t default_seed = 123456;
};


//////////////// Reading CNF FILE ///////////////////

// Put literals in ascending order of the variables
static bool abs_less(int x, int y) {
    return IABS(x) < IABS(y);
}


static int skip_line(FILE *infile) {
    int c;
    while ((c = getc(infile)) != EOF) {
	if (c == '\n')
	    return c;
    }
    return c;
}

// Skip over spaces, newlines, etc., until find something interesting
// Return last character encountered
static int find_token(FILE *infile) {
    int c;
    while ((c = getc(infile)) != EOF) {
	if (!isspace(c)) {
	    ungetc(c, infile);
	    break;
	}
    }
    return c;
}

// Read string token:
// Skip over spaces.
// Read contiguous non-space characters and store in dest.
// Set len to number of characters read.
// Return false if EOF encountered without getting string
static bool find_string_token(FILE *infile, char *dest, int maxlen, int *len) {
    int c;
    int rlen = 0;
    while ((c = getc(infile)) != EOF && rlen < maxlen-1) {
	if (isspace(c)) {
	    if (rlen > 0) {
		// Past token
		ungetc(c, infile);
		break;
	    }
	} else {
	    *(dest+rlen) = c;
	    rlen++;
	}
    }
    *(dest+rlen) = '\0';
    *len = rlen;
    return (c != EOF);
}

// Process comment, looking additional data variables & weights
// Return last character
static void process_comment(FILE *infile, std::unordered_set<int> *data_variables,
			    std::unordered_set<int> *tseitin_variables,
			    std::unordered_map<int,q25_ptr> *input_weights) {
    char buf[50];
    int len;
    if (find_string_token(infile, buf, 50, &len) && len == 1 && strncmp(buf, "p", 1) == 0
	&& find_string_token(infile, buf, 50, &len)) {
	bool show = true;
	if (len == 4 && ((show = (strncmp(buf, "show", 4) == 0))
			 || strncmp(buf, "tseitin", 7) == 0)) {
	    int var = -1;
	    std::unordered_set<int> *buf = show ? data_variables : tseitin_variables;
	    while (var != 0) {
		if (fscanf(infile, "%d", &var) != 1) {
		    err(false, "Couldn't read %s variable\n", show ? "data" : "Tseitin");
		    break;
		} else if (var != 0)
		    buf->insert(var);
	    }
	}
	else if (len == 6 && strncmp(buf, "weight", 6) == 0) {
	    int lit = 0;
	    if (fscanf(infile, "%d", &lit) != 1) {
		err(false, "Couldn't read weight literal (skipping)\n");
		skip_line(infile);
		return;
	    }
	    find_token(infile);
	    q25_ptr wt = q25_read(infile);
	    if (!q25_is_valid(wt)) {
		err(false, "Couldn't read weight for literal %d (skipping)\n", lit);
		skip_line(infile);
		return;
	    }
	    (*input_weights)[lit] = wt;
	    int zero;
	    if (fscanf(infile, "%d", &zero) != 1 || zero != 0) {
		err(false, "Couldn't read terminating zero in weight declaration for literal %d (accepting weight)\n", lit);
	    }

	}
    }
    skip_line(infile);
}		

Cnf::Cnf() {
    variable_type = NULL;
    data_variables = NULL;
    tseitin_variables = NULL;
    active_clauses = NULL;
    literal_clauses = NULL;
    input_weights = NULL;
    sat_elapsed = 0.0;
    initialize(0);
    promotion_try_count = promotion_success_count = 0;
}

Cnf::Cnf(int input_count) {
    variable_type = NULL;
    data_variables = NULL;
    tseitin_variables = NULL;
    active_clauses = NULL;
    literal_clauses = NULL;
    input_weights = NULL;
    sat_elapsed = 0.0;
    initialize(input_count);
    promotion_try_count = promotion_success_count = 0;
}

Cnf::~Cnf() {
    deallocate();
}

void Cnf::initialize(int input_count) {
    nvar = input_count;
    if (variable_type)
	delete variable_type;
    variable_type = new var_t[nvar];
    for (int v = 1; v <= nvar; v++)
	set_variable_type(v, VAR_UNUSED);
    clause_offset.clear();
    literal_sequence.clear();
    if (data_variables) {
	for (int v: *data_variables)
	    set_variable_type(v, VAR_DATA);
    } else
	data_variables = new std::unordered_set<int>;
    if (tseitin_variables) {
	for (int v: *tseitin_variables)
	    set_variable_type(v, VAR_TSEITIN_DETECT);
    } else
	tseitin_variables = new std::unordered_set<int>;
    if (active_clauses) {
	active_clauses->clear();
    } else
	active_clauses = new std::set<int>;
    if (literal_clauses) {
	literal_clauses->clear();
    } else
	literal_clauses = new std::unordered_map<int,std::unordered_set<int>>;
    if (!input_weights)
	input_weights = new std::unordered_map<int, q25_ptr>;
    new_clause();
    has_conflict = false;
    action_stack.clear();
    new_context();
    unit_literals.clear();
    bcp_unit_literals.clear();
    uquantified_variables.clear();
}

// Must explicitly deallocate sets
void Cnf::deallocate() {
    delete variable_type;
    delete data_variables;
    delete tseitin_variables;
    delete active_clauses;
    delete literal_clauses;
    delete input_weights;
    for (auto iter : defining_variables)
	delete iter.second;
    for (auto iter : defining_clauses)
	delete iter.second;
}

int Cnf::new_clause() {
    int cid = clause_offset.size();
    clause_offset.push_back(literal_sequence.size());
    if (cid > 0)
	active_clauses->insert(cid);
    return cid;
}

void Cnf::add_literal(int lit) {
    literal_sequence.push_back(lit);
    clause_offset.back() ++;
    int cid = clause_offset.size()-1;
    (*literal_clauses)[lit].insert(cid);
    int var = IABS(lit);
    if (get_variable_type(var) == VAR_UNUSED)
	set_variable_type(var, VAR_NONTSEITIN);
}

void Cnf::finish() {
    report(3, "CNF representation with %d inputs and %d clauses constructed\n",
	   variable_count(), maximum_clause_id());
}

 bool Cnf::import_file(FILE *infile, bool process_comments, bool skip_clauses) { 
    int expectedNclause = 0;
    bool read_failed = false;
    bool got_header = false;
    int c;
    bool eof = false;
    // Look for CNF header
    while ((c = getc(infile)) != EOF) {
	if (isspace(c)) 
	    continue;
	if (c == 'c') {
	    if (process_comments)
		process_comment(infile, data_variables, tseitin_variables, input_weights);
	    else
		skip_line(infile);
	    continue;
	}
	if (c == EOF) {
	    err(false, "Not valid CNF file.  No header line found\n");
	    return false;
	}
	if (c == 'p') {
	    char field[20];
	    if (fscanf(infile, "%s", field) != 1) {
		err(false, "Not valid CNF file.  Invalid header line\n");
		return false;
	    }
	    if (strcmp(field, "cnf") != 0) {
		err(false, "Not valid CNF file.  Header line shows type is '%s'\n", field);
		return false;
	    }
	    if (fscanf(infile, "%d %d", &nvar, &expectedNclause) != 2) {
		err(false, "Invalid CNF header\n");
		return false;
	    } 
	    initialize(nvar);
	    c = skip_line(infile);
	    got_header = true;
	    break;
	}
	if (c == EOF) {
	    err(false, "Invalid CNF.  EOF encountered before reading any clauses\n");
	    return false;
	}
    }
    if (!got_header) {
	err(false, "Not valid CNF.  No header line found\n");
	return false;
    }
    int clause_count = 0;
    while (clause_count < expectedNclause) {
	// Setup next clause
	if (!skip_clauses)
	    new_clause();
	bool starting_clause = true;
	while (true) {
	    int lit;
	    int c = find_token(infile);
	    if (c == EOF) {
		err(false, "Unexpected end of file\n");
		return false;
	    } else if (c == 'c' && starting_clause) {
		c = getc(infile);
		if (process_comments)
		    process_comment(infile, data_variables, tseitin_variables, input_weights);
		else
		    skip_line(infile);
		continue;
	    }
	    else if (fscanf(infile, "%d", &lit) != 1) {
		err(false, "Couldn't find literal or 0\n");
		return false;
	    }
	    if (lit == 0) {
		clause_count++;
		break;
	    }
	    else if (!skip_clauses)
		add_literal(lit);
	    starting_clause = false;
	}
    }
    while ((c = getc(infile)) != EOF) {
	if (isspace(c)) 
	    continue;
	if (c == 'c') {
	    if (process_comments)
		process_comment(infile, data_variables, tseitin_variables, input_weights);
	    else
		skip_line(infile);

	}
    }
    // If no data variables declared, assume all input variables are data variables
    if (data_variables->size() == 0) {
	for (int v = 1; v <= variable_count(); v++)
	    data_variables->insert(v);
    }
    for (int v : *data_variables)
	set_variable_type(v, VAR_DATA);
    incr_count_by(COUNT_INPUT_CLAUSE, maximum_clause_id());
    return true;
}

int Cnf::clause_length(int cid) {
    if (cid < 1 || cid > maximum_clause_id())
	err(true, "Invalid clause ID: %d\n", cid);
    return clause_offset[cid] - clause_offset[cid-1];
}

int Cnf::maximum_clause_id() {
    return clause_offset.size() - 1;
}

int Cnf::nonunit_clause_count() {
    return active_clauses->size();
}

int Cnf::current_clause_count() {
    return active_clauses->size() + bcp_unit_literals.size();
}

int Cnf::get_literal(int cid, int lid) {
    int len = clause_length(cid);
    int offset = clause_offset[cid-1];
    int lit = 0;
    if (lid >= 0 && lid < len)
	lit = literal_sequence[offset+lid];
    else
	err(true, "Invalid literal index %d for clause #%d.  Clause length = %d\n", lid, cid, len);
    return lit;
}

void Cnf::swap_literals(int cid, int i, int j) {
    int offset = clause_offset[cid-1];
    int tlit = literal_sequence[offset+i];
    literal_sequence[offset+i] = literal_sequence[offset+j];
    literal_sequence[offset+j] = tlit;
}

bool Cnf::show(FILE *outfile) {
    for (int lit : bcp_unit_literals)
	fprintf(outfile, "  UNIT: %d\n", lit);
    for (int cid : *active_clauses) {
	if (skip_clause(cid))
	    continue;
	int len = clause_length(cid);
	fprintf(outfile, "  %d:", cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (!skip_literal(lit))
		fprintf(outfile, " %d", get_literal(cid, lid));
	}
	fprintf(outfile, "\n");
    }
    return true;
}



// Helper routines for hashing
static int get_mapped_literal(int lit, int &var_count, Sequencer &seq,
		       std::unordered_map<int, int> &literal_map, std::vector<uint32_t> &literal_hash) {
    if (literal_map.find(lit) == literal_map.end()) {
	int nlit = ++var_count;
	literal_map[lit] = nlit;
	literal_map[-lit] = -nlit;
	literal_hash.push_back(seq.next());
	literal_hash.push_back(seq.next());
	return nlit;
    } else
	return literal_map[lit];
}

static uint32_t get_literal_hash(int lit, std::vector<uint32_t> &literal_hash) {
    int var = IABS(lit);
    int phase = lit > 0 ? 1 : 0;
    int idx = 2*(var-1) + phase;
    return literal_hash[idx];
}

#define HMOD 2147462143ULL
#define VWT 5281ULL
#define HWT 7919ULL

static uint32_t next_hash(uint32_t sofar, uint32_t val) {
    return ((uint64_t) sofar * VWT + (uint64_t) val * HWT) % HMOD;
}

std::vector<int> *Cnf::canonize(int prefix) {
    std::vector<int> *result = new std::vector<int>;
    std::unordered_map<int, int> literal_map;
    std::vector<uint32_t> literal_hash;
    Sequencer seq(123456);
    uint32_t zero_hash = seq.next();
    uint32_t formula_hash = seq.next();

    result->push_back(0);  // Placeholder for hash
    result->push_back(prefix); // Where result will be stored
    int var_count = 0;
#if 0
    for (int lit : bcp_unit_literals) {
	int nlit = get_mapped_literal(lit, var_count, seq, literal_map, literal_hash);
	result->push_back(nlit);
	formula_hash = next_hash(formula_hash, get_literal_hash(nlit, literal_hash));
	result->push_back(0);
	formula_hash = next_hash(formula_hash, zero_hash);
    }
#endif
    for (int cid : *active_clauses) {
	if (skip_clause(cid))
	    continue;
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (skip_literal(lit))
		continue;
	    int nlit = get_mapped_literal(lit, var_count, seq, literal_map, literal_hash);
	    result->push_back(nlit);
	    formula_hash = next_hash(formula_hash, get_literal_hash(nlit, literal_hash));
	}
	result->push_back(0);
	formula_hash = next_hash(formula_hash, zero_hash);
    }
    (*result)[0] = (int) formula_hash;
    if (verblevel >= 6) {
	report(6, "Canonizing formula gives:");
	for (int val : *result) {
	    printf(" %d", val);
	}
	printf("\n");
    }
    return result;
}

// Helper routines for formula cache lookup
static bool same_formulas(std::vector<int> *f1, std::vector<int> *f2) {
    if (f1->size() != f2->size())
	return false;
    // Compare hashes
    if ((*f1)[0] != (*f2)[0])
	return false;
    for (int idx = 2; idx < f1->size(); idx++) {
	if ((*f1)[idx] != (*f2)[idx]) {
	    return false;
	}
    }
    return true;
}

static bool get_sat(std::vector<int> *f) {
    if (f == NULL)
	return true;
    return (bool) (*f)[1];
}

static void set_sat(std::vector<int> *f, bool sat) {
    if (f == NULL)
	return;
    (*f)[1] = (int) sat;
}
    
// Look up formula.  If found, return stored version and free argument
std::vector<int> *Cnf::find_formula(std::vector<int> *f) {
    incr_count(COUNT_CACHE_LOOKUP);
    uint32_t h = (uint32_t) (*f)[0];
    auto bucket = formula_cache.equal_range(h);
    for (auto iter = bucket.first; iter != bucket.second; iter++) {
	std::vector<int> *g = iter->second;
	if (same_formulas(f, g)) {
	    incr_count(COUNT_CACHE_HIT);
	    return g;
	}
    }
    // Miss
    return NULL;
}

// Look up formula.  If found, return stored version and free argument
void Cnf::store_formula(std::vector<int> *f) {
    if (f == NULL)
	return;
    uint32_t h = (uint32_t) (*f)[0];
    formula_cache.insert({h,f});
}


bool Cnf::write(FILE *outfile, bool show_data_variables, bool show_tseitin_variables, bool show_weights) {
    int nvar = variable_count();
    // See if all projection variables have been eliminated
    bool degenerate = get_variable_type_count(VAR_NONTSEITIN) 
	+ get_variable_type_count(VAR_TSEITIN_DETECT) 
	+ get_variable_type_count(VAR_TSEITIN_PROMOTE) == 0;
    // Figure out which unit literals are data variables
    std::vector<int> data_literals;
    int removed_literals = 0;
    for (int lit : bcp_unit_literals) {
	int var = IABS(lit);
	if (is_data_variable(var))
	    data_literals.push_back(lit);
	else
	    removed_literals++;
    }
    if (show_data_variables && show_weights && input_weights->size() > 0)
	fprintf(outfile, degenerate ? "c t wmc\n" : "c t pwmc\n");
    else if (show_data_variables)
	fprintf(outfile, degenerate ? "c t mc\n" : "c t pmc\n");
    fprintf(outfile, "p cnf %d %d\n", nvar, current_clause_count()-removed_literals);
    if (!degenerate && show_data_variables && data_variables) {
	fprintf(outfile, "c p show");
	for (int v : *data_variables) 
	    fprintf(outfile, " %d", v);
	fprintf(outfile, " 0\n");
    }
    if (!degenerate && show_tseitin_variables && tseitin_variables && tseitin_variables->size() > 0) {
	fprintf(outfile, "c p forget");
	for (int v : *tseitin_variables) 
	    fprintf(outfile, " %d", v);
	fprintf(outfile, " 0\n");
    }
    if (show_weights && input_weights) {
	for (auto iter : *input_weights) {
	    int lit = iter.first;
	    q25_ptr weight = iter.second;
	    fprintf(outfile, "c p weight %d ", lit);
	    q25_write(weight, outfile);
	    fprintf(outfile, " 0\n");
	}
    }
    for (int lit : data_literals)
	fprintf(outfile, "%d 0\n", lit);
    for (int cid : *active_clauses) {
	if (skip_clause(cid))
	    // Hopefully won't need to do this
	    fprintf(outfile, "1 -1 0\n");
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (!skip_literal(lit))
		fprintf(outfile, "%d ", lit); 
	}
	fprintf(outfile, "0\n");
    }
    return true;
}

bool Cnf::is_satisfiable(bool bcp_only) {
    if (verblevel >= 5) {
	printf("Calling is_satisfiable for clauses:\n");
	show(stdout);
    }
    // Try for cache pre BCP
    std::vector<int> *f_prebcp = NULL;
    std::vector<int> *g = NULL;
#if 0    
    std::vector<int> *f_prebcp = canonize((int) true);
    std::vector<int> *g = find_formula(f_prebcp);
#endif
    if (g != NULL) {
	bool sat = get_sat(g);
 	report(5, "Cache lookup pre BCP yielded result %s\n", sat ? "SAT" : "UNSAT");
	if (!sat)
	    incr_count(COUNT_TSEITIN_CACHE_VAR);
	delete f_prebcp;
	return sat;
    }
    bcp(false);

    if (has_conflict) {
	set_sat(f_prebcp, false);
	store_formula(f_prebcp);
	incr_count(COUNT_TSEITIN_BCP_VAR);
	return false;
    }

    std::vector<int> *f_postbcp = canonize((int) true);
    // Try for cache post BCP
    g = find_formula(f_postbcp);
    if (g != NULL) {
	bool sat = get_sat(g);
	set_sat(f_prebcp, sat);
	store_formula(f_prebcp);
	delete f_postbcp;
	report(5, "Cache lookup post BCP yielded result %s\n", sat ? "SAT" : "UNSAT");
	if (!sat)
	    incr_count(COUNT_TSEITIN_CACHE_VAR);
	return sat;
    }

    if (bcp_only) {
	report(5, "BCP failed to find conflict.  Assuming to be SAT\n");
	return true;
    }

    double start = tod();
    int clause_count = 0;
    Glucose::Solver solver;
    solver.verbosity = 0;
    auto plit = new Glucose::Lit[nvar];
    auto nlit = new Glucose::Lit[nvar];
    for (int v = 1; v <= nvar; v++) {
	Glucose::Var gvar = solver.newVar(true, true);
	plit[v-1] = Glucose::mkLit(gvar, true);
	nlit[v-1] = Glucose::mkLit(gvar, false);
    }
    Glucose::vec<Glucose::Lit> gclause;
    for (int lit : bcp_unit_literals) {
	int var = IABS(lit);
	gclause.clear();
	Glucose::Lit glit = lit > 0 ? plit[var-1] : nlit[var-1];
	gclause.push(glit);
	solver.addClause(gclause);
	clause_count++;
    }
    for (int cid : *active_clauses) {
	if (skip_clause(cid))
	    continue;
	gclause.clear();
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (skip_literal(lit))
		continue;
	    int var = IABS(lit);
	    Glucose::Lit glit = lit > 0 ? plit[var-1] : nlit[var-1];
	    gclause.push(glit);
	}
	solver.addClause(gclause);
	clause_count++;
    }
    bool result = solver.solve();
    double elapsed = tod() - start;
    incr_timer(TIME_SAT, elapsed);
    incr_count(COUNT_SAT_CALL);
    incr_histo(HISTO_SAT_CLAUSES, clause_count);
    sat_elapsed += elapsed;
    delete[] plit;
    delete[] nlit;
    report(5, "Calling SAT solver on problem with %d variables and %d clauses yields %s\n",
	   nvar, clause_count, result ? "SAT" : "UNSAT");
    set_sat(f_prebcp, result);
    store_formula(f_prebcp);
    set_sat(f_postbcp, result);
    store_formula(f_postbcp);
    if (!result)
	incr_count(COUNT_TSEITIN_SAT_VAR);
    return result;
}

void Cnf::new_context() {
    action_stack.push_back({ACTION_START_CONTEXT, 0});    
}

void Cnf::pop_context() {
    while (true) {
	action_record ar = action_stack.back();
	active_record avr;
	action_stack.pop_back();
	switch (ar.action) {
	case ACTION_START_CONTEXT:      // Start of new context
	    return;
	case ACTION_CONFLICT:           // Found conflict
	    has_conflict = false;
	    break;
	case ACTION_DEACTIVATE_CLAUSE:  // Deactivated clause
	    activate_clause(ar.ele);
	    break;
	case ACTION_BCP:                // Derived unit literal by BCP
	    bcp_unit_literals.erase(ar.ele);
	    unit_literals.erase(ar.ele);
	    break;
	case ACTION_ASSERT:             // Asserted literal externally
	    unit_literals.erase(ar.ele);
	    break;
	case ACTION_ASSERT_FROM_BCP:    // Converted BCP unit literal into asserted literal
	    bcp_unit_literals.insert(ar.ele);
	    break;
	case ACTION_UQUANTIFY:           // Variable was universally quantfied
	    uquantified_variables.erase(ar.ele);
	    break;
	case ACTION_ACTIVE_CLAUSES:     // Changed set of active clauses
	    delete active_clauses;
	    delete literal_clauses;
	    avr = active_stack.back();
	    active_stack.pop_back();
	    active_clauses = avr.active_clauses;
	    literal_clauses = avr.literal_clauses;
	    break;
	default:
	    err(true, "Unknown action on action stack.  Value = %d\n", ar.action);
	}
    }
}

void Cnf::assign_literal(int lit, bool bcp) {
    int var = IABS(lit);
    if (var == 0 || var > nvar) {
	err(true, "Can't assign literal %d\n", lit);
    }
    bool was_unit = unit_literals.find(lit) != unit_literals.end();
    bool was_bcp_unit = bcp_unit_literals.find(lit) != bcp_unit_literals.end();
    
    if (unit_literals.find(-lit) != unit_literals.end()) {
	// Conflict
	trigger_conflict();
	return;
    }
    if (bcp) {
	if (was_unit) {
	    err(false, "Attempt to set literal %d by BCP that is already unit\n", lit);
	} else {
	    unit_literals.insert(lit);
	    bcp_unit_literals.insert(lit);
	    action_stack.push_back({ACTION_BCP, lit});
	}
    } else {
	if (was_unit && !was_bcp_unit) {
	    err(false, "Attempt to assert literal %d that is already unit\n", lit);
	} else {
	    if (was_bcp_unit) {
		bcp_unit_literals.erase(lit);
		action_stack.push_back({ACTION_ASSERT_FROM_BCP, lit});
	    } else {
		unit_literals.insert(lit);
		action_stack.push_back({ACTION_ASSERT, lit});
	    }
	}
    }
}

void Cnf::uquantify_variable(int var) {
    uquantified_variables.insert(var);
    action_stack.push_back({ACTION_UQUANTIFY, var});
}

void Cnf::activate_clause(int cid) {
    int len = clause_length(cid);
    for (int lid = 0; lid < len; lid++) {
	int lit = get_literal(cid, lid);
	(*literal_clauses)[lit].insert(cid);
    }
    active_clauses->insert(cid);
}

void Cnf::push_active(std::set<int> *nactive_clauses) {
    active_stack.push_back({active_clauses,literal_clauses});
    action_stack.push_back({ACTION_ACTIVE_CLAUSES,0});
    active_clauses = nactive_clauses;
    literal_clauses = new std::unordered_map<int,std::unordered_set<int>>;
    for (int cid : *active_clauses) {
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (!skip_literal(lit))
		(*literal_clauses)[lit].insert(cid);
	}
    }
}


// Mark clause for deactivation once iterator completes
// Clause is no longer considered part of clausal state
void Cnf::deactivate_clause(int cid) {
    int len = clause_length(cid);
    for (int lid = 0; lid < len; lid++) {
	int lit = get_literal(cid, lid);
	(*literal_clauses)[lit].erase(cid);
    }
    active_clauses->erase(cid);
    action_stack.push_back({ACTION_DEACTIVATE_CLAUSE, cid});
}

void Cnf::deactivate_clauses(std::vector<int> &remove) {
    for (int cid : remove)
	deactivate_clause(cid);
}

bool Cnf::skip_clause(int cid) {
    int len = clause_length(cid);
    for (int lid = 0; lid < len; lid++) {
	int lit = get_literal(cid, lid);
	if (unit_literals.find(lit) != unit_literals.end())
	    return true;
    }
    return false;
}

bool Cnf::skip_literal(int lit) {
    if  (unit_literals.find(-lit) != unit_literals.end())
	return true;
    int var = IABS(lit);
    if (uquantified_variables.find(var) != uquantified_variables.end())
	return true;
    return false;
}

void Cnf::trigger_conflict() {
    has_conflict = true;
    action_stack.push_back({ACTION_CONFLICT, 0});
}

// Return TAUTOLOGY, CONFLICT, propagated unit, or zero
int Cnf::propagate_clause(int cid) {
    int len = clause_length(cid);
    int result = CONFLICT;
    for (int lid = 0; lid < len; lid++) {
	int lit = get_literal(cid, lid);
	if (unit_literals.find(lit) != unit_literals.end()) {
	    result = TAUTOLOGY;
	    break;
	}
	if (skip_literal(lit))
	    continue;
	if (result == CONFLICT)
	    result = lit;
	else
	    result = 0;
    }
    return result;
}


int Cnf::bcp(bool preprocess) {
    double start = tod();
    unique_queue<int> clause_queue(active_clauses);
    int count = 0;
    while (!has_conflict && !clause_queue.empty()) {
	int cid = clause_queue.get_and_pop();
	if (active_clauses->find(cid) == active_clauses->end())
	    continue;
	int rval = propagate_clause(cid);
	if (rval == CONFLICT)
	    trigger_conflict();
	else if (rval == 0)
	    continue;
	else if (rval == TAUTOLOGY)
	    deactivate_clause(cid);
	else {
	    int lit = rval;
	    int var = IABS(lit);
	    if (preprocess)
		set_variable_type(var, VAR_ELIM);
	    std::vector<int> remove;
	    assign_literal(lit, true);
	    deactivate_clause(cid);
	    for (int ocid : (*literal_clauses)[lit]) {
		if (active_clauses->find(ocid) != active_clauses->end())
		    remove.push_back(ocid);
	    }
	    deactivate_clauses(remove);
	    if (literal_clauses->find(-lit) != literal_clauses->end()) {
		for (int ocid : (*literal_clauses)[-lit])
		    if (active_clauses->find(ocid) != active_clauses->end())
			clause_queue.push(ocid);
	    }
	    count++;
	}
    }
    incr_timer(TIME_BCP, tod() - start);
    return count;
}

void Cnf::set_variable_type(int var, var_t type) {
    if (var <= 0 || var > nvar)
	err(true, "Attempted to set type of variable %d to %d\n", var, (int) type);
    variable_type[var-1] = type;
}

var_t Cnf::get_variable_type(int var) {
    if (var <= 0 || var > nvar)
	err(true, "Attempted to get type of variable %d\n", var);
    return variable_type[var-1];
}


int Cnf::get_variable_type_count(var_t type) {
    int count = 0;
    for (int v = 1; v <= nvar; v++) {
	if (get_variable_type(v) == type)
	    count++;
    }
    return count;
}


///// Support for bounded variable elimination

int Cnf::resolve(int var, int cid1, int cid2) {
    // Merge two sets of literals
    std::vector<int> mlits;
    int len1 = clause_length(cid1);
    for (int lid1 = 0; lid1 < len1; lid1++) {
	int lit1 = get_literal(cid1, lid1);
	int var1 = IABS(lit1);
	if (var1 == var)
	    continue;
	if (skip_literal(lit1))
	    continue;
	mlits.push_back(lit1);
    }
    int len2 = clause_length(cid2);
    for (int lid2 = 0; lid2 < len2; lid2++) {
	int lit2 = get_literal(cid2, lid2);
	int var2 = IABS(lit2);
	if (var2 == var)
	    continue;
	if (skip_literal(lit2))
	    continue;
	mlits.push_back(lit2);
    }
    std::sort(mlits.begin(), mlits.end(), abs_less);
    int last_lit = 0;
    // Generate literals for resolvent
    std::vector<int> nlits;
    for (int lit : mlits) {
	if (lit == last_lit)
	    continue;
	if (lit == -last_lit) {
	    // Tautology
	    report(5, "Resolving clauses %d and %d (variable %d) yields tautology\n", cid1, cid2, var);
	    return 0;
	}
	nlits.push_back(lit);
	last_lit = lit;
    }
    int cid = new_clause();
    for (int lit : nlits)
	add_literal(lit);
    report(5, "Resolving clauses %d and %d (variable %d) yields clause %d\n", cid1, cid2, var, cid);
    return cid;
}

#if RANDOM_BVE
int Cnf::bve(bool preprocess, int maxdegree, bool preserve_literals) {
    // Literal preservation NOT implemented
    double start = tod();
    // Limit on number of added claues.  Based on number when have balanced elimination
    int maxadded = maxdegree*maxdegree - 2*maxdegree;
    // Projection variables
    std::unordered_set<int> proj_variables;
    // Variables with sufficiently low degree, queued by degree+random number
    std::map<int64_t,int> candidate_variables;
    std::unordered_map<int,int64_t> inverse_map;
    int seed = 123456;
    Sequencer seq(seed);
    int eliminated_count = 0;
    for (int cid: *active_clauses) {
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    int var = IABS(lit);
	    if (skip_literal(lit))
		continue;
	    if (is_data_variable(var))
		continue;
	    if (proj_variables.find(var) != proj_variables.end())
		continue;
	    proj_variables.insert(var);
	    int degree = IMIN((*literal_clauses)[lit].size(), (*literal_clauses)[-lit].size());;
 	    if (degree <= maxdegree) {
		int64_t key = pack(degree, seq.next());
		candidate_variables[key] = var;
		inverse_map[var] = key;
	    }
	    report(5, "Projection variable %d.  Degree = %d\n", var, degree);
	}
    }
    // Iteratively eliminate variables where one literal is of low degree
    while (candidate_variables.size() > 0) {
	auto find = candidate_variables.begin();
	int64_t key = find->first;
	int64_t var = find->second;
	candidate_variables.erase(key);
	inverse_map.erase(var);
	int dpos = (*literal_clauses)[var].size();
	int dneg = (*literal_clauses)[-var].size();
	int degree = IMIN(dpos,dneg);
	// Literal with lower degree
	int lit = dneg < dpos ? -var : var;
	int deprecated_clause_count =  dpos + dneg;
	int max_delta_clause_count = dpos * dneg - deprecated_clause_count;
	if (max_delta_clause_count > maxadded)
	    // Skip.  Might generate too many clauses
	    continue;
	// Perform BVE on var
	int new_clause_count = 0;
	eliminated_count++;
	if (preprocess)
	    set_variable_type(var, VAR_ELIM);
	std::unordered_set<int>change_variables;
	std::vector<int>deprecate_clauses;
	for (int cid1 : (*literal_clauses)[lit]) {
	    deprecate_clauses.push_back(cid1);
	    int len1 = clause_length(cid1);
	    for (int lid1 = 0; lid1 < len1; lid1++) {
		int lit1 = get_literal(cid1, lid1);
		if (skip_literal(lit1))
		    continue;
		if (lit1 == lit)
		    continue;
		int var1 = IABS(lit1);
		if (is_data_variable(var1))
		    continue;
		change_variables.insert(var1);
	    }
	}
	for (int cid2 : (*literal_clauses)[-lit]) {
	    deprecate_clauses.push_back(cid2);
	    int len2 = clause_length(cid2);
	    for (int lid2 = 0; lid2 < len2; lid2++) {
		int lit2 = get_literal(cid2, lid2);
		if (skip_literal(lit2))
		    continue;
		if (lit2 == -lit)
		    continue;
		int var2 = IABS(lit2);
		if (is_data_variable(var2))
		    continue;
		change_variables.insert(var2);
	    }
	}
	for (int cid1 : (*literal_clauses)[lit]) {
	    for (int cid2 : (*literal_clauses)[-lit]) {
		int ncid = resolve(var, cid1, cid2);
		if (ncid > 0)
		    new_clause_count++;
	    }
	}
	deactivate_clauses(deprecate_clauses);
	for (int ovar : change_variables) {
	    auto find = inverse_map.find(ovar);
	    if (find != inverse_map.end()) {
		int64_t key = find->second;
		inverse_map.erase(ovar);
		candidate_variables.erase(key);
	    }
	    int odegree = IMIN((*literal_clauses)[ovar].size(), (*literal_clauses)[-ovar].size());
	    if (odegree <= maxdegree) {
		int64_t okey = pack(odegree, seq.next());
		candidate_variables[okey] = ovar;
		inverse_map[ovar] = okey;
		report(5, "Projection variable %d.  Degree = %d\n", ovar, odegree);
	    }
	}
	if (degree == 0 && bcp_unit_literals.find(-lit) == bcp_unit_literals.end())
	    // Pure literal
	    assign_literal(-lit, true);
	report(3, "BVE on variable %d deprecated %d clauses and added %d new ones\n", var, deprecated_clause_count, new_clause_count);
	if (preprocess) {
	    incr_count_by(COUNT_BVE_ELIM_CLAUSE, deprecated_clause_count);
	    incr_count_by(COUNT_BVE_NEW_CLAUSE, new_clause_count);
	}
    }
    incr_timer(TIME_BVE, tod() - start);
    return eliminated_count;
}
#endif // RANDOM_BVE

#if !RANDOM_BVE
int Cnf::bve(bool preprocess, int maxdegree, bool preserve_literals) {
    double start = tod();
    // Limit on number of added claues.  Based on number when have balanced elimination
    int maxadded = maxdegree*maxdegree - 2*maxdegree;
    // Projection variables
    std::unordered_set<int> proj_variables;
    // Variables with sufficiently low degree
#if STACK_BVE
    std::vector<int> degree_variables[maxdegree+1];
#else
    std::unordered_set<int> degree_variables[maxdegree+1];
#endif
    // Eliminated variables
    std::unordered_set<int> eliminated_variables;
    for (int cid: *active_clauses) {
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    int var = IABS(lit);
	    if (skip_literal(lit))
		continue;
	    if (is_data_variable(var))
		continue;
	    if (proj_variables.find(var) != proj_variables.end())
		continue;
	    proj_variables.insert(var);
	    int degree = IMIN((*literal_clauses)[lit].size(), (*literal_clauses)[-lit].size());;
 	    if (degree <= maxdegree)
#if STACK_BVE
		degree_variables[degree].push_back(var);
#else
		degree_variables[degree].insert(var);
#endif
	    report(5, "Projection variable %d.  Degree = %d\n", var, degree);
	}
    }
    // Iteratively eliminate variables where one literal is of low degree
    while (true) {
	int var = 0;
	int lit = 0;  // Literal with lower degree
	int degree = 0;
	// Find variable contained in fewest number of clauses for some phase
	for (int d = 0; var == 0 && d <= maxdegree; d++) {
#if STACK_BVE
	    while (degree_variables[d].size() > 0) {
		int dvar = degree_variables[d].back();
		degree_variables[d].pop_back();
		int dpos = (*literal_clauses)[dvar].size();
		int dneg = (*literal_clauses)[-dvar].size();
		if (eliminated_variables.find(dvar) == eliminated_variables.end() 
		    && (dpos == d || dneg == d)) {
		    var = dvar;
		    lit = dpos <= dneg ? var : -var;
		    degree = d;
		    break;
		}
	    }

#else	    
	    std::vector<int> dequeue_variables;
	    for (int dvar : degree_variables[d]) {
		dequeue_variables.push_back(dvar);
		int dpos = (*literal_clauses)[dvar].size();
		int dneg = (*literal_clauses)[-dvar].size();
		if (eliminated_variables.find(dvar) == eliminated_variables.end() 
		    && (dpos == d || dneg == d)) {
		    var = dvar;
		    lit = dpos <= dneg ? var : -var;
		    degree = d;
		    break;
		}
	    }
	    // House cleaning.  Variables found or wrongly classified
	    for (int dvar : dequeue_variables)
		degree_variables[d].erase(dvar);
#endif
	}
	if (var == 0)
	    break;
	int dpos = (*literal_clauses)[var].size();
	int dneg = (*literal_clauses)[-var].size();
	int deprecated_clause_count =  dpos + dneg;
	int max_delta_clause_count = dpos * dneg - deprecated_clause_count;
	if (max_delta_clause_count > maxadded)
	    // Skip.  Might generate too many clauses
	    continue;
	if (preserve_literals && IMIN(dpos, dneg) == 1) {
	    // Check literal expansion
	    // For each phase, number of literals not including those of eliminated variable
	    int literal_count[2];
	    for (int phase = 0; phase <= 1; phase++) {
		literal_count[phase] = 0;
		int lit = phase ? var : -var;
		for (int cid : (*literal_clauses)[lit]) {
		    int len = clause_length(cid);
		    for (int lid = 0; lid < len; lid++) {
			int clit = get_literal(cid, lid);
			if (clit != lit && !skip_literal(clit))
			    literal_count[phase]++;
		    }
		}
	    }
	    report(2, "Literal expansion for variable %d.  Degrees = %d/%d, Current literals = %d.  BVE would give %d\n",
		   var, dpos, dneg, literal_count[0] + literal_count[1] + dpos + dneg, literal_count[0] * literal_count[1]);
	    if (literal_count[0] * literal_count[1] > literal_count[0] + literal_count[1] + dpos + dneg)
		continue;
	}
	// Perform BVE on var
	int new_clause_count = 0;
	eliminated_variables.insert(var);
	if (preprocess)
	    set_variable_type(var, VAR_ELIM);
	std::unordered_set<int>change_variables;
	std::vector<int>deprecate_clauses;
	for (int cid1 : (*literal_clauses)[lit]) {
	    deprecate_clauses.push_back(cid1);
	    int len1 = clause_length(cid1);
	    for (int lid1 = 0; lid1 < len1; lid1++) {
		int lit1 = get_literal(cid1, lid1);
		if (skip_literal(lit1))
		    continue;
		if (lit1 == lit)
		    continue;
		int var1 = IABS(lit1);
		if (is_data_variable(var1))
		    continue;
		change_variables.insert(var1);
	    }
	}
	for (int cid2 : (*literal_clauses)[-lit]) {
	    deprecate_clauses.push_back(cid2);
	    int len2 = clause_length(cid2);
	    for (int lid2 = 0; lid2 < len2; lid2++) {
		int lit2 = get_literal(cid2, lid2);
		if (skip_literal(lit2))
		    continue;
		if (lit2 == -lit)
		    continue;
		int var2 = IABS(lit2);
		if (is_data_variable(var2))
		    continue;
		change_variables.insert(var2);
	    }
	}
	for (int cid1 : (*literal_clauses)[lit]) {
	    for (int cid2 : (*literal_clauses)[-lit]) {
		int ncid = resolve(var, cid1, cid2);
		if (ncid > 0)
		    new_clause_count++;
	    }
	}
	deactivate_clauses(deprecate_clauses);
	for (int ovar : change_variables) {
	    int odegree = IMIN((*literal_clauses)[ovar].size(), (*literal_clauses)[-ovar].size());
	    if (odegree <= maxdegree) {
#if STACK_BVE
		degree_variables[odegree].push_back(ovar);
#else	      
		degree_variables[odegree].insert(ovar);
#endif
		report(5, "Projection variable %d.  Degree = %d\n", ovar, odegree);
	    }
	}
	if (degree == 0 && bcp_unit_literals.find(-lit) == bcp_unit_literals.end())
	    // Pure literal
	    assign_literal(-lit, true);
	report(3, "BVE on variable %d deprecated %d clauses and added %d new ones\n", var, deprecated_clause_count, new_clause_count);
	if (preprocess) {
	    incr_count_by(COUNT_BVE_ELIM_CLAUSE, deprecated_clause_count);
	    incr_count_by(COUNT_BVE_NEW_CLAUSE, new_clause_count);
	}
    }
    incr_timer(TIME_BVE, tod()-start);
    return (int) eliminated_variables.size();
}
#endif // !RANDOM_BVE




///// Support for Tseitin promotion

// Determine next index
static bool increment_indices(std::vector<int> &lengths, std::vector<int> &indices) {
    bool incremented = false;
    for (int i = 0; i < lengths.size(); i++) {
	if (indices[i] < lengths[i]-1) {
	    indices[i]++;
	    incremented = true;
	    break;
	} else
	    indices[i] = 0;
    }
    return incremented;
}

// Generate set of blocked clauses covering this literal & clauses
void Cnf::blocked_clause_expand(int lit, std::vector<int> &clause_list) {
    std::vector<int> clause_lengths;
    std::vector<int> clause_indices;
    for (int cid : clause_list) {
	// Stick all uninteresting literals at the end
	int len = clause_length(cid);
	int lid = 0;
	while (lid < len) {
	    int clit = get_literal(cid, lid);
	    if (clit == lit || skip_literal(clit))
		swap_literals(cid, lid, --len);
	    else
		lid++;
	}
	clause_lengths.push_back(len);
	clause_indices.push_back(0);
    }
    int running = true;
    int first_cid = 0;
    int last_cid = 0;
    while (running) {
	int ncid = new_clause();
	if (first_cid == 0)
	    first_cid = ncid;
	last_cid = ncid;
	add_literal(-lit);
	for (int i = 0; i < clause_list.size(); i++) {
	    int cid = clause_list[i];
	    int idx = clause_indices[i];
	    int clit = get_literal(cid, idx);
	    add_literal(-clit);
	}
	incr_count(COUNT_PROMOTE_CLAUSE);
	running = increment_indices(clause_lengths, clause_indices);
    }
    report(4, "Added blocked clauses #%d .. %d to promote variable %d\n", first_cid, last_cid, IABS(lit));
}


// Attempt to classify variable as Tseitin variable
// Fill list of fanout nodes in event of success, so that these can then be tested
bool Cnf::tseitin_variable_test(int var, int classify_flag, int sat_depth, std::set<int> &fanout_vars) {
    // Construct sets of clauses that contain only data & known Tseitin variables
    // For the variable
    std::set<int> *dt_var_clause_set = new std::set<int>;
    // For each phase
    std::vector<int> dt_lit_clause_list[2];
    // The other data and Tseitin literals that occur in these clauses
    std::unordered_set<int> dt_otherlit_set[2];
    // The other data and Tseitin variables that occur in these clauses
    std::unordered_set<int> dt_othervar_set;
    // Literals that occur in binary clause with either phase of variable
    std::unordered_set<int> dt_binarylit_set[2];
    // Clear the fanouts
    fanout_vars.clear();
    // For Xor detection.  0 = unassigned, -1 = not all same, >= 1 = uniform length
    int uniform_length = 0;
    bool result = false;

    for (int phase = 0; phase <= 1; phase ++) {
	int lit = (2*phase - 1) * var;
	for (int cid : (*literal_clauses)[lit]) {
	    if (skip_clause(cid))
		continue;
	    int len = clause_length(cid);
	    bool include = true;
	    std::vector<int> other_lits;
	    int lcount = 0;
	    int clause_other_lit = 0;
	    int clause_this_lit = 0;
	    for (int lid = 0; lid < len; lid++) {
		int clit = get_literal(cid, lid);
		if (skip_literal(clit))
		    continue;
		lcount++;
		int cvar = IABS(clit);
		if (cvar == var) {
		    clause_this_lit = clit;
		    continue;
		}
		if (data_variables->find(cvar) != data_variables->end() 
		    || tseitin_variables->find(cvar) != tseitin_variables->end()) {
		    clause_other_lit = clit;
		    other_lits.push_back(clit);
		} else {
		    include = false;
		    fanout_vars.insert(cvar);
		}
	    }
	    if (!include) 
		continue;

	    if (lcount == 1) {
		report(3, "Found unit variable %d.  Fanout size = %d\n", var, (int) fanout_vars.size());
		incr_count(COUNT_TSEITIN_UNIT_VAR);
		result = true;
		goto done;
	    }

	    if (lcount == uniform_length || uniform_length == 0)
		uniform_length = lcount;
	    else
		uniform_length = -1;

	    if (lcount == 2 && lit == clause_this_lit)
		dt_binarylit_set[phase].insert(clause_other_lit);

	    dt_var_clause_set->insert(cid);
	    dt_lit_clause_list[phase].push_back(cid);
	    for (int olit : other_lits) {
		dt_otherlit_set[phase].insert(olit);
		int ovar = IABS(olit);
		dt_othervar_set.insert(ovar);
	    }
	}
    }

    // Prepare for deep SAT analysis
    if (ALLOW_DETECT_SAT(classify_flag) && sat_depth >= 2) {
	defining_variables[var] = new std::unordered_set<int> (dt_othervar_set);
	defining_clauses[var] = new std::set<int> (*dt_var_clause_set);
    }


    if (ALLOW_DETECT_XOR(classify_flag)) {
	// See if Xor
	int ncount = dt_lit_clause_list[0].size();
	int pcount = dt_lit_clause_list[1].size();
	report(5, "Attempting XOR detection for variable %d.  ncount = %d, pcount = %d, uniform_length = %d, other var count = %d\n",
	       var, ncount, pcount, uniform_length, (int) dt_othervar_set.size());
	if (uniform_length >= 2 &&
	    ncount == pcount &&
	    ncount == 1 << (uniform_length-2) &&
	    dt_othervar_set.size() == uniform_length-1) {
	    // Assign power-of-two weights to the positive literals
	    std::unordered_map<int,int> literal_weights;
	    int wt = 1;
	    literal_weights[var] = wt;
	    literal_weights[-var] = 0;
	    wt *= 2;
	    for (int phase = 0; phase <= 1; phase++) {
		for (int olit : dt_otherlit_set[phase]) {
		    if (olit > 0) {
			if (literal_weights.find(olit) == literal_weights.end()) {
			    literal_weights[olit] = wt;
			    wt *= 2;
			}
		    } else
			literal_weights[olit] = 0;
		    
		}
		
	    }
	    // Compute weight for each clause.  Should be distinct
	    int clause_weight[ncount+pcount];
	    // Should only see single parity
	    bool found_odd_parity = false;
	    bool found_even_parity = false;
	    int idx = 0;
	    int ok = true;
	    for (int phase = 0; ok && phase <= 1; phase++) {
		int lit = 2 * phase - 1;
		for (int cid : dt_lit_clause_list[phase]) {
		    int cweight = 0;
		    bool odd_parity = false;
		    int len = clause_length(cid);
		    for (int lid = 0; lid < len; lid++) {
			int clit = get_literal(cid, lid);
			if (skip_literal(clit))
			    continue;
			if (clit > 0) 
			    odd_parity = !odd_parity;
			cweight += literal_weights[clit];
		    }
		    if (odd_parity) 
			found_odd_parity = true;
		    else
			found_even_parity = true;
		    if (found_even_parity && found_odd_parity) {
			ok = false;
			break;
		    }
		    for (int oidx = 0; ok && oidx < idx; oidx++) {
			if (cweight == clause_weight[oidx]) {
			    // Not distinct
			    ok = false;
			    break;
			}
		    }
		    clause_weight[idx++] = cweight;
		}
	    }
	    if (ok) {
		report(3, "Found Xor/Xnor structure for variable %d.  Fanout size = %d\n",
		       var, (int) fanout_vars.size());
		incr_count(COUNT_TSEITIN_XOR_VAR);
		result = true;
		goto done;
	    }
	}
	report(5, "Xor detection for variable %d failed\n", var);
    }
    if (ALLOW_DETECT_AND(classify_flag)) {
	// See if variant of And/Or operation
	report(5, "Attempting And/Or detection for variable %d.  binary_lits = %d/%d, clauses = %d/%d\n",
	       var, (int) dt_binarylit_set[0].size(), (int) dt_binarylit_set[1].size(),
	       (int) dt_lit_clause_list[0].size(), (int) dt_lit_clause_list[1].size());
	for (int phase = 0; phase <= 1; phase++) {
	    int lit = (2*phase - 1) * var;
	    report(5, "  Found %d binary clauses containing literal %d\n", (int) dt_binarylit_set[1-phase].size(), -lit);
	    if (dt_binarylit_set[1-phase].size() == 0)
		continue;
	    if (verblevel >= 5) {
		report(5, "  Found %d clauses containing literal %d:", (int) dt_lit_clause_list[phase].size(), lit);
		for (int cid : dt_lit_clause_list[phase])
		    printf(" %d", cid);
		printf("\n");
	    }
	    if (dt_lit_clause_list[phase].size() != 1) 
		continue;

	    for (int cid : dt_lit_clause_list[phase]) {
		report(5, "  Checking coverage of clause %d containing literal %d\n", cid, lit);
		bool covered = true;
		int len = clause_length(cid);
		for (int lid = 0; lid < len; lid++) {
		    int clit = get_literal(cid, lid);
		    if (skip_literal(clit))
			continue;
		    int cvar = IABS(clit);
		    if (cvar == var)
			continue;
		    if (dt_binarylit_set[1-phase].find(-clit) == dt_binarylit_set[1-phase].end()) {
			report(5, "  Literal %d in clause %d not covered by binary clause containing literal %d\n",
			       clit, cid, -lit);
			covered = false;
			break;
		    }
		}
		if (covered) {
		    report(3, "Found And/Or structure for variable %d.  Fanout size = %d\n",
			   var, (int) fanout_vars.size());
		    incr_count(COUNT_TSEITIN_ANDOR_VAR);
		    result = true;
		    goto done;
		}
	    }
	}
	report(5, "  And/Or detection for variable %d failed\n", var);
    }

    if (ALLOW_DETECT_BCP(classify_flag) || (ALLOW_DETECT_SAT(classify_flag) && sat_depth >= 1)) {

	for (int depth = 1; depth <= sat_depth; depth++) {
	    max_count(COUNT_TSEITIN_SAT_MAX_DEPTH_TRIED, depth);

	    if (depth > 1 && !ALLOW_DETECT_SAT(classify_flag))
		break;

	    if (depth == 1) 
		report(5, "Starting with %d dependent variables for variable %d at depth %d\n", 
		       (int) dt_othervar_set.size(), var, depth);
	    else {
		size_t vcount = dt_othervar_set.size();
		for (int dvar : dt_othervar_set) {
		    if (is_data_variable(dvar))
			continue;
		    if (defining_variables.find(dvar) == defining_variables.end()) {
			err(false, "Tseitin variable test for variable %d.  Can't find defining clauses for variable %d\n",
			    var, dvar);
			continue;
		    }
		    for (int ovar : *defining_variables[dvar]) {
			if (dt_othervar_set.find(ovar) == dt_othervar_set.end()) {
			    dt_othervar_set.insert(ovar);
			    if (defining_clauses.find(ovar) == defining_clauses.end())
				continue;
			    for (int cid : *defining_clauses[ovar])
				dt_var_clause_set->insert(cid);
			}
		    }
		}
		if (vcount == dt_othervar_set.size()) {
		    report(5, "Stopped Tseitin test for variable %d at depth %d.  No more layers to add\n", var, depth);
		    break;
		} else {
		    report(5, "Added %z more variables to dependency set for variable %d at depth %d\n", 
			   dt_othervar_set.size() - vcount, var, depth);
		}
	    }

	    report(5, "Attempting Tseitin detection for variable %d through quantification.  SAT depth = %d\n", var, depth);
	    new_context();
	    // Set gets consumed when context is popped.  Make fresh copy each time
	    std::set<int> *var_clause_set = new std::set<int> (*dt_var_clause_set);
	    push_active(var_clause_set);
	    uquantify_variable(var);
	    bool sat = is_satisfiable(!ALLOW_DETECT_SAT(classify_flag));
	    if (!sat)
		report(3, "Detected Tseitin variable %d through quantification.  SAT depth = %d.  Fanout size = %d\n",
		       var, depth, (int) fanout_vars.size());
	    if (verblevel >= 5) {
		report(5, "Tseitin test gives %s for variable %d with depth %d on clauses:", sat ? "failure" : "success", var, depth);
		for (int cid : *dt_var_clause_set)
		    printf(" %d", cid);
		printf("\n");
	    }
	    pop_context();
	    if (!sat) {
		max_count(COUNT_TSEITIN_SAT_MAX_DEPTH_SUCCESS, depth);
		result = true;
		goto done;
	    }
	}
    }
	

    if (ALLOW_PROMOTE(classify_flag)) {
	report(5, "Attempting Tseitin promotion for variable %d\n", var);
	promotion_try_count++;
	// See if can promote variable to Tseitin variable
	for (int phase = 0; phase <=  1; phase++) {
	    int lit = (2*phase - 1) * var;
	    // Make sure no other clauses contain this literal
	    if (dt_lit_clause_list[phase].size() < (*literal_clauses)[lit].size())
		continue;
	    // Make sure all other literals in these clauses are pure
	    bool pure = true;
	    for (int olit : dt_otherlit_set[phase]) {
		if (olit < 0)
		    // Only look to one side of negagtion
		    continue;
		if (dt_otherlit_set[phase].find(-olit) != dt_otherlit_set[phase].end()) {
		    pure = false;
		    break;
		}
	    }
	    // Go for it!
	    promotion_success_count++;
	    if (pure) {
		blocked_clause_expand(lit, dt_lit_clause_list[phase]);
		set_variable_type(var, VAR_TSEITIN_PROMOTE);
		report(3, "Promoted variable %d.  Fanout size = %d\n", var, (int) fanout_vars.size());
		result = true;
		goto done;
	    }
	}
    }
    fanout_vars.clear();
    result = false;
 done:
    delete dt_var_clause_set;
    return result;
}

void Cnf::classify_variables(int classify_flag, int sat_depth, double overall_total, double sat_total) {
    double start = tod();
    tseitin_variables->clear();
    // Set and queue of (potentially Tseitin) projection variables
    // Original list by definition occurrence, but then add fanouts of newly discovered/created variables
    unique_queue<int> pvar_queue;
    // Fanouts of identified or promoted Tseitin variables.  Use ordered set to guarantee reproducibility
    std::set<int> fanout_vars;
    // Variables that didn't get detected/promoted
    std::unordered_set<int> non_tseitin_vars;

    // Build mappings and ordered list of projection variables
    for (int cid : *active_clauses) {
	if (skip_clause(cid))
	    continue;
	int len = clause_length(cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    if (skip_literal(lit))
		continue;
	    int var = IABS(lit);
	    if (data_variables->find(var) != data_variables->end())
		continue;
	    if (pvar_queue.push(var))
		non_tseitin_vars.insert(var);
	}
    }
    while (!pvar_queue.empty()) {
	if (tod()-start >= overall_total) {
	    report(1, "c Exceeded overall time limit.  No further classification attempted\n");
	    break;
	}
	int var = pvar_queue.get_and_pop();
	if (ALLOW_DETECT_SAT(classify_flag) && sat_elapsed >= sat_total) {
	    classify_flag &= ~CLASSIFY_FLAG_DETECT_SAT;
	    err(false, "Exceeded SAT solver time limit %.1f.  Continuing other modes\n", sat_total);
	}
	if (tseitin_variable_test(var, classify_flag, sat_depth, fanout_vars)) {
	    if (get_variable_type(var) != VAR_TSEITIN_PROMOTE)
 		set_variable_type(var, VAR_TSEITIN_DETECT);
	    tseitin_variables->insert(var);
	    non_tseitin_vars.erase(var);
	}
	for (int fvar : fanout_vars) {
	    if (pvar_queue.push(fvar)) 
		report(3, "Added fanout variable %d for Tseitin variable %d\n", fvar, var);
	}
	incr_count(COUNT_TSEITIN_TEST);
    }
    report(3, "c Failed to detect/promote %d variables\n", (int) non_tseitin_vars.size());
    if (verblevel >= 5) {
	printf("c Non-Tseitin vars:");
	for (int ntvar : non_tseitin_vars)
	    printf(" %d", ntvar);
	printf("\n");
    }
    incr_timer(TIME_CLASSIFY, tod()-start);
}

