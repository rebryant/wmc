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
#include "cnf_info.hh"


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

#define BSIZE 1024

// Process comment, looking additional data variables & weights
// Return last character
static void process_comment(FILE *infile, std::unordered_set<int> *data_variables,
			    std::unordered_set<int> *forget_variables,
			    std::unordered_map<int,const char*> *input_weights) {
    char buf[BSIZE];
    char wbuf[BSIZE];
    int len;
    if (find_string_token(infile, buf, BSIZE, &len) && len == 1 && strncmp(buf, "p", 1) == 0
	&& find_string_token(infile, buf, BSIZE, &len)) {
	bool show = true;
	if (len == 4 && ((show = (strncmp(buf, "show", 4) == 0))
			 || strncmp(buf, "forget", 6) == 0)) {
	    int var = -1;
	    std::unordered_set<int> *buf = show ? data_variables : forget_variables;
	    while (var != 0) {
		if (fscanf(infile, "%d", &var) != 1) {
		    err(false, "Couldn't read %s variable\n", show ? "data" : "Forget");
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
	    if (find_string_token(infile, wbuf, BSIZE, &len))
		(*input_weights)[lit] = archive_string(wbuf);
	    else {
		err(false, "Couldn't read weight for literal %d (skipping)\n", lit);
		skip_line(infile);
		return;
	    }
	    int zero;
	    if (fscanf(infile, "%d", &zero) != 1 || zero != 0) {
		err(false, "Couldn't read terminating zero in weight declaration for literal %d (accepting weight)\n", lit);
	    }
	}
    }
    skip_line(infile);
}		

Cnf::Cnf() {
    data_variables = NULL;
    forget_variables = NULL;
    input_weights = NULL;
    initialize(0);
}

Cnf::Cnf(int input_count) {
    data_variables = NULL;
    forget_variables = NULL;
    input_weights = NULL;
    initialize(input_count);
}

Cnf::~Cnf() {
    deallocate();
}

void Cnf::initialize(int input_count) {
    nvar = input_count;
    clause_offset.clear();
    literal_sequence.clear();
    if (!data_variables)
	data_variables = new std::unordered_set<int>;
    if (!forget_variables)
	forget_variables = new std::unordered_set<int>;
    if (!input_weights)
	input_weights = new std::unordered_map<int, const char *>;
    new_clause();
}

// Must explicitly deallocate sets
void Cnf::deallocate() {
    delete forget_variables;
}

int Cnf::new_clause() {
    int cid = clause_offset.size();
    clause_offset.push_back(literal_sequence.size());
    return cid;
}

void Cnf::add_literal(int lit) {
    literal_sequence.push_back(lit);
    clause_offset.back() ++;
    int cid = clause_offset.size()-1;
    int var = IABS(lit);
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
		process_comment(infile, data_variables, forget_variables, input_weights);
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
		    process_comment(infile, data_variables, forget_variables, input_weights);
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
		process_comment(infile, data_variables, forget_variables, input_weights);
	    else
		skip_line(infile);

	}
    }
    // If no data variables declared, assume all input variables are data variables
    if (data_variables->size() == 0) {
	for (int v = 1; v <= variable_count(); v++)
	    data_variables->insert(v);
    }
    incr_count_by(COUNT_INPUT_CLAUSE, maximum_clause_id());
    incr_count_by(COUNT_DATA_VARIABLES, data_variables->size());
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

bool Cnf::show(FILE *outfile) {
    for (int cid = 1; cid <= clause_count(); cid++) {
	int len = clause_length(cid);
	fprintf(outfile, "  %d:", cid);
	for (int lid = 0; lid < len; lid++) {
	    int lit = get_literal(cid, lid);
	    fprintf(outfile, " %d", get_literal(cid, lid));
	}
	fprintf(outfile, "\n");
    }
    return true;
}
