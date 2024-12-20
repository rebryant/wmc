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


#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/param.h>

#include "report.h"
//#include "path.h"

int verblevel = 1;

FILE *errfile = NULL;
FILE *verbfile = NULL;

static char* prog_path = NULL;

static panic_function_t panic_function = NULL;

static const char *logfile_name = NULL;

static const char *datafile_name = "datafile.csv";

static double start_time = 0.0;

const char *archive_string(const char *tstring) {
    char *rstring = (char *) malloc(strlen(tstring)+1);
    strcpy(rstring, tstring);
    return (const char *) rstring;
}

#if 0
const char *find_program_path(const char *progname) {
    char buf[MAXPATHLEN];
    char patbuf[MAXPATHLEN];
    bool found = false;

    // Follow PATH environment variable

    if (!found) {
	char *origpath = NULL;
	if ((origpath = getenv("PATH")) == NULL)
	    return NULL;
	strncpy(patbuf, origpath, MAXPATHLEN);
    }
    char *path = patbuf;
    while (!found) {
	char *cp = index(path, ':');
	if (cp) {
	    *cp = '\0';
	}
	snprintf(buf, MAXPATHLEN, "%s/%s", path, progname);
	found = access(buf, X_OK) == 0;
	if (!found) {
	    if (cp)
		path = ++cp;
	    else
		break;
	}
    }

    // Look in directory of program executable

    if (!found) {
	snprintf(buf, MAXPATHLEN, "%s/%s", IPATH, progname);
	found = access(buf, X_OK) == 0;
    }

    return found ? archive_string(buf) : NULL;
}
#endif

//  Logging information
// Establish a log file
void set_logname(const char *fname) {
    if (fname == NULL) {
	logfile_name = NULL;
	return;
    }
    logfile_name = archive_string(fname);
    // Clear out whatever was there
    FILE *logfile = fopen(logfile_name, "w");
    if (logfile)
	fclose(logfile);
}


void set_verblevel(int level) {
    verblevel = level;
}

void set_panic(panic_function_t fun) {
    panic_function = fun;
}

void err(bool fatal, const char *fmt, ...) {
    if (!errfile)
	errfile = stdout;
    va_list ap;
    va_start(ap, fmt);
    if (fatal)
	fprintf(errfile, "c ERROR: ");
    else
	fprintf(errfile, "c WARNING: ");
    vfprintf(errfile, fmt, ap);
    fflush(errfile);
    va_end(ap);
    if (logfile_name) {
	FILE *logfile = fopen(logfile_name, "a");
	if (logfile) {
	    va_start(ap, fmt);
	    if (fatal)
		fprintf(logfile, "c ERROR: ");
	    else
		fprintf(logfile, "c WARNING: ");
	    vfprintf(logfile, fmt, ap);
	    va_end(ap);
	    fclose(logfile);
	}
    }
    if (fatal) {
	if (panic_function)
	    panic_function();
	exit(1);
    }
}

void report(int level, const char *fmt, ...) {
    if (!verbfile)
	verbfile = stdout;
    va_list ap;
    if (level <= verblevel) {
	fprintf(verbfile, "c ");
	va_start(ap, fmt);
	vfprintf(verbfile, fmt, ap);
	fflush(verbfile);
	va_end(ap);
	if (logfile_name) {
	    FILE *logfile = fopen(logfile_name, "a");
	    if (logfile) {
		fprintf(logfile, "c ");
		va_start(ap, fmt);
		vfprintf(logfile, fmt, ap);
		va_end(ap);
		fclose(logfile);
	    }
	}
    }
}

void lprintf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stdout, fmt, ap);
    fflush(stdout);
    va_end(ap);
    if (logfile_name) {
	FILE *logfile = fopen(logfile_name, "a");
	if (logfile) {
	    va_start(ap, fmt);
	    vfprintf(logfile, fmt, ap);
	    va_end(ap);
	    fclose(logfile);
	}
    }
}

void log_data(const char *fmt, ...) {
    va_list ap;
    if (datafile_name == NULL)
	return;
    FILE *datafile = fopen(datafile_name, "a");
    if (!datafile)
	return;
    va_start(ap, fmt);
    vfprintf(datafile, fmt, ap);
    va_end(ap);
    fclose(datafile);
}


double tod() {
    struct timeval tv;
    if (gettimeofday(&tv, NULL) == 0)
	return (double) tv.tv_sec + 1e-6 * tv.tv_usec;
    else
	return 0.0;
}

void start_timer() {
    start_time = tod();
}

double get_elapsed() {
    return tod() - start_time;
}

const char *b2a(bool b) {
    return b ? "True" : "False";
}
