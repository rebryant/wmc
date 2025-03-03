#!/usr/bin/python3

#####################################################################################
# Copyright (c) 2023 Randal E. Bryant, Carnegie Mellon University
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
########################################################################################

import csv
import sys


def usage(name):
    sys.stderr.write("Usage: %s [-v] [-w] [CSVPATH]\n" % name)
    sys.stderr.write("  -v      Verbose\n")
    sys.stderr.write("  -w      Print wide format\n")
    sys.stderr.write("  CSVPATH Directory with CSV files\n")

directory = "."
method_count = 5

method_mpf, method_mpfi, method_mpfi2, method_mpq, method_mpfi2_direct = range(method_count)
method_name = ["MPF", "MPFI", "MPFI2", "MPQ", "MPFI2-D"]
file_name = "tabulate-all.csv"

def err(fatal, msg):
    if fatal:
        sys.stderr.write("ERROR: %s\n" % msg)
        sys.exit(1)
    else:
        sys.stderr.write("WARNING: %s\n" % msg)

# Hack to account for times taken by failed runs
fail_time_dict = {"mc2024_track2-random_023" : 1592.0, "mc2024_track2-random_178" : 1634.2 }

def fail_time(instance):
    formula = instance[:-4]
    if formula not in fail_time_dict:
        sys.stderr.write("Couldn't find entry for instance %s (formula %s)\n" % (instance, formula))
        return 0.0
    return fail_time_dict[formula]

class Instance:
    name = None
    method = None
    mpf_seconds = None
    mpfi_seconds = None
    mpfi2_seconds = None
    mpq_seconds = None
    mpq_ok = False
    
    def __init__(self, dict):
        self.name = dict["bench"]
        smethod = dict["mode"]
        try:
            self.method = int(float(smethod)) - 1
        except:
            err(True, "Instance %s.  Couldn't extract method from '%s'" % (self.name, smethod))
        if self.method < 0 or self.method > method_mpq:
            err(True, "Instance %s.  Invalid method number %d" % (self.name, self.method))
        smpf = dict["mpf"]
        self.mpf_seconds = None
        try:
            self.mpf_seconds = float(smpf)
        except:
            pass
        smpfi = dict["mpfi"]
        self.mpfi_seconds = None
        try:
            self.mpfi_seconds = float(smpfi)
        except:
            pass
        smpfi2 = dict["mpfi2"]
        self.mpfi2_seconds = None
        try:
            self.mpfi2_seconds = float(smpfi2)
        except:
            pass
        smpq = dict["mpq"]
        self.mpq_seconds = None
        try:
            self.mpq_seconds = float(smpq)
            self.mpq_ok = True
        except:
            self.mpq_seconds = fail_time(self.name)
            self.mpq_ok = False
        if self.method == method_mpf and self.mpf_seconds is None: 
            err(True, "Instance %s.  Need MPF seconds" % str(self))
        elif self.method == method_mpfi and self.mpfi_seconds is None:
            err(True, "Instance %s.  Need MPFI seconds" % str(self))
        elif self.method == method_mpfi2 and (self.mpfi_seconds is None or self.mpfi2_seconds is None):
            err(True, "Instance %s.  Need MPFI and MPFI2 seconds" % str(self))
        elif self.method == method_mpq and (self.mpfi_seconds is None or self.mpfi2_seconds is None or self.mpq_seconds is None):
            err(True, "Instance %s.  Need MPFI, MPFI2, and MPQ seconds" % str(self))

    def __str__(self):
        return "I.%s (%s) MPF:%s MPFI:%s MPFI2:%s MPQ(%s):%s" % (self.name, method_name[self.method], str(self.mpf_seconds), str(self.mpfi_seconds), str(self.mpfi2_seconds), "S" if self.mpq_ok else "F", str(self.mpq_seconds))

instances = []
average_time = [0] * method_count

def load(verbose):
    global instances
    fname = directory + "/" + file_name
    try:
        infile = open(fname, "r")
    except:
        err(True, "Couldn't open file '%s'" % fname)
    creader = csv.DictReader(infile)
    for dict in creader:
        instance = Instance(dict)
        if verbose:
            sys.stderr.write("Created instance %s\n" % str(instance))
        instances.append(instance)

# Here, method refers to final strategy:
# MPQ: Use MPQ
# MPF: Use MPF + MPQ
# MPFI: Use MPFI + MPQ
# MPFI: Use MPF + MPFI + MPFI2 + MPQ

def tabulate(target_method, label):
    global average_time
    ok_count = [0] * method_count
    fail_count = [0] * method_count
    time = [0] * method_count
    total_ok_count = 0
    total_fail_count = 0
    total_time = 0
    for instance in instances:
        if target_method == method_mpq:
            time[method_mpq] += instance.mpq_seconds
            total_time += instance.mpq_seconds
            if instance.mpq_ok:
                ok_count[method_mpq] += 1
                total_ok_count += 1
            else:
                fail_count[method_mpq] += 1
                total_fail_count += 1
        elif target_method == method_mpf:
            if instance.method == method_mpf:
                ok_count[method_mpf] += 1
                time[method_mpf] += instance.mpf_seconds
                total_ok_count += 1
                total_time += instance.mpf_seconds
            else:
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
                if instance.mpq_ok:
                    ok_count[method_mpq] += 1
                    total_ok_count += 1
                else:
                    fail_count[method_mpq] += 1
                    total_fail_count += 1
        elif target_method == method_mpfi:
            if instance.method == method_mpf:
                ok_count[method_mpf] += 1
                time[method_mpf] += instance.mpf_seconds
                total_ok_count += 1
                total_time += instance.mpf_seconds
            elif instance.method == method_mpfi:
                ok_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_ok_count += 1
                total_time += instance.mpfi_seconds
            elif instance.method == method_mpfi2:
                fail_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
                if instance.mpq_seconds is not None:
                    ok_count[method_mpq] += 1
                    total_ok_count += 1
                else:
                    fail_count[method_mpq] += 1
                    total_fail_count += 1
            else:
                fail_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
                if instance.mpq_seconds is not None:
                    ok_count[method_mpq] += 1
                    total_ok_count += 1
                else:
                    fail_count[method_mpq] += 1
                    total_fail_count += 1
        elif target_method == method_mpfi2:
            if instance.method == method_mpf:
                ok_count[method_mpf] += 1
                time[method_mpf] += instance.mpf_seconds
                total_ok_count += 1
                total_time += instance.mpf_seconds
            elif instance.method == method_mpfi:
                ok_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_ok_count += 1
                total_time += instance.mpfi_seconds
            elif instance.method == method_mpfi2:
                fail_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
                ok_count[method_mpfi2] += 1
                time[method_mpfi2] += instance.mpfi2_seconds
                total_ok_count += 1
                total_time += instance.mpfi2_seconds
            else:
                fail_count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
                fail_count[method_mpfi2] += 1
                time[method_mpfi2] += instance.mpfi2_seconds
                total_time += instance.mpfi2_seconds
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
                if instance.mpq_seconds is not None:
                    ok_count[method_mpq] += 1
                    total_ok_count += 1
                else:
                    fail_count[method_mpq] += 1
                    total_fail_count += 1
        elif target_method == method_mpfi2_direct:
            if instance.method == method_mpf:
                ok_count[method_mpf] += 1
                time[method_mpf] += instance.mpf_seconds
                total_ok_count += 1
                total_time += instance.mpf_seconds
            elif instance.method in [method_mpfi, method_mpfi2]:
                if instance.mpfi2_seconds is None:
                    err(True, "Oops.  Can't process instance %s with MPFI2-direct" % (str(instance)))
                ok_count[method_mpfi2] += 1
                time[method_mpfi2] += instance.mpfi2_seconds
                total_ok_count += 1
                total_time += instance.mpfi2_seconds
            else:
                fail_count[method_mpfi2] += 1
                time[method_mpfi2] += instance.mpfi2_seconds
                total_time += instance.mpfi2_seconds
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
                if instance.mpq_seconds is not None:
                    ok_count[method_mpq] += 1
                    total_ok_count += 1
                else:
                    fail_count[method_mpq] += 1
                    total_fail_count += 1
    fields_top = [label, "Runs"]
    fields_bottom = ["", "Hours"]
    fields_both = [label]
    for comp_method in range(method_count):
        if comp_method == method_mpfi2_direct:
            continue
        mtime = time[comp_method]/3600.0
        rstring = "" if mtime == 0.0 else str(ok_count[comp_method]) + "+" + str(fail_count[comp_method])
        tstring = "" if mtime == 0.0 else "%.2f" % mtime
        fields_top += [rstring]
        fields_bottom +=  [tstring]
        fields_both += [rstring, tstring]
    total_count = total_ok_count + total_fail_count
    avg = total_time/total_count
    average_time[target_method] = avg
    ttime = total_time/3600.0
    rstring = str(total_ok_count) + "+" + str(total_fail_count)
    fields_top += [rstring]
    tstring = "%.2f" % ttime
    fields_bottom += [tstring]
    fields_both += [rstring, tstring]
#    fields_top += ["%.1f" % avg]
#    fields_bottom += [""]
    return (fields_top, fields_bottom, fields_both)

def table_entry(fields):
    print(" & ".join(fields) + ' \\\\')

def table_line():
    print("\\midrule")

def command(name, value):
    print("\\newcommand[\\global]{\\%s}{%s}" % (name, "%.2f" % value))


def run(name, args):
    global directory
    wide_format = False
    verbose = False
    i = argi = 0
    for i in range(len(args)):
        if args[argi][0] == '-':
            if len(args[argi]) != 2:
                sys.stderr.write("Invalid flag '%s'\n" % args[argi])
                usage(name)
                return
            if args[argi][1] == 'h':
                usage(name)
                return
            elif args[argi][1] == 'w':
                wide_format = True
            elif args[argi][1] == 'w':
                verbose = True
            else:
                sys.stderr.write("Invalid flag '%s'\n" % args[argi])
                usage(name)
                return
        else:
            directory = args[argi]
        argi += 1
    load(verbose)
    ft, fb, f = tabulate(method_mpq, "MPQ only")
    if wide_format:
        table_entry(f)
    else:
        table_entry(ft)
        table_entry(fb)
        table_line()
    ft, fb, f = tabulate(method_mpf, "MPF|MPQ")
    if wide_format:
        table_entry(f)
    else:
        table_entry(ft)
        table_entry(fb)
        table_line()
    ft, fb, f = tabulate(method_mpfi, "MPF|MPFI-128+MPQ")
    if wide_format:
        table_entry(f)
    else:
        table_entry(ft)
        table_entry(fb)
        table_line()
    ft, fb, f = tabulate(method_mpfi2, "\\textcolor{red}{MPF|MPFI$\\times$2+MPQ}")
    if wide_format:
        table_entry(f)
    else:
        table_entry(ft)
        table_entry(fb)
        table_line()
    ft, fb, f = tabulate(method_mpfi2_direct, "MPF|MPFI-256+MPQ")
    if wide_format:
        table_entry(f)
    else:
        table_entry(ft)
        table_entry(fb)

#    compare_mpq_hybrid = average_time[method_mpq]/average_time[method_mpfi]
#    compare_mpq_mpf = average_time[method_mpq]/average_time[method_mpf]
#    compare_mpf_hybrid = average_time[method_mpf]/average_time[method_mpfi]
#    command("avgMpqHybrid", compare_mpq_hybrid)
#    command("avgMpqMpf", compare_mpq_mpf)
#    command("avgMpfHybrid", compare_mpf_hybrid)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
