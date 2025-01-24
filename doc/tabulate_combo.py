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



directory = "."

method_count = 3

method_mpf, method_mpfi, method_mpq = range(method_count)
method_name = ["MPF", "MPFI", "MPQ"]
file_name = ["combo-mpf+mpq.csv", "combo-mpfi+mpq.csv", "combo-mpq+mpq.csv"]

class Instance:
    name = None
    method = None
    mpf_seconds = 0.0
    mpfi_seconds = 0.0
    mpq_seconds = 0.0
    
    def __init__(self, name, method, total_seconds, mpq_seconds):
        self.name = name
        self.mpf_seconds = 0
        self.mpfi_seconds = 0
        self.method = method
        self.mpq_seconds = mpq_seconds
        if method == method_mpf:
            self.mpf_seconds = total_seconds
        elif method == method_mpfi:
            self.mpfi_seconds = total_seconds
        else:
            self.mpfi_seconds = total_seconds - mpq_seconds

    def __str__(self):
        return "I.%s (%s) MPF:%.3f MPFI:%.3f MPQ:%.3f" % (self.name, method_name[self.method], self.mpf_seconds, self.mpfi_seconds, self.mpq_seconds)

instances = []
average_time = [0] * method_count

def load(verbose):
    global instances
    for method in range(3):
        icount = 0
        fname = directory + "/" + file_name[method]
        try:
            file = open(fname)
        except:
            print("Couldn't open file '%s'" % fname)
            sys.exit(1)
        creader = csv.reader(file)
        line = 0
        for entry in creader:
            line+=1
            if len(entry) != 3:
                print("File %s, Line %s.  Bad entry: Not enough fields" % (fname, line))
                sys.exit(1)
            name = entry[0]
            try:
                total_seconds = float(entry[1])
                mpq_seconds = float(entry[2])
            except:
                print("File %s, Line %s (%s).  Bad entry: Couldn't parse numbers" % (fname, line, name))
                sys.exit(1)
            instance = Instance(name, method, total_seconds, mpq_seconds)
            if (verbose):
                print("Created %s" % str(instance))
            instances.append(instance)
            icount += 1
        file.close()
        if verbose:
            print("Created %d instances from %s" % (icount, fname))

# Here, method refers to final strategy:
# MPF: Use MPF + MPQ
# MPFI: Use MPF + MPFI + MPQ
# MPQ: Use MPQ
def tabulate(target_method, label):
    global average_time
    count = [0] * 3
    time = [0] * 3
    total_count = 0
    total_time = 0
    for instance in instances:
        total_count += 1
        if instance.method == method_mpf:
            if target_method in [method_mpf, method_mpfi]:
                count[method_mpf] += 1
                time[method_mpf] += instance.mpf_seconds
                total_time += instance.mpf_seconds
            else:
                count[method_mpq] += 1
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
        elif instance.method == method_mpfi:
            if target_method in [method_mpf, method_mpq]:
                count[method_mpq] += 1
                time[method_mpq] += instance.mpq_seconds
                total_time += instance.mpq_seconds
            else:
                count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
        else:
            # Method mpq
            count[method_mpq] += 1
            time[method_mpq] += instance.mpq_seconds
            total_time += instance.mpq_seconds
            if target_method == method_mpfi:
                count[method_mpfi] += 1
                time[method_mpfi] += instance.mpfi_seconds
                total_time += instance.mpfi_seconds
    fields = [label]
    for comp_method in range(3):
        mtime = time[comp_method]/3600.0
        fields += [str(count[comp_method]), "%.2f" % mtime]
    avg = total_time/total_count
    average_time[target_method] = avg
    ttime = total_time/3600.0
    fields += [str(total_count), "%.2f" % ttime, "%.1f" % avg]
    return fields

def table_entry(fields):
    print(" & ".join(fields) + ' \\\\')

def command(name, value):
    print("\\newcommand[\\global]{\\%s}{%s}" % (name, "%.2f" % value))


def run(name, args):
    global directory
    if len(args) > 0:
        if args[0] == '-h':
            print("Usage: %s [CSV_DIRECTORY]" % name)
            sys.exit(0)
        directory = args[0]
    load(False)
    fields = tabulate(method_mpq, "MPQ only")
    table_entry(fields)
    fields = tabulate(method_mpf, "MPF+MPQ")
    table_entry(fields)
    fields = tabulate(method_mpfi, "Hybrid")
    table_entry(fields)
    compare_mpq_hybrid = average_time[method_mpq]/average_time[method_mpfi]
    compare_mpq_mpf = average_time[method_mpq]/average_time[method_mpf]
    compare_mpf_hybrid = average_time[method_mpf]/average_time[method_mpfi]
    command("avgMpqHybrid", compare_mpq_hybrid)
    command("avgMpqMpf", compare_mpq_mpf)
    command("avgMpfHybrid", compare_mpf_hybrid)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
              
        
            

