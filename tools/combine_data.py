#!/usr/bin/python3

import sys
import csv

# Merge multiple data sources with names of form XXXX_DDD into single, minimum value



def process(infile):
    entries = {}
    roots = []
    creader = csv.reader(infile)
    for fields in creader:
        name = fields[0]
        svalue = fields[1]
        parts = name.split("_")
        root = "_".join(parts[:-1])
        if root in entries:
            ovalue = float(entries[root])
            nvalue = float(svalue)
            if (nvalue < ovalue):
                entries[root] = svalue
        else:
            roots.append(root)
            entries[root] = svalue
    for root in roots:
        print("%s,%s" % (root, entries[root]))

def run(name, args):
    if len(args) == 0 or args[0] == '-h' or len(args) > 1:
        sys.stderr.write("Usage: %s FILE.CSV\n" % name)
        return
    try:
        infile = open(args[0], "r")
    except:
        sys.stderr.write("Couldn't open file '%s'\n" % args[0])
    process(infile)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
        
