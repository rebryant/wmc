#!/usr/bin/python3

import sys

id = "000"


def trim(s):
    while len(s) > 0 and s[-1] in "\r\n":
        s = s[:-1]
    return s


def processFile(root):
    iname = root + ".cnf"
    oname = root + "_" + id + ".cnf"
    try:
        infile = open(iname, "r")
    except:
        sys.stderr.write("Couldn't open input file '%s'\n" % iname)
        return False
    try:
        outfile = open(oname, "w")
    except:
        sys.stderr.write("Couldn't open output file '%s'\n" % oname)
        return False
    sys.stderr.write("Writing file '%s'\n" % oname)
    for line in infile:
        if len(line) == 0:
            continue
        if line[0] == 'c':
            outfile.write(line)
        elif line[0] == 'p':
            fields = line.split()
            if len(fields) != 4 or fields[1] != "cnf":
                sys.stderr.write("Unexpected line '%s' in file %s\n" % (trim(line), iname))
                return False
            fields[3] = "0"
            nline = " ".join(fields) + "\n"
            outfile.write(nline)
        else:
            continue
    infile.close()
    outfile.close()
    return True

def getRoot(fname):
    fields = fname.split(".")
    if len(fields) > 1:
        fields = fields[:-1]
    return ".".join(fields)

def run(name, args):
    elimit = 10
    if len(args) == 0 or args[0] == '-h':
        sys.stderr.write("Usage %s FILE1 ... FILEn\n" % name)
        return
    for fname in args:
        root = getRoot(fname)
        if not processFile(root):
            elimit -= 1
            if elimit <= 0:
                sys.stderr.write("Too many errors.  Aborting\n")
                return
        
if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)

