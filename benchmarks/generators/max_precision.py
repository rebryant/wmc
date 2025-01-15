#!/usr/bin/python3

# Generate CNF file that maximizes required precision

def usage(name):
    sys.stderr.write("Usage: %s [-h] [-d DIGITS] -n COUNT -N LOGCOUNT -r ROOT\n" % name)

import sys
import getopt
import os

import readwrite


def weightString(scaled, sigdigs):
    wstring = ""
    for i in range(sigdigs):
        digit = scaled % 10
        scaled = scaled // 10
        wstring = str(digit) + wstring
    wstring = str(scaled) + "." + wstring
    return wstring

# Attempt at encoding as CNF.  Not quite right.
def generateCnf(root, n, sigdigs):
    cnfName = root + ".cnf"
    cnf = readwrite.CnfWriter(n+4, fname=cnfName, verbLevel=2)
    cnf.doHeaderComment("Problem requiring very high precision")
    cnf.doHeaderComment("Variables: %d+1, Max weight 10e%d, Minimum weight 10e-%d" % (n, sigdigs, sigdigs))
    xvars = list(range(2, n+2))
    nxvars = [-x for x in xvars]
    z  = 1
    wdict = {}
    wdict[z]  =  "1.0"
    wdict[-z] =  "-1.0"
    for x in xvars:
        wdict[x] = weightString(int(100**sigdigs), sigdigs)
        wdict[-x] = weightString(1, sigdigs)
    cnf.addWeights(wdict)
    cnf.doComment("  z --> /\\ x_i")
    for x in xvars:
        cnf.doClause([z, x])
    cnf.doComment("  !z --> (/\\ x_i  \\/  /\\ !x_i)")
    x1 = xvars[0]
    xrest = xvars[1:]
    for x in xrest:
        cnf.doClause([-z, -x1, x])
        cnf.doClause([-z,   x1, -x])
    cnf.finish()

def generateWeightedCnf(root, n, sigdigs):
    cnfName = root + ".cnf"
    cnf = readwrite.CnfWriter(n+1, fname=cnfName, verbLevel=2)
    cnf.doHeaderComment("Problem requiring very high precision")
    cnf.doHeaderComment("Variables: %d+1, Max weight 10e%d, Minimum weight 10e-%d" % (n, sigdigs, sigdigs))
    xvars = list(range(2, n+2))
    nxvars = [-x for x in xvars]
    wdict = {}
    z  = 1
    wdict[z]   =  "1.0"
    wdict[-z]  = "-1.0"
    for x in xvars:
        wdict[x] = weightString(int(100**sigdigs), sigdigs)
        wdict[-x] = weightString(1, sigdigs)
    cnf.addWeights(wdict)
    cnf.finish()


def writeList(file, ls):
    slist = [str(v) for v in ls]
    s = " ".join(slist)
    s += '\n'
    file.write(s)

def generateNnf(root, n):
    nnfName = root + ".nnf"
    try:
        nfile = open(nnfName, "w")
    except:
        sys.stderr.write("Can't open output file '%s'" % nnfName)
        sys.exit(1)
    z = 1
    xlist = list(range(2, n+2))
    nxlist = [-x for x in xlist]
    writeList(nfile, ['o', 1, 0])
    writeList(nfile, ['o', 2, 0])
    writeList(nfile, ['o', 3, 0])
    writeList(nfile, ['t', 4, 0])
    writeList(nfile, [3, 4] + xlist + [0])
    writeList(nfile, [3, 4] + nxlist + [0])
    writeList(nfile, [2, 3, z, 0])
    writeList(nfile, [2, 4, -z] + xlist + [0])
    writeList(nfile, [1, 2, 0])
    nfile.close()
    

def run(name, args):
    n = None
    sigdigs = 9
    root = None
    optList, args = getopt.getopt(args, "hn:N:d:r:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-n':
            n = int(val)
        elif opt == '-N':
            pwr = int(val)
            n = int(10**pwr)
        elif opt == '-d':
            sigdigs = int(val)
        elif opt == '-r':
            root = val
    if n is None:
        sys.stderr.write("Need problem size\n")
        usage(name)
        return
    if root is None:
        sys.stderr.write("Need root name\n")
        usage(name)
        return
    generateCnf(root, n, sigdigs)
    generateNnf(root, n)
    
if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
    
