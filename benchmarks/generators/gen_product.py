#!/usr/bin/python3

# Generate CNF file based on product of variables
# Weight choices:
#   Default: Use list of weights in file
#   Uniform: Use specified weight for all vars

def usage(name):
    sys.stderr.write("Usage: %s [-h] [-u WEIGHT] [-s SEED] [-r u|e] [-c OUT.nnf] (-n COUNT | -p COUNT_LOG10)\n" % name)
    sys.stderr.write("  -u  WEIGHT  Used specified weight for all variables\n")
    sys.stderr.write("  -s  SEED    Set seed for randomly generated weights\n")
    sys.stderr.write("  -r  u|e     Use randomly generated weights (u: uniform, e:exponential)\n")
    sys.stderr.write("  -c OUT.nnf  Generate .nnf directly\n")
    sys.stderr.write("  -n COUNT    Specify variable count\n")
    sys.stderr.write("  -p CNT_L10  Specify LOG_10(count)\n")
    

import sys
import getopt
import os
import random
import math

import readwrite

weightFile = "weight-1000.txt"

weightList = []

defaultWeight = None

seed = 12345
randomMethod = None

def writeList(outFile, ls):
    sls = [str(v) for v in ls]
    outFile.write(" ".join(sls))
    outFile.write("\n")

def genNnf(outFile, n):
    writeList(outFile, ['o', 1, 0])
    writeList(outFile, ['t', 2, 0])
    ls = [1, 2] + [x+1 for x in range(n)] + [0]
    writeList(outFile, ls)


# For generating random weights with different distributions.  Copied from tools/reweight.py
def weightString(scaled, sigDigitCount):
    wstring = ""
    for i in range(sigDigitCount):
        digit = scaled % 10
        scaled = scaled // 10
        wstring = str(digit) + wstring
    wstring = str(scaled) + "." + wstring
    return wstring

class Weighter:
    seed = 123456
    wdict = None

    def __init__(self, seed = None):
        if seed is not None:
            self.seed = seed
        self.wdict = {}

    def reseed(self, seed):
        self.seed = seed

    def uprob(self, var, digits=9):
        smin = 1
        scaledOne = 10**digits
        smax = scaledOne - 1
        scaled = random.randint(smin, smax)
        nscaled = scaledOne - scaled
        self.wdict[var] = weightString(scaled, digits)
        self.wdict[-var] = weightString(nscaled, digits)

    def eprob(self, var, digits=9):
        smin = 1
        smax = (10**digits)/2
        lmin = math.log10(smin)
        lmax = math.log10(smax)
        lval = random.uniform(lmin, lmax)
        sval = int(10**lval)
        sval = max(smin, sval)
        sval = min(smax, sval)
        nsval = 10**digits - sval
        self.wdict[var] = weightString(sval, digits)
        self.wdict[-var] = weightString(nsval, digits)

    def uniformProbabilities(self, vlist, digits=9):
        random.seed(self.seed)
        svlist = sorted(vlist)
        for v in svlist:
            self.uprob(v, digits)
        
    def exponentialProbabilities(self, vlist, digits=9):
        random.seed(self.seed)
        svlist = sorted(vlist)
        for v in svlist:
            self.eprob(v, digits)

    def allProbabilities(self, n, digits=9, exponential=False):
        vlist = [i+1 for i in range(n)]
        if exponential:
            self.exponentialProbabilities(vlist, digits)
        else:
            self.uniformProbabilities(vlist, digits)

    def subsetProbabilities(self, vset, digits, exponential=False):
        vlist = list(sorted(vset))
        if exponential:
            self.exponentialProbabilities(vlist, digits)
        else:
            self.uniformProbabilities(vlist, digits)

    def updateCnf(self, cnf):
        cnf.addWeights(self.wdict)


def trim(s):
    while len(s) > 0 and s[0] in ' \t':
        s = s[1:]
    while len(s) > 0 and not s[0].isalnum():
        s = s[1:]
    while len(s) > 0 and s[-1] in '\r\n':
        s = s[:-1]
    return s

def initialize(progPath):
    global weightList
    weightList = []
    path = os.path.dirname(os.path.abspath(progPath)) + "/" + weightFile
    try:
        wfile = open(path, "r")
    except:
        sys.stderr.write("Couldn't open weight file '%s'\n", path)
        sys.exit(1)
    for line in wfile:
        line = trim(line)
        if len(line) > 0:
            weightList.append(line)
#    sys.stderr.write("Found %d weights\n" % len(weightList))
    wfile.close()
    if len(weightList) == 0:
        sys.stderr.write("No weights found in weight file '%s'" % weightFile)
        sys.exit(1)

def getWeight(v):
    if defaultWeight is None:
        idx = (v-1) % len(weightList) 
        return weightList[idx]
    else:
        return defaultWeight

def generate(name, n):
    cnf = readwrite.CnfWriter(n, fname=None, verbLevel=2)
    if defaultWeight is None:
        if randomMethod is None:
            cnf.doHeaderComment("Product of %d positive literals.  Weights extracted from file %s" % (n, weightFile))
        else:
            if randomMethod == 'u':
                descr = "uniformly"
                exponential = False
            elif randomMethod == 'e':
                descr = "exponentially"
                exponential = True
            else:
                sys.stderr.write("Unknown distribution type '%s'\n" % randomMethod)
                usage(name)
                sys.exit(1)
            cnf.doHeaderComment("Product of %d positive literals.  Weights generated %s with seed %d" % (n, descr, seed))
            w = Weighter(seed)
            w.allProbabilities(n, exponential=exponential)
            w.updateCnf(cnf)
    else:
        cnf.doHeaderComment("Product of %d positive literals.  All literal weights = %s" % (n, defaultWeight))
        for v in range(1, n+1):
            cnf.addWeight(v, getWeight(v))
    for v in range(1, n+1):
        cnf.doClause([v])
    cnf.finish()

def run(name, args):
    global defaultWeight
    global seed, randomMethod
    n = None
    needWeights = True
    outName = None
    optList, args = getopt.getopt(args, "hn:p:u:s:r:c:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-n':
            n = int(val)
        elif opt == '-p':
            n = 10**int(val)
        elif opt == '-s':
            seed = int(val)
        elif opt == '-r':
            randomMethod = val
            needWeights = False
        elif opt == '-u':
            defaultWeight = val
            needWeights = False
        elif opt == '-c':
            outName = val
    if n is None:
        sys.stderr.write("Need variable count\n")
        usage(name)
        return
    if (outName is not None):
        try:
            outFile = open(outName, "w")
        except:
            print("Couldn't open output file %s" % outName)
            return
        genNnf(outFile, n)
        outFile.close()
    else:
        if needWeights:
            initialize(name)
        generate(name, n)
    
if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
    
