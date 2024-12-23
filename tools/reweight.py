#!/usr/bin/python3

# Generate weights and create declaration in CNF file

import sys
import random
import getopt
import math

import readwrite

def usage(name):
    print("Usage: %s [-h] [-e] [-s SEED] [-n COUNT] [-d DIGITS] IN1 IN2 ..." % name)
    print("  -h        Print this message")
    print("  -s SEED   Seed for first file")
    print("  -n COUNT  Generate COUNT variants for each input file")
    print("  -e        Use exponential, rather than uniform, distribution")
    print("  -d DIGITS Number of significant digits")
    
def getRoot(path):
    fields = path.split("/")
    if len(fields) == 0:
        return ""
    fname = fields[-1]
    parts = fname.split(".")
    parts = parts[:-1]
    return ".".join(parts)

def replaceExtension(path, ext):
    fields = path.split(".")
    fields[-1] = ext
    return ".".join(fields)
    

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

def redoWeights(root, inputCnf, seed=123456, digits=9, exponential = False):
    w = Weighter(seed)
    outName = root + "_" + str(seed) + ".cnf"
    cwriter = readwrite.CnfWriter(inputCnf.nvar, outName, verbLevel = 2)
    descr = "Exponentially" if exponential else "Uniformly"
    cwriter.doHeaderComment("%s-distributed %d-digit probabilities generated with seed %d" % (descr, digits, seed))
    if inputCnf.showVariables is not None and len(inputCnf.showVariables) < inputCnf.nvar:
        cwriter.addShow(inputCnf.showVariables)
        w.subsetProbabilities(inputCnf.showVariables, digits, exponential)
    else:
        w.allProbabilities(inputCnf.nvar, digits, exponential)
    cwriter.addWeights(w.wdict)
    cwriter.finish()
    print("Generated file '%s'" % outName)
        
def process(inPath, seed, digits, count, exponential):
    try:
        inCnfName = replaceExtension(inPath, "cnf")
        inputCnf = readwrite.CnfReader(inCnfName, 2, False)
    except Exception as ex:
        print("Couldn't read input file '%s'.  (%s)" % (inCnfName, str(ex)))
        return
    root = getRoot(inPath)
    for i in range(count):
        redoWeights(root, inputCnf, seed+i, digits, exponential)

def run(name, args):
    seed = 123456
    digits = 9
    count = 1
    exponential = False
    optList, args = getopt.getopt(args, "hes:d:n:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == "-e":
            exponential = True
        elif opt == "-s":
            seed = int(val)
        elif opt == "-d":
            digits = int(val)
        elif opt == "-n":
            count = int(val)
    for inPath in args:
        process(inPath, seed, digits, count, exponential)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
