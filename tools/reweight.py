#!/usr/bin/python3

# Generate weights and create declaration in CNF file

import sys
import random
import getopt

import readwrite

def usage(name):
    print("Usage: %s [-h] [-v VLEVEL] -i INCNF [-s SEED] [-n COUNT] [-d DIGITS]" % name)
    
def getRoot(path):
    fields = path.split("/")
    if len(fields) == 0:
        return ""
    fname = fields[-1]
    parts = fname.split(".")
    parts = parts[:-1]
    return ".".join(parts)

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

    def uniformProbabilities(self, vlist, digits=9):
        random.seed(self.seed)
        svlist = sorted(vlist)
        for v in svlist:
            self.uprob(v, digits)
        
    def uniformAllProbabilities(self, n, digits=9):
        vlist = [i+1 for i in range(n)]
        self.uniformProbabilities(vlist, digits)

    def uniformSubsetProbabilities(self, vset, digits):
        vlist = list(vset)
        self.uniformProbabilities(vlist, digits)
            
    def updateCnf(self, cnf):
        cnf.addWeights(self.wdict)

def redoWeights(root, inputCnf, seed=123456, digits=9):
    w = Weighter(seed)
    outName = root + "_" + str(seed) + ".cnf"
    cwriter = readwrite.CnfWriter(inputCnf.nvar, outName, verbLevel = 2)
    cwriter.doHeaderComment("Uniform %d-digit probabilities generated with seed %d" % (digits, seed))
    if inputCnf.showVariables is not None and len(inputCnf.showVariables) < inputCnf.nvar:
        cwriter.addShow(inputCnf.showVariables)
        w.uniformSubsetProbabilities(inputCnf.showVariables, digits)
    else:
        w.uniformAllProbabilities(inputCnf.nvar, digits)
    print("Got weights %s" % (w.wdict))
    cwriter.addWeights(w.wdict)
    cwriter.finish()
    print("Generated file '%s'" % outName)
        
def process(inPath, seed, digits, count):
    try:
        inputCnf = readwrite.CnfReader(inPath, 2, False)
    except Exception as ex:
        print("Couldn't read input file '%s'.  (%s)" % (inPath, str(ex)))
        return
    root = getRoot(inPath)
    for i in range(count):
        redoWeights(root, inputCnf, seed+i, digits)

def run(name, args):
    seed = 123456
    digits = 9
    count = 1
    inPath = None
    optList, args = getopt.getopt(args, "hi:s:d:n:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == "-i":
            inPath = val
        elif opt == "-s":
            seed = int(val)
        elif opt == "-d":
            digits = int(val)
        elif opt == "-n":
            count = int(val)
    if inPath is None:
        print("Need input file name")
        usage(name)
        return
    process(inPath, seed, digits, count)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
