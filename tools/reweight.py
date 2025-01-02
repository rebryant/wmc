#!/usr/bin/python3

# Generate weights and create declaration in CNF file

import sys
import random
import getopt
import math

import readwrite

def usage(name):
    print("Usage: %s [-h] [-D u|e] [-R RANGE] [-N u|n|i] [-s SEED] [-n COUNT] [-d DIGITS] IN1 IN2 ..." % name)
    print("  -h          Print this message")
    print("  -s SEED     Seed for first file")
    print("  -n COUNT    Generate COUNT variants for each input file")
    print("  -D DIST     Specify distribution: uniform (u), single exponential (e)")
    print("  -R RANGE    Specify range of values.  Use open/closed interval notation with MIN,MAX")
    print("  -N NMETHOD  What show be relation between W(x) and W(-x): sum-to-one (u), negation (n) or indepedent (i)")
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
    
class IntegerGenerator:
    seed = 123456
    minValue = 1
    maxValue = 1000*1000*1000-1
    isExponential = False

    def __init__(self, minValue, maxValue, isExponential=False):
        self.minValue = minValue
        self.maxValue = maxValue
        self.isExponential = isExponential
        random.seed(self.seed)
        if self.isExponential and self.minValue <= 0:
            raise Exception("Invalid generator.  Can't have exponential distribution over nonnegative range")

    def reseed(self, seed):
        self.seed = seed
        random.seed(seed)

    def generate(self):
        if self.isExponential:
            lmin = math.log10(self.minValue)
            lmax = math.log10(self.maxValue)
            lval = random.uniform(lmin, lmax)
            sval = int(10**lval)
            sval = max(self.minValue, sval)
            sval = min(self.maxValue, sval)
            return sval
        else:
            return random.randint(self.minValue, self.maxValue)

    def randomBool(self):
        return random.randint(0,1) == 1

class Negator:
    codes = ['i', 'c', 'u']
    names = ["independent", "complementary", "sum-to-one"]
    independent, complementary, sumOne = range(3)
    method = None
    digits = 9

    def __init__(self, method, digits):
        self.digits = digits
        if method not in self.codes:
            raise Exception("Invalid negation code '%s'" % method)
        for n in range(3):
            if method == self.codes[n]:
                self.method = n
                break

    def name(self):
        return self.names[self.method]

    def negate(self, value):
        if self.method == self.complementary:
            return -value
        if self.method == self.sumOne:
            return int(10**self.digits) - value
        return None

class Range:
    minValue = 0
    maxValue = 1
    leftOpen = True
    rightOpen = True
    rangeString = ""

    def __init__(self, rangeArg):
        if len(rangeArg) < 5:
            raise Exception("Invalid range argument '%s'. Too short" % rangeArg)

        if rangeArg[0] == 'c':
            self.leftOpen = False
        elif rangeArg[0] == 'o':
            self.leftOpen = True
        else:
            raise Exception("Invalid range argument '%s'. Must start with 'o' or 'c'" % rangeArg)            
        rangeArg = rangeArg[1:]

        if rangeArg[-1] == 'c':
            self.rightOpen = False
        elif rangeArg[-1] == 'o':
            self.rightOpen = True
        else:
            raise Exception("Invalid range argument '%s'. Must end with 'o' or 'c'" % rangeArg)            
        rangeArg = rangeArg[:-1]

        fields = rangeArg.split(',')
        if len(fields) != 2:
            raise Exception("Invalid range argument '%s'. Must have range of form MIN,MAX" % rangeArg)
        try:
            self.minValue = float(fields[0])
            self.maxValue = float(fields[1])
        except:
            raise("Invalid range argument '%s'.  Must have numeric range of form MIN,MAX" % rangeArg)
        lsymbol = "(" if self.leftOpen else "["
        rsymbol = ")" if self.rightOpen else "]"
        self.rangeString =  lsymbol + str(self.minValue) + "," + str(self.maxValue) + rsymbol


class Weighter:
    digits = 9
    ranger = None
    negator = None
    generator = None
    wdict = None
    docString = ""
    bilateral = False

    def __init__(self, digits, rangeArg, negationMethod, isExponential):
        self.digits = digits
        self.ranger = Range(rangeArg)
        self.bilateral = False
        if isExponential and self.ranger.minValue == -self.ranger.maxValue:
            self.bilateral = True
            self.ranger.minValue = 0
            self.ranger.leftOpen = True
        self.negator = Negator(negationMethod, digits)
        scale = 10**self.digits
        minInteger = int(self.ranger.minValue * scale)
        if self.ranger.leftOpen:
            minInteger += 1
        maxInteger = int(self.ranger.maxValue * scale)
        if self.ranger.rightOpen:
            maxInteger -= 1
        self.generator = IntegerGenerator(minInteger, maxInteger, isExponential)
        self.wdict = {}
        nname = self.negator.name()
        distribution = "exponential" if isExponential else "uniform"
        if self.bilateral:
            distribution += " (bilateral)"
        rstring = self.ranger.rangeString
        self.docString = "Digits: %d, Range: %s, Distribution: %s, Negation: %s, Seed: " % (digits, rstring, distribution, nname)

    def reseed(self, seed):
        self.generator.reseed(seed)

    def document(self):
        return self.docString + str(self.generator.seed)

    def stringify(self, ival):
        negative = ival < 0
        ival = abs(ival)
        scale = 10**self.digits
        upper = ival // scale
        lower = ival % scale
        supper = str(upper)
        slower = str(lower)
        while len(slower) < self.digits:
            slower = '0' + slower
        return ("-" if negative else "") + supper + "." + slower

    def assignVariable(self, var):
        ival = self.generator.generate()
        if self.bilateral and self.generator.randomBool():
            ival = -ival
        sval = self.stringify(ival)
        self.wdict[var] = sval
        nival = self.negator.negate(ival)
        if nival is None:
            nival = self.generator.generate()
            if self.bilateral and self.generator.randomBool():
                nival = -nival
        snval = self.stringify(nival)
        self.wdict[-var] = snval

    def assignList(self, vlist):
        svlist = sorted(vlist)
        for v in svlist:
            self.assignVariable(v)

def redoWeights(root, inputCnf, seed, weighter):
    weighter.reseed(seed)
    outName = root + "_" + str(seed) + ".cnf"
    cwriter = readwrite.CnfWriter(inputCnf.nvar, outName, verbLevel = 2)
    cwriter.doHeaderComment(weighter.document())

    
    if inputCnf.showVariables is not None and len(inputCnf.showVariables) < inputCnf.nvar:
        cwriter.addShow(inputCnf.showVariables)
        vlist = list(inputCnf.showVariables)
    else:
        vlist = range(1, inputCnf.nvar+1)
    weighter.assignList(vlist)
    cwriter.addWeights(weighter.wdict)
    cwriter.finish()
    print("Generated file '%s'" % outName)
        
def process(inPath, seed, weighter, count):
    try:
        inCnfName = replaceExtension(inPath, "cnf")
        inputCnf = readwrite.CnfReader(inCnfName, 2, False)
    except Exception as ex:
        print("Couldn't read input file '%s'.  (%s)" % (inCnfName, str(ex)))
        return
    root = getRoot(inPath)
    for i in range(count):
        redoWeights(root, inputCnf, seed+i, weighter)

def run(name, args):
    seed = 123456
    digits = 9
    count = 1
    isExponential = False
    rangeArg = "o0,1o"
    negationMethod = "u"
    optList, args = getopt.getopt(args, "hs:n:D:R:N:d:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == "-s":
            seed = int(val)
        elif opt == "-D":
            if val == 'u':
                isExponential = False
            elif val == 'e':
                isExponential = True
            else:
                print("Invalid distribution '%s'" % val)
                usage(name)
                return
        elif opt == '-R':
            rangeArg = val
        elif opt == '-N':
            negationMethod = val
        elif opt == "-d":
            digits = int(val)
        elif opt == "-n":
            count = int(val)
    print("Range = %s.  Negation = %s" % (rangeArg, negationMethod))
    try:
        weighter = Weighter(digits, rangeArg, negationMethod, isExponential)
    except Exception as ex:
        print("Invalid parameters: %s" % str(ex))
        return
    for inPath in args:
        process(inPath, seed, weighter, count)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
