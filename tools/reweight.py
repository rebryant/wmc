#!/usr/bin/python3

# Generate weights and create declaration in CNF file

import sys
import random
import getopt
import math

import readwrite

def usage(name):
    print("Usage: %s [-h] [-D u|e|b|s] [-R RANGE] [-C u|n|r|e|i|I] [-N n|r] [-p PCOUNT] [-s SEED] [-n COUNT] [-d DIGITS] IN1 IN2 ..." % name)
    print("  -h          Print this message")
    print("  -s SEED     Seed for first file")
    print("  -n COUNT    Generate COUNT variants for each input file")
    print("  -p PCOUNT   How many variables to select for partial derivative computation")
    print("  -D DIST     Specify distribution: uniform (u), single exponential (e), boundary values (b), constant (seed+min) (s)")
    print("  -R RANGE    Specify range of values.  Use open/closed interval notation with MIN,MAX ('o' for open, 'c' for closed)")
    print("  -C CMETHOD  What should be relation between W(x) and W(-x):")
    print("     sum-to-one (u), negated (n), reciprocal (r), equal (e)independent (i), or independent with a nonzero sum (I)")
    print("  -N NMETHOD  How should negative weights be generated: none (n), random (r)")         
    print("  -d DIGITS   Number of significant digits")
    
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
    
class ValueGenerator:
    codes = ['u', 'e', 'b', 's']
    names = ["uniform", "exponential", "boundary", "sweep"]
    uniform, exponential, boundary, sweep = range(4)
    seed = 123456
    minValue = 1
    maxValue = 1000*1000*1000-1
    method = None

    def __init__(self, minValue, maxValue, code):
        self.minValue = minValue
        self.maxValue = maxValue

        if code not in self.codes:
            raise Exception("Invalid value generator code '%s'" % method)
        for n in range(len(self.codes)):
            if code == self.codes[n]:
                self.method = n
                break
        random.seed(self.seed)
        if self.method == self.exponential and self.minValue <= 0:
            raise Exception("Invalid generator.  Can't have exponential distribution with zero or negative values")

    def name(self):
        return self.names[self.method]

    def reseed(self, seed):
        self.seed = seed
        random.seed(seed)

    def generate(self):
        if self.method == self.exponential:
            lmin = math.log10(self.minValue)
            lmax = math.log10(self.maxValue)
            lval = random.uniform(lmin, lmax)
            sval = int(10**lval)
            sval = max(self.minValue, sval)
            sval = min(self.maxValue, sval)
            return sval
        elif self.method == self.uniform:
            return random.randint(self.minValue, self.maxValue)
        elif self.method == self.boundary:
            return self.minValue if self.randomBool() else self.maxValue
        elif self.method == self.sweep:
            return self.minValue + self.seed
        else:
            raise Exception("Invalid value generation method %d" % (self.method))

    def randomBool(self):
        return random.randint(0,1) == 1

class ComplementHandler:
    codes = ['n', 'u', 'r', 'e', 'i', 'I']
    names = ["negative", "sum-to-one", "reciprocal", "equal", "independent", "independent-nonzero-sum"]
    negative, sumOne, reciprocal, equal, independent, independentNonzero = range(6)
    method = None
    digits = 9
    minValue = None
    maxValue = None

    def __init__(self, code, digits, minValue, maxValue):
        self.digits = digits
        self.minValue = minValue
        self.maxValue = maxValue
        if code not in self.codes:
            raise Exception("Invalid complement code '%s'" % code)
        for n in range(len(self.codes)):
            if code == self.codes[n]:
                self.method = n
                break

    def name(self):
        return self.names[self.method]

    def allowZeroSum(self):
        return self.method != self.independentNonzero

    def complement(self, value):
        if self.method == self.negative:
            return -value
        if self.method == self.sumOne:
            return int(10**self.digits) - value
        if self.method == self.equal:
            return value
        if self.method == self.reciprocal:
            if value == 0:
                return 0
            negative = value < 0
            value = abs(value)
            fval = float(value) * 10**-self.digits
            rval = (1.0/fval)
            ival = int(10**self.digits * rval)
            if ival < self.minValue:
                ival = self.minValue
            if ival > self.maxValue:
                ival = self.maxValue
            if negative:
                ival = -ival
            return ival
        # Independents
        return None

class NegationHandler:
    codes = ['n', 'r']
    names = ["none", "random"]
    none, random = range(2)
    method = None
    boolGenerator = None


    # Require method for generating random Booleans
    # Provide by generators randomBool method
    def __init__(self, code, boolGenerator):
        self.boolGenerator = boolGenerator
        if code not in self.codes:
            raise Exception("Invalid negation code '%s'" % code)
        for n in range(len(self.codes)):
            if code == self.codes[n]:
                self.method = n
                break

    def name(self):
        return self.names[self.method]

    def modify(self, value):
        if self.method == self.none:
            return value
        elif self.method == self.random:
            return -value if self.boolGenerator() else value
        else:
            raise Exception("Invalid negation generation method %d" % self.method)

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
    generator = None
    complementor = None
    negator = None

    wdict = None
    docString = ""

    def __init__(self, digits, rangeArg, distributionCode, complementCode, negationCode):
        self.digits = digits
        self.ranger = Range(rangeArg)
        scale = 10**self.digits
        minInteger = int(self.ranger.minValue * scale)
        if self.ranger.leftOpen:
            minInteger += 1
        maxInteger = int(self.ranger.maxValue * scale)
        if self.ranger.rightOpen:
            maxInteger -= 1
        self.generator = ValueGenerator(minInteger, maxInteger, distributionCode)
        self.negator = NegationHandler(negationCode, self.generator.randomBool)
        self.complementor = ComplementHandler(complementCode, digits, minInteger, maxInteger)
        self.wdict = {}
        cname = self.complementor.name()
        nname = self.negator.name()
        dname = self.generator.name()
        rstring = self.ranger.rangeString
        self.docString = "Digits: %d, Range: %s, Distribution: %s, Complements: %s, Negation: %s, Seed: " % (digits, rstring, dname, cname, nname)

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
        while len(slower) > 1 and slower[-1] == '0':
            slower = slower[:-1]
        return ("-" if negative else "") + supper + "." + slower

    def assignVariable(self, var):
        ival = self.generator.generate()
        ival = self.negator.modify(ival)
        sval = self.stringify(ival)
        self.wdict[var] = sval
        nival = self.complementor.complement(ival)
        if nival is None:
            nival = self.generator.generate()
            nival = self.negator.modify(nival)
            while not self.complementor.allowZeroSum() and nival == -ival:
                nival = self.generator.generate()
                nival = self.negator.modify(nival)
        snval = self.stringify(nival)
        self.wdict[-var] = snval

    def assignList(self, vlist):
        svlist = sorted(vlist)
        for v in svlist:
            self.assignVariable(v)

def redoWeights(root, inputCnf, seed, weighter, partialCount):
    weighter.reseed(seed)
    outName = root + "_" + str(seed) + ".cnf"
    cwriter = readwrite.CnfWriter(inputCnf.nvar, outName, verbLevel = 2)
    cwriter.doHeaderComment(weighter.document())
    
    if inputCnf.showVariables is not None and len(inputCnf.showVariables) < inputCnf.nvar:
        cwriter.addShow(inputCnf.showVariables)
        vlist = list(inputCnf.showVariables)
    else:
        vlist = range(1, inputCnf.nvar+1)
    if partialCount > 0:
        if partialCount > len(vlist):
            print("Cannot select %d partial derivative variables.  Only %d input variables" % (partialCount, len(vlist)))
            return
        dlist = sorted(random.sample(vlist, partialCount))
        cwriter.addPartial(dlist)
        weighter.reseed(seed)
    weighter.assignList(vlist)
    cwriter.addWeights(weighter.wdict)
    cwriter.finish()
    print("Generated file '%s'" % outName)
        
def process(inPath, seed, weighter, partialCount, count):
    try:
        inCnfName = replaceExtension(inPath, "cnf")
        inputCnf = readwrite.CnfReader(inCnfName, 2, False)
    except Exception as ex:
        print("Couldn't read input file '%s'.  (%s)" % (inCnfName, str(ex)))
        return
    root = getRoot(inPath)
    for i in range(count):
        redoWeights(root, inputCnf, seed+i, weighter, partialCount)

def run(name, args):
    seed = 123456
    digits = 9
    count = 1
    rangeArg = "o0,1o"
    complementCode = "u"
    distributionCode = "u"
    negationCode = "n"
    partialCount = 0
    
    optList, args = getopt.getopt(args, "hs:n:D:R:C:N:p:d:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == "-s":
            seed = int(val)
        elif opt == "-D":
            distributionCode = val
        elif opt == '-R':
            rangeArg = val
        elif opt == '-C':
            complementCode = val
        elif opt == '-N':
            negationCode = val
        elif opt == '-p':
            partialCount = int(val)
        elif opt == "-d":
            digits = int(val)
        elif opt == "-n":
            count = int(val)
    try:
        weighter = Weighter(digits, rangeArg, distributionCode, complementCode, negationCode)
    except Exception as ex:
        print("Invalid parameters: %s" % str(ex))
        return
    for inPath in args:
        process(inPath, seed, weighter, partialCount, count)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
