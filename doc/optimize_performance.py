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
import math
import getopt

def usage(name):
    sys.stderr.write("Usage: %s [-v] [-l] [-c] [-d PATH] [-x XPREC] [-C (p|n|b)] [-m (c|t|e|a|w)] [-o OUT]\n" % name)
    sys.stderr.write("  -v       Verbose\n")
    sys.stderr.write("  -l       Use reduced set of digit precision values\n")    
    sys.stderr.write("  -c       Show output as CSV\n")
    sys.stderr.write("  -d PATH  Directory with CSV files\n")
    sys.stderr.write("  -x XPREC Extra digits when computing minimum MPFI precision\n")
    sys.stderr.write("  -C COLL  Specify which collections to include: positive (p), negative (n), or both (b)")
    sys.stderr.write("  -m MODE  Graphing mode: count (c) time (t) effort (e) avg. effort (a) work (w)\n")
    sys.stderr.write("  -o OUT   Specify output file\n")

directory = "."
verbose = False


extraDigits = 2

dataPoints = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]

limitedPoints = [1, 10, 15, 20, 30, 35, 40, 70]

constantFactor = 7

collection = "b"

csvOutput = False

def guaranteedPrecision(nvar, bitPrecision):
    return bitPrecision * math.log10(2) - math.log10(nvar) - math.log10(constantFactor)


# Hack to account for times taken by failed runs
failTimeDict = {"mc2024_track2-random_023" : 1592.0, "mc2024_track2-random_178" : 1634.2 }

def failTime(bench):
    formula = bench[:-4]
    if formula not in failTimeDict:
        sys.stderr.write("Couldn't find entry for instance %s (formula %s)\n" % (instance, formula))
        return 0.0
    return failTimeDict[formula]


precisionCount = 5
tabulationCount = 8

minPrecision = 1
maxPrecision = 70


class PrecisionType:
    erd, m64, m128, m256, mpq = range(precisionCount)
    nonnegNames = ["erd", "mpf64", "mpf128", "mpf256", "mpq"]
    negposNames = ["", "mpfi64", "mpfi128", "mpfi256", "mpq"]
    bits = [53, 63, 127, 255, 10000]
    
    nn_erd, nn_m64, nn_m128, nn_m256, np_m64, np_m128, np_m256, all_mpq = range(tabulationCount)
    tabulateNames = ["erd", "mpf-64", "mpf-128", "mpf-256", "mpfi-64", "mpfi-128", "mpfi-256", "mpq"]
    tabulateColors = ["erd", "mpflow", "mpfmed", "mpfhigh", "mpfilow", "mpfimed", "mpfihigh", "mpq"]


    nn_map = {}
    np_map = {}

    all_Types = []
    nn_Types = []
    np_Types = []

    def __init__(self):
        self.nn_map = { self.erd   :self.nn_erd,
                        self.m64   :self.nn_m64,
                        self.m128  :self.nn_m128,
                        self.m256  :self.nn_m256,
                        self.mpq   :self.all_mpq}
        self.np_map = { self.m64   : self.np_m64,
                        self.m128  : self.np_m128,
                        self.m256  : self.np_m256,
                        self.mpq   :self.all_mpq}
        self.nn_types = [self.nn_erd, self.nn_m64, self.nn_m128, self.nn_m256]
        self.np_types = [self.np_m64, self.np_m128, self.np_m256, self.all_mpq]
        self.all_types = self.nn_types + self.np_types
        
    # Accumulate information about solutions
    # histo is array of tabulationCount entries
    def newHistogram(self):
        return [0] * tabulationCount

    def scaleHistogram(self, histo, scale):
        return [h/scale  for h in histo]

    # Times for individual methods
    def accumulateEffort(self, isNonnegative, timeList, histo):
        mint = 0 if isNonnegative else 1
        for t in range(mint, precisionCount):
            tidx = self.nn_map[t] if isNonnegative else self.np_map[t]
            histo[tidx] += timeList[t]

    # Final method
    def accumulateCount(self, isNonnegative, ptype, histo):
        tidx = self.nn_map[ptype] if isNonnegative else self.np_map[ptype]
        histo[tidx] += 1

    # Total time by final method
    def accumulateTime(self, isNonnegative, ptype, total, histo):
        tidx = self.nn_map[ptype] if isNonnegative else self.np_map[ptype]
        histo[tidx] += total / 3600.0

    def divideWork(self, ptype, timeList):
        redundant = 0.0
        for t in range(ptype):
            redundant += timeList[t]
        real = timeList[ptype]
        return (real / 3600.0, redundant / 3600.0)
        


    def guaranteedPrecision(self, ptype, nvar):
        return guaranteedPrecision(nvar, self.bits[ptype])

ptyper = None

class Mode:
    percent, count, effort, aeffort, time, work  = range(6)
    modeCharacters = ['p', 'c', 'e', 'a', 't', 'w']

    def parseMode(self, ch):
        for m in range(len(self.modeCharacters)):
            if self.modeCharacters[m] == ch:
                return m
        raise "Can't parse mode character '%s'" % ch

moder = None


# Characterize optimal strategy for one instance and one target precision
class Solution:
    instance = None
    targetPrecision = 0
    achievedPrecision = 0
    ptype = None
    time = [0] * precisionCount
    
    def __init__(self, instance, targetPrecision):
        self.instance = instance
        self.targetPrecision = targetPrecision
        self.ptype = -1
        self.achievedPrecision = 0
        self.time = [0] * precisionCount

    def elapsed(self):
        return sum(self.time)

    def isNonnegative(self):
        return self.instance.isNonnegative

    def accumulateTime(self, histo):
        ptyper.accumulateTime(self.isNonnegative(), self.ptype, self.elapsed(), histo)

    def accumulateCount(self, histo):
        ptyper.accumulateCount(self.isNonnegative(), self.ptype, histo)

    def accumulateEffort(self, histo):
        ptyper.accumulateEffort(self.isNonnegative(), self.time, histo)

    def divideWork(self):
        return ptyper.divideWork(self.ptype, self.time)

    def __str__(self):
        s = self.instance.name()
        s += "Targ:%.1f Ach:%.1f Sec=%.3f [" % (self.targetPrecision, self.achievedPrecision, self.elapsed())
        mint = 0 if self.isNonnegative else 1
        stimes = ["%.3f" % self.time[t] for t in range(mint, precisionCount)]
        stimes[self.ptype] += "*"
        s += ':'.join(stimes) + "]"
        return s


class Instance:
    bench = ""
    nvar = 0
    isNonnegative = True
    precisions = [None] * precisionCount
    seconds =  [None] * precisionCount

    def __init__(self, isNonnegative, csvDict):
        self.isNonnegative = isNonnegative
        self.bench = csvDict["bench"]
        self.nvar = int(csvDict["var"])
        self.precisions = [None] * precisionCount
        self.seconds =  [None] * precisionCount
        if (isNonnegative):
            for t in range(precisionCount):
                tkey = ptyper.nonnegNames[t]
                if t == ptyper.mpq:
                    self.precisions[t] = 1e6
                else:
                    self.precisions[t] = ptyper.guaranteedPrecision(t, self.nvar)
                sval = csvDict[tkey + "-sec"]
                if len(sval) == 0:
                    sval = str(failTime(self.bench))
                self.seconds[t] = max(float(sval), 0.001)
        else:
            for t in range(1, precisionCount):
                tkey = ptyper.negposNames[t]
                if t == ptyper.mpq:
                    self.precisions[t] = 1e6
                else:
                    sval = csvDict[tkey + "-prec"]
                    if len(sval) == 0:
                        sval = "0.0"
                    self.precisions[t] = float(sval)
                sval = csvDict[tkey + "-sec"]
                if len(sval) == 0:
                    sval = str(failTime(self.bench))
                self.seconds[t] = max(float(sval), 0.001)

    def name(self):
        return self.bench + " V=%d T=%s " % (self.nvar, "NN" if self.isNonnegative else "NP")

    def __str__(self):
        s = self.name()
        if self.isNonnegative:
            for t in range(precisionCount):
                s += "[%s, D=%.3f, S=%.3f]" % (ptyper.nonnegNames[t], self.precisions[t], self.seconds[t])
        else:
            for t in range(1, precisionCount):
                s += "[%s, D=%.3f, S=%.3f]" % (ptyper.negposNames[t], self.precisions[t], self.seconds[t])
        return s
            
    # Would it be worthwhile to attempt this precision level
    def include(self, ptype, digitPrecision):
        if digitPrecision >= 70:
            return ptyper.guaranteedPrecision(ptype, self.nvar) > digitPrecision 
        else:
            return ptyper.guaranteedPrecision(ptype, self.nvar) > digitPrecision + extraDigits

    def solve(self, digitPrecision):
        s = Solution(self, digitPrecision)
        if self.isNonnegative:
            for t in range(precisionCount):
                if digitPrecision <= self.precisions[t]:
                    s.time[t] = self.seconds[t]
                    s.ptype = t
                    s.achievedPrecision = self.precisions[t]
                    break
        else:
            for t in range(1, precisionCount):
                if self.include(t, digitPrecision):
                    s.time[t] = self.seconds[t]
                    if self.precisions[t] >= digitPrecision:
                        s.ptype = t
                        s.achievedPrecision = self.precisions[t]
                        break
        return s

class SolutionSet:
    targetPrecision = 0
    solutionList = []

    def __init__(self, targetPrecision):
        self.targetPrecision = 0
        self.solutionList = []

    def addSolution(self, s):
        self.solutionList.append(s)

    def tabulateTime(self):
        histo = ptyper.newHistogram()
        for s in self.solutionList:
            s.accumulateTime(histo)
        return histo

    def tabulateCount(self):
        histo = ptyper.newHistogram()
        for s in self.solutionList:
            s.accumulateCount(histo)
        return histo

    def tabulateEffort(self):
        histo = ptyper.newHistogram()
        for s in self.solutionList:
            s.accumulateEffort(histo)
        return histo

    def tabulateWork(self, collection):
        (real, redundant) = (0.0, 0.0)
        for s in self.solutionList:
            if collection == 'n' and s.isNonnegative():
                continue
            if collection == 'p' and not s.isNonnegative():
                continue
            sreal, sredundant = s.divideWork()
            real += sreal
            redundant += sredundant
        return [real, redundant]

    def tabulate(self, mode, scale):
        if mode in [moder.percent, moder.count]:
            histo = self.tabulateCount()
        elif mode in [moder.effort, moder.aeffort]:
            histo = self.tabulateEffort()
        elif mode == moder.time:
            histo = self.tabulateTime()
        return ptyper.scaleHistogram(histo, scale)

        

class InstanceSet:
    instanceList = []

    def __init__(self):
        self.instanceList = []

    def load(self, csvFile, isNonnegative):
        creader = csv.DictReader(csvFile)
        count = 0
        for dict in creader:
            i = Instance(isNonnegative, dict)
            if verbose:
                print("Loaded instance %s" % str(i))
            self.instanceList.append(i)
            count += 1
        return count
    
    def solve(self, digitPrecision):
        solutions = SolutionSet(digitPrecision)
        for i in self.instanceList:
            s = i.solve(digitPrecision)
            solutions.addSolution(s)
            if verbose:
                print("Got solution %s" % str(s))
        return solutions

class SolutionRange:

    instances = None
    solutionSetList = []

    def __init__(self, instances):
        self.instances = instances
        self.solutionSetList = []
        for d in dataPoints:
            sset = instances.solve(d)
            self.solutionSetList.append(sset)
        
    def format(self, mode, outfile, types, count):
        scale = count if mode in [moder.percent, moder.aeffort] else 3600
        if mode == moder.percent:
            scale = scale * 0.01
        histoList = [ss.tabulate(mode, scale) for ss in self.solutionSetList]
        for t in types:
            outfile.write("\\addplot+[ybar, %s] plot coordinates {" % ptyper.tabulateColors[t])
            for i in range(len(self.solutionSetList)):
                d = dataPoints[i]
                v = histoList[i][t]
                outfile.write("(%d,%.3f)" % (d, v))
            outfile.write("};\n")

    def csvFormat(self, mode, outfile, types, count):
        scale = count if mode in [moder.percent, moder.aeffort] else 3600
        if mode == moder.percent:
            scale = scale * 0.01
        slist = ["Precision"] + [str(dataPoints[i]) for i in range(len(self.solutionSetList))]
        outfile.write(",".join(slist) + '\n')
        histoList = [ss.tabulate(mode, scale) for ss in self.solutionSetList]
        sums = [0 for i in range(len(self.solutionSetList))]
        for t in types:
            slist = [ptyper.tabulateNames[t]] + ["%.4f" % histoList[i][t] for i in range(len(self.solutionSetList))]
            sums = [sums[i] + histoList[i][t] for i in range(len(self.solutionSetList))]
            outfile.write(",".join(slist) + '\n')
        slist = ["Sum"] + ["%.4f" % sums[i] for i in range(len(self.solutionSetList))]
        outfile.write(",".join(slist) + '\n')        

    def csvFormatWork(self, outfile, collection, count):
        names = ["Real", "Redundant"]
        slist = ["Precision"] + [str(dataPoints[i]) for i in range(len(self.solutionSetList))]
        outfile.write(",".join(slist) + '\n')
        pairList = [ss.tabulateWork(collection) for ss in self.solutionSetList]
        sums = [0 for i in range(len(self.solutionSetList))]
        for t in range(2):
            slist = [names[t]] + ["%.4f" % pairList[i][t] for i in range(len(self.solutionSetList))]
            sums = [sums[i] + pairList[i][t] for i in range(len(self.solutionSetList))]
            outfile.write(",".join(slist) + '\n')
        slist = ["Sum"] + ["%.4f" % sums[i] for i in range(len(self.solutionSetList))]
        outfile.write(",".join(slist) + '\n')        
        fracs = [pairList[i][1]/sums[i] for i in range(len(self.solutionSetList))]
        flist = ["RFrac"] + ["%.4f" % fracs[i] for i in range(len(self.solutionSetList))]
        outfile.write(",".join(flist) + '\n')

def run(name, args):
    fname = ["nonneg-tabulate.csv", "negpos-tabulate.csv"]
    counts = [0, 0]
    nonnegative = [True, False]
    global verbose
    global directory
    global ptyper
    global moder
    global dataPoints
    global extraDigits
    global csvOutput
    global collection
    ptyper = PrecisionType()
    moder = Mode()
    modeCharacter = 'c'
    outName = None
    optList, args = getopt.getopt(args, "hvlcd:m:o:x:C:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-l':
            dataPoints = limitedPoints
        elif opt == '-v':
            verbose = True
        elif opt == '-c':
            csvOutput = True
        elif opt == '-d':
            directory = val
        elif opt == '-m':
            modeCharacter = val
        elif opt == '-o':
            outName = val
        elif opt == '-x':
            extraDigits = float(val)
        elif opt == '-C':
            collection = val
    instances = InstanceSet()
    
    for i in range(2):
        name = directory + "/" + fname[i]
        try:
            infile = open(name, 'r')
        except:
            print("Can't open file '%s'" % name)
            continue
        counts[i] = instances.load(infile, nonnegative[i])
        infile.close()

    types = ptyper.all_types
    count = counts[0] + counts[1]
    if collection == 'n':
        types = ptyper.np_types
        count = counts[1]
    if collection == 'p':
        types = ptyper.nn_types
        count = counts[0]

    srange = SolutionRange(instances)
    outfile = sys.stdout
    if outName is not None:
        try:
            outfile = open(outName, 'w')
        except:
            print("Couldn't open output file '%s'" % outName)
            return 1

    mode = moder.parseMode(modeCharacter)

    if mode == moder.work and not csvOutput:
        sys.stderr.write("Can only tabulate work when generating CSV\n")
        return 1

    if csvOutput:
        if mode == moder.work:
            srange.csvFormatWork(outfile, collection)
        else:
            srange.csvFormat(mode, outfile, types, count)
    else:
        srange.format(mode, outfile, types, count)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
    
    
