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
    sys.stderr.write("Usage: %s [-v] [-d PATH] [-m (c|t|e)] [-o OUT]\n" % name)
    sys.stderr.write("  -v      Verbose\n")
    sys.stderr.write("  -d PATH Directory with CSV files\n")
    sys.stderr.write("  -m MODE Graphing mode: count (c) time (t) effort (e)\n")
    sys.stderr.write("  -o OUT  Specify output file\n")

directory = "."
verbose = False

constantFactor = 7

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
    dbl, m64, m128, m256, mpq = range(precisionCount)
    nonnegNames = ["dbl", "mpf64", "mpf128", "mpf256", "mpq"]
    negposNames = ["", "mpfi64", "mpfi128", "mpfi256", "mpq"]
    bits = [53, 63, 127, 255, 10000]
    
    nn_dbl, nn_m64, nn_m128, nn_m256, np_m64, np_m128, np_m256, all_mpq = range(tabulationCount)
    tabulateNames = ["double", "mpf-64", "mpf-128", "mpf-256", "mpfi-64" "mpfi-128", "mpfi-256", "mpq"]
    tabulateColors = ["dbl", "mpflow", "mpfmed", "mpfhigh", "mpfilow", "mpfimed", "mpfihigh", "mpq"]


    nn_map = {}
    np_map = {}

    def __init__(self):
        self.nn_map = { self.dbl   :self.nn_dbl,
                        self.m64   :self.nn_m64,
                        self.m128  :self.nn_m128,
                        self.m256  :self.nn_m256,
                        self.mpq   :self.all_mpq}
        self.np_map = { self.m64   : self.np_m64,
                        self.m128  : self.np_m128,
                        self.m256  : self.np_m256,
                        self.mpq   :self.all_mpq}

    # Accumulate information about solutions
    # histo is array of tabulationCount entries
    def newHistogram(self):
        return [0] * tabulationCount

    # Times for individual methods
    def accumulateEffort(self, isNonnegative, timeList, histo):
        mint = 0 if isNonnegative else 1
        for t in range(mint, precisionCount):
            tidx = self.nn_map[t] if isNonnegative else self.np_map[t]
            histo[tidx] += timeList[t] / 3600.0

    # Final method
    def accumulateCount(self, isNonnegative, ptype, histo):
        tidx = self.nn_map[ptype] if isNonnegative else self.np_map[ptype]
        histo[tidx] += 1

    # Total time by final method
    def accumulateTime(self, isNonnegative, ptype, total, histo):
        tidx = self.nn_map[ptype] if isNonnegative else self.np_map[ptype]
        histo[tidx] += total / 3600.0


    def guaranteedPrecision(self, ptype, nvar):
        return guaranteedPrecision(nvar, self.bits[ptype])

ptyper = None


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
                if t == ptyper.dbl:
                    sval = csvDict[tkey + "-prec"]
                    if len(sval) == 0:
                        sval = "0.0"
                    if float(sval) > 0:
                        self.precisions[t] = ptyper.guaranteedPrecision(t, self.nvar)                        
                    else:
                        self.precisions[t] = 0.0
                elif t == ptyper.mpq:
                    self.precisions[t] = 1e6
                else:
                    self.precisions[t] = ptyper.guaranteedPrecision(t, self.nvar)
                sval = csvDict[tkey + "-sec"]
                if len(sval) == 0:
                    sval = "0.0" if t == ptyper.dbl else str(failTime(self.bench))
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
        return ptyper.guaranteedPrecision(ptype, self.nvar) > digitPrecision

    def solve(self, digitPrecision):
        s = Solution(self, digitPrecision)
        if self.isNonnegative:
            for t in range(precisionCount):
                if digitPrecision <= self.precisions[t]:
                    s.time[t] = self.seconds[t]
                    s.ptype = t
                    s.achievedPrecision = self.precisions[t]
                    break
                elif t == ptyper.dbl and self.include(t, digitPrecision):
                    s.time[t] = self.seconds[t]

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
    count, effort, time = range(3)
    modeCharacters = ['c', 'e', 't']

    def __init__(self, targetPrecision):
        self.targetPrecision = 0
        self.solutionList = []

    def getMode(self, c):
        for i in range(len(self.modeCharacters)):
            if c == self.modeCharacters[i]:
                return i
        return 0

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
    
    def tabulate(self, mode):
        if mode == self.count:
            return self.tabulateCount()
        elif mode == self.effort:
            return self.tabulateEffort()
        else:
            return self.tabulateTime()

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
    solutionPoints = [1] + [5*i for i in range(1, 15)]

    def __init__(self, instances):
        self.instances = instances
        self.solutionSetList = []
        for d in self.solutionPoints:
            sset = instances.solve(d)
            self.solutionSetList.append(sset)
        
    def format(self, modeCharacter, outfile):
        mode = self.solutionSetList[0].getMode(modeCharacter)
        histoList = [ss.tabulate(mode) for ss in self.solutionSetList]
        for t in range(tabulationCount):
            outfile.write("\\addplot+[ybar, %s] plot coordinates {" % ptyper.tabulateColors[t])
            for i in range(len(self.solutionSetList)):
                d = self.solutionPoints[i]
                v = histoList[i][t]
                outfile.write("(%d,%.3f)" % (d, v))
            outfile.write("};\n")

def run(name, args):
    fname = ["nonneg-tabulate.csv", "negpos-tabulate.csv"]
    nonnegative = [True, False]
    global verbose
    global directory
    global ptyper
    ptyper = PrecisionType()
    modeCharacter = 'c'
    outName = None
    optList, args = getopt.getopt(args, "hvd:m:o:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-v':
            verbose = True
        elif opt == '-d':
            directory = val
        elif opt == '-m':
            modeCharacter = val
        elif opt == '-o':
            outName = val
    instances = InstanceSet()
    for i in range(2):
        name = directory + "/" + fname[i]
        try:
            infile = open(name, 'r')
        except:
            print("Can't open file '%s'" % name)
            continue
        count = instances.load(infile, nonnegative[i])
        infile.close()
        print("Loaded %d instances from %s" % (count, name))
    srange = SolutionRange(instances)
    outfile = sys.stdout
    if outName is not None:
        try:
            outfile = open(outName, 'w')
        except:
            print("Couldn't open output file '%s'" % outName)
            return 1
    srange.format(modeCharacter, outfile)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)
    
    
