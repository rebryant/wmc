#####################################################################################
# Copyright (c) 2022 Randal E. Bryant, Carnegie Mellon University
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

import os
import sys

class ReadWriteException(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return "ReadWrite Exception: " + str(self.value)

# Code for reading and generating CNF, order, schedule, and cpog proof files

def trim(s):
    while len(s) > 0 and s[-1] in '\r\n':
        s = s[:-1]
    return s

def addPrefix(path, prefix):
    pfields = path.split("/")
    pfields[-1] = prefix + pfields[-1]
    return "/".join(pfields)

def changeExtension(path, extension):
    fields = path.split('.')
    fields[-1] = extension
    return ".".join(fields)
    

tautologyId = 1000 * 1000 * 1000

# Clean up clause.
# Remove duplicates
# Sort in reverse order of variable number
# Don't allow clause to have opposite literals (returns tautologyId)
def cleanClause(literalList):
    slist = sorted(literalList, key = lambda v: -abs(v))
    if len(slist) == 0:
        return tuple(slist)
    if slist[0] == tautologyId:
        return tautologyId
    if slist[0] == -tautologyId:
        slist = slist[1:]
        if slist[0] == tautologyId:
            return tautologyId
    if len(slist) == 1:
        return tuple(slist)
    nlist = [slist[0]]
    for i in range(1, len(slist)):
        if slist[i-1] == slist[i]:
            continue
        if slist[i-1] == -slist[i]:
            return tautologyId
        nlist.append(slist[i])
    return tuple(nlist)

# Fix up set of input clauses
# Flag error if any tautologies
def cleanClauses(clist, check = True):
    nlist = []
    id = 0
    for clause in clist:
        id += 1
        nclause = cleanClause(clause)
        if nclause == tautologyId and check:
            raise ReadWriteException("Tautologous clause #%d: %s" % (id, str(clause)))
        nlist.append(nclause)
    return nlist


def regularClause(clause):
    return clause is not None

def showClause(clause):
    if clause is None:
        return "NONE"
    return str(clause)


# Read header of CNF file
# Save number of variables and clauses
class CnfHeaderReader():
    nvar = 0
    nclause = 0
    file = None

    def __init__(self, fname):
        lineNumber = 0
        self.nclause = 0
        self.nvar = 0

        try:
            file = open(fname, 'r')
        except:
            raise ReadWriteException("Could not open file '%s'" % fname)
        for line in file:
            lineNumber += 1
            line = trim(line)
            if len(line) == 0:
                continue
            fields = line.split()
            if len(fields) == 0:
                continue
            elif line[0] == 'c':
                pass
            elif line[0] == 'p':
                fields = line[1:].split()
                if len(fields) != 3 or fields[0] != 'cnf':
                    raise ReadWriteException("Line %d.  Bad header line '%s'.  Not cnf" % (lineNumber, line))
                try:
                    self.nvar = int(fields[1])
                    self.nclause = int(fields[2])
                except Exception:
                    raise ReadWriteException("Line %d.  Bad header line '%s'.  Invalid number of variables or clauses" % (lineNumber, line))
                break
            else:
                raise ReadWriteException("Line %d.  No header line.  Not cnf" % (lineNumber))
        file.close()


# Read CNF file.
# Save list of clauses, each is a list of literals (zero at end removed)
# Also saves comment lines
class CnfReader():
    file = None
    commentLines = []
    clauses = []
    nvar = 0
    verbLevel = 1
    showVariables = None
    tseitinVariables = None
    # Dict indexed by literal & with strings for weights
    weights = None
   
    def __init__(self, fname = None, verbLevel = 1, check = True):
        self.verbLevel = verbLevel
        self.showVariables = None
        self.tseitinVariables = None
        if fname is None:
            opened = False
            self.file = sys.stdin
        else:
            opened = True
            try:
                self.file = open(fname, 'r')
            except Exception:
                raise ReadWriteException("Could not open file '%s'" % fname)
        self.clauses = []
        self.commentLines = []
        try:
            self.readCnf(check)
        except Exception as ex:
            if opened:
                self.file.close()
            raise ex
        if opened:
            self.file.close()
        if self.showVariables is None:
            self.showVariables = set(range(1, self.nvar+1))

    def processShow(self, fields):
        self.showVariables = set([])
        for s in fields[3:-1]:
            try:
                var = int(s)
            except:
                msg = "Couldn't parse '%s' as number" % s
                raise ReadWriteException(msg)
            if var < 1 or var > self.nvar:
                msg = "Invalid input variable %d" % var
                raise ReadWriteException(msg)
            self.showVariables.add(var)

    def processTseitin(self, fields):
        self.tseitinVariables = set([])
        for s in fields[3:-1]:
            try:
                var = int(s)
            except:
                msg = "Couldn't parse '%s' as number" % s
                raise ReadWriteException(msg)
            if var < 1 or var > self.nvar:
                msg = "Invalid input variable %d" % var
                raise ReadWriteException(msg)
            self.tseitinVariables.add(var)

    def processWeight(self, fields):
        if not self.weights:
            self.weights = {}
        try:
            lit = int(fields[3])
        except:
            msg = "Couldn't parse '%s' as literal" % fields[3]
            raise ReadWriteException(msg)
        var = abs(lit)
        if var < 1 or var > self.nvar:
            msg = "Invalid variable %d" % var
            raise ReadWriteException(msg)
        self.weights[lit] = fields[4]
            


    # See if there's anything interesting in the comment
    def processComment(self, line):
        if self.nvar == 0:
            fields = line.split()
            if len(fields) == 3 and fields[1] == 't' and fields[2] in ['pmc', 'pwmc']:
                self.showVariables = set([])
        else:
            fields = line.split()
            if self.showVariables is not None and len(fields) >= 3 and fields[1] == 'p' and fields[2] == 'show':
                self.processShow(fields)
            elif self.tseitinVariables is not None and len(fields) >= 3 and fields[1] == 'p' and fields[2] == 'forget':
                self.processTseitin(fields)
            elif len(fields) == 5 and fields[1] == 'p' and fields[2] == "weight":
                self.processWeight(fields)

    def readCnf(self, check):
        lineNumber = 0
        nclause = 0
        self.nvar = 0
        clauseCount = 0
        for line in self.file:
            lineNumber += 1
            line = trim(line)
            if len(line) == 0:
                continue
            fields = line.split()
            if len(fields) == 0:
                continue
            elif line[0] == 'c':
                self.processComment(line)
                if self.verbLevel > 1:
                    self.commentLines.append(line)
            elif line[0] == 'p':
                fields = line[1:].split()
                if len(fields) != 3 or fields[0] != 'cnf':
                    raise ReadWriteException("Line %d.  Bad header line '%s'.  Not cnf" % (lineNumber, line))
                try:
                    self.nvar = int(fields[1])
                    nclause = int(fields[2])
                except Exception:
                    raise ReadWriteException("Line %d.  Bad header line '%s'.  Invalid number of variables or clauses" % (lineNumber, line))
            else:
                if nclause == 0:
                    raise ReadWriteException("Line %d.  No header line.  Not cnf" % (lineNumber))
                try:
                    lits = [int(s) for s in line.split()]
                except:
                    raise ReadWriteException("Line %d.  Non-integer field" % lineNumber)
                # Last one should be 0
                if lits[-1] != 0:
                    raise ReadWriteException("Line %d.  Clause line should end with 0" % lineNumber)
                lits = lits[:-1]
                if check:
                    # Check formatting
                    vars = sorted([abs(l) for l in lits])
                    if len(vars) == 0:
                        raise ReadWriteException("Line %d.  Empty clause" % lineNumber)                    
                    if vars[-1] > self.nvar or vars[0] == 0:
                        raise ReadWriteException("Line %d.  Out-of-range literal" % lineNumber)
                    for i in range(len(vars) - 1):
                        if vars[i] == vars[i+1]:
                            raise ReadWriteException("Line %d.  Opposite or repeated literal" % lineNumber)
                self.clauses.append(lits)
                clauseCount += 1
        if clauseCount != nclause:
            raise ReadWriteException("Line %d: Got %d clauses.  Expected %d" % (lineNumber, clauseCount, nclause))
        self.file.close()


# Generic writer
class Writer:
    outfile = None
    suffix = None
    verbLevel = 1
    expectedVariableCount = None

    def __init__(self, count, fname, verbLevel = 1):
        self.expectedVariableCount = count
        self.verbLevel = verbLevel
        try:
            self.outfile = open(fname, 'w')
        except:
            print("Couldn't open file '%s'. Aborting" % fname)
            sys.exit(1)

    def trim(self, line):
        while len(line) > 0 and line[-1] == '\n':
            line = line[:-1]
        return line

    def vcount(self):
        return self.expectedVariableCount

    def show(self, line):
        line = self.trim(line)
        if self.verbLevel > 2:
            print(line)
        if self.outfile is not None:
            self.outfile.write(line + '\n')

    def finish(self):
        if self.outfile is None:
            return
        self.outfile.close()
        self.outfile = None

# Creating CNF
class CnfWriter(Writer):
    clauseCount = 0
    headerList = []
    outputList = []
    mcClass = "mc"

    # Track which variables actually occur
    vset = set([])

    def __init__(self, count, fname = None, verbLevel = 1):
        Writer.__init__(self, count, fname, verbLevel = verbLevel)
        self.clauseCount = 0
        self.headerList = []
        self.outputList = []
        self.vset = set([])

    # With CNF, must accumulate all of the clauses, since the file header
    # requires providing the number of clauses.

    def addShow(self, vars):
        svars = [str(var) for var in vars]
        self.doComment("p show %s 0" % " ".join(svars))
        if self.mcClass == "mc":
            self.mcClass = "pmc"
        elif self.mcClass == "wmc":
            self.mcClass = "pwmc"

    def addWeight(self, lit, weight):
        self.doComment("p weight %d %s 0" % (lit, str(weight)))
        if self.mcClass == "mc":
            self.mcClass = "wmc"
        elif self.mcClass == "pmc":
            self.mcClass = "pwmc"
        
    def addWeights(self, wdict):
        for lit in wdict.keys():
            self.addWeight(lit, wdict[lit])

    def doHeaderComment(self, line):
        self.headerList.append("c " + line)

    def doComment(self, line):
        self.outputList.append("c " + line)

    def doClause(self, literals):
        for lit in literals:
            var = abs(lit)
            if var <= 0 or var > self.expectedVariableCount:
                raise ReadWriteException("Variable %d out of range 1--%d" % (var, self.expectedVariableCount))
            self.vset.add(var)
        ilist = list(literals) + [0]
        self.outputList.append(" ".join([str(i) for i in ilist]))
        self.clauseCount += 1
        return self.clauseCount

    def variableCount(self):
        return len(self.vset)

    def finish(self, incremental = False):
        if self.outfile is None:
            return
        self.show("c t %s" % self.mcClass)
        for line in self.headerList:
            self.show(line)
        if incremental:
            self.show("p inccnf")
        else:
            self.show("p cnf %d %d" % (self.expectedVariableCount, self.clauseCount))
        for line in self.outputList:
            self.show(line)
        self.outfile.close()
        self.outfile = None


# Creating CNF incrementally.  Don't know number of variables in advance
class LazyCnfWriter:

    variableCount = 0
    # Set of tuples (T/F, item)
    # Boolean T for clause F for comment
    # item: list of literals for clause, string for comment
    items = []
    fname = ""
    verbLevel = 1
    clauseCount = 0

    def __init__(self, fname, verbLevel = 1):
        self.variableCount = 0
        self.items = []
        self.fname = fname
        self.verbLevel = verbLevel
        self.clauseCount = 0


    def newVariable(self):
        self.variableCount += 1
        return self.variableCount

    def vcount(self):
        return self.variableCount

    def newVariables(self, n):
        return [self.newVariable() for i in range(n)]
    
    def doComment(self, line):
        self.items.append((False, line))

    def doClause(self, lits):
        self.items.append((True, lits))
        self.clauseCount += 1

    def clauseList(self):
        clist = []
        for (isClause, value) in self.items:
            if isClause:
                clist.append(value)
        return clist

    def finish(self):
        writer = CnfWriter(self.variableCount, self.fname, self.verbLevel)
        for (isClause, value) in self.items:
            if isClause:
                writer.doClause(value)
            else:
                writer.doComment(value)
        writer.finish()
        print("c File '%s' has %d variables and %d clauses" % (self.fname, self.variableCount, writer.clauseCount))

class PogWriter(Writer):
    nextVariable = 0
    
    def __init__(self, variableCount, fname, verbLevel = 1):
        Writer.__init__(self, variableCount, fname, verbLevel = verbLevel)
        self.nextVariable = variableCount + 1

    def doOp(self, symbol, argList):
        var = self.nextVariable
        self.nextVariable += 1
        args = [symbol, var] + argList
        sargs = [str(a) for a in args]
        self.show(" ".join(sargs))

    def doComment(self, line):
        self.show("c " + line)

    def doAnd(self, argList):
        self.doOp('p', argList)

    def doOr(self, argList):
        self.doOp('s', argList)

    def doRoot(self, lit):
        self.show("r %d" % lit)
