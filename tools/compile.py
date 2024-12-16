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

import getopt
import sys
import os.path
import subprocess
import datetime
import time
import glob
import csv

import parallel

p = parallel.Printer()

def usage(name):
    p.print("Usage: %s [-h] [-f] [-D PATH] [-N THRDS] [-t TIME] [-l NFILE] [FILE.EXT ...]" % name)
    p.print("  -h       Print this message")
    p.print("  -f       Force regeneration of all files")
    p.print("  -D PATH  Directory for files")
    p.print("  -N THRDS Run N threads concurrently")
    p.print("  -t TIME  Limit time for each of the programs")
    p.print("  -l NFILE Specify file containing root names")
    p.print("  EXT can be any extension for wild-card matching (e.g., cnf, nnf)")

# Defaults

force = False
timeLimit = 1000
# Set to number if using multiple threads
threadCount = None

# Pathnames
def genProgramPath(progName, subdirectory = "bin"):
    ppath = os.path.abspath(__file__)
    parts = ppath.split('/')
    if len(parts) >= 2:
        parts[-2] = subdirectory
        parts[-1] = progName
    npath = "/".join(parts)
    if not os.path.exists(npath):
        raise Exception("Couldn't find program '%s'" % npath)
    return npath



commentChar = 'c'

def trim(s):
    while len(s) > 0 and not s[0].isalnum():
        s = s[1:]
    while len(s) > 0 and s[-1] in '\r\n':
        s = s[:-1]
    return s

def setTimeLimit(t):
    global timeLimit
    timeLimit = t

                
def runProgram(prefix, root, commandList, logFile, extraLogName = None, tlimit = None):
    result = ""
    cstring = " ".join(commandList)
    if tlimit is None:
        p.print("%s. %s: Running '%s' with no external time limit" % (root, prefix, cstring))
    else:
        p.print("%s. %s: Running '%s' with time limit of %d seconds" % (root, prefix, cstring, tlimit))
    logFile.write("%s LOG: Running %s\n" % (prefix, cstring))
    if tlimit is not None:
        logFile.write("%s LOG: Time limit %d seconds\n" % (prefix, tlimit))
    start = datetime.datetime.now()
    try:
        if tlimit is None:
            cp = subprocess.run(commandList, capture_output = True, text=True)
        else:
            cp = subprocess.run(commandList, capture_output = True, timeout=tlimit, text=True)            
    except subprocess.TimeoutExpired as ex:
        # Incorporate information recorded by external logging
        if (extraLogName is not None):
            try:
                xlog = open(extraLogName, "r")
                for line in xlog:
                    logFile.write(line)
                xlog.close()
            except:
                pass
        p.print("%s. %s Program timed out after %d seconds" % (root, prefix, tlimit))
        result += "%s ERROR: Timeout after %d seconds\n" % (prefix, tlimit)
        delta = datetime.datetime.now() - start
        seconds = delta.seconds + 1e-6 * delta.microseconds
        result += "%s LOG: Elapsed time = %.3f seconds\n" % (prefix, seconds)
        result += "%s OUTCOME: Timeout\n" % (prefix)
        logFile.write(result)
        logFile.close()
        return False
    ok = True
    if cp.returncode != 0:
        result += "%s ERROR: Return code = %d\n" % (prefix, cp.returncode)
        logFile.write("%s.  Output from stderr:" % (prefix))
        logFile.write(cp.stderr)
        ok = False
    outcome = "normal" if ok else "failed"
    delta = datetime.datetime.now() - start
    seconds = delta.seconds + 1e-6 * delta.microseconds
    result += "%s LOG: Elapsed time = %.3f seconds\n" % (prefix, seconds)
    result += "%s OUTCOME: %s\n" % (prefix, outcome)
    p.print("%s. %s: OUTCOME: %s" % (root, prefix, outcome))
    p.print("%s. %s: Elapsed time: %.3f seconds" % (root, prefix, seconds))
    logFile.write(cp.stdout)
    logFile.write(result)
    return ok

def genLogName(root, home):
    return root + ".d4v2_log"


def runD4(root, home):
    cnfExtension = ".cnf"
    cnfName = home + "/" + root + cnfExtension
    nnfName = home + "/" + root + ".nnf"
    logName = genLogName(root, home)
    cmd = [genProgramPath("d4v2")]
    cmd += ["-i", cnfName]
    cmd += ["-m", "ddnnf-compiler"]
    cmd += ["--dump-ddnnf", nnfName]
#    cmd += ["--quiet", "on"]
    try:
        logFile = open(logName, "w")
    except:
        p.print("%s ERROR:Couldn't open file '%s'" % (root, logName))
        return
    prefix = "D4v2"
    ok = runProgram(prefix, root, cmd, logFile, tlimit = timeLimit)
    return ok

def stripSuffix(fname):
    fields = fname.split(".")
    if len(fields) > 1:
        fields = fields[:-1]
    return ".".join(fields)

def runJob(home, root, force):
    logName = genLogName(root, home)
    if os.path.exists(logName) and not force:
        p.print("File %s exists.  Skipping" % logName)
        return
    ok = runD4(root, home)
    if not ok:
        p.print("Error encountered running on %s" % root)
        sys.exit(1)
    p.print("File %s written" % logName)

    
class Job:
    home = None
    root = None
    force = None

    def __init__(self, home, root, force):
        self.home = home
        self.root = root
        self.force = force

    def run(self):
        runJob(self.home, self.root, self.force)

def runBatch(home, fileList, force):
    roots = [stripSuffix(f) for f in fileList]
    roots = [r for r in roots if r is not None]
    p.print("Running on roots %s" % roots)
    if threadCount is None:
        for r in roots:
            runJob(home, r, force)
    else:
        s = parallel.Scheduler(threadCount)
        p.activate()
        for r in roots:
            j = Job(home, r, force)
            s.schedule(j)
        s.wait()


def run(name, args):
    global force
    global threadCount
    home = "."
    nameFile = None
    optList, args = getopt.getopt(args, "hfD:N:t:l:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-f':
            force = True
        elif opt == '-D':
            home = val
        elif opt == '-N':
            threadCount = int(val)
        elif opt == '-t':
            setTimeLimit(int(val))
        elif opt == '-l':
            nameFile = val
        else:
            p.print("Unknown option '%s'" % opt)
            usage(name)
            return
    fileList = args
    if nameFile is not None:
        try:
            nfile = open(nameFile, 'r')
        except:
            p.print("Couldn't open name file '%s'" % nameFile)
            usage(name)
            return
        for line in nfile:
            fname = trim(line)
            fileList.append(fname)
        nfile.close
            
    runBatch(home, fileList, force)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)

