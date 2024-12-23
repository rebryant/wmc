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
    p.print("Usage: %s [-f] [-h] [-s] [-D SPATH] [-N THRDS] P1.NNF P2.NNF ... " % name)
    p.print("  -h       Print this message")
    p.print("  -f       Force regeneration of all files")
    p.print("  -D SPATH  Directory for source NNF files")
    p.print("  -N THRDS Run N threads concurrently")
    p.print("  -t TIME  Limit time for each of the programs")

# Defaults
force = False
# Set to number if using multiple threads
threadCount = None
# Use smoothing?
smooth = False

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

def getRoot(fname):
    fields = fname.split("/")
    name = fields[-1]
    parts = name.split(".")
    if len(parts) > 1:
        parts = parts[:-1]
    return ".".join(parts)

def runCommand(cmd):
    p.print("Running '%s'" % " ".join(cmd))
    try:
        cp = subprocess.run(cmd)
    except Exception as ex:
        p.print("Exception encountered (%s)" % str(ex))
        if cp.stderr is not None:
            p.print("Output on stderr")
            p.print(cp.stderr)
        return False
    ok = cp.returncode == 0
    if not ok:
        p.print("Command failed.  Return code = %d" % cp.returncode)
        if cp.stderr is not None:
            p.print("Output on stderr")
            p.print(cp.stderr)
    return ok

def runJob(nnfPath):
    root = getRoot(nnfPath)
    cnfNames = glob.glob(root + "*.cnf")
    if not force:
        newNames = []
        for name in cnfNames:
            cname = getRoot(name) + (".scount" if smooth else ".count")
            if not os.path.exists(cname):
                newNames.append(name)
        cnfNames = newNames
    if len(cnfNames) == 0:
        p.print("No new CNF files for NNF file %s" % nnfPath)
    else:
        cmd = [genProgramPath("nnfcount")]
        if smooth:
            cmd += ['-s']
        cmd += [nnfPath] + cnfNames
        ok = runCommand(cmd)
        if not ok:
            p.print("Error encountered running on %s" % root)

    
class Job:
    nnfName = None
    def __init__(self, name):
        self.nnfName = name

    def run(self):
        runJob(self.nnfName)

def runBatch(nnfList):
    roots = [getRoot(name) for name in nnfList]
    p.print("Running on roots %s" % str(roots))
    if threadCount is None:
        for f in nnfList:
            runJob(f)
    else:
        s = parallel.Scheduler(threadCount)
        p.activate()
        for f in nnfList:
            j = Job(home, f)
            s.schedule(j)
        s.wait()


def run(name, args):
    global force
    global threadCount
    global smooth
    home = None
    optList, args = getopt.getopt(args, "hfsD:N:")
    for (opt, val) in optList:
        if opt == '-h':
            usage(name)
            return
        elif opt == '-f':
            force = True
        elif opt == '-s':
            smooth = True
        elif opt == '-D':
            home = val
        elif opt == '-N':
            threadCount = int(val)
        else:
            p.print("Unknown option '%s'" % opt)
            usage(name)
            return
    nnfList = args
    if home is not None:
        nnfList += glob.glob(home + "/*.nnf")
    runBatch(nnfList)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
    sys.exit(0)

