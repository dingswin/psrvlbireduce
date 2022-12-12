#!/usr/bin/env python
import os, sys

if not len(sys.argv) == 3:
    print "Usage: %s <observation> <rootdir>" % sys.argv[0]
    sys.exit()

observation = sys.argv[1]
rootdir = sys.argv[2]
directory = rootdir + '/' + observation + '/'
logfile = directory + observation + '.datacheck.log'
if not os.path.exists(logfile):
    print logfile + " does not exist!"
    sys.exit()
try:
    auxdir = os.environ['PSRVLBAUXDIR']
except KeyError:
    print 'Environment variable PSRVLBAUXDIR must be defined'
    exit(0)

lines = open(logfile).readlines()
outputfile = directory + observation + '.pipelinesummary'
output = open(outputfile, "w")

writtenlines = []
for line in lines:
    if (("Runlevel" in line and ':' in line) or "solution" in line or ("rms" in line and "CALIB" in line)) and not line in writtenlines:
        output.write(line)
        writtenlines.append(line)

output.close()
#diagdir = auxdir +'/finaldiagnostics/' + observation + '/'
#if not os.path.exists(diagdir):
#    print diagdir + " does not exist! create one.."
#    os.system("mkdir %s" % diagdir)
#    sys.exit()
