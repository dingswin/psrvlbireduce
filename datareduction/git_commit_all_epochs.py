#!/usr/bin/env python
################################################################################
## run under the targetdir
################################################################################
import numpy as np
import os, sys, glob
import mspsrpifun, howfun
from optparse import OptionParser

################################################################################
## Main code
################################################################################
usage = "usage: %prog [options] (run under targetdir)"
parser = OptionParser(usage)
parser.add_option("-x", "--exceptions", dest="exceptions", default=[''], help="not to git commit on one or several epochs")
parser.add_option("-c", "--comment", dest="comment", default='', help="comment for commits")
(options, junk) = parser.parse_args()
comment      = options.comment
exceptions   = options.exceptions

if type(exceptions) == str: #only works for one exception!
    exceptions = [exceptions]

expnos = []
vexfiles = glob.glob(r'*/*.vex')
vexfiles.sort()
print(vexfiles)
for vexfile in vexfiles:
    expno = vexfile.split('/')[-2].strip()
    if expno not in exceptions:
        expnos.append(expno)
print(expnos)

targetdir = os.getcwd()
print(targetdir)
for expno in expnos:
    os.chdir(targetdir + '/' + expno)
    os.system('git add .')
    os.system("git commit -m %s" % comment)
