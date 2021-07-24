#!/usr/bin/env python
"""
Usage
-----
1. run under the targetdir.
2. %prog N (exception_epoch1 exception_epoch2 ...), where N stands for the Nth latest commit,
    and the exceptions are provided as sys.argv.
"""
import os, sys, glob
#import mspsrpifun, howfun
#from optparse import OptionParser

################################################################################
## Main code
################################################################################
N = int(sys.argv[1])
exceptions = sys.argv[2:]
expnos = []
vexfiles = glob.glob(r'*/*.vex')
vexfiles.sort()
for vexfile in vexfiles:
    expno = vexfile.split('/')[-2].strip()
    if expno not in exceptions:
        expnos.append(expno)

tempfile = '.the_last_git_commits'
if os.path.exists(tempfile):
    os.remove(tempfile)
os.system('''run_command_to_all_epochs.py -c "git log --oneline|sed -n '%dp'" >> %s''' % (N, tempfile))
lines = open(tempfile).readlines()[:len(expnos)]
os.remove(tempfile)

commit_NOs = []
for line in lines:
    print(line)
    commit_NO = line.split(' ')[0].strip()
    commit_NOs.append(commit_NO)

targetdir = os.getcwd()
print("Now git checkout to the previous commit for all epochs...")
for i in range(len(expnos)):
    expdir = targetdir + '/' + expnos[i]
    os.chdir(expdir)
    os.system('git checkout %s' % commit_NOs[i])
print("\ngit checkout done!")
