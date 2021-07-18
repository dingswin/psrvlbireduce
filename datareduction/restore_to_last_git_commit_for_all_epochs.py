#!/usr/bin/env python
"""
Usage
-----
1. run under the targetdir.
2. %prog (exception_epoch1 exception_epoch2 ...), the exceptions are provided as sys.argv.
"""
import os, sys, glob
#import mspsrpifun, howfun
#from optparse import OptionParser

################################################################################
## Main code
################################################################################
exceptions = sys.argv[1:]
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
os.system("run_command_to_all_epochs.py -c 'git log --oneline|head -n 1' >> %s" % tempfile)
lines = open(tempfile).readlines()[:len(expnos)]
os.remove(tempfile)

commit_NOs = []
for line in lines:
    print(line)
    commit_NO = line.split(' ')[0].strip()
    commit_NOs.append(commit_NO)

targetdir = os.getcwd()
answer = raw_input("You must confirm before going further.\nDo you like to run 'git reset --hard ' to the last commit (see the above commits) for all epochs? (Y/N):\n")
if answer in ['N', 'n', 'no', 'No', 'NO']:
    print("Process terminated as required.")
    sys.exit()
else:
    print("Now git-reset (--hard) to the last commit for all epochs...")
for i in range(len(expnos)):
    expdir = targetdir + '/' + expnos[i]
    os.chdir(expdir)
    os.system('git reset %s --hard' % commit_NOs[i])
print("git-reset done!")
