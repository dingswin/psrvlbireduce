#!/usr/bin/env python
###########################################################################
###-- remove sn,bp,ps tables in "tables" folders and rerun final*.py for a set of observations
###-- tailored for MSPSRPI
###-- if -e bd223a, just rerun bd223a, if -s J1012+5307, rerun all epochs regarding the target 
############################################################################
import os,sys,glob,yaml,howfun,datetime
from optparse import OptionParser

auxdir    = os.environ['PSRVLBAUXDIR']
codedir   = os.environ['PSRVLBICODEDIR']
configdir = auxdir + '/configs/'

usage = "usage: %prog []\n-e --experiment\n-t --target\n-p --prepare\n-o --prepareonly\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Experiment name to rerun the final_astrometric.py")
parser.add_option("-t", "--target", dest="target", default="",
                  help="target to rerun the final_astrometric_reduce.py")
parser.add_option("-p", "--prepare", dest="prepare", default=False,
                  action="store_true",help="run prepare_astrometric_epoch.py alongside final_astrometric_reduce.py")
parser.add_option("-o", "--prepareonly", dest="prepareonly", default=False,
                  action="store_true",help="run prepare_astrometric_epoch.py only")
(options, junk) = parser.parse_args()
targetname      = options.target
experiment      = options.experiment
prepare         = options.prepare
prepareonly     = options.prepareonly
targetdir = auxdir + '/processing/' + targetname

def exp2expdir(string):
    expconfigfile = configdir + string + '.yaml'
    print expconfigfile
    if not os.path.exists(expconfigfile):
        parser.error("Experiment config file %s does not exist!" % expconfigfile)
    expconfig     = yaml.load(open(expconfigfile))
    expdir        = expconfig['rootdir'] + '/' + experiment + '/'
    return expdir

if ((experiment=='') and (targetname=='')) or ((experiment!='') and (targetname!='')):
    print usage
    sys.exit()
if experiment!='' and targetname=='':
    expdir = exp2expdir(experiment)
    if prepare:
        os.chdir(expdir)
        print "\nrun prepare_astrometric_epoch.py...\n"
        os.system("prepare_astrometric_epoch.py %s.vex" % (experiment))
    if prepareonly:
        sys.exit()
    print "deleting sn, bp and ps tables...\n"
    try:
        os.system('rm %s/tables/{*sn,*bp,*ps}' % (expdir))
    except OSError:
        print "Some of the tables are already removed\n"
    os.system('final_astrometric_reduce.py -e %s --clearcatalog' % (experiment))
    sys.exit()

# if experiment=='' and targetname!=''
targetdir = auxdir + '/processing/' + targetname
if not os.path.exists(targetdir):
    print("%s folder doesn't exist; aborting\n" % (targetdir))
    sys.exit()

sys.stdout = howfun.Logger(targetdir + "/rerun_finalastrometricreduce_runlog.txt")
current_time=datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
print "%s\nRunning script at %s (UTC)\n%s\n" % (70*'=',current_time,70*'=')

tablefolders=glob.glob(r'%s/*/tables' % (targetdir)) 
tablefolders.sort() 
print 'tablefolders=' 
print tablefolders
for tablefolder in tablefolders:
    experiment = tablefolder.split('/')[-2].strip()
    expdir     = exp2expdir(experiment)
    if prepare:
        os.chdir(expdir)
        print "\nrun prepare_astrometric_epoch.py...\n"
        os.system("prepare_astrometric_epoch.py %s.vex" % (experiment))
    if prepareonly:
        continue
    print "\ndeleting sn, bp and ps tables for %s...\n" % (experiment)
    try:
        os.system('rm %s/{*sn,*bp,*ps}' % (tablefolder))
    except OSError:
        print "Some of the tables are already removed\n"
    print "rerun final_astrometric_reduce.py for %s" % (experiment)
    #os.chdir(codedir)
    os.system('final_astrometric_reduce.py -e %s --clearcatalog' % (experiment))
current_time=datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
print "\nSuccessfully rerun throuhg all epochs for %s at %s (UTC)\n" % (targetname,current_time)
