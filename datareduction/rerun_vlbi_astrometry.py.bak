#!/usr/bin/env python
###########################################################################
###-- remove sn,bp,ps tables in "tables" folders and rerun vlbi_astrometry.py for a set of observations
###-- tailored for MSPSRPI
###-- if -e bd223a, just rerun bd223a, if -s J1012+5307, rerun all epochs regarding the target 
############################################################################
import os,sys,glob,yaml,howfun,datetime
from optparse import OptionParser

def exp2expdir(string):
    expconfigfile = configdir + string + '.yaml'
    print expconfigfile
    if not os.path.exists(expconfigfile):
        parser.error("Experiment config file %s does not exist!" % expconfigfile)
    expconfig     = yaml.load(open(expconfigfile))
    expdir        = expconfig['rootdir'] + '/' + experiment + '/'
    return expdir
def prepare_given_exp(experiment):
    expdir = exp2expdir(experiment)
    tablefolder = expdir + '/tables/'
    if prepare:
        os.chdir(expdir)
        print "\nrun prepare_astrometric_epoch.py...\n"
        os.system("prepare_astrometric_epoch.py %s.vex" % (experiment))
def run_pipeline_given_exp(experiment, runlevel, skipdiagnosticplots=False):
    """
    Note
    ----
    1. Use vlbi_astrometry.py instead of  rerunvlbi_astrometry.py from a high runlevel if you want to keep the solution tables.
        However, a keep_talbes_runlevel_threshold is set (currently 32) to avoid unwanted deletions of solution tables.
    2. Since the tables will be deleted for runfromlevel < keep_tables_runlevel_threshold, one either set a small 
        runfromlevel where no tables have been made, or from >= keep_tables_runlevel_threshold. Otherwise, there will be
        complaint.

    Input parameters
    ----------------
    experiment : str
        Experiment number.
    runlevel : str
        Containing 'runfromlevel' or 'runfromlevel,runtolevel' info.
        e.g. '2', or '2,5'.
    skipdiagnosticplots : bool (default : False)
        If True, skip making diagnostic plots.
    """
    expdir = exp2expdir(experiment)
    tablefolder = expdir + '/tables/'
    keep_tables_runlevel_threshold = 32
    try:
        runfromlevel = int(runlevel)
    except ValueError:
        runfromlevel = int(runlevel.split(',')[0])
    if runfromlevel < keep_tables_runlevel_threshold:
        print "deleting sn, bp and ps tables...\n"
        try:
            os.system('rm %s/{*sn,*bp,*ps}' % (tablefolder))
        except OSError:
            print "Some of the tables are already removed\n"
    #solutions to keep being used need to save to *.sn.save or *.bp.save. 
    solutions_to_use = glob.glob(r'%s/*.save' % (tablefolder))
    for solution in solutions_to_use:
        os.system('cp %s %s' % (solution, solution.replace('.save', '')))
    if skipdiagnosticplots:
        os.system('vlbi_astrometry.py -e %s -r %s -k --alwayssaved' % (experiment, runlevel))
    else:
        os.system('vlbi_astrometry.py -e %s -r %s --alwayssaved' % (experiment, runlevel))

auxdir    = os.environ['PSRVLBAUXDIR']
codedir   = os.environ['PSRVLBICODEDIR']
configdir = auxdir + '/configs/'

usage = "usage: %prog []\n-e --experiment\n-t --target\n-p --prepare\n-o --prepareonly\n-r --runlevel\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Experiment name to rerun the vlbi_astrometry.py")
parser.add_option("-t", "--target", dest="target", default="",
                  help="target to rerun the vlbi_astrometry.py")
parser.add_option("-p", "--prepare", dest="prepare", default=False,
                  action="store_true",help="run prepare_astrometric_epoch.py alongside vlbi_astrometry.py")
parser.add_option("-o", "--prepareonly", dest="prepareonly", default=False,
                  action="store_true",help="run prepare_astrometric_epoch.py only")
parser.add_option("-r", "--runlevel", dest="runlevel", default=1,
                  help="runlevel at which to start (and stop), e.g. -r 1,2 or -r 2. Since the tables will be deleted for\
                  runfromlevel < keep_tables_runlevel_threshold, one either set a small runfromlevel where no tables have\
                  been made, or from >= keep_tables_runlevel_threshold. Otherwise, there will be complaint.")
parser.add_option("-k", "--skipdiagnosticplots", dest="skipdiagnosticplots", default=False,
                  action="store_true", help="Do not make diagnostic plots")
(options, junk) = parser.parse_args()
targetname      = options.target
experiment      = options.experiment
prepare         = options.prepare
prepareonly     = options.prepareonly
runlevel        = options.runlevel
skipdiagnosticplots = options.skipdiagnosticplots
targetdir = auxdir + '/processing/' + targetname


if ((experiment=='') and (targetname=='')) or ((experiment!='') and (targetname!='')):
    print usage
    sys.exit()
if experiment!='' and targetname=='':
    prepare_given_exp(experiment)
    if prepareonly:
        sys.exit()
    run_pipeline_given_exp(experiment, runlevel, skipdiagnosticplots)
    sys.exit()

# if experiment=='' and targetname!=''
targetdir = auxdir + '/processing/' + targetname
if not os.path.exists(targetdir):
    print("%s folder doesn't exist; aborting\n" % (targetdir))
    sys.exit()

sys.stdout = howfun.Logger(targetdir + "/rerun_vlbi_astrometry_runlog.txt")
current_time=datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
print "%s\nRunning script at %s (UTC)\n%s\n" % (70*'=',current_time,70*'=')

vexfiles = glob.glob(r'%s/*/*.vex' % (targetdir)) 
vexfiles.sort() 
print('vexfiles=\n %s' % vexfiles) 
for vexfile in vexfiles:
    experiment = vexfile.split('/')[-2].strip()
    prepare_given_exp(experiment)
    if prepareonly:
        continue
    run_pipeline_given_exp(experiment, runlevel, skipdiagnosticplots)
    """
    expdir     = exp2expdir(experiment)
    tablefolder = expdir + '/tables/'
    if prepare:
        os.chdir(expdir)
        print "\nrun prepare_astrometric_epoch.py...\n"
        os.system("prepare_astrometric_epoch.py %s.vex" % (experiment))
    if prepareonly:
        continue
    print "rerun vlbi_astrometry.py for %s" % (experiment)
    #os.chdir(codedir)
    
    print "\ndeleting sn, bp and ps tables for %s...\n" % (experiment)
    try:
        os.system('rm %s/{*sn,*bp,*ps}' % (tablefolder))
    except OSError:
        print "Some of the tables are already removed\n"
    #solutions to keep being used need to save to *.sn.save or *.bp.save. 
    solutions_to_use = glob.glob(r'%s/*.save' % (tablefolder))
    for solution in solutions_to_use:
        os.system('cp %s %s' % (solution, solution.replace('.save', '')))
        
    os.system('vlbi_astrometry.py -e %s -r %s --alwayssaved' % (experiment, runlevel))
    """
current_time=datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
print "\nSuccessfully rerun through all epochs for %s at %s (UTC)\n" % (targetname,current_time)
