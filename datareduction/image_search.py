#!/usr/bin/env ParselTongue
################################################################################
## image_search.py: A ParselTongue script for searching for inbeam calibrators
## Adam Deller, 03 Mar 2010
################################################################################

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

################################################################################
# General python imports
################################################################################
import sys, os, string, math, warnings, subprocess, signal, glob
import interaction, vlbatasks
from astroobsresult import *
from optparse import OptionParser
from interaction import yesno, setalwaysuse
warnings.defaultaction = "always"

################################################################################
# Global variables and option parsing
################################################################################
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    aipsver = '31DEC20'
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("--rootdir", dest="rootdir",
                  default="",
                  help="The base directory for the experiment series")
parser.add_option("--sourcefile", dest="sourcefile", default="",
                  help="Name of the file with all the sources listed")
parser.add_option("--seqno", dest="seqno", default="1",
                  help="Sequence number for AIPS UV datasets")
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Experiment name")
parser.add_option("-k", "--klass", dest="klass", default="UVDATA",
                  help="AIPS data class used")
parser.add_option("-r", "--runlevel", dest="runlevel", default="1",
                  help="Runlevel to start[,stop] at")
parser.add_option("--refant", dest="refant", default="5",
                  help="Reference antenna number")
parser.add_option("--ampcalscan", dest="ampcalscan", default="1",
                  help="Take this scan of the amplitude calibrator " + \
                        "for bandpass calibration etc")
parser.add_option("--ampcalmins", dest="ampcalmins", default="2",
                  help="Solint (mins) for amp cal on bandpass " + \
                       "calibrator, neg for none")
parser.add_option("--skipampcalfring", dest="skipampcalfring", \
                  default=False, action="store_true",
                  help="Do absolutely nothing with the amplitude calibrator")
parser.add_option("--limitedtargetlist", dest="limitedtargetlist", default="",
                  help="Comma sep list of targets to reduce, blank=all")
parser.add_option("--notsysfix", dest="notsysfix", default=False,
                  action="store_true", help="Don't try to fix HN amplitudes")
parser.add_option("--targetonly", dest="targetonly", default=False,
                  action="store_true", help="Reduce target only")
parser.add_option("--calonly", dest="calonly", default=False,
                  action="store_true", help="Reduce calibrator only")
parser.add_option("--skiptecor", dest="skiptecor", default=False,
                  action="store_true", help="No iono corrections")
parser.add_option("--skippcal", dest="skippcal", default=False,
                  action="store_true", help="No pcal corrections")
parser.add_option("--skippbcor", dest="skippbcor", default=False,
                action="store_true", help="No primary beam corrections")
parser.add_option("--skipaccor", dest="skipaccor", default=False,
                action="store_true", help="No ACCOR corrections")
parser.add_option("--skipeops", dest="skipeops", default=False,
                action="store_true", help="No EOP corrections")
parser.add_option("--skipmissingpbcor", dest="skipmissingpbcor", default=False,
                action="store_true", help="Skip missing sources in PBCOR step")
parser.add_option("--douvsrt", dest="douvsrt", default=False,
                  action="store_true", help="Sort into TB order")
parser.add_option("--alwaysusesaved", dest="alwaysusesaved", default=False,
                  action="store_true", help="Always used saved SN/BP tables")
parser.add_option("--fullauto", dest="fullauto", default=False,
                  action="store_true", help="Do no editing when imaging")
parser.add_option("--startlocaltv", dest="startlocaltv", default=False,
                  action="store_true", help="Start a local TV")
parser.add_option("--searchselfcalsnr", dest="searchselfcalsnr", default="8",
                  help="Initial limit for including a source in test selfcal")
parser.add_option("--searchselfcaldur", dest="searchselfcaldur", default="2",
                  help="Solution interval (mins) for search selfcal")
parser.add_option("--searchselfcalseparateifs", dest="searchselfcalsumifs", default=True,
                  action="store_false", help="Do a per-IF selfcal solution")
parser.add_option("--inttimesecs", dest="inttimesecs", default="4",
                  help="The integration time of the data in seconds")
parser.add_option("--phaserefsolint", dest="phaserefsolint", default="2",
                  help="Phase reference source solution interval (minutes)")
parser.add_option("--phaserefcalibmins", dest="phaserefcalibmins", default="3",
                  help="Phase ref solution interval for CALIB (minutes)")
parser.add_option("--skipphsrefcalib", dest="skipphsrefcalib", default=False,
                  action="store_true", help="No CALIB on phase reference sources")
parser.add_option("--dophsrefampcalib", dest="dophsrefampcalib", default=False,
                  action="store_true", help="Solve for amplitude as well during " + \
                                            "CALIB on phase reference sources")
parser.add_option("--skiptarget0", dest="skiptarget0", default=False,
                  action="store_true", help="Don't try to process file 0 for targets")
parser.add_option("--skipsnedt", dest="skipsnedt", default=False, action="store_true",
                  help="Don't SNEDT, just make plots to view later")
parser.add_option("--dosinglecalib", dest="dosinglecalib", default=False, 
                  action="store_true", help="Sum IFs for phase ref CALIB")
parser.add_option("--calibsoltype", dest="calibsoltype", default="",
                  help="Sol type for CALIB eg L1R, R, blank...")
parser.add_option("--clearcatalog", dest="clearcatalog", default=False,
                  action="store_true", help="Clear catalog before starting")
parser.add_option("--zapconfusing", dest="zapconfusing", default=False,
                  action="store_true", help="Flag potentially confused baselines")
parser.add_option("--skipuntil", dest="skipuntil", default="",
                  help="Don't start imaging until you reach this source")
parser.add_option("--ntargetsearchpix", dest="ntargetsearchpix", default="8192",
                  help="The number of pixels in the search image")
parser.add_option("--targetsearchpixmas", dest="targetsearchpixmas", default="1.5",
                  help="The size of the pixels in the search image (mas)")
parser.add_option("--noconcat", dest="noconcat", default=False, 
                  action="store_true", 
                  help="Don't concatenate sources with same names")
parser.add_option("--elevationflag", dest="elevationflag", default="12",
                  help="Elevation flag limit (set to -99 to turn off)")
parser.add_option("--writefullrestargets", dest="writefullrestargets",
                  default=False, action="store_true",
                  help="Write out the targets with no freq averaging")
parser.add_option("--noimage", dest="noimage", default=False, action="store_true",
                  help="Don't actually do any imaging (want to concatenate 2 exps")
parser.add_option("--norenumber", dest="norenumber", default=False, action="store_true",
                  help="Don't try to renumber the antennas (set true for non-VLBA)")
parser.add_option("--combineobs", dest="combineobs", default="", 
                  help="Load and DBCON full-res data from these experiments" + \
                  "(e.g., --combineobs=bd161aa,bd161ca)")
parser.add_option("--userno", dest="userno", default="2575", 
                  help="AIPS userno to use")
parser.add_option("--finalavtime", dest="finalavtime", default="0",
                  help="Averaging time in seconds for final target datasets")
parser.add_option("--targetfinalcc", dest="targetfinalcc", default="400",
                  help="Max number of clean components in final target image")
parser.add_option("--targetfinalflux", dest="targetfinalflux", default="0.5", 
                  help="Min flux to keep cleaning in target images (mJy)")
parser.add_option("--targettaper", dest="targettaper", default="18",
                  help="Taper in Mlambda for final target images")
(options, junk) = parser.parse_args()
AIPS.userno     = int(options.userno)
sourcefile      = options.sourcefile
experiment      = options.experiment.upper()
klass           = options.klass.upper()
uvsequence      = int(options.seqno)
ampcalscan      = int(options.ampcalscan)
ampcalmins      = float(options.ampcalmins)
skipampcalfring = options.skipampcalfring
rootdir         = options.rootdir
clversion       = 1
snversion       = 1
limitedtargetlist = options.limitedtargetlist.split(',')
notsysfix       = options.notsysfix
runsplit        = options.runlevel.split(',')
runfromlevel    = int(runsplit[0])
runtolevel      = 999
ntargetsearchpix= int(options.ntargetsearchpix)
targetsearchpixmas= float(options.targetsearchpixmas)
refant          = int(options.refant)
skipuntil       = options.skipuntil
targetonly      = options.targetonly
calonly         = options.calonly
zapconfusing    = options.zapconfusing
skiptecor       = options.skiptecor
skippcal        = options.skippcal
skippbcor       = options.skippbcor
skipsnedt       = options.skipsnedt
skipaccor       = options.skipaccor
skipeops        = options.skipeops
skipmissingpbcor= options.skipmissingpbcor
skipphsrefcalib = options.skipphsrefcalib
dophsrefampcalib= options.dophsrefampcalib
no_renumber     = options.norenumber
fullauto        = options.fullauto
startlocaltv    = options.startlocaltv
phaserefsolint  = float(options.phaserefsolint)
phsrefcalibmins = float(options.phaserefcalibmins)
dosinglecalib   = options.dosinglecalib
inttimesecs     = float(options.inttimesecs)
clearcatalog    = options.clearcatalog
douvsrt         = options.douvsrt
doconcat        = not options.noconcat
elevationflag   = int(options.elevationflag)
skiptarget0     = options.skiptarget0
writefullrestargets=options.writefullrestargets
noimage         = options.noimage
combineobs      = options.combineobs.split(",")
searchselfcalsnr= float(options.searchselfcalsnr)
searchselfcaldur= float(options.searchselfcaldur)
searchselfcalsumifs = options.searchselfcalsumifs
finalavtime     = float(options.finalavtime)
targetfinalcc   = float(options.targetfinalcc)
targetfinalflux = float(options.targetfinalflux)/1000.0
targettaper     = float(options.targettaper)
clearlyinsnr    = 6.5
clearlyinflux   = 5 # mJy
UNFLAGGED_THRESHOLD = 0.33
setalwaysuse(options.alwaysusesaved)
if len(runsplit) > 1:
    runtolevel = int(runsplit[1])

################################################################################
# Check validity of inputs
################################################################################
if rootdir == "":
    parser.error("You must supply a root directory!")
if experiment == "":
    parser.error("You must supply an experiment name")
if runfromlevel > runtolevel:
    parser.error("Runtolevel (" + str(runtolevel) + ") must be greater " \
                 "than or equal to runfromlevel (" + str(runfromlevel) + ")")
directory   = rootdir + '/' + experiment.lower()
calmodeldir = rootdir + '/models/'
solutiondir = rootdir + '/solutions/'
tabledir    = directory + "/tables/"
logdir      = directory + "/logs/"
vexfile     = directory + '/' + experiment.lower() + '.vex'
if not os.path.exists(solutiondir):
    os.mkdir(solutiondir)
if not os.path.exists(tabledir):
    os.mkdir(tabledir)
if not os.path.exists(logdir):
    parser.error("The log directory " + logdir + " does not exist!")
if not os.path.exists(vexfile):
    parser.error("The vex file " + vexfile + " does not exist!")
if len(sourcefile) == 0:
    parser.error("You must supply a source file!")
if not sourcefile[0] == '/':
    sourcefile = directory + '/' + sourcefile
if skipuntil != "" and runfromlevel != runtolevel:
    parser.error("skipuntil can only be used with imaging step - set runlevel appropriately")
if len(combineobs) > 1 and (runfromlevel != 20 or runtolevel != 20):
    parser.error("combineobs can only be used with runlevel=20,20")

################################################################################
# Parse the source file
################################################################################
if not os.path.exists(sourcefile):
    parser.error("Source file " + sourcefile + " does not exist!")

sourcein = open(sourcefile)
sourcelines = sourcein.readlines()
sourcein.close()
atline = 0
targetuvfiles = []

if atline >= len(sourcelines):
    print "Problem parsing source file " + sourcefile + " - aborting!"
    sys.exit()
keyval = sourcelines[atline].split(':')
if len(keyval) == 2 and keyval[0] == "CALIBRATOR FILE":
    calibsuvfile = keyval[1].strip()
    splitcalibsfiles = calibsuvfile.split(',')
    count = 1
    calibsuvfile = ""
    for f in splitcalibsfiles:
        csplit = f.rsplit('/',1)
        newfile = csplit[0] + "/srccal-" + str(count) + ".fitsidi"
        os.system("rm -f " + newfile)
        os.system("ln -s " + f  + " " + newfile)
        count += 1
        calibsuvfile = calibsuvfile + newfile + ":"
    calibsuvfile = calibsuvfile[:-1]
else:
    print "Problem parsing source file " + sourcefile + " - bad cal file - aborting!"
    sys.exit()
atline += 1

if atline >= len(sourcelines):
    print "Problem parsing source file " + sourcefile + " - aborting!"
    sys.exit()
keyval = sourcelines[atline].split(':')
if len(keyval) == 2 and keyval[0] == "AMP CAL SOURCE":
    ampcalsrc = keyval[1].strip()
else:
    print "Problem parsing source file " + sourcefile + " - bad ampcal source - aborting!"
    sys.exit()
atline += 1

if atline >= len(sourcelines):
    print "Problem parsing source file " + sourcefile + " - aborting!"
    sys.exit()
keyval = sourcelines[atline].split(':')
if len(keyval) == 2 and keyval[0] == "MAX PHASE CENTRES":
    maxphasecentres = int(keyval[1].strip())
else:
    print "Problem parsing source file " + sourcefile + " - bad maxphscentres - aborting!"
    sys.exit()
atline += 1

for i in range(maxphasecentres):
    if atline >= len(sourcelines):
        print "Problem parsing source file " + sourcefile + " - aborting!"
        sys.exit()
    keyval = sourcelines[atline].split(':')
    if len(keyval) == 2 and keyval[0] == ("PHASE CENTRE %d FILE" % (i)):
        targetuvfiles.append(keyval[1].strip())
        splittargetfiles = targetuvfiles[-1].split(',')
        count = 1
        targetuvfiles[-1] = ""
        for f in splittargetfiles:
            tsplit = f.rsplit('/',1)
            newfile = tsplit[0] + "/src" + str(i) + "-" + str(count) + ".fitsidi"
            os.system("rm -f " + newfile)
            os.system("ln -s " + f  + " " + newfile)
            count += 1
            targetuvfiles[-1] = targetuvfiles[-1] + newfile + ":"
        targetuvfiles[-1] = targetuvfiles[-1][:-1]
    else:
        print "Problem parsing source file " + sourcefile + " - bad phscentrefile - aborting!"
        sys.exit()
    atline += 1

if atline >= len(sourcelines):
    print "Problem parsing source file " + sourcefile + " - aborting!"
    sys.exit()
keyval = sourcelines[atline].split(':')
if len(keyval) == 2 and keyval[0] == "NUM TARGETS":
    numtargets = int(keyval[1].strip())
else:
    print "Problem parsing source file " + sourcefile + " - bad numtargets - aborting!"
    sys.exit()
atline += 1

phscalsrc = []
numfields = []
numfieldsources = []
fieldsourcenames = []
pointingcentrenames = []
for i in range(numtargets):
    if atline >= len(sourcelines):
        print "Problem parsing source file " + sourcefile + " - aborting!"
        sys.exit()
    keyval = sourcelines[atline].split(':')
    if len(keyval) == 2 and keyval[0] == ("TARGET %d PHASE REF"%(i)):
        phscalsrc.append(keyval[1].strip())
    else:
        print "Problem parsing source file " + sourcefile + " - bad phaseref name - aborting!"
        sys.exit()
    atline += 1
    if atline >= len(sourcelines):
        print "Problem parsing source file " + sourcefile + " - aborting!"
        sys.exit()
    keyval = sourcelines[atline].split(':')
    if len(keyval) == 2 and keyval[0] == ("TARGET %d NUM FIELDS"%(i)):
        numfields.append(int(keyval[1].strip()))
        numfieldsources.append([])
        fieldsourcenames.append([])
        pointingcentrenames.append([])
    else:
        print "Problem parsing source file " + sourcefile + " - bad numfields - aborting!"
        sys.exit()
    atline += 1
    for j in range(numfields[i]):
        if atline >= len(sourcelines):
            print "Problem parsing source file " + sourcefile + " - aborting!"
            sys.exit()
        keyval = sourcelines[atline].split(':')
        if len(keyval) == 2 and keyval[0] == ("TARGET %d FIELD %d NUM SOURCES"%(i,j)):
            numfieldsources[i].append(int(keyval[1].strip()))
            fieldsourcenames[i].append([])
        else:
            print "Problem parsing source file " + sourcefile + " - bad numsources - aborting!"
            sys.exit()
        atline += 1
        for k in range(numfieldsources[i][j]):
            if atline >= len(sourcelines):
                print "Problem parsing source file " + sourcefile + " (hit the end of the file): aborting!"
                sys.exit()
            keyval = sourcelines[atline].split(':')
            if len(keyval) == 2 and keyval[0] == ("TARGET %d FIELD %d SOURCE %d NAME"%(i,j,k)):
                fieldsourcenames[i][j].append(keyval[1].strip()[0:16])
                if k==0:
                    pointingcentrenames[i].append('%s-%d' % (fieldsourcenames[i][j][0], j+1))
            else:
                print "Problem parsing source file " + sourcefile + " - bad targetname - aborting!"
                print "Expected something beginning with TARGET %d FIELD %d SOURCE %d NAME"%(i,j,k)
                print "Got %s" % sourcelines[atline]
                sys.exit()
            atline += 1
print phscalsrc

if not targetonly and runfromlevel == 1 and not calibsuvfile[0] == '/':
    calibsuvfile = directory + '/' + calibsuvfile
for i in range(maxphasecentres):
    if len(targetuvfiles[i]) > 0 and runfromlevel == 1 and not targetuvfiles[i][0] == '/':
        targetuvfiles[i] = directory + '/' + targetuvfiles[i]
if not targetonly:
    calibsuvdata = AIPSUVData(experiment + "_C", klass, 1, uvsequence)
else:
    calibsuvdata = None
targetuvdata = []
for i in range(maxphasecentres):
    targetuvdata.append(AIPSUVData(experiment + "_T" + str(i), klass, 1, uvsequence))

if startlocaltv:
    os.system("rm -f junkjunk.txt")
    os.system("ps -A | grep XAS > junkjunk.txt")
    xaslines = open("junkjunk.txt").readlines()
    os.system("rm -f junkjunk.txt")
    for line in xaslines:
        if len(line) > 0:
            pid = int(line.split()[0])
            os.kill(pid, signal.SIGKILL)
    #tv = AIPSTV("local")
    tv = AIPSTV()
    tv.start()

runlevel = 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Load the uv data ############################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Loading UVDATA from disk"
    if not targetonly:
        if options.clearcatalog and calibsuvdata.exists():
            calibsuvdata.zap()
        print calibsuvfile
        vlbatasks.fitld_vlba(calibsuvfile, calibsuvdata, [])
    if not calonly:
        for i in range(maxphasecentres):
            if skiptarget0 and i == 0: continue
            if options.clearcatalog and targetuvdata[i].exists():
                targetuvdata[i].zap()
            vlbatasks.fitld_vlba(targetuvfiles[i], targetuvdata[i], [])
else:
    print "Skipping UVDATA load"

## Run uvsrt if necessary ######################################################
runlevel += 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
if runfromlevel <= runlevel and runtolevel >= runlevel and douvsrt and \
   not calonly:
    print "Runlevel " + str(runlevel) + ": Running uvsrt on target file only"
    sortedtargetuvdata = []
    tosub = 0
    if skiptarget0:
        tosub = 1
    for i in range(maxphasecentres):
        if skiptarget0 and i == 0: continue
        sortedtargetuvdata.append(AIPSUVData(experiment + "_T" + str(i), klass, 
                                             1, uvsequence+1))
        if sortedtargetuvdata[-1].exists():
            if not clearcatalog and \
               not yesno("Delete existing UVSRT output data file (no aborts pipeline)?"):
                sys.exit()
            sortedtargetuvdata[-1].zap()
        vlbatasks.uvsrt(targetuvdata[i], sortedtargetuvdata[i-tosub])
        if not no_renumber:
            vlbatasks.renumber_vlba(sortedtargetuvdata[i-tosub])
        targetuvdata[i].zap()
    targetuvdata = sortedtargetuvdata
elif douvsrt: #Still need to update the names
    targetuvdata = []
    for i in range(maxphasecentres):
        if skiptarget0 and i == 0: continue
        unsorteddata = AIPSUVData(experiment + "_T" + str(i), klass, 1, uvsequence)
        if unsorteddata.exists():
            unsorteddata.zap()
        targetuvdata.append(AIPSUVData(experiment + "_T" + str(i), klass, 1, 
                                       uvsequence+1))

########## HOPEFULLY THIS IS NOT NECESSARY EVER ################################
#runlevel = runlevel + 1
## Run VLBAFIX #################################################################
#if runfromlevel <= runlevel and runtolevel >= runlevel:
#    vlbatasks.vlbafix(ampcaluvdata)
#    vlbatasks.vlbafix(phscaluvdata)
#    vlbatasks.vlbafix(inbeamuvdata)
#    vlbatasks.vlbafix(targetuvdata)
################################################################################

if skiptarget0:
    maxphasecentres -= 1

runlevel = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Load the user flags (if any) ################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Loading user flags"
    userflagfile = tabledir + "additionaledit.flag"
    if os.path.exists(userflagfile):
        if not targetonly:
            vlbatasks.userflag(calibsuvdata, 1, userflagfile)
        if not calonly:
            for i in range(maxphasecentres):
                vlbatasks.userflag(targetuvdata[i], 1, userflagfile)
    else:
        print "No user flag file - skipping"
else:
    print "Skipping user flagging"

runlevel = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Flag on elevation unless told not to #######################################
if runfromlevel <= runlevel and runtolevel >= runlevel and elevationflag > 0:
    print "Runlevel " + str(runlevel) + ": Flagging on elevation at " + \
          str(elevationflag) + " degrees"
    if not targetonly:
        vlbatasks.elevationflag(calibsuvdata, 1, elevationflag)
    if not calonly: 
        for i in range(maxphasecentres):
            vlbatasks.elevationflag(targetuvdata[i], 1, elevationflag)
else:
    print "Skipping elevation based flagging"

runlevel = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run TECOR ###################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skiptecor:
    print "Runlevel " + str(runlevel) + ": Running TECOR to correct ionosphere"
    if not targetonly:
        vlbatasks.correct_iono(calibsuvdata, logdir, clversion)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.correct_iono(targetuvdata[i], logdir, clversion)
else:
    print "Skipping ionospheric corrections"

runlevel  = runlevel + 1
if not skiptecor:
    clversion = clversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run CLCOR to correct EOPs ###################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipeops:
    print "Runlevel " + str(runlevel) + ": Running CLCOR to correct for EOPs"
    if not targetonly:
        vlbatasks.correct_eops(calibsuvdata, logdir, clversion)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.correct_eops(targetuvdata[i], logdir, clversion)
else:
    print "Skipping EOP corrections"

runlevel  = runlevel + 1
if not skipeops:
    clversion = clversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Do the amplitude calibration (either load existing table or do and then edit) and inspect
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Doing amplitude calibration"
    if targetonly and ((not skipaccor and not os.path.exists(tabledir + 'accor.sn')) or \
                       not os.path.exists(tabledir + 'apcal.sn')):
        print "For target-only, the SN files must exist already - aborting!"
        sys.exit(1)
    uvdatas = []
    if not targetonly:
        uvdatas.append(calibsuvdata)
    if not calonly:
        for uvdata in targetuvdata:
            uvdatas.append(uvdata)
    if not targetonly and (((not skipaccor and not os.path.exists(tabledir + 'accor.sn')) and not\
        os.path.exists(tabledir + 'apcal.sn')) or not \
       yesno("Do you wish to used saved SN table of amplitude calibration?")):
        for i in range(7):
            vlbatasks.deletetable(calibsuvdata, 'SN', i+snversion)
        dotsysfix = False
        for ant in calibsuvdata.antennas:
            if ant == 'HN' and not notsysfix:
                dotsysfix = True
        snplus = 0
        if not skipaccor:
            vlbatasks.accor(calibsuvdata)
            vlbatasks.snsmo(calibsuvdata, 'AMP', 10, 0.04, 0.0, 0.0, snversion, refant)
            snplus = 2
        vlbatasks.apcal(calibsuvdata, snversion+snplus)
        vlbatasks.snsmo(calibsuvdata, 'AMP', 10, 1.50, 0.0, 0.0, snversion+snplus, refant)
        if dotsysfix:
            vlbatasks.fixtsys(calibsuvdata, snversion+snplus+1)
        else:
            vlbatasks.tacop(calibsuvdata, 'SN', snversion+snplus+1, calibsuvdata, snversion+snplus+2)
        snoutver1 = snversion + 1
        snoutver2 = snversion + 4
        if not skipsnedt:
            vlbatasks.snedt(calibsuvdata, snversion+1)
            vlbatasks.snedt(calibsuvdata, snversion+4)
            snoutver1 = snversion+5
            snoutver2 = snversion+6
        if skipaccor:
            snoutver2 -= 2
        else:
            if os.path.exists(tabledir + 'accor.sn'):
                os.remove(tabledir + 'accor.sn')
        if os.path.exists(tabledir + 'apcal.sn'):
            os.remove(tabledir + 'apcal.sn')
        if not skipaccor:
            vlbatasks.writetable(calibsuvdata, 'SN', snoutver1, tabledir + 'accor.sn')
            vlbatasks.plottops(calibsuvdata, 'SN', snoutver1, 'AMP', 0, 2, 4, tabledir + 'accor.ps')
        vlbatasks.writetable(calibsuvdata, 'SN', snoutver2, tabledir + 'apcal.sn')
        vlbatasks.plottops(calibsuvdata, 'SN', snoutver2, 'AMP', 0, 2, 4, tabledir + 'apcal.ps')
        for i in range(7):
            vlbatasks.deletetable(calibsuvdata, 'SN', i+snversion)
    for uvdata in uvdatas:
        if skipaccor:
            vlbatasks.loadtable(uvdata, tabledir + 'apcal.sn', snversion)
            vlbatasks.applysntable(uvdata, snversion, '2PT', clversion, refant)
        else:
            vlbatasks.loadtable(uvdata, tabledir + 'accor.sn', snversion)
            vlbatasks.loadtable(uvdata, tabledir + 'apcal.sn', snversion+1)
            vlbatasks.applysntable(uvdata, snversion, '2PT', clversion, refant)
            vlbatasks.applysntable(uvdata, snversion+1, '2PT', clversion+1, refant)
else:
    print "Skipping amplitude calibration"

onsourcetimes = []
runlevel  = runlevel + 1
if skipaccor:
    clversion += 1
    snversion += 1
else:
    clversion = clversion + 2
    snversion = snversion + 2
scanlist = vlbatasks.getvexscaninfo(vexfile)
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Correct for primary beam attenuation ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skippbcor:
    print "Runlevel " + str(runlevel) + ": Doing primary beam correction"
    if not targetonly:
        pbsntable = tabledir + 'pbcor.cal.sn'
        if not os.path.exists(pbsntable) or not \
               yesno("Do you wish to used saved SN table of primary beam correction?"):
            vlbatasks.deletetable(calibsuvdata, 'SN', snversion)
            ost = vlbatasks.correct_primarybeam(calibsuvdata, snversion-1, 0, scanlist, fieldsourcenames, 
                                                True, True, skiptarget0, False, skipmissingpbcor)
            onsourcetimes.append(ost)
            if os.path.exists(pbsntable):
                os.remove(pbsntable)
            vlbatasks.writetable(calibsuvdata, 'SN', snversion, pbsntable)
            vlbatasks.deletetable(calibsuvdata, 'SN', snversion)
        vlbatasks.loadtable(calibsuvdata, pbsntable, snversion)
        vlbatasks.applysntable(calibsuvdata, snversion, '2PT', clversion, refant)
    if not calonly:
        count = 0
        print "There are " + str(len(targetuvdata)) + " uv databases"
        for uvdata in targetuvdata:
            pbsntable = tabledir + 'pbcor.' + str(count) + '.sn'
            #if not os.path.exists(pbsntable) or (not targetonly and not \
            #   yesno("Do you wish to used saved SN table of primary beam correction?")):
            #    vlbatasks.deletetable(uvdata, 'SN', snversion)
            #    vlbatasks.correct_primarybeam(uvdata, snversion-1, count, scanlist, fieldsourcenames, False)
            #    if os.path.exists(pbsntable):
            #        os.remove(pbsntable)
            #    vlbatasks.writetable(uvdata, 'SN', snversion, pbsntable)
            #    vlbatasks.deletetable(uvdata, 'SN', snversion)
            vlbatasks.deletetable(uvdata, 'SN', snversion)
            ost = vlbatasks.correct_primarybeam(uvdata, snversion-1, count, scanlist, fieldsourcenames, 
                                                False, True, skiptarget0, False, skipmissingpbcor)
            onsourcetimes.append(ost)
            if os.path.exists(pbsntable):
                os.remove(pbsntable)
            vlbatasks.writetable(uvdata, 'SN', snversion, pbsntable)
            vlbatasks.deletetable(uvdata, 'SN', snversion)
            vlbatasks.loadtable(uvdata, pbsntable, snversion)
            vlbatasks.applysntable(uvdata, snversion, 'SELN', clversion, refant)
            count += 1
else:
    print "Skipping primary beam correction"
    if not calonly and not noimage and len(combineobs) <= 1:
        count = 0
        for uvdata in targetuvdata:
            ost = vlbatasks.correct_primarybeam(uvdata, snversion-1, count, 
                                                scanlist, fieldsourcenames,
                                                False, True, skiptarget0, True,
                                                skipmissingpbcor)
            onsourcetimes.append(ost)
            count += 1


runlevel  = runlevel + 1
if not skippbcor:
    clversion = clversion + 1
    snversion = snversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Do PCAL correction and inspect ##############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skippcal:
    print "Runlevel " + str(runlevel) + ": Doing pulse cal calibration"
    if targetonly and (not os.path.exists(tabledir + 'pccor.sn')):
        print "For target-only, the PC file must exist already - aborting!"
        sys.exit(1)
    if  not (targetonly or (os.path.exists(tabledir + 'pccor.sn') and \
                            yesno("Do you wish to used saved SN table " + \
                                  "for pulse cal?"))):
        try:
            calibsuvdata.table('SN', snversion).zap()
        except IOError:
            print "No need to delete old SN table" 
        vlbatasks.pccor(calibsuvdata, ampcalsrc, snversion, ampcalscan, refant)
        if os.path.exists(tabledir + 'pccor.sn'):
            os.remove(tabledir + 'pccor.sn')
        vlbatasks.writetable(calibsuvdata, 'SN', snversion, tabledir + 'pccor.sn')
        calibsuvdata.table('SN', snversion).zap()
    if not targetonly:
        vlbatasks.loadtable(calibsuvdata, tabledir + 'pccor.sn', snversion)
        vlbatasks.applysntable(calibsuvdata, snversion, '2PT', clversion, refant)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.loadtable(targetuvdata[i], tabledir + 'pccor.sn', snversion)
            vlbatasks.applysntable(targetuvdata[i], snversion, '2PT', clversion, refant)
else:
    print "Skipping pulse cal corrections"

runlevel  = runlevel + 1
if not skippcal:
    snversion = snversion + 1
    clversion = clversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run FRING (bandpass calibrator) #############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    not skipampcalfring:
    if not os.path.exists(tabledir + 'ampcalfring.sn') or \
           not yesno("Do you wish to used saved SN table for ampcal FRING?"):
        print "Runlevel " + str(runlevel) + ": FRING'ing amplitude calibrator"
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+1)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+2)
        bandpassfringmins = 3
        vlbatasks.fring(calibsuvdata, snversion, clversion, bandpassfringmins, inttimesecs, ampcalsrc, 
                        refant)
        print "Runlevel " + str(runlevel) + ": editing FRING for ampcal calibrator"
        vlbatasks.snsmo(calibsuvdata, 'DELA', 20, 0.0, 0.0, 10.0, snversion, 
                        refant)
        snoutver = snversion+1
        if not skipsnedt:
            vlbatasks.snedt(calibsuvdata, snversion+1)
            snoutver = snversion+2
        if os.path.exists(tabledir + 'ampcalfring.sn'):
            os.remove(tabledir + 'ampcalfring.sn')
        vlbatasks.writetable(calibsuvdata, 'SN', snoutver, tabledir + \
                             'ampcalfring.sn')
        vlbatasks.plottops(calibsuvdata, 'SN', snoutver, 'DELA', 0, 2, 4, 
                           tabledir + 'ampcalfring.ps')
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+1)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+2)
else:
    print "Skipping FRING and SNSMO/SNEDT of ampcal FRING results"

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Copy the FRING SN table around and apply it #################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipampcalfring:
    print "Runlevel " + str(runlevel) + ": Loading ampcal FRING SN table & calibrating"
    if targetonly and (not os.path.exists(tabledir + 'ampcalfring.sn')):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly:
        vlbatasks.loadtable(calibsuvdata, tabledir + 'ampcalfring.sn', snversion)
        vlbatasks.applysntable(calibsuvdata, snversion, '2PT', clversion, refant)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.loadtable(targetuvdata[i], tabledir + 'ampcalfring.sn', snversion)
            vlbatasks.applysntable(targetuvdata[i], snversion, '2PT', clversion, refant)
else:
    print "Skipping calibration of FRING results"

runlevel  = runlevel + 1
if not skipampcalfring:
    snversion = snversion + 1
    clversion = clversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run BPASS ###################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    not skipampcalfring:
    print "Runlevel " + str(runlevel) + ": Generating bandpass corrections"
    if not os.path.exists(tabledir + 'bpass.bp') or not \
       yesno("Do you wish to used saved BP table for bandpass?"):
        vlbatasks.bpass(calibsuvdata, ampcalsrc, clversion, ampcalscan)
        if os.path.exists(tabledir + 'bpass.bp'):
            os.remove(tabledir + 'bpass.bp')
        vlbatasks.writetable(calibsuvdata, 'BP', 1, tabledir + 'bpass.bp')
        vlbatasks.deletetable(calibsuvdata, 'BP', 1)
else:
    print "Skipping bandpass corrections"

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Load BPASS ##################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipampcalfring:
    print "Runlevel " + str(runlevel) + ": Loading bandpass corrections"
    if not os.path.exists(tabledir + 'bpass.bp'):
        print "Error - bandpass table " + tabledir + 'bpass.bp does not exist'
        sys.exit(1)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.loadtable(targetuvdata[i], tabledir + 'bpass.bp', 1)
    if not targetonly:
        vlbatasks.loadtable(calibsuvdata, tabledir + 'bpass.bp', 1)
else:
    print "Skipping loading of bandpass corrections"

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Split and run CALIB on the bandpass calibrator ##############################
if runfromlevel <= runlevel and runtolevel >= runlevel and ampcalmins > 0.0 and \
    not skipampcalfring:
    print "Splitting bandpass calibrator and running amplitude CALIB"
    if ampcalsrc in phscalsrc:
        calibtablepath = tabledir + 'ampcalsrc.calib.sn'
    else:
        calibtablepath = tabledir + ampcalsrc + '.calib.sn'
    if targetonly and not os.path.exists(calibtablepath):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly and (not os.path.exists(calibtablepath) or not \
       yesno("Do you wish to used saved SN table for amp cal CALIB?")):
        ampcalsplitdata = AIPSUVData(ampcalsrc, 'CALIB', 1, 1)
        calibsoltype = 'L1R'
        calibsnr = 5
        if ampcalsplitdata.exists():
            ampcalsplitdata.zap()
        vlbatasks.split(calibsuvdata, clversion, 'CALIB', ampcalsrc)
        ampcalsplitdata.table('NX', 1).zap()
        ampcalmodeldata = None
        ampcalmodelfile = calmodeldir + ampcalsrc + '.clean.fits'
        if os.path.exists(ampcalmodelfile):
            ampcalmodeldata = AIPSImage(ampcalsrc, 'CLNMOD', 1, 1)
            if ampcalmodeldata.exists():
                ampcalmodeldata.zap()
            vlbatasks.fitld_image(ampcalmodelfile, ampcalmodeldata)
        else:
            s = raw_input("You have no model of amplitude calibrator " + ampcalsrc + \
                          " - this could lead to very bad amp cal solutions! " + \
                          "Type CONTINUE to continue, otherwise I'll split it off for you...")
            if not s == "CONTINUE":
                ampcalfile = directory + '/' + ampcalsrc + "_formodel_uv.fits"
                vlbatasks.writedata(ampcalsplitdata, ampcalfile, True)
                sys.exit()
        vlbatasks.singlesource_calib(ampcalsplitdata, ampcalmodeldata, 1, refant, True, ampcalmins,
                                     False, calibsoltype, calibsnr, False)
        vlbatasks.snsmo(ampcalsplitdata, 'BOTH', 20, 0.0, 30.0, 0.0, 1, refant)
        if os.path.exists(calibtablepath):
            os.remove(calibtablepath)
        snoutver = 2
        if not skipsnedt:
            vlbatasks.snedt(ampcalsplitdata, 2)
            snoutver += 1
        vlbatasks.writetable(ampcalsplitdata, 'SN', snoutver, calibtablepath)
        vlbatasks.plottops(calibsuvdata, 'SN', snoutver, 'DELA', 0, 2, 4, 
                           tabledir + 'ampcalcalib.ps')
        if not os.path.exists(ampcalmodelfile):
            renorm = vlbatasks.norm_snamplitudes(calibtablepath)
            ampcalsplitdata.zap()
            vlbatasks.split(calibsuvdata, clversion, 'CALIB', ampcalsrc)
            ampcalsplitdata.table('NX', 1).zap()
            vlbatasks.singlesource_calib(ampcalsplitdata, renorm, 1, refant, True, ampcalmins,
                                         False, calibsoltype, calibsnr, False)
            vlbatasks.snsmo(ampcalsplitdata, 'BOTH', 20, 0.0, 30.0, 0.0, 1, refant)
            os.remove(calibtablepath)
            vlbatasks.writetable(ampcalsplitdata, 'SN', 2, calibtablepath)
else:
    print "Skipping amp cal CALIB"

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Load CALIB solutions from the bandpass calibrator ###########################
if runfromlevel <= runlevel and runtolevel >= runlevel and ampcalmins > 0.0 and \
    not skipampcalfring:
    print "Runlevel " + str(runlevel) + ": Loading CALIB solutions from bandpass"
    if ampcalsrc in phscalsrc:
        calibtablepath = tabledir + 'ampcalsrc.calib.sn'
    else:
        calibtablepath = tabledir + ampcalsrc + '.calib.sn'
    if not os.path.exists(calibtablepath):
        print "Error - bandpass cal amplitude calibration table " + \
              calibtablepath + ' does not exist'
        sys.exit(1)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.loadtable(targetuvdata[i], calibtablepath, snversion)
            vlbatasks.applysntable(targetuvdata[i], snversion, '2PT', clversion, refant)
    if not targetonly:
        vlbatasks.loadtable(calibsuvdata, calibtablepath, snversion)
        vlbatasks.applysntable(calibsuvdata, snversion, '2PT', clversion, refant)
else:
    print "Skipping loading of CALIB solutions from bandpass"

runlevel  = runlevel + 1
doband = True
if ampcalmins > 0.0 and not skipampcalfring:
    snversion = snversion + 1
    clversion = clversion + 1
if skipampcalfring:
    doband = False
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run FRING (phase reference calibrator) ######################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    phaserefsolint > 0:
    if not os.path.exists(tabledir + 'phsreffring.sn') or \
           not yesno("Do you wish to used saved SN table for phsref FRING?"):
        print "Runlevel " + str(runlevel) + ": FRING'ing phase calibrators"
        try:
            calibsuvdata.table('SN', snversion).zap()
            calibsuvdata.table('SN', snversion+1).zap()
            calibsuvdata.table('SN', snversion+2).zap()
        except IOError:
            print "No need to delete old SN table(s)"
        phscalnames = []
        for phscal in phscalsrc:
            if phscal in phscalnames:
                continue
            phscalmodeldata = None
            phscalmodelfile = calmodeldir + phscal + '.clean.fits'
            if os.path.exists(phscalmodelfile):
                print "Using final model for " + phscal
                phscalmodeldata = AIPSImage(phscal, 'CLNMOD', 1, 1)
                if phscalmodeldata.exists():
                    phscalmodeldata.zap()
                vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
            else:
                phscalmodelfile = calmodeldir + 'initial/' + phscal + '.clean.fits'
                if os.path.exists(phscalmodelfile):
                    print "Using initial model for " + phscal
                    phscalmodeldata = AIPSImage(phscal, 'CLNMOD', 1, 1)
                    if phscalmodeldata.exists():
                        phscalmodeldata.zap()
                    vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
                else:
                    print "Currently no model for " + phscal

            phscalnames.append(phscal)
            print "About to fring source " + phscal
            try:
                doband = True
                if skipampcalfring:
                    doband = False
                fringsnr = 4.5
                vlbatasks.fring(calibsuvdata, snversion, clversion, phaserefsolint, 
                                inttimesecs, phscal, refant, doband, fringsnr, False, phscalmodeldata)
            except RuntimeError:
                if not (raw_input("Problem with source %s - enter 'skip' to continue, anything else will abort" % (phscal)) == "skip"):
                    sys.exit()
        #vlbatasks.snsmo(calibsuvdata, 'VLBI', 60, 0.0, 0.0, 6.0, snversion,
        #                refant, True, 4.0)
        if skiptarget0: # Proxy for mJIVE observations, where we have frequent sampling
            vlbatasks.snmwfclip(calibsuvdata, 12, 0.0, 0.0, 10.0, 16.0, snversion, refant)
        else:
            vlbatasks.snmwfclip(calibsuvdata, 30, 0.0, 0.0, 16.0, 16.0, snversion, refant)
        snoutver = snversion+1
        if not skipsnedt:
            vlbatasks.snedt(calibsuvdata, snversion+1)
            snoutver += 1
        if os.path.exists(tabledir + 'phsreffring.sn'):
            os.remove(tabledir + 'phsreffring.sn')
        vlbatasks.writetable(calibsuvdata, 'SN', snoutver, tabledir + \
                             'phsreffring.sn')
        vlbatasks.plottops(calibsuvdata, 'SN', snoutver, 'DELA', 0, 2, 4,
                           tabledir + 'phsreffring.delay.ps')
        vlbatasks.plottops(calibsuvdata, 'SN', snoutver, 'RATE', 0, 2, 4,
                           tabledir + 'phsreffring.delay.ps')
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+1)
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+2)
else: 
    print "Skipping FRING and SNSMO/SNEDT of FRING results" 

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Copy the phsref FRING SN table around and apply it ##########################
if runfromlevel <= runlevel and runtolevel >= runlevel and phaserefsolint > 0:
    print "Runlevel " + str(runlevel) + ": Loading phsref FRING SN table & calibrating"
    if targetonly and (not os.path.exists(tabledir + 'phsreffring.sn')):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly:
        try:
            calibsuvdata.table('CL', clversion+1).zap()
        except IOError:
            print "No need to delete old CL table"
        try:
            calibsuvdata.table('SN', snversion).zap()
        except IOError:
            print "No need to delete old SN table(s)"
        vlbatasks.loadtable(calibsuvdata, tabledir + 'phsreffring.sn', snversion)
        vlbatasks.applysntable(calibsuvdata, snversion, '2PT', clversion, refant)
    if not calonly:
        for i in range(maxphasecentres):
            try:
                targetuvdata[i].table('CL', clversion+1).zap()
            except IOError:
                print "No need to delete old CL table"
            try:
                targetuvdata[i].table('SN', snversion).zap()
            except IOError:
                print "No need to delete old SN table(s)"
            vlbatasks.loadtable(targetuvdata[i], tabledir + 'phsreffring.sn', snversion)
            vlbatasks.applysntable(targetuvdata[i], snversion, '2PT', clversion, refant)
else:
    print "Skipping calibration of FRING results"

runlevel  = runlevel + 1
if phaserefsolint > 0:
    snversion = snversion + 1
    clversion = clversion + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Run CALIB on the phase reference sources ####################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipphsrefcalib:
    if not targetonly and (not os.path.exists(tabledir + 'phsrefcalib.sn') or \
           not yesno("Do you wish to used saved SN table for phsref CALIB?")):
        print "Runlevel " + str(runlevel) + ": Running CALIB on phsref sources"
        phscalnames = []
        for phscal in phscalsrc:
            if phscal in phscalnames:
                continue
            phscalnames.append(phscal)
            for i in range(100): #Clear any old CALIB split catalog entries
                phscal_uv_data = AIPSUVData(phscal, 'CALIB', 1, i)
                if phscal_uv_data.exists():
                    phscal_uv_data.zap()
            phscal_uv_data = AIPSUVData(phscal, 'CALIB', 1, 1)
            doband = True
            if skipampcalfring:
                doband = False
            if doband:
                vlbatasks.split(calibsuvdata, clversion, 'CALIB', phscal)
            else:
                vlbatasks.splittoseqnoband(calibsuvdata, clversion, 'CALIB', phscal, 1)
            phscal_uv_data.table('NX', 1).zap()
            phscalmodeldata = None
            phscalmodelfile = calmodeldir + phscal + '.clean.fits'
            if os.path.exists(phscalmodelfile):
                print "Using final model for " + phscal
                phscalmodeldata = AIPSImage(phscal, 'CLNMOD', 1, 1)
                if phscalmodeldata.exists():
                    phscalmodeldata.zap()
                vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
            else:
                phscalmodelfile = calmodeldir + 'initial/' + phscal + '.clean.fits'
                if os.path.exists(phscalmodelfile):
                    print "Using initial model for " + phscal
                    phscalmodeldata = AIPSImage(phscal, 'CLNMOD', 1, 1)
                    if phscalmodeldata.exists():
                        phscalmodeldata.zap()
                    vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
                else:
                    print "Currently no model for " + phscal
            calibsoltype = 'L1R'
            calibtablepath = tabledir + phscal + '.calib.sn'
            calibsnr = 4.25
            vlbatasks.singlesource_calib(phscal_uv_data, phscalmodeldata,
                                         1, refant, dophsrefampcalib, phsrefcalibmins,
                                         dosinglecalib, calibsoltype, calibsnr, 
                                         dosinglecalib)
            #vlbatasks.snsmo(phscal_uv_data, 'BOTH', 20, 0.0,
            #                120.0, 0.0, 1, refant)
            if os.path.exists(calibtablepath):
                os.remove(calibtablepath)
            #vlbatasks.writetable(phscal_uv_data, 'SN', 2,
            #                     calibtablepath)
            plottype = "PHAS"
            if dophsrefampcalib:
                plottype = "AMP"
            vlbatasks.writetable(phscal_uv_data, 'SN', 1, calibtablepath)
            vlbatasks.plottops(phscal_uv_data, 'SN', 1, plottype, 0, 2, 4, 
                               tabledir + phscal + '.calib.ps')
            phscal_uv_data.zap()
else:
    print "Skipping CALIB on the phase reference source"

runlevel  = runlevel + 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Load all the CALIB solutions ###############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipphsrefcalib:
    print "Runlevel " + str(runlevel) + ": Loading CALIB solutions and applying"
    sncount = 0
    for phscal in phscalsrc:
        if not targetonly:
            vlbatasks.deletetable(calibsuvdata, 'SN', snversion+sncount)
        if not calonly:
            for i in range(maxphasecentres):
                vlbatasks.deletetable(targetuvdata[i], 'SN', snversion+sncount)
        sncount += 1
    sncount = 0
    phscalnames = []
    for phscal in phscalsrc:
        if phscal in phscalnames:
            continue
        phscalnames.append(phscal)
        calibtablepath = tabledir + phscal + '.calib.sn'
        if targetonly and (not os.path.exists(calibtablepath)):
            print "For target-only, the SN file must exist already - aborting!"
            sys.exit(1)
        if not targetonly:
            vlbatasks.loadtable(calibsuvdata, calibtablepath, 
                                snversion+sncount)
        if not calonly:
            for i in range(maxphasecentres):
                vlbatasks.loadtable(targetuvdata[i], calibtablepath, 
                                    snversion+sncount)
        sncount += 1
    if not targetonly:
        vlbatasks.deletetable(calibsuvdata, 'SN', snversion+sncount)
        vlbatasks.mergesntables(calibsuvdata, snversion, sncount, refant)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.deletetable(targetuvdata[i], 'SN', snversion+sncount)
            vlbatasks.mergesntables(targetuvdata[i], snversion, sncount, refant)
    if not targetonly:
        vlbatasks.deletetable(calibsuvdata, 'CL', clversion+1)
        vlbatasks.applysntable(calibsuvdata, snversion+sncount, '2PT', clversion, refant)
    if not calonly:
        for i in range(maxphasecentres):
            vlbatasks.deletetable(targetuvdata[i], 'CL', clversion+1)
            vlbatasks.applysntable(targetuvdata[i], snversion+sncount, '2PT', clversion, refant)
else:
    print "Skipping loading/application of phase reference source CALIB solutions"
    sncount = 0
    phscalnames = []
    for phscal in phscalsrc:
        if phscal in phscalnames:
            continue
        phscalnames.append(phscal)
        sncount += 1

runlevel  = runlevel + 1
if not skipphsrefcalib:
    snversion = snversion + sncount + 1
    clversion = clversion + 1
targetjmfitfilelist = []
targetfixedjmfitfilelist = []
phscaljmfitfilelist = []
selfcallist         = []
#aipsnameslist       = []
cleanclass          = 'ICL001'
finalclass          = 'SPLIT'
# Imaging parameters
targetsearchpix     = ntargetsearchpix
targetsearchpixsize = targetsearchpixmas
targetsearchcc      = 0
targetsearchflux    = 0
targetfinalpix      = 1024
targetfinalpixsize  = 0.75
calibpix            = 1024
calibpixsize        = 1
calibnumcc          = 1000
calibtaper          = 0
aipsclass           = 'WFSPT'
skipping            = False
if skipuntil != "":
    skipping = True
catfiles = {}
if len(combineobs) > 1:
    duplicatelist = open('/home/adeller/svn_codebase/nraosvn/mjive20/schedules/duplicates.txt').readlines()
    duplicates = []
    for line in duplicatelist:
        duplicates.append(line.split(':'))
    for o in combineobs:
        wfsptfiles = glob.glob(rootdir + '/fullresdatasets/*fullres_uv.fits')
        for f in wfsptfiles:
            src = f.split('/')[-1].split('_')[1]
            for d in duplicates:
                if src == d[0]:
                    src = d[1].strip()
            if not src in catfiles:
                catfiles[src] = [f]
            else:
                catfiles[src].append(f)
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Image each potential source using IMAGR (widefield) and JMFIT ###############
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Doing wide-field image search"
    phscalnames = []
    targetnames = []
    targetcount = []
    orgtargetindex = []
    # First do the phs ref sources (if requested)
    if not targetonly:
        phscalnames = []
        for i in range(50): #Clear any old FINAL split catalog entries
            for phscal in phscalsrc:
                splitdata = AIPSUVData(phscal, 'FINAL', 1, i)
                if splitdata.exists():
                    splitdata.zap()
                splitdata = AIPSUVData(phscal, 'WFSPT', 1, i)
                if splitdata.exists():
                    splitdata.zap()
        for i in range(numtargets):
            if not phscalsrc[i] in phscalnames:
                phscalnames.append(phscalsrc[i])
                splitdata = AIPSUVData(phscalsrc[i], 'FINAL', 1, 1)
                doband = True
                if skipampcalfring:
                    doband = False
                if doband:
                    vlbatasks.split(calibsuvdata, clversion, 'FINAL', phscalsrc[i])
                else:
                    vlbatasks.splittoseqnoband(calibsuvdata, clversion, 'FINAL', phscalsrc[i], 1)
                if zapconfusing:
                    vlbatasks.splittoseq(calibsuvdata, clversion, aipsclass, phscalsrc[i], 1, True, doband)
                phscaluvfile = directory + '/' + phscalsrc[i] + '_uv.fits'
                vlbatasks.writedata(splitdata, phscaluvfile, True)
                cleanradiusmas = 30
                cleanimage = AIPSImage(phscalsrc[i], 'ICL001', 1, 1)
                vlbatasks.widefieldimage(splitdata, phscalsrc[i], calibpix, calibpixsize, True, 0.0007, calibtaper, 0.0, 0.0, calibnumcc, cleanradiusmas)
                aipsname = phscalsrc[i]
                imagrimagefile = directory + '/' + experiment + '_' + aipsname + ".imagr.clean.fits"
                imagrjmfitfile = directory + '/' + experiment + '_' + aipsname + ".imagr.jmfit"
                vlbatasks.writedata(cleanimage, imagrimagefile, True)
                snr1, snr2 = vlbatasks.nonpulsarjmfit(imagrimagefile, imagrjmfitfile, aipsname)
                phscaljmfitfilelist.append(imagrjmfitfile)
                psfile = directory + '/images/' + aipsname + '.clean.ps'
                jpgfile  = directory + '/images/' + aipsname + '.clean.jpg'
                vlbatasks.make_contour_plot(cleanimage, psfile, calibpix/2,
                                            calibpix/2, jpgfile, snr1)
    if not calonly:
        for i in range(numtargets):
            selfcallist.append([])
        doband = True
        if skipampcalfring:
            doband = False
        if len(combineobs) > 1:
            directory = rootdir + '/combinedobs/'
            experiment = ""
            aipsnames = []
            for c in combineobs:
                experiment += c[-2:]
            for s in catfiles:
                targetnames.append(s)
                aipsnames.append(s)
                targetcount.append(1)
                orgtargetindex.append(0)
                inputuvdata = []
                outputuvdata = AIPSUVData(s, aipsclass, 1, 1)
                if outputuvdata.exists():
                    outputuvdata.zap()
                count = 1
                print catfiles[s]
                for f in catfiles[s]:
                    inputuvdata.append(AIPSUVData(s, 'JUNK', 1, count))
                    if inputuvdata[-1].exists():
                        inputuvdata[-1].zap()
                    vlbatasks.fitld_uvfits(f, inputuvdata[-1], "")
                    count += 1
                if len(inputuvdata) == 1:
                    inputuvdata[0].rename(s, aipsclass, 1)
                else:
                    vlbatasks.dbcon(inputuvdata, outputuvdata)
                    for junkdata in inputuvdata:
                        junkdata.zap()
        else:
            targetnames, aipsnames, targetcount, orgtargetindex = vlbatasks.splitanddbcon(targetuvdata, numtargets, 
            numfields, numfieldsources, fieldsourcenames, clversion, aipsclass, limitedtargetlist, doconcat, skipping,
            doband)
        # Identify any bright sources that may need flagging for some targets
        zapsources = []
        if zapconfusing:
            print "Checking for bright, nearby confusing sources..."
            # Assume a run has previously been done
            for i in range(numtargets):
                jmfitfile = directory + '/' + experiment + '_' + phscalsrc[i] + ".imagr.jmfit.stats"
                try:
                    result = ExtendedJmfitResult(phscalsrc[i], experiment, jmfitfile)
                    flux = result.intflux
                except ValueError:
                    result = BasicJmfitResult(phscalsrc[i], experiment, jmfitfile)
                    flux = result.peakflux
                if flux > 50:
                    zapsources.append(result)
                    print "Including " + result.object
            for t in targetnames:
                jmfitfile = directory + '/' + experiment + '_' + t + ".imagr.jmfit.stats"
                try:
                    result = ExtendedJmfitResult(t, experiment, jmfitfile)
                    flux = result.intflux
                except ValueError:
                    result = BasicJmfitResult(t, experiment, jmfitfile)
                    flux = result.peakflux
                if flux > 50 and result.snr > 10:
                    #Check that this source is not already included as the phase reference...
                    avoid = False
                    for z in zapsources:
                        arcsecdiff = astro_utils.posdiff(result.rarad, result.decrad, \
                                                     z.rarad, z.decrad)*180*60*60/math.pi
                        if arcsecdiff < 5:
                            avoid = True
                    if not avoid:
                        zapsources.append(result)
                        print "Including " + result.object
        flaggedpercentage = 0.0
        totallyflagged = 0
        totallyflaggedsources = []
        for i in range(len(targetnames)):
            imagrimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                             ".imagr.clean.fits"
            imagrjmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                             ".imagr.jmfit"
            limitsfile = directory + '/' + experiment + '_' + targetnames[i] + ".limits"
            imagrfixedjmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                                  ".imagr.fixed.jmfit"
            outfile = directory + '/' + targetnames[i] + '_shifted_uv.fits'
            if skipping:
                if targetnames[i].strip() == skipuntil.strip():
                    skipping = False
                else:
                    targetjmfitfilelist.append(imagrjmfitfile)
                    targetfixedjmfitfilelist.append(imagrfixedjmfitfile)
                    continue
            runto = 1
            aipsname = aipsnames[i]
            if not doconcat:
                runto = targetcount[i]
            # Do the flagging of confusing sources if needed
            print "About to do the zapconfusing stuff..."
            print zapconfusing
            if zapconfusing:
                print "Checking against the zap list (length " + str(len(zapsources)) + ")"
                currentzaplist = []
                jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + ".imagr.jmfit.stats"
                try:
                    result = ExtendedJmfitResult(targetnames[i], experiment, jmfitfile)
                except ValueError:
                    result = BasicJmfitResult(targetnames[i], experiment, jmfitfile)
                for zapsource in zapsources:
                    try:
                        flux = z.intflux
                    except AttributeError:
                        flux = z.peakflux
                    if targetnames[i] == zapsource.object: continue
                    arcsecdiff = astro_utils.posdiff(result.rarad, result.decrad, \
                                                     zapsource.rarad, zapsource.decrad)*180*60*60/math.pi
                    print "arcsecdiff to " + zapsource.object + " is " + str(arcsecdiff) + ", flux is " + str(flux) + " mJy"
                    if arcsecdiff > 1 and arcsecdiff < 90:
                        currentzaplist.append(zapsource)
                    elif arcsecdiff > 1 and arcsecdiff < 300 and flux > 1000:
                        currentzaplist.append(zapsource)
                    elif arcsecdiff > 1 and arcsecdiff < 1500 and flux > 250 and math.fabs(result.decrad) < 10*math.pi/180.0:
                        currentzaplist.append(zapsource)
                currentdata = AIPSUVData(aipsname, aipsclass, 1, 1)
                currentflagged = []
                for z in currentzaplist:
                    try:
                        flux = z.intflux
                    except AttributeError:
                        flux = z.peakflux
                    print "For " + aipsname + ", zapping bright source " + z.object + \
                          ", separated by " + str(arcsecdiff) + " arcseconds and with flux " + \
                          str(flux) + " mJy integrated flux"
                    brightsrcclass = aipsclass
                    brightsrcuvdata = AIPSUVData(z.object, brightsrcclass, 1, 1) 
                    sourceflagged = vlbatasks.smearingflag(brightsrcuvdata, currentdata, 50./float(flux))
                    currentflagged.append(sourceflagged)
                    flaggedpercentage += sourceflagged/len(targetnames)
                tempdata = AIPSUVData(aipsname,"JUNK",1,1)
                if tempdata.exists():
                    tempdata.zap()
                try:
                    gooddataestimate = 1.0
                    for cf in currentflagged:
                        gooddataestimate *= (1-cf)
                    print "For " + aipsname + ", gooddata estimate is " + str(gooddataestimate)
                    if gooddataestimate < UNFLAGGED_THRESHOLD:
                        raise RuntimeError("Too much flagged data")
                    vlbatasks.splat(currentdata,1,[0],tempdata)
                except RuntimeError:
                    print "Looks like " + aipsname + " was completely flagged due to an interfering source - continuing..."
                    totallyflagged += 1
                    totallyflaggedsources.append(aipsname)
                    if writefullrestargets:
                        os.system("rm -f " + rootdir + '/fullresdatasets/' + experiment + '_' + aipsname + '_fullres_uv.fits')
                    continue
                currentdata.zap()
                tempdata.rename(targetnames[i], aipsclass, 1)
            for j in range(runto):
                #Make a widefield image
                initimage = AIPSImage(aipsname, 'IIM001', 1, 1)
                beamimage = AIPSImage(aipsname, 'IBM001', 1, 1)
                cleanimage = AIPSImage(aipsname, cleanclass, 1, 1)
                if initimage.exists():
                    initimage.zap()
                if beamimage.exists():
                    beamimage.zap()
                if cleanimage.exists():
                    cleanimage.zap()
                imagrdata = AIPSUVData(aipsname, aipsclass, 1, j+1)
                cleanradiusmas = -1
                if writefullrestargets and not aipsname in totallyflaggedsources:
                    vlbatasks.writedata(imagrdata, rootdir + '/fullresdatasets/' + experiment + '_' + \
                                        aipsname + '_fullres_uv.fits', True)
                if noimage:
                     continue
                vlbatasks.widefieldimage(imagrdata, aipsname, targetsearchpix, targetsearchpixsize, 
                                         False, targetsearchflux, targettaper, 0.0, 0.0, targetsearchcc, 
                                         cleanradiusmas)
                guardband = 0.05
                stats = vlbatasks.getimagestats(initimage, guardband) # [peak, rms, peakx, peaky, rastr, decstr]
                if stats[0] == 0.0:
                    print "Skipping bad source " + aipsname
                    continue
                rashiftmas  = -(stats[2]-targetsearchpix/2)*targetsearchpixsize
                decshiftmas = (stats[3]-targetsearchpix/2)*targetsearchpixsize
                print "RA shift mas is " + str(rashiftmas)
                print "Dec shift mas is " + str(decshiftmas)
                print "Init peak was " + str(stats[0])
                print "init rms was " + str(stats[1])
                equivonsource = 0.0
                if len(onsourcetimes) > 0 and len(combineobs) <= 1:
                    for ost in onsourcetimes:
                        try:
                            e = ost[targetnames[i]]
                        except KeyError:
                            e = 0.0
                        equivonsource += e
                        print "Equiv on source time is now " + str(equivonsource)
                limitsout = open(limitsfile, "w")
                statssize = (1.0 - 2*guardband)*targetsearchpix/1000.0
                limitsout.write("Image size in arcsec for stats: %.2f x %0.2f\n" % \
                                (statssize, statssize))
                limitsout.write("Peak (mJy): %.3f\n" % (stats[0]*1000))
                limitsout.write("Measured RMS (mJy): %.3f\n" % (stats[1]*1000))
                if equivonsource > 0.0:
                    limitsout.write("Predicted RMS (mJy): %.3f\n" % (0.076*math.sqrt(3600.0/equivonsource)))
                else:
                    limitsout.write("Predicted RMS (mJy): -1.0\n")
                savedrms = stats[1]

                uvfixdata = AIPSUVData(aipsname, 'UVFIX', 1, 1)
                if uvfixdata.exists():
                    uvfixdata.zap()
                nchan = imagrdata.table('FQ', 1)[0].total_bandwidth[0]/\
				imagrdata.table('FQ', 1)[0].ch_width[0]
                vlbatasks.uvfix(imagrdata, uvfixdata, rashiftmas, decshiftmas)
                imagrdata.zap()
                splatdata = AIPSUVData(aipsname, finalclass, 1, 1)
                if splatdata.exists():
                    splatdata.zap()
                vlbatasks.splat(uvfixdata, nchan, [0], splatdata, -1, finalavtime/60.0)
                uvfixdata.zap()
                vlbatasks.writedata(splatdata, outfile, True)

                initimage.zap()
                AIPSImage(aipsname, 'IBM001', 1, 1).zap()
                cleanradiusmas     = 15
                vlbatasks.widefieldimage(splatdata, aipsname, targetfinalpix, targetfinalpixsize,
                                         True, targetfinalflux, targettaper, 0.0,
                                         0.0, targetfinalcc, cleanradiusmas)
                stats  = vlbatasks.getimagestats(cleanimage, 0.2) # [peak, rms, peakx, peaky, rastr, decstr]
                finalrapeakpixel  = stats[2]
                finaldecpeakpixel = stats[3]
                print "Final peak was " + str(stats[0])
                print "Final rms was " + str(stats[1])
                print "Which is probably a lot less than the original rms: " + str(savedrms)
                vlbatasks.writedata(cleanimage, imagrimagefile, True)
                peak1, snr1 = vlbatasks.nonpulsarjmfit(imagrimagefile, imagrjmfitfile, 
                                                       targetnames[i], finalrapeakpixel, 
                                                       finaldecpeakpixel, True)
                targetjmfitfilelist.append(imagrjmfitfile)
                peak2, snr2 = vlbatasks.nonpulsarjmfit(imagrimagefile, imagrfixedjmfitfile,
                                                       targetnames[i], finalrapeakpixel,
                                                       finaldecpeakpixel, False)
                targetfixedjmfitfilelist.append(imagrfixedjmfitfile)

                #Produce a greyscale plot with contours
                #peak   = float(stats[0])
                #rms    = float(stats[1])
                if snr1 > clearlyinsnr or peak1 > clearlyinflux or \
                   snr2 > clearlyinsnr or peak2 > clearlyinflux or \
                   (peak1 < 0.75*peak2 and peak2 > 0.75*clearlyinflux):
                    psfile = directory + '/images/' + targetnames[i] + '.clean.ps'
                    jpgfile  = directory + '/images/' + targetnames[i] + '.clean.jpg'
                    vlbatasks.make_contour_plot(cleanimage, psfile, finalrapeakpixel, 
                                                finaldecpeakpixel, jpgfile, snr1)
                if (snr1 > searchselfcalsnr or snr2 > searchselfcalsnr) and len(combineobs) <= 1:
                    print len(selfcallist)
                    print "i is " + str(i)
                    print len(orgtargetindex)
                    print orgtargetindex[i]
                    snr = snr1
                    if snr2 > snr1: snr = snr2
                    selfcallist[orgtargetindex[i]].append([aipsname, finalclass, cleanclass, snr])
    if not calonly:
        flagsummaryfile = directory + '/' + experiment.lower() + '.flagsummary'
        flagsumout = open(flagsummaryfile, "w")
        flagsummary = "%.3f %% of the data was flagged due to potential confusing sources" % (100*flaggedpercentage)
        flagsummary += ("\n%d / %d sources were completely flagged\n" % (totallyflagged, len(targetnames)))
        for s in totallyflaggedsources:
            flagsummary += s + '\n'
        print flagsummary
        flagsumout.write(flagsummary)
        flagsumout.close()

    # Write out the lists of jmfit files, selfcal files and aips names
    if not targetonly:
        jmfitlistout = open(directory + '/phscaljmfitfiles.list', 'w')
        for filename in phscaljmfitfilelist:
            jmfitlistout.write(filename + '\n')
        jmfitlistout.close()
    if not calonly:
        jmfitlistout = open(directory + '/targetjmfitfiles.list', 'w')
        for i in range(len(targetjmfitfilelist)):
            jmfitlistout.write(targetjmfitfilelist[i] + ' ' + targetfixedjmfitfilelist[i] + '\n')
        jmfitlistout.close()
        selfcallistout = open(directory + '/selfcal.list', 'w')
        count = 0
        for ts in selfcallist:
            for s in ts:
                print s[0]
                print s[1]
                print s[2]
                towrite = "%d %s %s %s %0.2f\n" % (count, s[0],s[1],s[2], s[3])
                selfcallistout.write(towrite)
            count += 1
        selfcallistout.close()
#        aipsnameout = open(directory + '/aipsnames.list', 'w')
#        for a in aipsnameslist:
#            aipsnameout.write("%d %d %s %s\n" % (a[0], a[1], a[2], a[3]))
#        aipsnameout.close()
else:
    print "Skipping wide-field imaging of targets"

if calonly:
    print "bye"
    sys.exit()

if os.path.exists(directory + '/selfcal.list'):
    selfcallistin = open(directory + '/selfcal.list')
    selfcallistlines = selfcallistin.readlines()
    selfcallistin.close()
else:
    print "Selfcallist does not exist - hope you are just running through " + \
          "to look at runlevels, otherwise will crash soon!"
    selfcallistlines = []
selfcallist = []
#aipsnamelistin = open(directory + '/aipsnames.list')
#aipsnamelistlines = aipsnamelistin.readlines()
#aipsnamelistin.close()
#aipsnameslist = []
for i in range(numtargets):
    selfcallist.append([])
#    aipsnameslist.append([])
for line in selfcallistlines:
    selfcallist[int(line.split()[0])].append(line.split()[1:])
#for line in aipsnamelistlines:
#    aipsnameslist[int(line.split()[0])].append(line.split()[1:])

runlevel += 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Do a joint selfcal on any decent sources found ##############################
if runfromlevel <= runlevel and runtolevel >= runlevel and not calonly:
    print "Runlevel " + str(runlevel) + ": Selfcal'ing on targets"
    count = 1
    for sublist in selfcallist:
        if len(sublist) > 0:
            normed = []
            for s in sublist:
                srcname  = s[0]
                uvklass  = s[1]
                clnklass = s[2]
                snr      = s[3]
                uvdata   = AIPSUVData(srcname, uvklass, 1, 1)
                if '/' in clnklass: # This is the file to load, not an AIPS class
                    manualcleanfile = clnklass
                    clnimage = AIPSImage(srcname, "MANIM", 1, 1)
                    if clnimage.exists():
                        clnimage.zap()
                    vlbatasks.fitld_image(manualcleanfile, clnimage)
                else:
                    clnimage = AIPSImage(srcname, clnklass, 1, 1)
                if not uvdata.exists():
                    print "The required source data with name %s and class %s does not exist!" % (srcname, uvklass)
                    print "You need to go back to runlevel %d and reimage" % (runlevel-1)
                    sys.exit()
                if not clnimage.exists():
                    print "The required image data with name %s and class %s does not exist!" % (srcname, clnklass)
                    print "You need to go back to runlevel %d and reimage" % (runlevel-1)
                    sys.exit()
                normdata = AIPSUVData(srcname, 'NORMUV', 1, 1)
                if normdata.exists():
                    normdata.zap()
                vlbatasks.normaliseUVData(uvdata, clnimage, normdata)
                normed.append(normdata)
            if len(normed) == 1:
                jointdata = normed[0]
            else:
                concatuvdata = AIPSUVData('CONCAT' + str(count), 'DBCON', 1, 1)
                if concatuvdata.exists():
                    concatuvdata.zap()
                for n in normed[1:]:
                    vlbatasks.match_headersource(normed[0], n)
                vlbatasks.dbcon(normed, concatuvdata)
                jointdata = AIPSUVData('CONCAT' + str(count), 'UVSRT', 1, 1)
                if jointdata.exists():
                    jointdata.zap()
                vlbatasks.uvsrt(concatuvdata, jointdata)
                for n in normed:
                    n.zap()
            vlbatasks.singlesource_calib(jointdata, 1.0, 1, refant, False, searchselfcaldur,
                                         True, 'L1R', 4.5, searchselfcalsumifs)
            sntable = tabledir + 'targetselfcal.%d.sn' % (count)
            if os.path.exists(sntable):
                os.remove(sntable)
            vlbatasks.writetable(jointdata, 'SN', 1, sntable)
        count += 1
    sncount = 0
    deleted = vlbatasks.deletetable(targetuvdata[0], 'SN', snversion)
    if deleted:
        print "Had to delete an old SN table - hope you are not squishing anything!"
    for i in range(numtargets):
        sntable = tabledir + 'targetselfcal.%d.sn' % (i+1)
        if os.path.exists(sntable):
            vlbatasks.loadtable(targetuvdata[0], sntable, snversion+i)
            sncount += 1
    vlbatasks.mergesntables(targetuvdata[0], snversion, sncount, refant)
    sntable = tabledir + 'targetselfcal.merged.sn'
    if os.path.exists(sntable):
        os.system("rm -f %s" % (sntable))
    vlbatasks.writetable(targetuvdata[0], 'SN', snversion+sncount, sntable)
    for i in range(numtargets+1):
        vlbatasks.deletetable(targetuvdata[0], 'SN', snversion+i)
else:
    print "Skipping joint target selfcal"

runlevel += 1
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Apply joint selfcal results #################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not calonly:
    print "Runlevel " + str(runlevel) + ": Applying selfcal solutions"
    for uvdata in targetuvdata:
        sntable = tabledir + 'targetselfcal.merged.sn'
        if not os.path.exists(sntable):
            print "Can't find merged selfcal file %s" % (sntable)
            sys.exit()
        deleted = vlbatasks.deletetable(uvdata, 'SN', snversion)
        if deleted:
            print "Had to delete an old SN table - hope you are not squishing anything!"
        vlbatasks.loadtable(uvdata, sntable, snversion)
        vlbatasks.applysntable(uvdata, snversion, '2PT', clversion, refant)
else:
    print "Skipping application of joint target selfcal"

runlevel += 1
snversion += 1
clversion += 1
selfcalklass = 'WFSCAL'
print "Runlevel is " + str(runlevel) + " and clversion is " + str(clversion)
## Image using the selfcal results #############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not calonly:
    print "Runlevel " + str(runlevel) + ": Redoing the imaging, with selfcal corrections"
    selfcaljmfitfilelist = []
    selfcalfixedjmfitfilelist = []
    if skipampcalfring:
        doband = False
    #for i in range(len(phscalsrc)):
    #    if doband:
    #        vlbatasks.split(calibsuvdata, clversion, 'FINAL', phscalsrc[i])
    #    else:
    #        vlbatasks.splittoseqnoband(calibsuvdata, clversion, 'FINAL', phscalsrc[i], 1)
    targetnames, aipsnames, targetcount, orgtargetindex = vlbatasks.splitanddbcon(targetuvdata,
        numtargets, numfields, numfieldsources, fieldsourcenames, clversion, selfcalklass, 
        limitedtargetlist, doconcat, False, doband)    
#    for c in range(50): #Clear any old WFSCAL split catalog entries
#        for j in range(numtargets):
#            for a in range(len(aipsnameslist[j])):
#                name = aipsnameslist[j][a][3]
#                if len(name) > 16:
#                    name = name[:12]
#                splitdata = AIPSUVData(name, selfcalklass, 1, c)
#                if splitdata.exists():
#                    splitdata.zap()
#    for j in range(numtargets):
#        for a in range(len(aipsnameslist[j])):
#            print "About to split " + aipsnameslist[j][a][3]
#            vlbatasks.splitmulti(targetuvdata[aipsnameslist[j][a][1]], clversion, selfcalklass,
#                                 aipsnameslist[j][a][3], 1)
#            if len(aipsnameslist[j][a][3]) >= 16: # Crappy old style name
#                splitdata = AIPSUVData(aipsnameslist[j][a][3][:12], selfcalklass, 1, 1)
#                splitdata.rename(aipsnameslist[j][a][2], selfcalklass, 1)
#    for i in range(numtargets):
#        for a in range(len(aipsnameslist[i])):
#            aipsname = aipsnamelist[i][a][2]
    count = -1
    for aipsname in aipsnames:
        count += 1
        #Make a widefield image
        initimage = AIPSImage(aipsname, 'IIM001', 1, 1)
        beamimage = AIPSImage(aipsname, 'IBM001', 1, 1)
        cleanimage = AIPSImage(aipsname, cleanclass, 1, 1)
        if initimage.exists():
            initimage.zap()
        if beamimage.exists():
            beamimage.zap()
        if cleanimage.exists():
            cleanimage.zap()
        imagrdata = AIPSUVData(aipsname, selfcalklass, 1, 1)
        cleanradiusmas = -1
        vlbatasks.widefieldimage(imagrdata, aipsname, targetsearchpix, targetsearchpixsize,
                                 False, targetsearchflux, targettaper, 0.0, 0.0, targetsearchcc,
                                 cleanradiusmas)
        stats = vlbatasks.getimagestats(initimage, 0.05) # [peak, rms, peakx, peaky, rastr, decstr]
        if stats[0] == 0.0:
            print "Skipping bad source " + aipsname
            continue
        rashiftmas  = -(stats[2]-targetsearchpix/2)*targetsearchpixsize
        decshiftmas = (stats[3]-targetsearchpix/2)*targetsearchpixsize

        uvfixdata = AIPSUVData(aipsname, 'UVFIX2', 1, 1)
        if uvfixdata.exists():
            uvfixdata.zap()
        nchan = imagrdata.table('FQ', 1)[0].total_bandwidth[0]/\
                imagrdata.table('FQ', 1)[0].ch_width[0]
        vlbatasks.uvfix(imagrdata, uvfixdata, rashiftmas, decshiftmas)
        imagrdata.zap()
        finalclass = 'SSPLIT'
        splatdata = AIPSUVData(aipsname, finalclass, 1, 1)
        outfile = directory + '/' + aipsname + '_shiftedselfcald_uv.fits'
        if splatdata.exists():
            splatdata.zap()
        vlbatasks.splat(uvfixdata, nchan, [0], splatdata, -1, finalavtime/60.0)
        uvfixdata.zap()
        vlbatasks.writedata(splatdata, outfile, True)

        initimage.zap()
        AIPSImage(aipsname, 'IBM001', 1, 1).zap()
        cleanradiusmas     = 15
        vlbatasks.widefieldimage(splatdata, aipsname, targetfinalpix, targetfinalpixsize,
                                 True, targetfinalflux, targettaper, 0.0,
                                 0.0, targetfinalcc, cleanradiusmas)
        stats  = vlbatasks.getimagestats(cleanimage, 0.2) # [peak, rms, peakx, peaky, rastr, decstr]
        finalrapeakpixel  = stats[2]
        finaldecpeakpixel = stats[3]
        imagrimagefile = directory + '/' + experiment + '_' + aipsname + \
                         ".imagr.clean.selfcal.fits"
        imagrjmfitfile = directory + '/' + experiment + '_' + aipsname + \
                         ".imagr.jmfit.selfcal"
        imagrfixedjmfitfile = directory + '/' + experiment + '_' + aipsname + \
                              ".imagr.fixed.jmfit.selfcal"
        vlbatasks.writedata(cleanimage, imagrimagefile, True)
        peak1, snr1 = vlbatasks.nonpulsarjmfit(imagrimagefile, imagrjmfitfile,
                                               aipsname, finalrapeakpixel,
                                               finaldecpeakpixel, True)
        selfcaljmfitfilelist.append(imagrjmfitfile)
        peak2, snr2 = vlbatasks.nonpulsarjmfit(imagrimagefile, imagrfixedjmfitfile,
                                               aipsname, finalrapeakpixel,
                                               finaldecpeakpixel, False)
        selfcalfixedjmfitfilelist.append(imagrfixedjmfitfile)

        #Produce a greyscale plot with contours
        if snr1 > clearlyinsnr or peak1 > clearlyinflux or \
           snr2 > clearlyinsnr or peak2 > clearlyinflux or \
          (peak1 < 0.75*peak2 and peak2 > 0.75*clearlyinflux):
            psfile = directory + '/images/' + targetnames[count] + '.clean.selfcal.ps'
            jpgfile  = directory + '/images/' + targetnames[count] + '.clean.selfcal.jpg'
            vlbatasks.make_greyscale_plot(cleanimage, psfile, finalrapeakpixel,
                                          finaldecpeakpixel, jpgfile=jpgfile)
    # Write out the list of selfcal jmfit files
    jmfitlistout = open(directory + '/targetselfcaldjmfitfiles.list', 'w')
    for i in range(len(selfcaljmfitfilelist)):
        jmfitlistout.write(selfcaljmfitfilelist[i] + ' ' + selfcalfixedjmfitfilelist[i] + '\n')
    jmfitlistout.close()
else:
    print "Skipping imaging with selfcal corrections"    

#if startlocaltv:
#    tv.kill()

