#!/usr/bin/env ParselTongue

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
import sys, os, string, math, yaml
import vlbatasks
from time import gmtime, strftime
from optparse import OptionParser

try:
    svndir = os.environ['PSRVLBAUXDIR']
except KeyError:
    print "PSRVLBAUXDIR is not defined - aborting!"
    sys.exit(1)
print "PSRVLBAUXDIR is ", svndir
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    aipsver = '31DEC18'

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Experiment name")
parser.add_option("--userno", dest="userno", default="204",
                  help="The AIPS userno to use")
(options, junk) = parser.parse_args()
AIPS.userno     = int(options.userno)
experiment      = options.experiment
configdir       = svndir + '/configs/'

if experiment == "":
    parser.error("You must supply an experiment name")

expconfigfile = configdir + experiment + '.yaml'
print expconfigfile
if not os.path.exists(expconfigfile):
    parser.error("Experiment config file %s does not exist!" % expconfigfile)
expconfig     = yaml.load(open(expconfigfile))
expdir        = expconfig['rootdir'] + '/' + experiment + '/'
numtargets    = len(expconfig['targets'])
targetconfigs = []
for i in range(numtargets):
    targetconfigfile = configdir + expconfig["targets"][i] + '.yaml'
    if not os.path.exists(targetconfigfile):
        parser.error("Target config file %s does not exist!" % targetconfigfile)
    targetconfigs.append(yaml.load(open(targetconfigfile)))

for i in range(numtargets):
    targetname = expconfig["targets"][i]
    targetconfig = targetconfigs[i]
    sourcefile = expdir + experiment + '.source'
    sourcefilelines = open(sourcefile).readlines()
    phsrefname = ""
    modeltypestr = 'final'
    for j in range(len(sourcefilelines)):
        if targetname in sourcefilelines[j]:
            phsrefname = sourcefilelines[j+1].split(':')[-1].strip()
    if phsrefname == "":
        print "Couldn't find the phase reference source for " + targetname
        sys.exit()
    if expconfig['useprelimmodels']:
        modeltypestr = 'preliminary'
    aips8phsrefname = phsrefname
    aips12phsrefname = phsrefname
    if len(phsrefname) > 8:
        aips8phsrefname = phsrefname[:8]
    if len(phsrefname) > 12:
        aips12phsrefname = phsrefname[:12]
    phsrefdatafile = expdir + experiment + '_' + phsrefname + '_pipeline_uv.fits'
    phsrefimagefile = svndir + 'sourcemodels/' + modeltypestr + '/' + phsrefname + '.clean.fits'
    shiftedfile = expdir + experiment + '_' + phsrefname + '_ibshiftdiv_uv.fits'
    jmfitfile = expdir + experiment + '_' + phsrefname + '_ibshiftdiv.jmfit'
    phsrefdata = AIPSUVData(aips8phsrefname, "TEMP", 1, 1)
    phsrefimage = AIPSImage(aips12phsrefname, "TMPCLN", 1, 1)
    splatdata = AIPSUVData(aips8phsrefname, "SPLAT", 1, 1)
    normdata = AIPSUVData(aips8phsrefname, "NORM", 1, 1)
    shiftedimage = AIPSImage(aips12phsrefname, "ICL001", 1, 1)
    shiftedbeam = AIPSImage(aips12phsrefname, "IBM001", 1, 1)
    if phsrefdata.exists():
        phsrefdata.zap()
    if splatdata.exists():
        splatdata.zap()
    if phsrefimage.exists():
        phsrefimage.zap()
    if normdata.exists():
        normdata.zap()
    if shiftedimage.exists():
        shiftedimage.zap()
    if shiftedbeam.exists():
        shiftedbeam.zap()
    vlbatasks.fitld_uvfits(phsrefdatafile, phsrefdata, "")
    vlbatasks.fitld_image(phsrefimagefile, phsrefimage)
    if "," in targetconfig['primaryinbeam']:
        caltable = expdir + 'tables/CONCAT0.icalib.p1.sn'
    else:
        caltable = expdir + 'tables/' + targetconfig['primaryinbeam'] + '.icalib.p1.sn'
    if not os.path.exists(caltable):
        print caltable + ' does not exist'
        sys.exit()
    vlbatasks.loadtable(phsrefdata, caltable, 1)
    vlbatasks.splat(phsrefdata, 1, [0], splatdata, 1, 0)
    vlbatasks.normaliseUVData(splatdata, phsrefimage, normdata)
    vlbatasks.writedata(normdata, shiftedfile, True)
    vlbatasks.widefieldimage(normdata, aips12phsrefname, 256, 0.75, True, 0.050,
                   0, 0, 0, 100, 20)
    vlbatasks.nonpulsarjmfit("", jmfitfile, aips12phsrefname, -1, -1, True,
                             False,shiftedimage,48)
