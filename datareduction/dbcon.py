#!/usr/bin/env ParselTongue
#Imports ########################################################
from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData
from optparse import OptionParser
import os, sys, glob

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
                  default="/export/home/marathon2/data/1023+0038/",
                  help="The base directory for the experiment series")
parser.add_option("-e", "--expseries", dest="expseries", default="",
                  help="Experiment series name")
parser.add_option("-c", "--expcodes", dest="expcodes", default="",
                  help="The codes for the experiments e.g. A,B,C")
#parser.add_option("-i", "--inbeamname", dest="inbeamname", default="inbeam",
#                  help="The name of the inbeam in the filenames")
parser.add_option("--middlestring", dest="middlestring", default="J1023+0031_formodel",
                  help="The string between eg BD141B_ and _uv.fits")
parser.add_option("--nofirststring", dest="nofirststring", default=False,
                  action="store_true", help="Skip the BDXXX_ prefix on filenames")
parser.add_option("--ifs", dest="ifs", default="",
                  help="beginIF,endIF - leave blank for all")
parser.add_option("--filelist", dest="filelist", default="",
                  help="List of files (alternative override to specify files)")
parser.add_option("--filelistdir", dest="filelistdir", default="",
                   help="Prepend all filelist files with this")
(options, junk) = parser.parse_args()
AIPS.userno     = 2575
rootdir         = options.rootdir
expseries       = options.expseries
expcodes        = options.expcodes.split(',')
nofirststring   = options.nofirststring
#inbeamname      = options.inbeamname
middlestring    = options.middlestring
ifs             = options.ifs
filelist        = options.filelist
filelistdir     = options.filelistdir

################################################################################
# Check validity of inputs
################################################################################
if filelist == "" and expseries == "":
    print "Doing a glob of all *.fits files in the current directory"
    doglob = True
    flist = glob.glob("*.fits")
    if len(flist) > 0:
        for f in flist:
            filelist = filelist + f + ','
    filelist = filelist[:-1]
    filelistdir = os.getcwd()
    if len(flist) < 2:
        parser.error("Only found " + str(len(flist)) + " fits files!")
if filelist == "" and len(expcodes) < 2:
    parser.error("You must supply at least two experiments!")

# Set up the fitld, fittp and dbcon objects ####################################
dbcon = AIPSTask('dbcon', version = aipsver)
dbcon.reweight[1:] = [0]
dbcon.dopos[1:][1:] = [0]
dbcon.doarray = 0
dbcon.fqtol = -1
fitld = AIPSTask('fitld', version = aipsver)
fitld.optype = ''
fitld.ncount = 0
fitld.dotable = 1
fitld.douvcomp = 1
fitld.doconcat = -1
fitld.digicor = -1
if ifs != "":
    fitld.bif = int(ifs.split(',')[0])
    fitld.eif = int(ifs.split(',')[1])
fittp = AIPSTask('fittp', version = aipsver)
#fittp.dostokes = -1
#fittp.donewtab = 1
fittp.format = 0

# Remove any existing dbcon file ###############################################
#outfile = rootdir + '/models/dbcon_' + inbeamname + '_uv.fits'
outfile = rootdir + '/models/dbcon_' + middlestring + '_uv.fits'
if filelist != "":
    outfile = os.getcwd() + "/dbcon_uv.fits"
if os.path.exists(outfile):
    os.system('rm -f ' + outfile)

# Find all the files we can ####################################################
count = 0
activeexps = []
if filelist != "":
    for f in filelist.split(','):
        if len(f) > 1 and f[0] != '/':
            f = filelistdir + '/' + f
        if os.path.exists(f):
            modeldata = AIPSUVData('DBCON', "%03d" % count, 1, 1)
            if modeldata.exists():
                modeldata.zap()
            os.system("rm -f " + filelistdir + "/templink.fits")
            os.system("ln -s %s %s/templink.fits" % (f, filelistdir))
            fitld.datain = filelistdir + "/templink.fits"
            fitld.outdata = modeldata
            fitld()
            os.system("rm -f " + filelistdir + "/templink.fits")
            activeexps.append("%03d" % count)
            count = count+1
        else:
            print "Couldn't find pipeline'd FITS file " + f
else:
    for expcode in expcodes:
        experiment = expseries + expcode
        firststring = experiment.upper() + '_'
        laststring = '_uv.fits'
        if nofirststring:
            firststring = ""
            laststring = 'uv.fits'
        infile = rootdir + '/' + experiment.lower() + '/' + firststring + \
                 middlestring + laststring
#                 '_' + inbeamname + '_pipeline_uv.fits'
        if os.path.exists(infile):
            modeldata = AIPSUVData('DBCON', expcode, 1, 1)
#           modeldata = AIPSUVData(inbeamname, experiment.upper(), 1, 1)
            if modeldata.exists():
                modeldata.zap()
            os.system("rm -f templink.fits")
            os.system("ln -s %s templink.fits" % infile)
            fitld.datain = "templink.fits"
            fitld.outdata = modeldata
            fitld()
            os.system("rm -f templink.fits")
            count = count+1
            activeexps.append(expcode)
        else:
            print "Couldn't find pipeline'd FITS file for experiment " + experiment + \
                  "(" + infile + ")"

# Do the dbcon #################################################################
if count > 1:
#    dbcon.indata = AIPSUVData(inbeamname, activeexps[0], 1, 1)
#    dbcon.in2data = AIPSUVData(inbeamname, activeexps[1], 1, 1)
    dbcon.indata  = AIPSUVData('DBCON', activeexps[0], 1, 1)
    dbcon.in2data = AIPSUVData('DBCON', activeexps[1], 1, 1)
    junkodata = AIPSUVData('JUNK', 'JUNK', 1, 1)
    dbcon.outdata = junkodata
    if junkodata.exists():
        junkodata.zap()
    dbcon()
    for i in range(count-2):
#        dbcon.indata = AIPSUVData(inbeamname, activeexps[i+2], 1, 1)
        dbcon.indata = AIPSUVData('DBCON', activeexps[i+2], 1, 1)
        junkidata =  AIPSUVData('JUNK', 'JUNK', 1, i+1)
        dbcon.in2data = junkidata
        junkodata = AIPSUVData('JUNK', 'JUNK', 1, i+2)
        dbcon.outdata =  junkodata
        if junkodata.exists():
            junkodata.zap()
        dbcon()
    #Write the results out
    fittp.indata = junkodata
    fittp.dataout = outfile
    fittp()
else:
#    print "Only found " + str(count) + " for " + inbeamname + " - aborting"
    print "Only found " + str(count) + " for " + middlestring + " - aborting"

# Zap all the loaded data and dbcon'd data #####################################
for e in activeexps:
#    modeldata = AIPSUVData(inbeamname, experiment.upper(), 1, 1)
    modeldata = AIPSUVData('DBCON', e, 1, 1)
    if modeldata.exists():
        modeldata.zap()
for i in range(count-1):
    modeldata = AIPSUVData('JUNK', 'JUNK', 1, i+1)
    modeldata.zap()

