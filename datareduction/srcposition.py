#!/usr/bin/env ParselTongue
import os,glob,sys,yaml,pickle
import howfun,mspsrpifun
#import numpy as np
#from optparse import OptionParser
#np.set_printoptions(suppress=True)

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

def modeltype(targetname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = mspsrpifun.prepare_path_source(targetname)
    sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
    sourcefile = sourcefiles[0]
    expno = sourcefile.split('/')[-1].split('.')[0].strip().lower()
    expconfigfile = configdir + '/' +  expno + '.yaml'
    if not os.path.exists(expconfigfile):
        print("%s does not exist; aborting\n" % expconfigfile)
        sys.exit()
    expconfig = yaml.load(open(expconfigfile))
    userno = expconfig['userno']
    useprelimmodels = expconfig['useprelimmodels']
    if useprelimmodels:
        modeltype = 'preliminary'
    else:
        modeltype = 'final'
    return modeltype, userno
def src_RA_Dec(targetname, srcname):
    [modeltype1, AIPS.userno] = modeltype(targetname)
    print modeltype(targetname)
    [auxdir, configdir, targetdir, phscalname, prIBCname] = mspsrpifun.prepare_path_source(targetname)
    srcmodel = auxdir + '/sourcemodels/' + modeltype1 + '/' + srcname + '.clean.fits'
    print srcmodel
    modeldata = AIPSImage(srcname, "CLEAN", 1, 1)
    if not modeldata.exists() and os.path.exists(srscmodel):
        vlbatasks.fitld_image(srcmodel, modeldata)
    if not modeldata.exists() and not os.path.exists(srcmodel):
        print("neither srcmodel nor loaded counterpart in AIPS is found; aborting\n")
        sys.exit()
    RA_deg = modeldata.header.crval[0]
    Dec_deg = modeldata.header.crval[1]
    print RA_deg/15, Dec_deg
    RA_dms = howfun.deg2dms(RA_deg/15)
    Dec_dms = howfun.deg2dms(Dec_deg)
    return RA_dms, Dec_dms
## main program ########################################################################################
if len(sys.argv) < 3:
    print("srcposition.py targetname srcname")
    sys.exit()
[junk, targetname, srcname] = sys.argv
if len(srcname) > 12:
    srcname = srcname[:12]
position = src_RA_Dec(targetname, srcname)
print position
## use pickle to save result, thus making the results accessible by outside scripts
[auxdir, configdir, targetdir, phscalname, prIBCname] = mspsrpifun.prepare_path_source(targetname)
tmpfile = targetdir + '/.variables.tmp'
writefile = open(tmpfile, 'w')
pickle.dump(position, writefile)
writefile.close()
