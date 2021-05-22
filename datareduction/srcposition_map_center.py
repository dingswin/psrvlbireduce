#!/usr/bin/env ParselTongue
import os,glob,sys,pickle,yaml
import howfun,mspsrpifun
import vlbatasks
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

def src_RA_Dec(targetname, srcname): #use different (for target/other calibrators) AIPS catalog entry to get phase center
    [modeltype1, AIPS.userno] = mspsrpifun.modeltype(targetname)
    print mspsrpifun.modeltype(targetname)
    [auxdir, configdir, targetdir, phscalname, prIBCname] = mspsrpifun.prepare_path_source(targetname)
    aa = mspsrpifun.measure_the_angular_broadened_size_of_the_target(targetname)
    targetname_AIPS = mspsrpifun.expno2sources(aa.expnos[0])[0]
    srcmodel = auxdir + '/sourcemodels/' + modeltype1 + '/' + srcname + '.clean.fits'
    print srcmodel
    if srcname != targetname:
        modeldata = AIPSImage(srcname, "CLEAN", 1, 1)
    else:

        modeldata = AIPSUVData(targetname_AIPS, "GFINL", 1, 1) #gated final
    print modeldata
    if not modeldata.exists() and not os.path.exists(srcmodel):
        print("neither srcmodel nor loaded counterpart in AIPS is found; aborting\n")
        sys.exit()
    if os.path.exists(srcmodel):
        if modeldata.exists():
            modeldata.zap()
        vlbatasks.fitld_image(srcmodel, modeldata)
        #fitld = AIPSTask('fitld', version = '31DEC18')
        #fitld.datain = srcmodel
        #fitld.outdata = modeldata
        #fitld.dotable = 1
        #fitld.douvcomp = -1
        #fitld()
        
    if srcname != targetname:
        RA_deg = modeldata.header.crval[0]
        Dec_deg = modeldata.header.crval[1]
    else:
        RA_deg = modeldata.header.crval[4]
        Dec_deg = modeldata.header.crval[5]
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
tmpfile = targetdir + '/.' + srcname + '.map.center.position.tmp'
writefile = open(tmpfile, 'w')
pickle.dump(position, writefile)
writefile.close()
