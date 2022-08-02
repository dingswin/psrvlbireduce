#!/usr/bin/env ParselTongue
###############################################################################
# purpose: get inbeamname from idifitsfile
###############################################################################
#import os,glob,sys,pickle,yaml
#import howfun,mspsrpifun
import vlbatasks, sys
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

usage = 'summarize*.py idifitsfile'
if len(sys.argv) != 2:
    print(usage)
    sys.exit()
uvdatafile = sys.argv[1]
AIPS.userno = 123

uvdata = AIPSUVData("JUNK","JUNK",1,1)
if uvdata.exists():
    uvdata.zap()
vlbatasks.fitld_uvfits(uvdatafile, uvdata, [])
sutable = uvdata.table('SU', 1)
print(sutable)
for row in sutable:
    print(row.source)
    #if 'IBC' in (row.source)
    #    inbeamname = row.source
    #    print('inbeamname is %s')
