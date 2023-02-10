#!/usr/bin/env ParselTongue
################################################################################
# AIPS imports
################################################################################
from __future__ import print_function
from builtins import input
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
#import matplotlib
#matplotlib.use('Agg')
import matplotlib
#matplotlib.use('Qt5Agg', warn=False, force=True) ##non-interactive way, which has to be claimed before import plt
matplotlib.use('Qt5Agg', force=True) ##non-interactive way, which has to be claimed before import plt
import matplotlib.pyplot as plt
from scipy.special import jn
import sys, os, subprocess, math, datetime, glob
import interaction, pylab, ephem, astro_utils
import numpy as np
from interaction import yesno
from test_imagr import *

if __name__ == "__main__":
    main()

def main():
    AIPS.userno = 223
    uvdataset = AIPSUVData('BD186D1_I2','UVDATA',1,1)
    srcname = 'IBC18600024'
    widefieldimage(uvdataset, srcname, 256,0.75,False,0.05,0,0,0,100,20)

def aipsversion():
    try:
        aipsver = os.environ['PSRVLBAIPSVER']
    except KeyError:
        try:
            aipsver = os.environ['AIPS_VERSION'].split('/')[-1]
        except KeyError:
            aipsver = '31DEC22'
    return aipsver


##### Image a single-source file, optionally cleaning using autobox ############
def widefieldimage(uvdataset, srcname, numcells, cellmas, doclean, stopflux,
                   taperml, rashiftmas, decshiftmas, numcc, cleanradius):
    """
    Note
    ----
    1. Auto-box is turned on when "cleanradius" <= 0. In this case, a negative cleanradius can assign the number of search boxes (with imagr.nboxes = -cleanradius). The auto-box functionality is required for the mapping of multi-component sources, where components are separated by ~100 mas.
    """
    aipsver = aipsversion() ## Set up AIPS version
    imagr = AIPSTask('imagr', version = aipsver)
    imagr.indata = uvdataset
    imagr.nfield = 1
    imagr.sources[1] = srcname
    imagr.stokes = 'I'
    imagr.outname = srcname
    imagr.flux = stopflux
    #imagr.outseq = 1 ## this is the problem!!
    imagr.flagver = 1
    imagr.docal = -1
    imagr.doband = -1
    imagr.do3dimag = 1
    imagr.rashift[1] = rashiftmas/1000.0
    imagr.decshift[1] = decshiftmas/1000.0
    try:
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth[0]/\
                uvdataset.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth/\
                uvdataset.table('FQ', 1)[0].ch_width
    imagr.bchan = 0
    imagr.echan = 0
    imagr.bif = 0
    imagr.eif = 0
    imagr.nchav = nchan
    imagr.chinc = 1
    imagr.cellsize[1] = cellmas/1000.0
    imagr.cellsize[2] = cellmas/1000.0
    imagr.imsize[1] = numcells
    imagr.imsize[2] = numcells
    imagr.uvtaper[1] = taperml*1000.0
    imagr.uvtaper[2] = taperml*1000.0
    imagr.uvrange[1] = 0
    imagr.uvrange[2] = taperml*1000.0*3
    imagr.robust = 5
    imagr.gain = 0.1
    imagr.minpatch = 51
    imagr.uvwtfn = 'N'
    imagr.nboxes = 0
    
    if doclean:
        imagr.niter = numcc
        if cleanradius > 0:
            imagr.nboxes = 1
            imagr.clbox[1][1] = -1
            imagr.clbox[1][2] = int(cleanradius/cellmas)
            imagr.clbox[1][3] = numcells/2
            imagr.clbox[1][4] = numcells/2
        else:
            imagr.nboxes= 5
            if cleanradius < 0:
                imagr.nboxes = -cleanradius
            imagr.im2parm[1:]  = [0]
            imagr.im2parm[1:6] = [1, 6.5, 9, 0.005, 0] ## use auto boxes
    else:
        imagr.niter = 0
    imagr.dotv = -1
    imagr()

