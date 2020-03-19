#!/usr/bin/env python
#####################################################################################
## compile stats files (for phscal and primary IBC) to get systematics
## tailored for MSPSRPI datasets
## by Hao Ding on 5 May 2019
#####################################################################################
import os,glob,sys,yaml,math
import howfun,mspsrpifun
import numpy as np
#import matplotlib.pyplot as plt
#from scipy.stats import norm
#from astropy.time import Time
from optparse import OptionParser
np.set_printoptions(suppress=True)
## optparse setup ####################################################
usage = "usage: %prog []\n-t --target\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-t", "--target", dest="target", default="",
                   help="target on which we will gather and parse statsfiles")
(options, junk) = parser.parse_args()
targetname      = options.target
## path setup ########################################################
auxdir    = os.environ['PSRVLBAUXDIR']
configdir = auxdir + '/configs/'
expconfigfile = configdir + targetname + '.yaml'
targetdir = auxdir + '/processing/' + targetname
if not os.path.exists(targetdir):
    print("%s doesn't exist; aborting\n" % targetdir)
    sys.exit()
## parse sourcenames  ###############################################
expconfig = yaml.load(open(expconfigfile))
prIBCname = expconfig['primaryinbeam']
phscalname = mspsrpifun.target2cals(targetname)[0]
print phscalname, prIBCname
## find and parse statsfiles, reach RAs/Decs #########################
def statsfiles2positions(statsfiles):
    statsfiles.sort()
    print statsfiles
    RAs = np.array([])
    Decs = np.array([])
    for statsfile in statsfiles:
        lines = open(statsfile).readlines()
        for line in lines:
            if 'Actual RA' in line:
                temp = line.split('RA:')[-1].strip()
                RAs = np.append(RAs, temp)
            if 'Actual Dec' in line:
                temp = line.split('Dec:')[-1].strip()
                Decs = np.append(Decs, temp)
    return RAs, Decs
phscalstatsfiles = glob.glob(r'%s/*/*_ibshiftdiv_difmap.jmfit.stats' % targetdir) #find statsfile for each epoch
[phscalRAs, phscalDecs] = statsfiles2positions(phscalstatsfiles)
print phscalRAs, phscalDecs
prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.divided.difmap.jmfit.stats' % targetdir)
[prIBC_RAs, prIBC_Decs] = statsfiles2positions(prIBCstatsfiles)
print prIBC_RAs, prIBC_Decs
## statistics about RAs/Decs ##########################################
def dms_positions2stat(RAs, Decs):
    [RA_average, RA_std_mas] = dms_array2stat(RAs)
    [Dec_average, Dec_std_mas] = dms_array2stat(Decs)
    Dec_av_rad = Dec_average*math.pi/180
    RA_std_mas = RA_std_mas*15*math.cos(Dec_av_rad)
    RA_average_dms = howfun.deg2dms(RA_average)
    Dec_average_dms = howfun.deg2dms(Dec_average)
    return [RA_average_dms, RA_std_mas], [Dec_average_dms, Dec_std_mas]
def dms_array2stat(array):
    if len(array) < 2:
        print("len(list)<2, thus unable to carry out statistics; abort")
        sys.exit()
    degrees = howfun.dms2deg(array)
    average = np.average(degrees)
    std_degree = np.std(degrees)
    std_mas = std_degree*3600*1000
    return average, std_mas
[phscalstatsRAs, phscalstatsDecs] = dms_positions2stat(phscalRAs, phscalDecs)
print phscalstatsRAs, phscalstatsDecs
[prIBCstatsRAs, prIBCstatsDecs] = dms_positions2stat(prIBC_RAs, prIBC_Decs)
print prIBCstatsRAs, prIBCstatsDecs
