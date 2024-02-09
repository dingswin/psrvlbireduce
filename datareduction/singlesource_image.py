#!/usr/bin/env ParselTongue
## written by Adam Deller

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
import sys, os, math
import vlbatasks
from optparse import OptionParser

################################################################################
# Option parsing and defaulted global variables
################################################################################
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--fitsfile", dest="fitsfile", default="",
                  help="Fits file to bisect")
parser.add_option("-s", "--sourcename", dest="sourcename", default="",
                  help="Name of the source (blank will try to guess from filename")
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Name of the experiment (blank will try to guess from filename")
parser.add_option("--endif", dest="endif", default=-1,
                  help="End IF number (default -1 does only combined")
parser.add_option("--pixelmas", dest="pixelmas", default="0.75",
                  help="Pixel size in milliarcseconds")
parser.add_option("--pixelwindow", dest="pixelwindow", default="20",
                  help="Number of pixels to fit")
parser.add_option("--imagesize", dest="imagesize", default="2048",
                  help="Size of the initial difmap image in pixels")
parser.add_option("--finalimagesize", dest="finalimagesize", default="256",
                  help="Size of the final difmap image in pixels")
parser.add_option("--weightstring", dest="weightstring", default="0,-1",
                  help="Difmap weight string to use (default 0,-2)")

(options, junk) = parser.parse_args()
aipsver         = '31DEC20'
fitsfile        = options.fitsfile
sourcename      = options.sourcename
experiment      = options.experiment
endif           = int(options.endif)
pixelmas        = float(options.pixelmas)
pixelwindow     = int(options.pixelwindow)
imagesize       = int(options.imagesize)
finalimagesize  = int(options.finalimagesize)
weightstring    = options.weightstring
AIPS.userno     = 2
beginif         = 1
prefix          = os.getcwd() + '/'

if fitsfile == "":
    parser.error("You must supply a filename with -f or --fitsfile")
if sourcename == "":
    sourcename = fitsfile.split('/')[-1].split('_')[1]
if experiment == "":
    experiment = fitsfile.split('/')[-1].split('_')[0]
if fitsfile.rfind('/') >= 0:
    prefix = fitsfile[:fitsfile.rfind('/')] + '/'
else:
    fitsfile = prefix + fitsfile

os.system("rm -f " + os.getcwd() + "/templink.fits")
os.system("ln -s %s templink.fits" % fitsfile)
# Load the file
uvdata = AIPSUVData("JUNK", "JUNK", 1, 1)
if uvdata.exists():
    uvdata.zap()
vlbatasks.fitld_uvfits(os.getcwd() + "/templink.fits", uvdata, [])
try:
    fullmjd  = vlbatasks.get_dataset_mjd_midpoint(uvdata)
except ValueError:
    print("Couldn't get MJD")
    fullmjd = -1
uvdata.zap()

fullauto = True
stokesi = True
npixels = imagesize
gaussiantarget = False
beginif = 1
uvaverstr = '20,False'
uvtaperstr = '0.99,1000'
fullimagefile = fitsfile[:fitsfile.rfind('.')] + ".clean"
fulljmfitfile = fitsfile[:fitsfile.rfind('.')] + ".clean"

#vlbatasks.difmap_maptarget(fitsfile, fullimagefile, fullauto, stokesi,
#                           pixelmas, npixels, weightstring, uvaverstr, gaussiantarget,
#                           beginif, endif, "", finalimagesize)
#vlbatasks.jmfit(fullimagefile, fulljmfitfile, sourcename, stokesi, endif, pixelwindow, fullmjd)
vlbatasks.difmap_maptarget(fitsfile, fullimagefile, fullauto, stokesi,
                           pixelmas, npixels, weightstring, uvaverstr, uvtaperstr, gaussiantarget,
                           beginif, endif)
vlbatasks.jmfit(fullimagefile, fulljmfitfile, sourcename, stokesi, endif)

