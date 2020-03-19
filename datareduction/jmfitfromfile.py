#!/usr/bin/env ParselTongue
################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

################################################################################
# General imports
################################################################################
import os,sys,math

################################################################################
# Check invocation
################################################################################
usage =  "usage: jmfitfromfile.py <image fits file> <output stats file> [pixwindow | pixwindowx,pixwindowy | blc,blc,trc,trc | blc,blc,blc,trc,trc,trc]]\n"
if not len(sys.argv) == 3 and not len(sys.argv) == 4:
    print usage
    sys.exit()

################################################################################
# Set variables
################################################################################
AIPS.userno = 123
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    aipsver = '31DEC18'
imagefile = sys.argv[1]
#if len(imagefile) > 10:
#    imagefile = imagefile[-10:]
jmfitfile = sys.argv[2]
pixwindowx = 50
pixwindowy = 50
blc=[None]
trc=[None]
imageslice = -1
if len(sys.argv) == 4:
    pixstring = sys.argv[3].split(',')
    pixwindowx = int(pixstring[0])
    pixwindowy = int(pixstring[0])
    if len(pixstring) > 1:
        if len(pixstring) == 6:
            blc.append(int(pixstring[0]))
            blc.append(int(pixstring[1]))
            blc.append(int(pixstring[2]))
            trc.append(int(pixstring[3]))
            trc.append(int(pixstring[4]))
            trc.append(int(pixstring[5]))
            pixwindowx = trc[1] - blc[1]
            pixwindowy = trc[2] - blc[2]
            imageslice = int(pixstring[2])
            if imageslice != int(pixstring[5]):
                print "Can't JMFIT across multiple slices!"
                sys.exit()
        elif len(pixstring) == 4:
            blc.append(int(pixstring[0]))
            blc.append(int(pixstring[1]))
            trc.append(int(pixstring[2]))
            trc.append(int(pixstring[3]))
            pixwindowx = trc[1] - blc[1]
            pixwindowy = trc[2] - blc[2]
        elif len(pixstring) == 2:
            pixwindowy = int(pixstring[1])
        else:
            print usage, "len(pixstring)=", len(pixstring)
            sys.exit()

if not imagefile[0] == '/':
    imagefile = os.getcwd() + '/' + imagefile
if not jmfitfile[0] == '/':
    jmfitfile = os.getcwd() + '/' + jmfitfile

################################################################################
# Run JMFIT
################################################################################
print "JMFIT'ing over a box of width " + str(pixwindowx) + " x " + str(pixwindowy)
wrong_imagedata = AIPSUVData('JUNK', 'IMG', 1, 1)
if wrong_imagedata.exists():
    wrong_imagedata.zap()
imagedata = AIPSImage('JUNK', 'IMG', 1, 1)
outdata = AIPSImage('JUNK', 'IMG', 1, 2)
targetstatout = open(jmfitfile, "w")
if imagedata.exists():
    imagedata.zap()
if outdata.exists():
    outdata.zap()

fitld = AIPSTask('fitld', version = aipsver)
fitld.datain = imagefile
fitld.outdata = imagedata
fitld.dotable = 1
fitld.douvcomp = -1
fitld()

obsdate = imagedata.header.date_obs
year = int(obsdate[0:4])
month = int(obsdate[5:7])
day = int(obsdate[8:10])
mjd = year*367 - int(7*(year + int((month + 9)/12))/4) + \
      int(275*month/9) + day - 678987

# Get clean beam size if needed
if imageslice >= 0:
    uktable = imagedata.table("UK",1)
    cleanbeammaj = cleanbeammin = cleanbeampa = -1
    for row in uktable:
        if row.chan == imageslice:
            cleanbeammaj = row.bmaj / math.fabs(imagedata.header.cdelt[1]*3600) # convert from beamsize in arcseconds to pixels
            cleanbeammin = row.bmin / math.fabs(imagedata.header.cdelt[1]*3600) # convert from beamsize in arcseconds to pixels
            cleanbeampa  = row.bpa
            print "Clean beam details:", cleanbeammaj, cleanbeammin, cleanbeampa, row
            break
    if cleanbeammaj < 0:
        print "Couldn't find clean beam details!"
        sys.exit()

if len(blc) == 1:
    imean = AIPSTask("imean")
    imean.indata = imagedata
    imean.go()
    lines = imean.message()
    xpix = 128
    ypix = 128
    for line in lines:
        if "Maximum" in line:
            splitline = line.split()
            xpix = float(splitline[4])
            ypix = float(splitline[5])
    blc = [None,0,0]
    trc = [None,0,0]
    blc[1:] = [xpix-pixwindowx/2, ypix-pixwindowy/2]
    trc[1:] = [xpix+pixwindowx/2, ypix+pixwindowy/2]
    print blc
    print trc
    if blc[1] < 0: blc[1] = 1
    if blc[2] < 0: blc[2] = 1
    if trc[1] > int(imagedata.header.naxis[0]): trc[1] = int(imagedata.header.naxis[0])-1
    if trc[2] > int(imagedata.header.naxis[1]): trc[2] = int(imagedata.header.naxis[1])-1

jmfit = AIPSTask('jmfit')
jmfit.indata = imagedata
jmfit.outdata = outdata
jmfit.niter = 4000
jmfit.ngauss = 1
if imageslice >= 0:
    jmfit.gwidth[1][1] = cleanbeammaj
    jmfit.gwidth[1][2] = cleanbeammin
    jmfit.gwidth[1][3] = cleanbeampa
jmfit.fwidth[1][1] = 0
jmfit.fwidth[1][2] = 0
jmfit.fwidth[1][3] = 0
jmfit.domax[1:] = [1]
jmfit.dopos[1][1] = 1
jmfit.dopos[1][2] = 1
jmfit.dowidth[1][1] = 1
jmfit.dowidth[1][2] = 1
jmfit.dowidth[1][3] = 1
jmfit.blc = blc
jmfit.trc = trc
jmfit.go()
jmfitmessage = jmfit.message()
msgindex = 0
exciselinenos = []
for i in range(len(jmfitmessage)-1):
    if jmfitmessage[len(jmfitmessage)-(i+1)][:5] != 'JMFIT':
        exciselinenos.append(len(jmfitmessage)-(i+1))
for e in exciselinenos:
    jmfitmessage = jmfitmessage[:e] + jmfitmessage[e+1:]
while jmfitmessage[msgindex].find('X-ref') == -1 and \
          msgindex < len(jmfitmessage):
    msgindex = msgindex + 1
if msgindex >= len(jmfitmessage):
    print "Did not find solution for " + imagefile
    sys.exit(1)
centrerasplit = jmfitmessage[msgindex].split()
centredecsplit = jmfitmessage[msgindex+1].split()
rapixsize = math.fabs(float(centrerasplit[-1])*1000.0)
decpixsize = math.fabs(float(centredecsplit[-1])*1000.0)
pixsize = rapixsize
if rapixsize != decpixsize:
    print "*****WARNING*****: Different pixel sizes in RA and dec"
    print "*****WARNING*****: Beam size measurement will be approximate"
    pixsize = (rapixsize + decpixsize)/2.0
msgindex = msgindex + 2
while not ("Beam=" in jmfitmessage[msgindex] or "NO CLEAN BEAM" in jmfitmessage[msgindex]) and msgindex < len(jmfitmessage):
    msgindex = msgindex + 1
if msgindex >= len(jmfitmessage):
    print "Did not find solution for " + imagefile
    sys.exit(1)
if "NO CLEAN BEAM" in jmfitmessage[msgindex]:
    print "There was no clean beam (probably you're JMFIT'ing from a cube).  Ignore all the deconvolved params!"
    beammaj = -1
    beammin = -1
    beampa  = -1
    nocleanbeam = True
else:
    beamsplit = jmfitmessage[msgindex].split()
    if not len(beamsplit) == 10:
         print "*****WARNING*****: Poorly formatted beam string!"
    beammaj = float(beamsplit[2])*pixsize
    beammin = float(beamsplit[4])*pixsize
    beampa  = float(beamsplit[8])
    nocleanbeam = False
msgindex += 1
while jmfitmessage[msgindex].find('********* Solution from JMFIT') == -1 \
          and msgindex < len(jmfitmessage):
    msgindex = msgindex + 1
if msgindex >= len(jmfitmessage):
    print "Did not find solution for " + imagefile
    sys.exit(1)
fluxsplit = jmfitmessage[msgindex+3].split()
if nocleanbeam:
    msgindex -= 1
sourcerasplit = jmfitmessage[msgindex+7].split()
sourcedecsplit = jmfitmessage[msgindex+8].split()
fitmajsplit = jmfitmessage[msgindex+12].split()
fitminsplit = jmfitmessage[msgindex+13].split()
fitpasplit = jmfitmessage[msgindex+14].split()
try:
    fitmaj = float(fitmajsplit[4])
    fitmin = float(fitminsplit[4])
except ValueError:
    print "Bad value for one or more of fit major axis " +  fitmajsplit[4] + \
          ", minor axis " + fitminsplit[4] + ", both will be set to zero"
    fitmaj = 0.0
    fitmin = 0.0
targetstatout.write("Pulsar XXXXXXXX: MJD " + str(mjd) + ":\n")
targetstatout.write("Centre RA:            " + centrerasplit[-7] + ":" + \
                    centrerasplit[-6] + ":" + centrerasplit[-5] + "\n")
targetstatout.write("Centre Dec:           " + centredecsplit[-7] + ":" + \
                    centredecsplit[-6] + ":" + centredecsplit[-5] + "\n")
if fluxsplit[4] == '+/-':
    flux = float(fluxsplit[3][1:])*1000.0
    rms = float(fluxsplit[5])*1000.0
else:
    flux = float(fluxsplit[4])*1000.0
    rms = float(fluxsplit[6])*1000.0
targetstatout.write("Flux (mJy):           " + str(flux) + "\n")
targetstatout.write("S/N:                  " + str(flux/rms) + "\n")
centrerahours = float(centrerasplit[-7]) +  float(centrerasplit[-6])/60.0 \
                + float(centrerasplit[-5])/3600.0
centredecdegs = float(centredecsplit[-7]) +  float(centredecsplit[-6])/60.0 \
                + float(centredecsplit[-5])/3600.0
srcrahours = float(sourcerasplit[2]) +  float(sourcerasplit[3])/60.0 + \
             float(sourcerasplit[4])/3600.0
srcdecdegs = float(sourcedecsplit[2]) +  float(sourcedecsplit[3])/60.0 + \
             float(sourcedecsplit[4])/3600.0
raoff = (srcrahours-centrerahours)*15*60*60*1000*\
        math.cos(centredecdegs*math.pi/180.0)
decoff = (srcdecdegs-centredecdegs)*60*60*1000
targetstatout.write("RA offset (mas):      " + str(raoff) + "\n")
targetstatout.write("Dec offset (mas):     " + str(decoff) + "\n")
targetstatout.write("Actual RA:            " + sourcerasplit[2] + ":" + \
                    sourcerasplit[3] + ":" + sourcerasplit[4] + "\n")
targetstatout.write("Actual Dec:           " + sourcedecsplit[2] + ":" + \
                    sourcedecsplit[3] + ":" + sourcedecsplit[4] + "\n")
targetstatout.write("Fit:                  " + str(fitmin*1000) + "x" + \
                    str(fitmaj*1000) + " at " + fitpasplit[4] + \
                    " degrees; beam " + str(beammin) + "x" + \
                    str(beammaj) + " at " + str(beampa) + " degrees\n")
raerr = float(sourcerasplit[6])*1000
decerr = float(sourcedecsplit[6])*1000
targetstatout.write("Est. RA error (mas):  " + \
                    str(raerr*math.cos(centredecdegs*math.pi/180.0)*15.0) \
                    + "\n")
targetstatout.write("Est. RA error (hms):  " + str(raerr) + "\n")
targetstatout.write("Est. Dec error (mas): " + str(decerr) + "\n")
targetstatout.close()

