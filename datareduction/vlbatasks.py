################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
#import matplotlib
#matplotlib.use('Agg')
from scipy.special import jn
import sys, os, subprocess, math, datetime, glob
import interaction, pylab, ephem, astro_utils
import numpy as np
from interaction import yesno

# Set up AIPS version
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    try:
        aipsver = os.environ['AIPS_VERSION'].split('/')[-1]
    except KeyError:
        aipsver = '31DEC20'

speedOfLight = 299792458.

vlbadiameter = 25.47415721 # metres
gbtdiameter = 100 # metres
vlbabeams = {}
# Squint(R-L) Az, Squint(R-L) EL, RCP Az beamwidth, RCP EL beamwidth, LCP Az beamwidth, LCP El beamwidth, reffreq
vlbabeams['SC'] = [-1.53, -0.56, 29.58, 29.57, 29.61, 29.36, 1438.0e6]
vlbabeams['HN'] = [-1.35, -0.67, 29.17, 29.21, 29.26, 29.09, 1438.0e6]
vlbabeams['NL'] = [-1.56, -0.62, 29.44, 29.55, 29.51, 29.47, 1438.0e6]
vlbabeams['FD'] = [-1.56, -0.58, 29.27, 29.24, 29.41, 29.27, 1438.0e6]
vlbabeams['LA'] = [-1.62, -0.67, 29.12, 29.03, 29.20, 29.00, 1438.0e6]
vlbabeams['PT'] = [-1.59, -0.63, 29.27, 29.18, 29.28, 29.14, 1438.0e6]
vlbabeams['KP'] = [-1.61, -0.59, 29.30, 29.32, 29.19, 29.04, 1438.0e6]
vlbabeams['OV'] = [-1.76, -0.95, 29.57, 29.64, 29.67, 29.62, 1438.0e6]
vlbabeams['BR'] = [-1.44, -0.50, 29.19, 29.36, 29.55, 29.44, 1438.0e6]
vlbabeams['MK'] = [-1.56, -0.60, 29.20, 29.26, 29.24, 29.25, 1438.0e6]
vlbabeams['GB'] = [0, 0, 7.3, 7.3, 7.3, 7.3, 1438.0e6] # Total guess, based on scaling from VLBA
vlbalocations = {}
# Latitude, longitude, elevation(m)
vlbalocations['SC'] = ['17:45:23.68', '-64:35:01.07',  16  ]
vlbalocations['HN'] = ['42:56:00.99', '-71:59:11.69',  296 ]
vlbalocations['NL'] = ['41:46:17.13', '-91:34:26.88',  222 ]
vlbalocations['FD'] = ['30:38:06.11', '-103:56:41.34', 1606]
vlbalocations['LA'] = ['35:46:30.45', '-106:14:44.15', 1962]
vlbalocations['PT'] = ['34:18:03.61', '-108:07:09.06', 2365]
vlbalocations['KP'] = ['31:57:22.70', '-111:36:44.72', 1902]
vlbalocations['OV'] = ['37:13:53.95', '-118:16:37.37', 1196]
vlbalocations['BR'] = ['48:07:52.42', '-119:40:59.80', 250 ]
vlbalocations['MK'] = ['19:48:04.97', '-155:27:19.81', 3763]
vlbalocations['GB'] = ['38:26:16.16', '-100:09:51.20', 844]

################################################################################
# Container classes
################################################################################
class VexScan:
    def __init__(self, scanname, starttime, stoptime, source, modename, vlaapdet):
        self.scanname  = scanname
        self.starttime = starttime
        self.stoptime  = stoptime
        self.source    = source
        self.modename  = modename
        self.startmjd  = getVexMJD(starttime)
        self.stopmjd   = getVexMJD(stoptime)
        self.vlaautophasedetermine = vlaapdet

    def incStopTime(self, incsecs):
        refmjd = int(self.stopmjd)
        doy = self.getStartDOY()
        self.stopmjd += incsecs/86400.0
        fracmjd = self.stopmjd-refmjd
        if fracmjd > 1:
            doy += 1
            fracmjd -= 1
        hh = int(fracmjd*24.0)
        mm = int(fracmjd*1440.0 - 60*hh)
        ss = int(fracmjd*86400.0 - (3600*hh + 60*mm))
        self.stoptime = "%s%03dd%02dh%02dm%02ds" % (self.stoptime[:5], doy, hh, mm, ss)

    def setStopTime(self, stoptime):
        self.stoptime = stoptime
        self.stopmjd = getVexMJD(stoptime)

    def getStartDOY(self):
        return int(self.starttime[5:8])

    def getStartHour(self):
        return int(self.starttime[9:11])

    def getStartMinute(self):
        return int(self.starttime[12:14])

    def getStartSecond(self):
        return int(self.starttime[15:17])

    def getStopDOY(self):
        return int(self.stoptime[5:8])

    def getStopHour(self):
        return int(self.stoptime[9:11])

    def getStopMinute(self):
        return int(self.stoptime[12:14])

    def getStopSecond(self):
        return int(self.stoptime[15:17])

class Telescope:
    def __init__(self, name, x, y, z, axisoff):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.axisoff = axisoff

class Source:
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra
        self.dec = dec

    def getSepArcmin(self, rarad, decrad):
        sinsqdecdiff = math.sin((decrad-self.dec)/2.0)
        sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
        sinsqradiff  = math.sin((rarad-self.ra)/2.0)
        sinsqradiff  = sinsqradiff*sinsqradiff

        return 180.0*60.0*2*math.asin(math.sqrt(sinsqdecdiff +
                           math.cos(decrad)*math.cos(self.dec)*sinsqradiff))/math.pi

class Scan:
    def __init__(self, imlines, calclines, ilinenum, clinenum, numantennas,
                 polyorder, polyinterval, isdifx15, sources, startmjd, startsec,
                 compatibilityMode):
        self.startmjd = startmjd
        self.startsec = startsec
        self.polyorder = polyorder
        self.polyinterval = polyinterval
        self.polymjd = []
        self.polysec = []
        self.delays  = []
        self.u       = []
        self.v       = []
        self.w       = []
        self.source = []

        if compatibilityMode: #"im" is really delay, "calc" is really uvw.
                              #Initialise accordingly
            self.initialiseCompatibility(imlines, calclines, ilinenum, clinenum, 
                                         numantennas, polyorder, polyinterval, 
                                         sources, startmjd, startsec)
            return

        # Do the calc side first
        self.scanstart = int(calclines[clinenum+1].split(':')[1].strip())
        if isdifx15:
            self.scandur = int(calclines[clinenum+0].split(':')[1].strip())
            name = calclines[clinenum+3].split(':')[1].strip()
            ra = float(calclines[clinenum+4].split(':')[1])
            dec = float(calclines[clinenum+5].split(':')[1])
            self.source.append(Source(name, ra, dec))
            clinenum += 8
            self.numphasecentres = 0
        else:
            self.scandur = int(calclines[clinenum+2].split(':')[1])
            self.numphasecentres = int(calclines[clinenum+7].split(':')[1])
            self.source.append(sources[int(calclines[clinenum+6].split(':')[1])])
            for i in range(self.numphasecentres):
                self.source.append(sources[int(calclines[clinenum+8+i].split(':')[1])])
            clinenum += self.numphasecentres+8

        self.startsec += self.scanstart

        # Then then im side
        if isdifx15:
            self.numpoly = int(imlines[ilinenum+1].split(':')[1])
            ilinenum += 2
        else:
            self.numpoly = int(imlines[ilinenum+self.numphasecentres+2].split(':')[1])
            ilinenum += self.numphasecentres + 3
        for i in range(self.numpoly):
            self.polymjd.append(int(imlines[ilinenum].split(':')[1]))
            ilinenum += 1
            self.polysec.append(int(imlines[ilinenum].split(':')[1]))
            ilinenum += 1
            self.delays.append([])
            self.u.append([])
            self.v.append([])
            self.w.append([])
            for j in range(self.numphasecentres+1):
                self.delays[i].append([])
                self.u[i].append([])
                self.v[i].append([])
                self.w[i].append([])
                for k in range(numantennas):
                    self.delays[i][j].append([])
                    self.u[i][j].append([])
                    self.v[i][j].append([])
                    self.w[i][j].append([])
                    offsettou = 3
                    while not "U (m)" in imlines[ilinenum+offsettou]:
                        offsettou += 1
                    d = (imlines[ilinenum].split(':')[1]).split()
                    u = (imlines[ilinenum+offsettou].split(':')[1]).split()
                    v = (imlines[ilinenum+offsettou+1].split(':')[1]).split()
                    w = (imlines[ilinenum+offsettou+2].split(':')[1]).split()
                    for l in range(polyorder+1):
                        self.delays[i][j][k].append(float(d[l]))
                        self.u[i][j][k].append(float(u[l]))
                        self.v[i][j][k].append(float(v[l]))
                        self.w[i][j][k].append(float(w[l]))
                    ilinenum += 3 + offsettou

        # Finally store the end line numbers
        self.endclinenum = clinenum
        self.endilinenum = ilinenum

    def initialiseCompatibility(self, delaylines, uvwlines, dlinenum, ulinenum,
                                numantennas, polyorder, polyinterval,
                                sources, startmjd, startsec):
        self.scanstart = int(delaylines[dlinenum+1].split(':')[1].strip())
        self.scandur   = int(delaylines[dlinenum+0].split(':')[1].strip())
        name = uvwlines[ulinenum+2].split(':')[1].strip()
        ra   = float(uvwlines[ulinenum+3].split(':')[1])
        dec  = float(uvwlines[ulinenum+4].split(':')[1])
        self.source.append(Source(name, ra, dec))
        dlinenum += 4
        ulinenum += 6
        self.numphasecentres = 0
        self.startsec += self.scanstart
        self.numpoly = self.scandur

        for i in range(self.numpoly):
            lastdelays    = delaylines[i + dlinenum - 1][20:].split()
            currentdelays = delaylines[i + dlinenum][20:].split()
            nextdelays    = delaylines[i + dlinenum + 1][20:].split()
            lastuvws      = uvwlines[i + ulinenum - 1][20:].split()
            currentuvws   = uvwlines[i + ulinenum][20:].split()
            nextuvws      = uvwlines[i + ulinenum + 1][20:].split()
            self.polymjd.append(self.startmjd)
            self.polysec.append(self.startsec + i)
            self.delays.append([])
            self.u.append([])
            self.v.append([])
            self.w.append([])
            for j in range(self.numphasecentres+1):
                self.delays[i].append([])
                self.u[i].append([])
                self.v[i].append([])
                self.w[i].append([])
                for k in range(numantennas):
                    d0 = float(lastdelays[k])
                    d1 = float(currentdelays[k])
                    d2 = float(nextdelays[k])
                    u0 = float(lastuvws[k*3+0])
                    u1 = float(currentuvws[k*3+0])
                    u2 = float(nextuvws[k*3+0])
                    v0 = float(lastuvws[k*3+1])
                    v1 = float(currentuvws[k*3+1])
                    v2 = float(nextuvws[k*3+1])
                    w0 = float(lastuvws[k*3+2])
                    w1 = float(currentuvws[k*3+2])
                    w2 = float(nextuvws[k*3+2])
                    self.delays[i][j].append([])
                    self.u[i][j].append([])
                    self.v[i][j].append([])
                    self.w[i][j].append([])
                    self.delays[i][j][k].append(d1)
                    self.delays[i][j][k].append(d2/2.0 - d0/2.0)
                    self.delays[i][j][k].append(d2/2.0 + d0/2.0 - d1)
                    #if i == 0 or i == self.numpoly-1:
                    #    print self.delays[i][j][k]
                    self.u[i][j][k].append(u1)
                    self.u[i][j][k].append(u2/2.0 - u0/2.0)
                    self.u[i][j][k].append(u2/2.0 + u0/2.0 - u1)
                    self.v[i][j][k].append(v1)
                    self.v[i][j][k].append(v2/2.0 - v0/2.0)
                    self.v[i][j][k].append(v2/2.0 + v0/2.0 - v1)
                    self.w[i][j][k].append(w1)
                    self.w[i][j][k].append(w2/2.0 - w0/2.0)
                    self.w[i][j][k].append(w2/2.0 + w0/2.0 - w1)

        # Finally store the end line numbers
        self.enddlinenum = dlinenum + self.scandur + 2
        self.endulinenum = ulinenum + self.scandur + 2

    def getSourceIndex(self, sourcename):
        for i in range(self.numphasecentres+1):
            #print "Comparing *" + self.source[i].name.strip() + "* with *" + sourcename.strip() + "*"
            #print self.source[i].name.strip() == sourcename.strip()
            if self.source[i].name.strip() == sourcename.strip():
                return i
        if self.source[i].name.strip() == "SCAN_GAP":
            return -1
        else:
            print "Couldn't find source " + sourcename.strip() + " - my pointing centre was " + self.source[0].name.strip()
            return -999

    def containsTime(self, mjd, sec):
#        print "my mjd and sec is %d, %d, being asked for %d, %f" % (self.startmjd, self.startsec, mjd, sec)
        offset = (mjd-self.startmjd)*86400 + sec - self.startsec
        if offset >= 0 and offset < self.scandur + 0.1: #Crappy rpfits time stamps!
            return True
        return False

    def getPolyEntryAndOffset(self, mjd, second):
        polyentry = -1
        for i in range(self.numpoly):
            polyoffset = (mjd - self.polymjd[i])*86400 + second - self.polysec[i]
            if polyoffset >= 0 and polyoffset <= self.polyinterval:
                polyentry = i
                break
        if polyentry < 0:
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan - aborting!"
            polyoffset = 0
        return polyentry, polyoffset

    def getUVW(self, mjd, second, antenna1index, antenna2index, sourceindex):
        toreturn = [0.0, 0.0, 0.0]
        if sourceindex > self.numphasecentres:
            print "Trying to get source index " + str(sourceindex) + " when this " + \
                  "scan only has " + str(self.numphasecentres) + " phase centres"
            return toreturn
        polyentry, polyoffset = self.getPolyEntryAndOffset(mjd, second)
        if polyentry < 0:
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan"
            return toreturn
        xval = 1.0
        for i in range(self.polyorder+1):
            toreturn[0] += self.u[polyentry][sourceindex][antenna1index][i]*xval
            toreturn[1] += self.v[polyentry][sourceindex][antenna1index][i]*xval
            toreturn[2] += self.w[polyentry][sourceindex][antenna1index][i]*xval
            xval *= polyoffset
        xval = 1.0
        for i in range(self.polyorder+1):
            toreturn[0] -= self.u[polyentry][sourceindex][antenna2index][i]*xval
            toreturn[1] -= self.v[polyentry][sourceindex][antenna2index][i]*xval
            toreturn[2] -= self.w[polyentry][sourceindex][antenna2index][i]*xval
            xval *= polyoffset
        if math.fabs(toreturn[0]) < 1e-6 and math.fabs(toreturn[1]) < 1e-6 and math.fabs(toreturn[2]) < 1e-6:
            toreturn = [0.0, 0.0, 0.0]
        return toreturn

    def getDelayAndRate(self, mjd, second, antennaindex, sourceindex):
        if sourceindex > self.numphasecentres:
            print "Trying to get source index " + str(sourceindex) + " when this " + \
                  "scan only has " + str(self.numphasecentres) + " phase centres"
            return 0.0, 0.0
        polyentry, polyoffset = self.getPolyEntryAndOffset(mjd, second)
        if polyentry < 0:
            print "Couldn't find a poly entry for time " + str(mjd) + ", sec " + \
                  str(second) + " in this scan"
            return 0.0, 0.0
        correction = 0.0
        delay = 0.0
        rate = 0.0
        xval = 1.0
        #print "Polyentry is " + str(polyentry) + ", polyoffset is " + str(polyoffset)
        #print self.delays[polyentry][sourceindex][antennaindex]
        for i in range(self.polyorder+1):
            delay += self.delays[polyentry][sourceindex][antennaindex][i]*xval
            rate += i*self.delays[polyentry][sourceindex][antennaindex][i]*xval/polyoffset
            xval *= polyoffset
        delay += delay*rate/1e6
        #print "Calculated delay is %10.6f, calculated rate is %10.6f" % (delay, rate)
        return delay/1000000.0, rate/1000000.0 #microsec -> seconds and us/sec -> sec/sec

class Model:
    def __init__(self, calcfile, delayfile=None, uvwfile=None):
        if calcfile == None:
            if not os.path.exists(delayfile):
                print "Delay file (compatibility mode) does not exist- aborting!"
                sys.exit()
            if not os.path.exists(uvwfile):
                print "UVW file (compatibility mode) does not exist- aborting!"
                sys.exit()
            self.initialiseCompatibility(delayfile, uvwfile)
            return
        if not os.path.exists(calcfile):
            print "Calc file " + calcfile + " does not exist - aborting"
            sys.exit()
        self.compatibilityMode = False
        calcin = open(calcfile)
        self.filename = calcfile
        calclines = calcin.readlines()
        calcin.close()
        clinenum = 0
        clinenum, val = getToLine("DIFX VERSION", calclines, clinenum)
        difxverstr = val.strip()
        clinenum, val = getToLine("IM FILENAME", calclines, clinenum)
        imfilename = val.strip()
        self.isDifx15 = False
        if "1.5" in difxverstr:
            self.isDifx15 = True
        imin = open(imfilename)
        imlines = imin.readlines()
        imin.close()
        ilinenum = 0
        ilinenum, yearstr = getToLine("START YEAR", imlines, ilinenum)
        ilinenum, monthstr = getToLine("START MONTH", imlines, ilinenum)
        ilinenum, daystr = getToLine("START DAY", imlines, ilinenum)
        ilinenum, hourstr = getToLine("START HOUR", imlines, ilinenum)
        ilinenum, minstr = getToLine("START MINUTE", imlines, ilinenum)
        ilinenum, secstr = getToLine("START SECOND", imlines, ilinenum)
        self.startmjd = ymd2mjd(int(yearstr), int(monthstr), int(daystr))
        self.startsec = int(hourstr)*3600 + int(minstr)*60 + int(secstr)
        ilinenum, polyorderstr = getToLine("POLYNOMIAL ORDER", imlines, ilinenum)
        self.polyorder = int(polyorderstr)
        ilinenum, polyintervalstr = getToLine("INTERVAL (SECS)", imlines, ilinenum)
        self.polyinterval = int(polyintervalstr)
        ilinenum, numtelstr = getToLine("NUM TELESCOPES", imlines, ilinenum)
        self.numantennas = int(numtelstr)
        clinenum = 0
        clinenum, numtelstr = getToLine("NUM TELESCOPES", calclines, clinenum)
        if self.numantennas != int(numtelstr):
            print "Number of telescopes doesn't match between calc and IM files!"
            sys.exit()
        self.antennas = []
        self.antennamap = {}
        for i in range(self.numantennas):
            clinenum, val = getToLine("TELESCOPE %d NAME" % i, calclines, clinenum)
            antennaname = val.strip()
            clinenum, val = getToLine("TELESCOPE %d OFFSET (m)" % i, calclines, clinenum)
            axisoffset = float(val)
            clinenum, val = getToLine("TELESCOPE %d X (m)" % i, calclines, clinenum)
            x = float(val.strip())
            clinenum, val = getToLine("TELESCOPE %d Y (m)" % i, calclines, clinenum)
            y = float(val.strip())
            clinenum, val = getToLine("TELESCOPE %d Z (m)" % i, calclines, clinenum)
            z = float(val.strip())
            self.antennas.append(Telescope(antennaname, x, y, z, axisoffset))
            self.antennamap[antennaname] = i
        self.sources = []
        if not self.isDifx15:
            clinenum, val = getToLine("NUM SOURCES", calclines, clinenum)
            self.numsources = int(val)
            for i in range(self.numsources):
                clinenum, val = getToLine("SOURCE %d NAME" % i, calclines, clinenum)
                name = val.strip()
                clinenum, val  = getToLine("SOURCE %d RA" % i, calclines, clinenum)
                ra = float(val.strip())
                clinenum, val  = getToLine("SOURCE %d DEC" % i, calclines, clinenum)
                dec = float(val.strip())
                self.sources.append(Source(name, ra, dec))
        clinenum, val = getToLine("NUM SCANS", calclines, clinenum)
        self.numscans = int(val)
        ilinenum, val = getToLine("NUM SCANS", imlines, ilinenum)
        if not self.numscans == int(val):
            print "Error - .im and .calc file disagree on number of scans!"
            sys.exit(1)
        clinenum += 1
        ilinenum += 1
        self.scans = []
        for j in range(self.numscans):
            toadd = Scan(imlines, calclines, ilinenum, clinenum,
                         self.numantennas, self.polyorder, self.polyinterval,
                         self.isDifx15, self.sources, self.startmjd, self.startsec,
                         self.compatibilityMode)
            self.scans.append(toadd)
            ilinenum = toadd.endilinenum
            clinenum = toadd.endclinenum
        self.lastscannum = 0
        self.changed = True

    def initialiseCompatibility(self, delayfile, uvwfile):
        self.isDifx15 = False
        self.compatibilityMode = True
        delayin = open(delayfile)
        delaylines = delayin.readlines()
        delayin.close()
        uvwin = open(uvwfile)
        uvwlines = uvwin.readlines()
        uvwin.close()
        dlinenum = 0
        ulinenum = 0
        dlinenum, yearstr = getToLine("START YEAR", delaylines, dlinenum)
        dlinenum, monthstr = getToLine("START MONTH", delaylines, dlinenum)
        dlinenum, daystr = getToLine("START DAY", delaylines, dlinenum)
        dlinenum, hourstr = getToLine("START HOUR", delaylines, dlinenum)
        dlinenum, minstr = getToLine("START MINUTE", delaylines, dlinenum)
        dlinenum, secstr = getToLine("START SECOND", delaylines, dlinenum)
        self.startmjd = ymd2mjd(int(yearstr), int(monthstr), int(daystr))
        self.startsec = int(hourstr)*3600 + int(minstr)*60 + int(secstr)
        self.polyorder = 2
        self.polyinterval = 1
        ulinenum, numtelstr = getToLine("NUM TELESCOPES", uvwlines, ulinenum)
        self.numantennas = int(numtelstr)
        self.antennas = []
        self.antennamap = {}
        for i in range(self.numantennas):
            ulinenum, val = getToLine("TELESCOPE %d NAME" % i, uvwlines, ulinenum)
            antennaname = val.strip()
            ulinenum, val = getToLine("TELESCOPE %d X (m)" % i, uvwlines, ulinenum)
            x = float(val.strip())
            ulinenum, val = getToLine("TELESCOPE %d Y (m)" % i, uvwlines, ulinenum)
            y = float(val.strip())
            ulinenum, val = getToLine("TELESCOPE %d Z (m)" % i, uvwlines, ulinenum)
            z = float(val.strip())
            self.antennas.append(Telescope(antennaname, x, y, z, 0.0)) #Lucky the axis offset isn't needed
            self.antennamap[antennaname] = i
        self.sources = []
        dlinenum, val = getToLine("NUM SCANS", delaylines, dlinenum)
        self.numscans = int(val)
        ulinenum, val = getToLine("NUM SCANS", uvwlines, ulinenum)
        if not self.numscans == int(val):
            print "Error - .delay and .uvw file disagree on number of scans!"
            sys.exit(1)
        dlinenum += 1
        ulinenum += 1
        self.scans = []
        for j in range(self.numscans):
            toadd = Scan(delaylines, uvwlines, dlinenum, ulinenum,
                         self.numantennas, self.polyorder, self.polyinterval,
                         self.isDifx15, self.sources, self.startmjd, self.startsec,
                         self.compatibilityMode)
            self.scans.append(toadd)
            dlinenum = toadd.enddlinenum
            ulinenum = toadd.endulinenum
        self.lastscannum = 0
        self.changed = True

    def getDelayAndRate(self, mjd, second, antennaname, sourcename):
        while self.lastscannum < self.numscans and \
              not self.scans[self.lastscannum].containsTime(mjd, second):
            self.lastscannum += 1
            self.changed = True
        if self.lastscannum == self.numscans:
            self.lastscannum = 0
            while self.lastscannum < self.numscans and \
                  not self.scans[self.lastscannum].containsTime(mjd, second):
                self.lastscannum += 1
                self.changed = True
            if self.lastscannum == self.numscans:
                print "Couldn't find scan for time %d/%f" % (mjd, second)
                return -9e99,-9e99
        if self.changed:
            #if self.compatibilityMode:
            #    print "Running in compatibility mode"
            #print "Looking for a source at time " + str(mjd) + "/" + str(second)
            self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
            if self.sourceindex < -1: #Couldn't find a source and should have been able to
                print "Was looking for a source at time " + str(mjd) + "/" + str(second)
                print "Current scan runs from " + str(self.scans[self.lastscannum].startmjd) + " + " + str(self.scans[self.lastscannum].scanstart) + \
                      " to " +  str(self.scans[self.lastscannum].scanstart + self.scans[self.lastscannum].scandur)
                self.lastscannum += 1
                if self.scans[self.lastscannum].containsTime(mjd, second):
                    self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
                    if self.sourceindex < -1: #Still couldn't find a source and should have been able to
                        print "Still can't find right source - aborting"
                        sys.exit()
                else:
                    print "The next scan ran from " + str(self.scans[self.lastscannum].startmjd) + " + " + str(self.scans[self.lastscannum].scanstart) + \
                      " to " +  str(self.scans[self.lastscannum].scanstart + self.scans[self.lastscannum].scandur)
                    sys.exit()
            self.changed = False
        try:
            antennaindex = self.antennamap[antennaname.strip()]
        except KeyError:
            print "Could not find antenna " + antennaname.strip() + \
                  " in the antenna map from file " + self.filename
            sys.exit()
        if self.sourceindex < 0:
            return -9e99,-9e99
        return self.scans[self.lastscannum].getDelayAndRate(mjd, second, antennaindex,
                                                       self.sourceindex)
    def getUVW(self, antenna1name, antenna2name, sourcename, mjd, second):
        while self.lastscannum < self.numscans and \
              not self.scans[self.lastscannum].containsTime(mjd, second):
            self.lastscannum += 1
            self.changed = True
        if self.lastscannum == self.numscans:
            self.lastscannum = 0
            while self.lastscannum < self.numscans and \
                  not self.scans[self.lastscannum].containsTime(mjd, second):
                self.lastscannum += 1
                self.changed = True
            if self.lastscannum == self.numscans:
                print "Couldn't find scan for time %d/%f" % (mjd, second)
                return [-9e99,-9e99,-9e99]
        if self.changed:
            self.sourceindex = self.scans[self.lastscannum].getSourceIndex(sourcename)
            self.changed = False
        antenna1index = self.antennamap[antenna1name]
        antenna2index = self.antennamap[antenna2name]
        if self.sourceindex < 0:
            return [-9e99,-9e99,-9e99]
        return self.scans[self.lastscannum].getUVW(mjd, second, antenna1index, antenna2index,
                                              self.sourceindex)

#    def getScanStartStop(self, scannum):
#        return self.scans[scannum].scanstart, \
#               self.scans[scannum].scanstart + self.scans[scannum].scandur
def wizCorrectScintAdvanced(uvdatafile, outputsntable, solmins, rashiftmas, 
                            decshiftmas, maxuvradiuslambda, examplesntable, 
                            averageifs=False):
    uvdata = AIPSUVData("JUNK","JUNK",1,1)
    if uvdata.exists():
        uvdata.zap()
    uvdata2 = AIPSUVData("JUNK","JUNK",1,2)
    if uvdata2.exists():
        uvdata2.zap()
    if not os.path.exists(examplesntable):
        print "Couldn't find example dataset " + examplesntable + " - aborting"
        sys.exit()
    if os.path.exists(outputsntable):
        if yesno("Can I delete " + outputsntable + " (no will abort)? "):
            os.remove(outputsntable)
        else:
            sys.exit()

    # Load the uvdata and the example table, grab some info needed
    fitld_uvfits(uvdatafile, uvdata, [])
    loadtable(uvdata, examplesntable, 1)
    wizuvdata      = WizAIPSUVData(uvdata)
    examplesntable = wizuvdata.table('SN', 1)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol
    numantennas    = len(uvdata.antennas)
    examplerow     = examplesntable[0]
    annumbers      = []
    maxanno        = 0
    antable        = uvdata.table('AN', 1)
    for row in antable:
        if row.nosta > maxanno:
            maxanno = row.nosta
    for i in range(maxanno):
        annumbers.append(i+1)

    nflagged = 0
    ntotal = 0
    # Flag anything that falls outside the uv radius set
    for vis in wizuvdata:
        uvrad = math.sqrt(vis.uvw[0]**2 + vis.uvw[1]**2)
        #if ntotal % 100 == 0:
        #    print "UV radius is " + str(uvrad) + ", limit is " + str(uvradius)
        if uvrad > maxuvradiuslambda:
            for i in range(num_if):
                for j in range(num_pol):
                    vis.visibility[i][0][j][2] = -vis.visibility[i][0][j][2]
            vis.update()
            nflagged += 1
        ntotal += 1
    print "Flagged " + str(nflagged) + " out of " + str(ntotal) + " visibilities due to uvrange"

    # Add new SN table to the loaded uv dataset
    newsntable = wizuvdata.attach_table('SN', 2, no_term=num_snterms)
    newsntable.keywords['NO_IF'] = num_if
    newsntable.keywords['NO_POL'] = num_pol
    newsntable.keywords['NO_TERM'] = num_snterms
    newsntable.keywords['NO_ANT'] = maxanno
    newimag1 = []
    newimag2 = []
    newdelay1 = []
    newrate1  = []
    newdelay2 = []
    newrate2  = []
    for j in range(num_if):
        newimag1.append(0.0)
        newimag2.append(0.0)
        newdelay1.append(0.0)
        newdelay2.append(0.0)
        newrate1.append(0.0)
        newrate2.append(0.0)
    examplerow.imag1 = newimag1
    examplerow.imag2 = newimag2
    examplerow.delay_1 = newdelay1
    examplerow.delay_2 = newdelay2
    examplerow.rate_1 = newrate1
    examplerow.rate_2 = newrate2
    examplerow.time_interval = solmins/1400.0

    # Run TBAVG to generate the dataset we'll use to get the scint values
    tbavg = AIPSTask("tbavg", version = aipsver)
    tbavg.indata = uvdata
    tbavg.outdata = uvdata2
    tbavg.solint = solmins*60
    tbavg.shift[1] = rashiftmas/1000.0
    tbavg.shift[2] = decshiftmas/1000.0
    tbavg()

    # Go through the TBAVG'd dataset and get the amplitudes, storing in the new SN table
    wizuvdata2 = WizAIPSUVData(uvdata2)
    ampsum = 0
    ampcount = 0
    for vis in wizuvdata2:
        for i in range(num_if):
            for j in range(num_pol):
                amp = math.sqrt(vis.visibility[i][0][j][0]**2 + vis.visibility[i][0][j][1]**2)
                ampsum += amp
                ampcount += 1

    avamp = ampsum / ampcount
    for vis in wizuvdata2:
        examplerow.time = vis.time
        newreal1 = []
        newreal2 = []
        newweight1 = []
        newweight2 = []
        asum = []
        acount = []
        for i in range(num_if):
            asum.append(0.0)
            acount.append(0)
        for i in range(num_if):
            dest = i
            if averageifs:
                dest = 0
            for j in range(num_pol):
                if vis.visibility[i][0][j][2] > 0.0:
                    amp = math.sqrt(vis.visibility[i][0][j][0]**2 + \
                                    vis.visibility[i][0][j][1]**2)
                    asum[dest] += amp
                    asum[dest] += 1
        for i in range(num_if):
            source = i
            if averageifs:
                source = 0
            correction = 1.0
            corrweight = 1.0
            if acount[source] > 0:
                amp = asum[source] / acount[source]
                correction = math.sqrt(avamp/amp)
                corrweight = 1.0/(correction*correction)
            newreal1.append(correction)
            newweight1.append(corrweight)
            if num_pol > 1:
                newreal2.append(correction)
                newweight2.append(corrweight)
        examplerow.real1 = newreal1
        examplerow.weight_1 = newweight1
        if num_pol > 1:
            examplerow.real2 = newreal2
            examplerow.weight_2 = newweight2
        for annum in annumbers:
            examplerow.antenna_no = annum
            newsntable.append(examplerow)
    newsntable.close()
    writetable(uvdata, 'SN', 2, outputsntable)

def wizCorrectScint(parentuvdata, examplesnversion, outputsnversion, splituvdata, 
                    solmins, outputsntable):
    # Get an example SN table, prefill all the entries that can be prefilled
    wizuvdata      = WizAIPSUVData(parentuvdata)
    examplesntable = wizuvdata.table('SN', examplesnversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol
    numantennas    = len(parentuvdata.antennas)
    examplerow     = examplesntable[0]
    annumbers      = []
    maxanno        = 0
    antable = parentuvdata.table('AN', 1)
    for row in antable:
        if row.nosta > maxanno:
            maxanno = row.nosta
    for i in range(maxanno):
        annumbers.append(i+1)
    tables = parentuvdata.tables
    for table in tables:
        if table[1][-2:] == 'SN' and table[0] >= outputsnversion:
            if yesno("Can i delete " + table[1][-2:] + " table version " + \
                     str(table[0]) + "? (No will abort pipeline)"):
                parentuvdata.table('SN', table[0]).zap()
            else:
                sys.exit()
    wizuvdata  = WizAIPSUVData(parentuvdata)
    newsntable = wizuvdata.attach_table('SN', outputsnversion, no_term=num_snterms)
    newsntable.keywords['NO_IF'] = num_if
    newsntable.keywords['NO_POL'] = num_pol
    newsntable.keywords['NO_TERM'] = num_snterms
    newsntable.keywords['NO_ANT'] = maxanno
    newimag1 = []
    newimag2 = []
    newdelay1 = []
    newrate1  = []
    newdelay2 = []
    newrate2  = []
    for j in range(num_if):
        newimag1.append(0.0)
        newimag2.append(0.0)
        newdelay1.append(0.0)
        newdelay2.append(0.0)
        newrate1.append(0.0)
        newrate2.append(0.0)
    examplerow.imag1 = newimag1
    examplerow.imag2 = newimag2
    examplerow.delay_1 = newdelay1
    examplerow.delay_2 = newdelay2
    examplerow.rate_1 = newrate1
    examplerow.rate_2 = newrate2
    examplerow.time_interval = solmins/1400.0

    # Run SPLAT to generate the dataset that will be used
    outputuvdata = AIPSUVData(splituvdata.name,"SCNTT",1,1)
    if outputuvdata.exists():
        outputuvdata.zap()
    splat =  AIPSTask('splat', version = aipsver)
    splat.indata = splituvdata
    splat.outdata = outputuvdata
    splat.solint = solmins
    splat()

    # Loop through the data, summing all amplitudes of all baselines
    wizuvdata2 = WizAIPSUVData(outputuvdata)
    ampsum = 0
    count = 0
    for visibility in wizuvdata2:
        for i in range(num_if):
            for j in range(num_pol):
                if visibility.visibility[i][0][j][2] > 0.0:
                    amp = math.sqrt(visibility.visibility[i][0][j][0]**2 + \
                                    visibility.visibility[i][0][j][1]**2)
                    ampsum += amp
                    count += 1
    avamp = ampsum/count
    lasttime = -999
    ifampsum = []
    ifcount  = []
    for i in range(num_if):
        ifampsum.append(0)
        ifcount.append(0)
    for visibility in wizuvdata2:
        if visibility.time > lasttime + examplerow.time_interval/2.0:
            if lasttime > 0:
                examplerow.time = lasttime
                newreal1 = []
                newreal2 = []
                newweight1 = []
                newweight2 = []
                for i in range(num_if):
                    correction = 1.0
                    corrweight = 1.0
                    if ifcount[i] > 0:
                        ifampsum[i] /= ifcount[i]
                        correction = math.sqrt(avamp/ifampsum[i])
                        corrweight = 1.0/(correction*correction)
                    newreal1.append(correction)
                    newreal2.append(correction)
                    newweight1.append(corrweight)
                    newweight2.append(corrweight)
                examplerow.real1 = newreal1
                examplerow.real2 = newreal2
                examplerow.weight_1 = newweight1
                examplerow.weight_2 = newweight2
                for annum in annumbers:
                    examplerow.antenna_no = annum
                    newsntable.append(examplerow)
            lasttime = visibility.time
            for i in range(num_if):
                ifampsum[i] = 0
                ifcount[i]  = 0
        for i in range(num_if):
            for j in range(num_pol):
                if visibility.visibility[i][0][j][2] > 0.0:
                    amp = math.sqrt(visibility.visibility[i][0][j][0]**2 + \
                                    visibility.visibility[i][0][j][1]**2)
                    ifampsum[i] += amp
                    ifcount[i] += 1
    newsntable.close()
    if os.path.exists(outputsntable):
        os.system("rm -f " + outputsntable)
    writetable(parentuvdata, 'SN', outputsnversion, outputsntable)
    #parentuvdata.zap_table('SN', outputsnversion)

################################################################################
# Functions
################################################################################

####### PARSE A VEX FILE AND GET THE REFERENCE DATE ############################
def getrefdate(vexfile):
    refdate = ""
    vexlines = open(vexfile).readlines()
    for line in vexlines:
        if "exper_nominal_start" in line:
            splitline = line.split('=')
            vexmjd = getVexMJD(splitline[-1][:-1])
            year, month, day, hour, minute, second = astro_utils.mjd2ymdhms(vexmjd)
            refdate = "%04d%02d%02d" % (year, month, day)
            break
    return refdate

####### PARSE A VEX FILE INTO A LIST OF SCANS ##################################
def getvexscaninfo(vexfile):
    scanlist = []
    srclist  = []
    vexin    = open(vexfile)
    vexlines = vexin.readlines()
    vexin.close()

    vlaautophasedetermine = False
    atline = 0
    numlines = len(vexlines)
    while atline < numlines and not '$SOURCE;' in vexlines[atline]:
        atline += 1
    if atline == numlines:
        print "Couldn't find source block in vex file"
        sys.exit()

    atline += 1
    while atline < numlines and not '$' in vexlines[atline]:
        splitline = vexlines[atline].split()
        if len(splitline) > 1 and splitline[0] == 'def':
            atline += 1
            while atline < numlines and not 'source_name' in vexlines[atline]:
                atline += 1
            if atline == numlines:
                print "Got lost in vex source block"
                sys.exit()
            srcname = vexlines[atline].split()[-1][:-1]
            atline += 1
            while atline < numlines and not 'J2000' in vexlines[atline]:
                atline += 1
            if atline == numlines:
                print "Got lost in vex source block"
                sys.exit()
            splitline = vexlines[atline].split()
            ra  = splitline[2][:-1]
            dec = splitline[5][:-1]
            hh = float(ra[:2])
            mm = float(ra[3:5])
            ss = float(ra[6:-1])
            rarad = (hh/24.0 + mm/1440.0 + ss/86400.0)*2.0*math.pi
            if dec[0] == '-':
                multiplier = -1.0
                dd = float(dec[1:3])
                mm = float(dec[4:6])
                ss = float(dec[7:-1])
            else:
                multiplier = 1.0
                dd = float(dec[:2])
                mm = float(dec[3:5])
                ss = float(dec[6:-1])
            decrad = multiplier*(dd + mm/60.0 + ss/3600.0)*math.pi/180.0
            srclist.append(Source(srcname, rarad, decrad))
        atline += 1
    print "Found " + str(len(srclist)) + " sources"

    while atline < numlines and not '$SCHED;' in vexlines[atline]:
        atline += 1
    if atline == numlines:
        print "Couldn't find scan block in vex file"
        sys.exit()

    lastsec = 0
    while atline < numlines:
        while atline < numlines and (not 'start' in vexlines[atline] or '*' in vexlines[atline]):
            if 'sec:' in vexlines[atline]:
                lastsec = int(vexlines[atline].split('sec:')[1].strip())
            if "VLA:AUTOPHASE_DETERMINE" in vexlines[atline]:
                vlaautophasedetermine = True
            if "VLA:AUTOPHASE_APPLY" in vexlines[atline] or "VLA:AUTOPHASE_OFF" in vexlines[atline]:
                vlaautophasedetermine = False
            atline += 1
        if atline == numlines:
            break
        scanname = "unknown"
        if "scan" in vexlines[atline-1]:
            scanname = vexlines[atline-1].split()[1][:-1]
        splitline = vexlines[atline].split()
        starttime = splitline[0].split('=')[1][:-1]
        if len(scanlist) > 0:
            #print "Updating scanlist[-1] with a scan duration of %d seconds" % lastsec
            scanlist[-1].incStopTime(lastsec)
        stoptime = starttime
        scansourcename = splitline[2].split('=')[1][:-1]
        scanmodename = splitline[1].split('=')[1][:-1]
        source = None
        for src in srclist:
            if src.name == scansourcename:
                source = src
                break
        if source == None:
            print "Couldn't find source " + scansourcename + " in source list!"
            sys.exit()
        scanlist.append(VexScan(scanname, starttime, stoptime, source, scanmodename, vlaautophasedetermine))
        atline += 1
    if len(scanlist) > 0:
        #print "Updating scanlist[-1] with a scan duration of %d seconds" % lastsec
        scanlist[-1].incStopTime(lastsec)
    if "bd179" in vexfile:
        for source in srclist:
            if source.name[-2:] == "AP":
                source.name = source.name + "T"
    print "Found " + str(len(scanlist)) + " scans in vex file " + vexfile
    return scanlist

####### COMPUTE THE MEDIAN ABSOLUTE DIFFERENCE OF AN ARRAY #####################
def mad(array):
    med=n.median(array)
    mad=n.median(n.abs(array-med))
    return med, mad

####### FUNCTION TO CONVERT SECONDS TO HH:MM:SS.SS FORMAT, RETURNS A STRING ####
def time2hms(seconds):
    h=int(seconds/3600)
    m=int(seconds % 3600)/60
    s=seconds-(h*3600)-(m*60)
    h=`h`
    m=`m`
    s="%4.2f" % s
    hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
    return hms

####### CONVERT A FRACTIONAL DAY INTO AIPS FORMAT ##############################
def time2aips(fracday):
    nsec=int(round(fracday*24*3600))
    dd=nsec/(24*3600)
    hh=(nsec-dd*24*3600)/3600
    mm=(nsec-dd*24*3600-hh*3600)/60
    ss=nsec-dd*24*3600-hh*3600-mm*60
    return '%02i,%02i,%02i,%02i' % (dd, hh, mm, ss)

####### FIND A LINE STARTING WITH A GIVEN KEY ##################################
def getToLine(key, imlines, linenum):
    while linenum < len(imlines) and not imlines[linenum].split(':')[0] == key:
        linenum += 1
    if linenum == len(imlines):
        print "Couldn't find key " + key + " - aborting!"
        sys.exit(1)
    return linenum, imlines[linenum].split(':')[1]

####### CONVERT A YEAR MONTH DAY DATE TO MJD ###################################
def ymd2mjd(year, month, day):
    return year*367 - int(7*(year + int((month + 9)/12))/4) + \
           int(275*month/9) + day - 678987

####### CONVERT AN AIPS HEADER DATE STRING TO MJD ##############################
def getAIPSMJD(aipsdate):
    year = int(aipsdate[0:4])
    month = int(aipsdate[5:7])
    day = int(aipsdate[8:10])
    return ymd2mjd(year, month, day)

####### CONVERT A VEX HEADER DATE STRING TO MJD ################################
def getVexMJD(vexdate):
    year = int(vexdate[0:4])
    doy = int(vexdate[5:8])
    hh = int(vexdate[9:11])
    mm = int(vexdate[12:14])
    ss = int(vexdate[15:17])
    jan1mjd = ymd2mjd(year, 1, 1)
    return jan1mjd + (doy - 1) + float(hh)/24.0 + float(mm)/1440.0 + \
           float(ss)/86400.0

####### NORMALISE THE AMPLITUDES IN AN SN TABLE SO AVERAGE IS 1 ################
def norm_snamplitudes(sntable):
    snin = open(sntable)
    snlines = snin.readlines()
    snin.close()

    lineno = 0
    numlines = len(snlines)
    wtindex = [-1,-1]
    reindex     = [-1,-1]
    imindex     = [-1,-1]
    numweightfound = 0
    numrefound = 0
    numimfound = 0
    while lineno < numlines and (wtindex[0] < 0 or wtindex[1] < 0 or \
          reindex[0] < 0 or reindex[1] < 0 or imindex[0] < 0 or imindex[1] < 0):
        splitline = snlines[lineno].split()
        if len(splitline) > 3 and splitline[2] == "'WEIGHT":
            wtindex[numweightfound] = int(splitline[-1])
            numweightfound += 1
        if len(splitline) > 3 and splitline[2][0:5] == "'REAL":
            reindex[numrefound] = int(splitline[-1])
            numrefound += 1
        if len(splitline) > 3 and splitline[2][0:5] == "'IMAG":
            imindex[numimfound] = int(splitline[-1])
            numimfound += 1
        lineno += 1
    if lineno == numlines:
        print "Couldn't find all weight/im/re row_noss - aborting!"
        sys.exit()
    wtcolindex = [-1,-1]
    recolindex = [-1,-1]
    imcolindex = [-1,-1]
    wtsum = 0.0
    re = [[], []]
    im = [[], []]
    wt = [[], []]
    count = 0
    while lineno < numlines:
        splitline = snlines[lineno].split()
        if len(splitline) > 0 and splitline[0] == "COL.":
            break
        lineno += 1
    if lineno == numlines:
        print "Couldn't find start of table - aborting!"
        sys.exit()
    wtdone = [False, False]
    redone = [False, False]
    imdone = [False, False]
    while lineno < numlines:
        splitline = snlines[lineno].split()
        if len(splitline) > 0 and splitline[0] == "COL.":
            wtcolindex = [-1,-1]
            recolindex = [-1,-1]
            imcolindex = [-1,-1]
            for i in range(len(splitline)):
                for j in range(2):
                    if splitline[i] == str(wtindex[j]):
                        wtcolindex[j] = i-1
                        wtdone[j] = True
                    if splitline[i] == str(reindex[j]):
                        recolindex[j] = i-1
                        redone[j] = True
                    if splitline[i] == str(imindex[j]):
                        imcolindex[j] = i-1
                        imdone[j] = True
            lineno += 4
            splitline = snlines[lineno].split()
        if splitline[0] == "***END*PASS***":
            lineno += 1
            continue
        for i in range(2):
#            if wtcolindex[i] > 0:
#                wtsum += math.sqrt(float(splitline[wtcolindex[i]]))
#                count += 1
            if recolindex[i] > 0 and not 'INDE' in splitline[recolindex[i]]:
                re[i].append(float(splitline[recolindex[i]]))
            if imcolindex[i] > 0 and not 'INDE' in splitline[imcolindex[i]]:
                im[i].append(float(splitline[imcolindex[i]]))
        lineno += 1
    if not redone[0] or not redone[1] or not imdone[0] or not imdone[1]:
        print "Couldn't find reals and imags in the table"
        sys.exit()
    if len(re[0]) != len(im[0]) or len(re[1]) != len(im[1]) or len(re[0]) != len(re[1]):
        print "Ended up with different length real and imag"
        sys.exit()
    count = 0
    ampsqsum = 0
    runto = len(re[0])
    if len(re[1]) < len(re[0]):
        runto = len(re[1])
    for i in range(runto):
        ampsqsum += re[0][i]*re[0][i] + im[0][i]*im[0][i]
        ampsqsum += re[1][i]*re[1][i] + im[1][i]*im[1][i]
        count += 2
    return count/ampsqsum

    '''
    correction = wtsum/count
    print "Correction is %f" % (correction)
    lineno = 0
    wtdone = [False, False]
    redone     = [False, False]
    imdone     = [False, False]
    wtcolindex = [-1,-1]
    recolindex = [-1,-1]
    imcolindex = [-1,-1]
    while lineno < numlines:
        splitline = snlines[lineno].split()
        if len(splitline) > 0 and splitline[0] == "COL.":
            break
        lineno += 1
    while lineno < numlines:
        splitline = snlines[lineno].split()
        if len(splitline) > 0 and splitline[0] == "COL.":
            wtcolindex = [-1,-1]
            recolindex = [-1,-1]
            imcolindex = [-1,-1]
            for i in range(len(splitline)):
                for j in range(2):
                    if splitline[i] == str(wtindex[j]):
                        wtcolindex[j] = i-1
                        wtdone[j] = True
                    if splitline[i] == str(reindex[j]):
                        recolindex[j] = i-1
                        redone[j] = True
                    if splitline[i] == str(imindex[j]):
                        imcolindex[j] = i-1
                        imdone[j] = True
            lineno += 4
            splitline = snlines[lineno].split()
        if splitline[0] == "***END*PASS***":
            lineno += 1
            continue
        for i in range(2):
            if wtcolindex[i] > 0:
                replacement = "%15.6E" % (float(splitline[wtcolindex[i]])/(correction*correction))
                toreplace = "%15.6E" % (float(splitline[wtcolindex[i]]))
                snlines[lineno] = snlines[lineno].replace(toreplace, replacement)
            if recolindex[i] > 0:
                replacement = "%15.6E" % (float(splitline[recolindex[i]])*correction)
                toreplace = "%15.6E" % (float(splitline[recolindex[i]]))
                print "replacing %s with %s" % (toreplace, replacement)
                snlines[lineno] = snlines[lineno].replace(toreplace, replacement)
            if imcolindex[i] > 0:
                replacement = "%15.6E" % (float(splitline[imcolindex[i]])*correction)
                toreplace = "%15.6E" % (float(splitline[imcolindex[i]]))
                snlines[lineno] = snlines[lineno].replace(toreplace, replacement)
        lineno += 1
    if not wtdone[0] or not wtdone[1] or not redone[0] or \
       not redone[1] or not imdone[0] or not imdone[1]:
        print "Couldn't find both weights/re/im in the table"
        sys.exit()
    output = open(sntable, "w")
    for line in snlines:
        output.write(line)
    output.close()'''

####### MAKE DIFFERENTIAL CL TABLE BASED ON TWO CALC/IM FILESETS ###############
def fixmodel(uvdata, inputclversion, refant, calcfile1, calcfile2start, startmjd,
             calcfile2end="", startday=-1.0, endday=-1.0):
    if len(calcfile1.split(',')) == 2: #Compatibility mode
        model1  = Model(None,calcfile1.split(',')[0],calcfile1.split(',')[1])
    else:
        model1  = Model(calcfile1)
    model2start = Model(calcfile2start)
    fixpropermotion = False
    if calcfile2end != "":
        model2end = Model(calcfile2end)
        pmduration = endday - startday
        fixpropermotion = True
    wizuvdata = WizAIPSUVData(uvdata)
    cltablefound = False
    for table in uvdata.tables:
        if table[1] == 'AIPS CL' and table[0] == inputclversion:
            cltablefound = True
        if table[1] == 'AIPS CL' and table[0] == inputclversion+1:
            print "Output CL table version number " + str(inputclversion+1) + \
                  "already exists - aborting!"
            sys.exit()
    if not cltablefound:
        print "Couldn't find input CL table " + str(inputclversion) - " aborting!"
        sys.exit()

    #Make a list of sources and antennas
    maxsuid = 0
    maxanid = 0
    sutable = uvdata.table('SU', 1)
    antable = uvdata.table('AN', 1)
    for row in sutable:
        if row.id__no > maxsuid:
            maxsuid = row.id__no
    for row in antable:
        if row.nosta > maxanid:
            maxanid = row.nosta
    sunames = []
    annames = []
    for i in range(maxsuid):
        sunames.append("")
    for i in range(maxanid):
        annames.append("")
    for row in sutable:
        sunames[row.id__no-1] = row.source.strip()
    for row in antable:
        annames[row.nosta-1] = row.anname.strip()
        print "Row " + str(row.nosta-1) + " of the annames has been set to " + row.anname.strip()

    #Work out the reference frequencies
    reffreqs = []
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            for iffreq in row.if_freq:
                freqentry = float(iffreq) + float(uvdata.header.crval[2])
                reffreqs.append(float(iffreq) + float(uvdata.header.crval[2]))
        except (AttributeError, TypeError):
                freqentry = float(row.if_freq) + float(uvdata.header.crval[2])
                reffreqs.append(float(row.if_freq) + float(uvdata.header.crval[2]))
    refwavelength = 299792458.0/reffreqs[0]
    print "Reference wavelength is " + str(refwavelength)

    # Make new CL table
    inputcltable = wizuvdata.table('CL', inputclversion)
    newcltable = wizuvdata.attach_table('CL', inputclversion+1,
                                        no_term=inputcltable.keywords['NO_TERM'])
    newcltable.keywords['NO_ANT'] = inputcltable.keywords['NO_ANT']
    newcltable.keywords['NO_POL'] = inputcltable.keywords['NO_POL']
    newcltable.keywords['NO_IF']  = inputcltable.keywords['NO_IF']
    numpol = inputcltable.keywords['NO_POL']
    lasttime = -999.99
    for row in inputcltable:
        if row.time != lasttime:
            reforgdelay, reforgrate = model1.getDelayAndRate(startmjd, row.time*86400.0,
                                                             annames[refant - 1],
                                                             sunames[row.source_id - 1])
            refdstart, refrstart = model2start.getDelayAndRate(startmjd, row.time*86400.0,
                                                               annames[refant - 1],
                                                               sunames[row.source_id - 1])
            if fixpropermotion:
                refdend, refrend = model2end.getDelayAndRate(startmjd, row.time*86400.0,
                                                             annames[refant - 1],
                                                             sunames[row.source_id - 1])
                refdstart += (refdend - refdstart) * (row.time - startday)/(endday - startday)
                refrstart += (refrend - refrstart) * (row.time - startday)/(endday - startday)
            refdeltad = refdstart - reforgdelay
            refdeltar = refrstart - reforgrate
            lasttime = row.time
        #print "About to do time " + str(row.time*86400.0) + ", antenna " + annames[row.antenna_no - 1]
        dorg, rorg = model1.getDelayAndRate(startmjd, row.time*86400.0, 
                                            annames[row.antenna_no - 1],
                                            sunames[row.source_id - 1])
        dstart, rstart = model2start.getDelayAndRate(startmjd, row.time*86400.0,
                                            annames[row.antenna_no - 1],
                                            sunames[row.source_id - 1])
        if fixpropermotion:
            dend, rend = model2end.getDelayAndRate(startmjd, row.time*86400.0,
                                            annames[row.antenna_no - 1],
                                            sunames[row.source_id - 1])
            dstart += (dend - dstart) * (row.time - startday)/(endday - startday)
            rstart += (rend - rstart) * (row.time - startday)/(endday - startday)
        deltad    = (dstart - dorg) - refdeltad
        deltar    = (rstart - rorg) - refdeltar
        if dstart < -1e99 or dend < -1e99:
            deltad = 0.0
            deltar = 0.0
        newdelay1 = []
        newrate1  = []
        newreal1  = []
        newimag1  = []
        newdelay2 = []
        newrate2  = []
        newreal2  = []
        newimag2  = []
        for i in range(len(reffreqs)):
            newdelay1.append(deltad)
            newrate1.append(deltar)
            phase = deltad*reffreqs[i]*2.0*math.pi
            newreal1.append(math.cos(phase))
            newimag1.append(math.sin(phase))
            if numpol > 1:
                newdelay2.append(deltad)
                newrate2.append(deltar)
                newreal2.append(math.cos(phase))
                newimag2.append(math.sin(phase))
        row.delay_1 = newdelay1
        row.rate_1 = newrate1
        row.real1 = newreal1
        row.imag1 = newimag1
        if numpol > 1:
            row.delay_2 = newdelay2
            row.rate_2 = newrate2
            row.real2 = newreal2
            row.imag2 = newimag2
        newcltable.append(row)
    newcltable.close()
    print "Finished creating new CL table"

    #Fix the antenna table
    antable = wizuvdata.table('AN', 1)
    antennas_updated = []
    antenna_row_noss  = []
    for row in antable:
        antennas_updated.append(False)
    for ant in model2start.antennas:
        count = 0
        for row in antable:
            if row.anname.strip() == ant.name.strip():
                row.stabxyz = [ant.x,ant.y,ant.z]
                row.update()
                antennas_updated[count] = True
            count += 1
    count = 0
    for row in antable:
        if not antennas_updated[count]:
            print "Antenna " + row.anname.strip() + " not updated!"
            sys.exit()
    antable.close()

    #Fix the source table
    sutable = wizuvdata.table('SU', 1)
    for row in sutable:
        found = False
        for scan in model2start.scans:
            if found:
                break
            for s in scan.source:
                if s.name.strip() == row.source.strip():
                    found = True
                    row.raapp  = s.ra*180.0/math.pi  - row.raepo  + row.raapp
                    row.decapp = s.dec*180.0/math.pi - row.decepo + row.decapp
                    row.raepo  = s.ra*180.0/math.pi
                    row.decepo = s.dec*180.0/math.pi
                    row.update()
                    break
        if not found:
            print "Warning - source " + row.source.strip() + " not found in model file!"
    sutable.close()

    #Update the actual uvw values
    for vis in wizuvdata:
        t1 = annames[vis.baseline[0]-1]
        t2 = annames[vis.baseline[1]-1]
        uvw = model2start.getUVW(t1, t2, sunames[int(vis.source)-1], startmjd, 
                                 vis.time*86400.0)
        if uvw[0] < -1e99:
            continue
        uvw[:] = [-x/refwavelength for x in uvw] 
        vis.uvw = uvw
        vis.update()

####### ESTIMATE IMGAE RMS BASED ON DATA WEIGHTS ###############################
def estimateimagerms(uvdata, inttime, bwhz):
    wizuvdata = WizAIPSUVData(uvdata)
    weighttotal = 0.0
    numifs = 0
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            print row
            for iffreq in row.if_freq:
                numifs += 1
        except (AttributeError, TypeError):
            numifs += 1
    print "Number of IFs is " + str(numifs)
    numstokes = len(uvdata.stokes)
    print "Number of stokes products is " + str(numstokes)
    for visibility in wizuvdata:
        for i in range(numifs):
            for j in range(numstokes):
                weight = visibility.visibility[i][0][j][2]
                if weight > 0:
                    #print weight
                    weighttotal += weight
    return math.sqrt(1.0/(weighttotal*inttime*bwhz))

####### INSERT MISSING TSYS ENTRIES ############################################
def copytsys(uvdata, inputsnversion, inantname, outantname, scalefactor, inif=-1, outif=-1):
    numantennas    = len(uvdata.antennas)
    wizuvdata      = WizAIPSUVData(uvdata)
    examplesntable = wizuvdata.table('SN', inputsnversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol

    #Make a list of antennas
    inid = -99
    outid = -99
    maxanid = 0
    antable = uvdata.table('AN', 1)
    for row in antable:
        if row.nosta > maxanid:
            maxanid = row.nosta
    annames = []
    for i in range(maxanid):
        annames.append("")
    for row in antable:
        if row.anname.strip() == inantname.strip():
            inid = row.nosta
        if row.anname.strip() == outantname.strip():
            outid = row.nosta
        annames[row.nosta-1] = row.anname.strip()
    if inid < 0:
        print "Couldn't find input antenna " + inantname
        sys.exit()
    if outid < 0:
        print "Couldn't find output antenna " + outantname
        sys.exit()                                                  

    newsntable = wizuvdata.attach_table('SN', inputsnversion+1, no_term=num_snterms)
    #newsntable = wizuvdata.attach_table('SN', inputsnversion+1)
    newsntable.keywords['NO_IF'] = num_if
    newsntable.keywords['NO_POL'] = num_pol
    newsntable.keywords['NO_TERM'] = num_snterms
    newsntable.keywords['NO_ANT'] = numantennas

    if inid != outid:
        for row in examplesntable:
            if row.antenna_no == outid:
                print "Skipping existing entry for " + outantname
                continue
            if row.antenna_no == inid:
                newsntable.append(row) # First copy this row as is
                newreal1 = []
                newweight1 = []
                if num_pol > 1:
                    newreal2 = []
                    newweight2 = []
                for i in range(num_if):
                    newreal1.append(row.real1[i]*scalefactor)
                    newweight1.append(row.weight_1[i]/(scalefactor*scalefactor))
                    if num_pol > 1:
                        newreal2.append(row.real2[i]*scalefactor)
                        newweight2.append(row.weight_2[i]/(scalefactor*scalefactor))
                row.real1 = newreal1
                row.weight_1 = newweight1
                if num_pol > 1:
                    row.real2 = newreal2
                    row.weight_2 = newweight2
                row.antenna_no = outid
                newsntable.append(row) # Now copy in also for the output antenna
            else:
                newsntable.append(row)
    else:
        if inif < 0:
            print "To copy to the same antenna, you must set inif and outif!"
            sys.exit()
        for row in examplesntable:
            if row.antenna_no == inid:
                newreal1 = []
                newweight1 = []
                if num_pol > 1:
                    newreal2 = []
                    newweight2 = []
                for i in range(num_if):
                    ifindex = i
                    if i+1 == outif or outif < 0:
                        ifindex = inif-1
                    newreal1.append(row.real1[ifindex]*scalefactor)
                    newweight1.append(row.weight_1[ifindex]/(scalefactor*scalefactor))
                    if num_pol > 1:
                        newreal2.append(row.real2[ifindex]*scalefactor)
                        newweight2.append(row.weight_2[ifindex]/(scalefactor*scalefactor))
                row.real1 = newreal1
                row.weight_1 = newweight1
                if num_pol > 1:
                    row.real2 = newreal2
                    row.weight_2 = newweight2
                newsntable.append(row) # Copy in the updated row
            else:
                newsntable.append(row)
    newsntable.close()

####### FIX BOGUS HN TSYS ENTRIES ##############################################
def fixtsys(uvdata, snversion):
    numantennas    = len(uvdata.antennas)
    wizuvdata      = WizAIPSUVData(uvdata)
    examplesntable = wizuvdata.table('SN', snversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol

    #Make a list of antennas
    maxanid = 0
    antable = uvdata.table('AN', 1)
    for row in antable:
        if row.nosta > maxanid:
            maxanid = row.nosta
    annames = []
    for i in range(maxanid):
        annames.append("")
    for row in antable:
        annames[row.nosta-1] = row.anname.strip()

    examplehnrow = None
    hnid = -1
    for row in examplesntable:
        if annames[row.antenna_no - 1] == "HN":
            examplehnrow = row
            hnid = row.antenna_no
            break
    if examplehnrow == None:
        print "Couldn't find an example HN row!"
        sys.exit()
    newreal1 = []
    newreal2 = []
    newimag1 = []
    newimag2 = []
    newweight1 = []
    newweight2 = []
    for j in range(num_if):
        newreal1.append(20.0)
        newreal2.append(20.0)
        newimag1.append(0.0)
        newimag2.append(0.0)
        newweight1.append(1.0/400.0)
        newweight2.append(1.0/400.0)
    examplehnrow.real1 = newreal1
    examplehnrow.real2 = newreal2
    examplehnrow.imag1 = newimag1
    examplehnrow.imag2 = newimag2
    examplehnrow.weight_1 = newweight1
    examplehnrow.weight_2 = newweight2

    newsntable = wizuvdata.attach_table('SN', snversion+1, no_term=num_snterms)
    #newsntable = wizuvdata.attach_table('SN', snversion+1)
    newsntable.keywords['NO_IF'] = num_if
    newsntable.keywords['NO_POL'] = num_pol
    newsntable.keywords['NO_TERM'] = num_snterms
    newsntable.keywords['NO_ANT'] = numantennas
 
    examplehnrow.time = 99.9
    for row in examplesntable:
        if row.antenna_no == hnid:
            examplehnrow.time = row.time
            newreal1 = []
            newreal2 = []
            newweight1 = []
            newweight2 = []
            for i in range(num_if):
                if row.real1[i] > 25:
                    newreal1.append(examplehnrow.real1[i])
                    newweight1.append(examplehnrow.weight_1[i])
                else:
                    newreal1.append(row.real1[i])
                    newweight1.append(row.weight_1[i])
                if row.real2[i] > 25:
                    newreal2.append(examplehnrow.real2[i])
                    newweight2.append(examplehnrow.weight_2[i])
                else:
                    newreal2.append(row.real2[i])
                    newweight2.append(row.weight_2[i])
            row.real1 = newreal1
            row.weight_1 = newweight1
            row.real2 = newreal2
            row.weight_2 = newweight2
            examplehnrow.real1 = newreal1
            examplehnrow.weight_1 = newweight1
            examplehnrow.real2 = newreal2
            examplehnrow.weight_2 = newweight2
        else:
            #Check for wacky values
            badindex = []
            goodval = -1.0
            for i in range(num_if):
                if row.real1[i] < 0.0 or row.real1[i] > 400.0:
                    badindex.append(i)
                else:
                    goodval = row.real1[i]
            if len(badindex) > 0 and goodval > 0.0:
                for i in badindex:
                    row.real1[i] = goodval
            badindex = []
            goodval = -1.0
            for i in range(num_if):
                if row.real2[i] < 0.0 or row.real2[i] > 400.0:
                    badindex.append(i)
                else:
                    goodval = row.real2[i]
            if len(badindex) > 0 and goodval > 0.0:
                for i in badindex:
                    row.real2[i] = goodval
        newsntable.append(row)
        if row.time - examplehnrow.time > 1.0/1440.0:
            examplehnrow.time = row.time
            row.antenna_no = hnid
            newreal1 = []
            newreal2 = []
            newweight1 = []
            newweight2 = []
            for i in range(num_if):
                newreal1.append(examplehnrow.real1[i])
                newreal2.append(examplehnrow.real2[i])
                newweight1.append(examplehnrow.weight_1[i])
                newweight2.append(examplehnrow.weight_2[i])
            row.real1 = newreal1
            row.weight_1 = newweight1
            row.real2 = newreal2
            row.weight_2 = newweight2
            newsntable.append(row)
    newsntable.close()

####### CONVERTS A DAY WITH FRACTIONAL COMPONENT TO 4-VAL AIPS LIST TIME #######
def fracdayToAIPSTime(fracday):
    day = int(fracday)
    hour = int(24*(fracday-day))
    minute = int(60*(24*(fracday-day) - hour))
    second = int(60*(60*(24*(fracday-day) - hour) - minute))
    return [day, hour, minute, second]

####### GIVES THE GREAT CIRCLE ANGULAR DISTANCE BETWEEN TWO POINTS #############
def haversine(az1, el1, az2, el2):
    """
    Determine the great circle angular distance between two points
    in spherical polar coordinates

    This should give the distance from the pointing centre (az1, el1) of a 2nd point (az2, el2)
    """
    a = math.sin((el1 - el2)/2)**2.0 + (math.sin((az1 - az2)/2)**2.0 * \
        math.cos(el1) * math.cos(el2))
    return 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0-a))

####### GIVES THE INITIAL BEARING FROM NORTH OF POINT 1 FROM POINT 2 ###########
def azbearing(az1, el1, az2, el2):
    """
    Initial bearing from North of point 1 from point 2

    For an az/el mounted telescope, North is always 'up' so this should
    give the bearing within the primary beam of (az2, el2) from a 
    pointing centre (az1, el1).

    The exact form is such that if the two points share the same RA, and the decl(2) > decl(1),
    the angle returned is the Parallactic Angle

    NB not tested for southern hemisphere antennas.

    """
    return -math.atan2(math.sin(az2-az1)*math.cos(el1), \
           (math.cos(el2)*math.sin(el1)) - (math.sin(el2)*math.cos(el1)*math.cos(az2-az1)))

####### GIVES THE SEPARATION BETWEEN TWO POINTINGS IN AZ AND EL ################
def get_beam_position(pointing_az, pointing_al, target_az, target_al):
    """
    return cartesian offset (N and E) of pointing from target
    inputs in radians
    return position of target w.r.t. pointing in cartesian coordinates
           az, al both in radians
    """
    # work out distance from pointing to target 
    # (doesn't really need recalculating every time but this saves 
    # passing it to this function)
    r = haversine(target_az, target_al, pointing_az, pointing_al)
    # Work out angle from North in the primary beam
    theta = azbearing(target_az, target_al, pointing_az, pointing_al)
    return r*math.sin(theta), r*math.cos(theta)

####### GIVES THE ANTENNA RESPONSE FOR A GIVEN POINTING ########################
def ant_response(now, antcode, cephem, tephem, wavelength, verbose):
    """
    Calculates the beam response of an antenna given its
    pointing centre, the target coordinates, bmaj in Az and El, and the
    amount and direction of beam squint
    """
    # Get the antenna and beam location, and the wavelength and scale compared to reference wavelength
    ant = vlbalocations[antcode]
    beam = vlbabeams[antcode]
    f = 299792458.0/wavelength
    beamcorr = beam[-1]/f

    # get date and location of antenna and calculate the separation between
    # pointing centre and target
    date=ephem.Date(now)
    lat=astro_utils.stringToRad(ant[0], False)
    lon=astro_utils.stringToRad(ant[1], False)
    obs=ephem.Observer()
    obs.long, obs.lat, obs.elevation=ant[1], ant[0], int(ant[2])
    obs.date=date

    # pyephem calculates Az/El using an atmospheric model, set pressure to zero to switch this off
    # obs.pressure=0.0
    # this is the "international standard atmosphere"
    # see http://de.wikipedia.org/wiki/Barometrische_H%C3%B6henformel#Internationale_H.C3.B6henformel
    obs.pressure=1013.25*(1-0.0065*obs.elevation/288.15)**5.255
    cephem.compute(obs)
    tephem.compute(obs)

    # calculate how much the beams differ in rad
    # FIXME: This is inverted, but seems to work??
    DAz=-beam[0]*beamcorr*(math.pi/180.0)/2.0/60.0
    DEl=-beam[1]*beamcorr*(math.pi/180.0)/2.0/60.0


    # Get the separation between pointing center and target
    sep=get_beam_position(float(cephem.az), float(cephem.alt), \
                          float(tephem.az), float(tephem.alt))

    # Add squint offsets, assuming everything is in a plane as angles are small
    Az_LCP=sep[0]+DAz
    Az_RCP=sep[0]-DAz
    El_LCP=sep[1]+DEl
    El_RCP=sep[1]-DEl

    # calculate the total offsets, in radians
    sep_LCP=math.sqrt(Az_LCP**2.0 + El_LCP**2.0)
    sep_RCP=math.sqrt(Az_RCP**2.0 + El_RCP**2.0)

    # get the primary beam responses as a function of that separation
    if antcode == "GB":
        LCP = airyresponse(sep_LCP, gbtdiameter, wavelength)
        RCP = airyresponse(sep_RCP, gbtdiameter, wavelength)
    else:
        LCP = airyresponse(sep_LCP, vlbadiameter, wavelength)
        RCP = airyresponse(sep_RCP, vlbadiameter, wavelength)

#    # compute coords where the beams are actually pointed at
#    # assuming everything is in a plane as angles are small
#    Az_LCP=cephem.az-DAz
#    Az_RCP=cephem.az+DAz
#    El_LCP=cephem.alt-DEl
#    El_RCP=cephem.alt+DEl
#
#    # now calculate the separation between the target and the two beams
#    sep_LCP   =ephem.separation([Az_LCP, El_LCP],    [tephem.az, tephem.alt])
#    #sep_LCP_az=ephem.separation([Az_LCP, tephem.alt],[tephem.az, tephem.alt])
#    #sep_LCP_el=ephem.separation([tephem.az, El_LCP], [tephem.az, tephem.alt])
#    sep_RCP   =ephem.separation([Az_RCP, El_RCP],    [tephem.az, tephem.alt])
#    #sep_RCP_az=ephem.separation([Az_RCP, tephem.alt],[tephem.az, tephem.alt])
#    #sep_RCP_el=ephem.separation([tephem.az, El_RCP], [tephem.az, tephem.alt])
#
#    # get the primary beam responses as a function of that separation
#    LCP = airyresponse(sep_LCP, vlbadiameter, wavelength)
#    RCP = airyresponse(sep_RCP, vlbadiameter, wavelength)
#    #print "Antenna %s R/L ratio: %.3f" % (antcode, RCP/LCP)

    # make lots of noise if desired
    if verbose:
        junk, azlcpstr = astro_utils.posradians2string(Az_LCP, Az_LCP)
        junk, ellcpstr = astro_utils.posradians2string(El_LCP, El_LCP)
        junk, azrcpstr = astro_utils.posradians2string(Az_RCP, Az_RCP)
        junk, elrcpstr = astro_utils.posradians2string(El_RCP, El_RCP)
        print "\nTime                : %s" % date
        print "Station             : %s  Long: %s  Lat: %s  Elev: %im  Press: %.1f" % (antcode, obs.long, obs.lat, obs.elevation, obs.pressure)
        print "Pointing centre     : ra=%s  dec=%s  az=%s  el=%s" % (cephem.ra, cephem.dec, cephem.az, cephem.alt)
        print "LCP pointing centre :                                az=%s  el=%s" % (Az_LCP, El_LCP)
        print "RCP pointing centre :                                az=%s  el=%s" % (Az_RCP, El_RCP)
        print "Target              : ra=%s  dec=%s  az=%s  el=%s" % (tephem.ra, tephem.dec, tephem.az, tephem.alt)
        print "Scaling factor      : %.4f" % beamcorr
        print "Beam offsets        : %.4f'  %.4f'" % (beam[0], beam[1])
        print "Scaled beam offsets : %.4f'  %.4f'" % (beam[0]*beamcorr, beam[1]*beamcorr)
        print "Beam sizes          : %.4f'  %.4f'  %.4f'  %.4f'" % (beam[2], beam[3], beam[4], beam[5])
        print "Scaled beam sizes   : %.4f'  %.4f'  %.4f'  %.4f'" % (beam[2]*beamcorr, beam[3]*beamcorr, beam[4]*beamcorr, beam[5]*beamcorr)
        print "LCP-Target          : %.4f'" % (sep_LCP*180*60/math.pi)
        print "RCP-Target          : %.4f'" % (sep_RCP*180*60/math.pi)
        print "LCP response        : %.5f" % LCP
        print "RCP response        : %.5f" % RCP

    # Return the results
    return (RCP, LCP)

####### THEORETICAL RESPONSE OF A FILLED DISH  #################################
def airyresponse(theta, D, lam):
    u = math.sin(theta)
    if u==0:
        return 1.0
    else:
        return (2*jn(1,(math.pi*u*D/lam))/(math.pi*u*D/lam))**2

####### MAKE AN AMPLITUDE CORRECTION FOR PRIMARY BEAM EFFECTS ##################
def correct_primarybeam(uvdata, examplesnversion, phasecentrenum, scanlist, fieldsourcenames, 
                        iscal, issearch=True, isonepointing=False, onlygettimes=False, 
                        skipmissingsources=False):
    numantennas    = len(uvdata.antennas)
    wizuvdata      = WizAIPSUVData(uvdata)
    examplesntable = wizuvdata.table('SN', examplesnversion)
    num_if         = examplesntable.keywords['NO_IF']
    num_pol        = examplesntable.keywords['NO_POL']
    num_snterms    = num_if*num_pol
    #antenna_diam   = 25.0 # metres
    # Changed for best-fit to Airy disk response
    #antenna_diam   = 25.474 # metres

    antable = uvdata.table('AN', 1)
    annumbers = []
    annames = []
    maxannumber = 0
    for row in antable:
        annumbers.append(row.nosta)
        annames.append( row.anname.strip())
        if row.nosta > maxannumber:
            maxannumber = row.nosta

    freqs    = []
    #beamhwhm = []
    lambdas  = []
    equivtimeonsource = {}
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            #print row
            for iffreq, ifbw, ifsideband in zip(row.if_freq, row.total_bandwidth, row.sideband):
                freqentry = float(iffreq) + float(uvdata.header.crval[2]) + float(ifbw)*float(ifsideband)/2.0
                freqs.append(freqentry)
                #print freqentry
        except (AttributeError, TypeError):
            freqentry = float(row.if_freq) + float(uvdata.header.crval[2]) + float(row.total_bandwidth)*float(row.sideband)/2.0
            freqs.append(freqentry)
    for f in freqs:
        lambdas.append(299792458.0/f)
        #beamhwhm.append((1.22*299792458/(2*f*antenna_diam))*180.8*60.0/math.pi)

    if not onlygettimes:
        newsntable     = wizuvdata.attach_table('SN', examplesnversion+1, no_term=num_snterms)
        newsntable.keywords['NO_IF'] = num_if
        newsntable.keywords['NO_POL'] = num_pol
        newsntable.keywords['NO_TERM'] = num_snterms
        newsntable.keywords['NO_ANT'] = maxannumber
    row = examplesntable[0]
    gainlength = num_if
    newdelay1 = []
    newdelay2 = []
    newrate1 = []
    newrate2 = []
    newimag1 = []
    newimag2 = []
    for i in range(gainlength):
        newdelay1.append(0.0)
        newdelay2.append(0.0)
        newrate1.append(0.0)
        newrate2.append(0.0)
        newimag1.append(0.0)
        newimag2.append(0.0)
    row.delay_1 = newdelay1
    row.delay_2 = newdelay2
    row.rate_1 = newrate1
    row.rate_2 = newrate2
    row.imag1 = newimag1
    row.imag2 = newimag2
    row.i_far_rot = 0
    row.time_interval = 0.0
    row.source_id = 1

    aipsstartmjd = getAIPSMJD(uvdata.header.date_obs)
    sutable = uvdata.table('SU', 1)
    numsusources = 0
    for srow in sutable:
        numsusources += 1
    for scan in scanlist:
        meancorr = 0.0
        meancorrcount = 0
        #pbcor = []
        #for a in uvdata.antennas:
        #    pbcor.append([])
        sutablerow = None
        scansource = scan.source
        print "Looking for match to vexname ", scansource.name
        if iscal:
            #for i in range(num_if):
            #    pbcor.append(1.0)
            localname = scan.source.name
            if numsusources == 1: #Probably a bd161 experiment, must be this source
                for srow in sutable:
                    sutablerow = srow
                    break
            if (localname[-2] == '-' or "PT" in localname[-4:]) and numsusources > 1:
                continue #Target source, not important
            for srow in sutable:
                if srow.source.strip() == localname:
                    sutablerow = srow
                    break
            if sutablerow == None:
                print "Couldn't find source " + localname + " in SU table!"
                if skipmissingsources:
                    continue
                else:
                    sys.exit()
        else:
            if issearch:
                if not scansource.name[-2] == '-' and not "-PT" in scansource.name[-4:-1] and not isonepointing: #A calibrator
                    print "Skipping source " + scansource.name
                    continue
                #print "Looking at source " + scansource.name
                #print "Ra is %f, dec is %f" % (scansource.ra, scansource.dec)
                ssname = scansource.name[:-2]
                if "-PT" in scansource.name[-4:-1]:
                    ssname = scansource.name[:-4]
                tcount = 0
                # Figure out what source this is, calculate PB attenuation
                if len(fieldsourcenames) > 1:
                    for s in fieldsourcenames:
                        #print ssname + " " + str(s[0])
                        if ssname == s[0][0] or (len(ssname) == 5 and ssname in s[0][0]):
                            break
                        elif ssname[-1] == "A" and ssname[:-1] in s[0][0] and s[0][0][-1] == "A":
                            break
                        tcount += 1
                    if tcount == len(fieldsourcenames):
                        print "Couldn't find source " + scansource.name + " (was looking for " + ssname + ")"
                        if skipmissingsources:
                            continue
                        else:
                            sys.exit()
                else:
                    print "There is only one target, ssname is " + ssname
                fcount = int(scansource.name[-1]) - 1
                print fcount
                if phasecentrenum >= len(fieldsourcenames[tcount][fcount]):
                    print "No phase centre %d for pointing %s" % (phasecentrenum, scansource.name)
                    continue
                localname = fieldsourcenames[tcount][fcount][phasecentrenum]
            else:
                try:
                    localname = fieldsourcenames[scan.source.name]
                    #print "scanname is %s, localname is %s" % (scan.source.name, localname)
                except KeyError:
                    print "No inbeam for source " + scan.source.name + " for this pass"
                    continue
            for srow in sutable:
                print srow.source.strip()
                if srow.source.strip() == localname:
                    sutablerow = srow
                    break
            if sutablerow == None:
                print "Couldn't find source " + localname + " in SU table!"
                if issearch:
                    sys.exit()
                else:
                    continue
            #print sutablerow
        # create centre and target objects in pyephem
        crastr, cdecstr = astro_utils.posradians2string(scansource.ra, scansource.dec)
        cephem = ephem.readdb("%s,f|M|F7,%s, %s,2.02,2000" % (scansource.name, crastr, cdecstr))
        trastr, tdecstr = astro_utils.posradians2string(math.pi*sutablerow.raepo/180.0, math.pi*sutablerow.decepo/180.0)
        tephem = ephem.readdb("%s,f|M|F7,%s, %s,2.02,2000" % (scansource.name, trastr, tdecstr))
        # Now changed to include squint...
        #separcmin = scansource.getSepArcmin(math.pi*sutablerow.raepo/180.0, math.pi*sutablerow.decepo/180.0)
        #print "For source %s, with tablerow RA %.6f and vex RA %.6f, separcmin is %f" % (scansource.name, scansource.ra, math.pi*sutablerow.raepo/180.0, separcmin)
        #print "Sep in arcmin is %f for source %s since positions were(%f,%f), (%f,%f)" % \
            #      (separcmin, sutablerow.source, sutablerow.raepo*math.pi/180.0, 
            #       sutablerow.decepo*math.pi/180.0, scansource.ra, scansource.dec)
        #for i in range(num_if):
        #    #x = separcmin/beamhwhm[i]
        #    #pbcor.append(math.sqrt(math.exp(-(x*x)/1.442695)))
        #    # Now changed to use Airy disk response
        #    #theta = math.pi*separcmin/(60*180)
        #    #pbcor.append(math.sqrt(airyresponse(theta, antenna_diam, lambdas[i])))
        #    # Now changed again to allow for squint
        #    for i in range(len(uvdata.antennas)):
        #        pbcor[i].append(0)
        #print "Separcmin was " + str(separcmin)
        #print "X was " + str(x)
        #print pbcor[-1]

        #Do three points - 5 sec after start of scan, one in middle of scan, one 5 sec from end of scan
        scanpoints = []
        scanpoints.append(scan.startmjd - aipsstartmjd + 5.0/86400.0)
        scanpoints.append((scan.startmjd+scan.stopmjd)/2.0 - aipsstartmjd)
        scanpoints.append(scan.stopmjd - aipsstartmjd - 5.0/86400.0)
        #print "%f %f" % (scan.startmjd, scan.stopmjd)
        for scanpoint in scanpoints:
            year, month, day, hour, minute, second = astro_utils.mjd2ymdhms(scanpoint + aipsstartmjd)
            now = datetime.datetime(year, month, day, hour, minute, int(second))
            #count = 0        
            #for a in uvdata.antennas:
            #    pbcor[count][count].append(ant_response(now, a, cephem, tephem, lambdas[i], False))
            #    count += 1

            row.time = scanpoint
            row.source_id = sutablerow.id__no
            # Put in an entry for every antenna
            for i in range(numantennas):
                row.antenna_no = annumbers[i]
                
                newreal1 = []
                newreal2 = []
                newweight1 = []
                newweight2 = []
                for j in range(gainlength):
                    rcpcorr, lcpcorr = ant_response(now, annames[i], cephem, tephem, lambdas[j], False)
                    meancorr += rcpcorr + lcpcorr
                    meancorrcount += 2
                    newreal1.append(1.0/math.sqrt(rcpcorr))
                    newreal2.append(1.0/math.sqrt(lcpcorr))
                    newweight1.append(rcpcorr)
                    newweight2.append(lcpcorr)
                row.real1 = newreal1
                row.real2 = newreal2
                row.weight_1 = newweight1
                row.weight_2 = newweight2
                print "Doing time %f for antenna %d" % (scanpoint, i+1)
                if not onlygettimes:
                    newsntable.append(row)
        try:
            test = equivtimeonsource[localname]
        except KeyError:
            equivtimeonsource[localname] = 0.0
        scanlen = (scan.stopmjd - scan.startmjd)*86400.0
        meancorr /= meancorrcount
        equivtimeonsource[localname] += scanlen*meancorr*meancorr
        print "Adding %f seconds (weighted to %f) to %s onsourcetime, it is now %f, meancorr was %f" % \
              (scanlen, scanlen*meancorr*meancorr, localname, equivtimeonsource[localname], meancorr)
    if not onlygettimes:
        newsntable.close()
    return equivtimeonsource

####### RUN SOURCE EXTRACTOR TO ESTIMATE RMS IN A CLEAN MAP ####################
def run_sextractor(filename):
    sextractor_bin='/usr/bin/sextractor'
    conf_file='~/packages/src/python_scripts/s_extractor.conf'
    command='%s %s -c %s' % (sextractor_bin, filename, conf_file)
    print command
    os.system(command)

####### RUN BLOBCAT TO ESTIMATE INTEGRATED FLUX ################################
def blobcat_intflux(cleanfits, statsfile, doplot=False, dsnr=6.5, fsnr=5.0, edgepixels=100):
    run_sextractor(cleanfits)
    plotstring = ""
    if doplot:
        plotstring = "--plot "
    command='blobcat.py %s --rmsmap rms.fits %s --fSNR=%.3f --dSNR=%.3f --ppe=0.01 --pasbe=0.15 --cpeRA=0.0001 --cpeDec=0.0001 --edgemin=%d' % \
            (cleanfits, plotstring, fsnr, dsnr, edgepixels)
    print command
    os.system(command)
    blobcatoutfile = cleanfits[:-5] + '_blobs.txt'
    cleanfitsdir = blobcatoutfile[:blobcatoutfile.rfind('/')-1]
    blobcatdir = cleanfitsdir[:cleanfitsdir.rfind('/')-1] + '/blobcat/'
    resultslines = open(blobcatoutfile).readlines()
    if os.path.exists(blobcatdir):
        os.system("mv %s %s" % (blobcatoutfile, blobcatdir))
    intfluxes  = []
    interrs    = []
    peakfluxes = []
    peakerrs   = []
    peaksnrs   = []
    ras        = []
    decs       = []
    raerrs     = []
    decerrs    = []
    intfluxsum = 0
    interrsum = 0
    snrsquaresum = 0
    for line in resultslines:
        if len(line) < 2 or line[0] == '#':
            continue
        splitline = line.split()
        intfluxes.append(float(splitline[37])*1000) 
        interrs.append(float(splitline[38])*1000) 
        intfluxsum += intfluxes[-1]
        interrsum += interrs[-1]
        peakfluxes.append(float(splitline[32])*1000)
        peakerrs.append(float(splitline[33])*1000)
        peaksnrs.append(float(splitline[27]))
        snrsquaresum += peaksnrs[-1]*peaksnrs[-1]
        ras.append(float(splitline[4]))
        decs.append(float(splitline[5]))
        raerrs.append(float(splitline[6]))
        decerrs.append(float(splitline[7]))
    if not os.path.exists(statsfile): 
        # We want to start with a JMFIT file and replace some segments
        print "Stats file " + statsfile + " does not exist!"
        sys.exit()
    if len(ras) == 0:
        print "No blobcat detection"
        return 0, 0 # Note exit here - we leave stats file unchanged
    statslines = open(statsfile).readlines()
    statsout = open(statsfile, 'w')
    statsout.write(statslines[0])
    statsout.write("Peak flux (mJy/beam): " + str(peakfluxes[0]) + '\n')
    statsout.write("Blob int flux (mJy):  %.3f +/- %.3f\n" % (intfluxsum, interrsum))
    statsout.write(statslines[3])
    statsout.write("S/N:                  %.3f\n" % math.sqrt(snrsquaresum))
    print "Replacing S/N of " + statslines[4].split()[-1] + " with " + str(math.sqrt(snrsquaresum))
    #rarad = ras[0]*math.pi/180
    #decrad = decs[0]*math.pi/180
    #rastr, decstr = astro_utils.posradians2string(rarad, decrad)
    #statsout.write("Actual RA:            %s\n" % rastr)
    #statsout.write("Actual Dec:           %s\n" % decstr)
    statsout.write(statslines[5])
    statsout.write(statslines[6])
    jmfit1compfitstring = statslines[7][22:-1]
    jmfitdeconvstring = statslines[8][22:-1]
    statsout.write("JMFIT 1-gauss fit:    %s\n" % jmfit1compfitstring)
    statsout.write("Deconvolved 1-g fit:  %s\n" % jmfitdeconvstring)
    raerrmas = raerrs[0]*3.6e6
    raerrhms = raerrmas/(15*math.cos(decs[0]*math.pi/180.0))
    decerrmas = decerrs[0]*3.6e6
    statsout.write("Est. RA error (mas):  %.6f\n" % raerrmas)
    statsout.write("Est. RA error (hms):  %.6f\n" % raerrhms)
    statsout.write("Est. Dec error (mas): %.6f\n" % decerrmas)
    statsout.close()
    return peakfluxes[0], math.sqrt(snrsquaresum)

####### FLAG ONE UVDATASET BASED ON EXPECTED TIME/BANDWIDTH DECORRELATION ######
def smearingflag(brightsrcuvdata, targetuvdata, correlationlimit, fgver=1):
    # First check that the files match acceptably
    if len(brightsrcuvdata.antennas) != len(targetuvdata.antennas):
        print "Antenna tables do not match! Aborting."
        return -1.0
    maxantennano = 0
    antable = brightsrcuvdata.table("AIPS AN", 1)
    for row in antable:
       print row
       if row.nosta > maxantennano:
           maxantennano = row.nosta
    antennas = []
    for i in range(maxantennano):
        antennas.append("??")
    for row in antable:
        antennas[row.nosta-1] = row.anname.strip()
    for a1, a2 in zip(brightsrcuvdata.antennas, targetuvdata.antennas):
        if not a1==a2:
            print "Antenna tables do not match! Aborting."
            return -1.0
    fqtable1 = brightsrcuvdata.table("AIPS FQ", 1)
    for row in fqtable1:
        print brightsrcuvdata.header
        try:
            matchfreq = row.if_freq[0] + float(brightsrcuvdata.header.crval[2]) - brightsrcuvdata.header.crpix[2]*brightsrcuvdata.header.cdelt[2]
            reffreq = row.if_freq[0] + float(brightsrcuvdata.header.crval[2])
        except TypeError:
            matchfreq = row.if_freq + float(brightsrcuvdata.header.crval[2]) - brightsrcuvdata.header.crpix[2]*brightsrcuvdata.header.cdelt[2]
            reffreq = row.if_freq + float(brightsrcuvdata.header.crval[2])
        break
    refwavelength1 = speedOfLight / reffreq
    fqtable2 = targetuvdata.table("AIPS FQ", 1)
    for row in fqtable2:
        try:
            matchfreq2 = row.if_freq[0] + float(targetuvdata.header.crval[2]) - targetuvdata.header.crpix[2]*targetuvdata.header.cdelt[2]
            reffreq2 = row.if_freq[0] + float(targetuvdata.header.crval[2])
        except TypeError:
            matchfreq2 = row.if_freq + float(targetuvdata.header.crval[2]) - targetuvdata.header.crpix[2]*targetuvdata.header.cdelt[2]
            reffreq2 = row.if_freq + float(targetuvdata.header.crval[2])
        break
    refwavelength2 = speedOfLight / reffreq2
    if math.fabs(matchfreq2 - matchfreq) > 10:
        print "Mismatched frequencies? %.1f, %.1f; Aborting." % (reffreq, reffreq2)
        return -1.0
    refbw = getchannelbandwidthhz(brightsrcuvdata)
    print "Refbw is " + str(refbw)
    if math.fabs(getchannelbandwidthhz(targetuvdata) - refbw) > 10:
        print "Mismatched bandwidths? Aborting."
        return -1.0
    # Now loop through checking the decorrelation at 5 minute intervals
    # SHOULD ALSO CHECK THAT IT IS A SINGLE SOURCE FILE...
    wizuvdata1 = WizAIPSUVData(brightsrcuvdata)
    wizuvdata2 = WizAIPSUVData(targetuvdata)
    times1 = {}
    wvals1 = {}
    row_noss1 = {}
    times2 = {}
    wvals2 = {}
    row_noss2 = {}
    flags = {}
    mintime = 9999
    maxtime = -9999
    for a in range(1, maxantennano):
        for b in range(a+1,maxantennano+1):
            baselinenum = 256*a + b
            baselinestr = antennas[a-1] + '-' + antennas[b-1]
            times1[baselinestr] = []
            wvals1[baselinestr] = []
            times2[baselinestr] = []
            wvals2[baselinestr] = []
            flags[baselinestr] = []
            row_noss1[baselinestr] = 0
            row_noss2[baselinestr] = 0
    lasttime = -9999
    counts = {}
    for visibility in wizuvdata1:
        baselinenum = 256*visibility.baseline[0] + visibility.baseline[1]
        baselinestr = antennas[visibility.baseline[0]-1] + '-' + antennas[visibility.baseline[1]-1]
        w = visibility.uvw[2]
        if visibility.time > lasttime:
            if lasttime > 0:
                inttime = visibility.time - lasttime
                if str(inttime) in counts:
                    counts[str(inttime)] += 1
                else:
                    counts[str(inttime)] = 1
            lasttime = visibility.time
        times1[baselinestr].append(visibility.time)
        wvals1[baselinestr].append(w)
    maxcount = 0
    for inttime in counts:
        if counts[inttime] > maxcount:
            maxcount = counts[inttime]
            bestinttime = float(inttime)*86400
    for visibility in wizuvdata2:
        baselinenum = 256*visibility.baseline[0] + visibility.baseline[1]
        baselinestr = antennas[visibility.baseline[0]-1] + '-' + antennas[visibility.baseline[1]-1]
        w = visibility.uvw[2]
        times2[baselinestr].append(visibility.time)
        wvals2[baselinestr].append(w)
    skipbaselines = []
    for a in range(1, maxantennano):
        for b in range(a+1,maxantennano+1):
            baselinenum = 256*a + b
            baselinestr = antennas[a-1] + '-' + antennas[b-1]
            if len(times1[baselinestr]) == 0 or len(times2[baselinestr]) == 0:
                print "Skipping baseline " + baselinestr + " (not enough times)"
                skipbaselines.append(baselinestr)
            else:
                print "%s: times1 spans %f to %f, times2 spans %f to %f" % (baselinestr, times1[baselinestr][0], times1[baselinestr][-1], times2[baselinestr][0], times2[baselinestr][-1])
                if times1[baselinestr][0] < mintime:
                    mintime = times1[baselinestr][0]
                if times2[baselinestr][0] < mintime:
                    mintime = times2[baselinestr][0]
                if times1[baselinestr][-1] > maxtime:
                    maxtime = times1[baselinestr][-1]
                if times2[baselinestr][-1] > maxtime:
                    maxtime = times2[baselinestr][-1]
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = targetuvdata
    uvflg.outfgver = fgver
    uvflg.aparm[1:] = [0]
    uvflg.reason = "Confusing source"
    uvflg.opcode = "FLAG"
    flaggedtime = 0.0
    unflaggedtime = 0.0
    for a in range(1, maxantennano):
        for b in range(a+1,maxantennano+1):
            baselinenum = 256*a + b
            baselinestr = antennas[a-1] + '-' + antennas[b-1]
            if baselinestr in skipbaselines:
                continue
            currentindex1 = 0
            currentindex2 = 0
            lastdelaydiff = -9e99
            lasttime = -9e99
            checktime = mintime + 1.0/1440.0
            #print "Starting at time " + str(checktime) + ", maxtime is " + str(maxtime)
            while checktime < maxtime - 1.0/1440.0:
                print "time " + str(checktime) + ", maxtime is " + str(maxtime)
                while currentindex1 < len(times1[baselinestr]) and times1[baselinestr][currentindex1] < checktime - 0.05/86400.0:
                    #print str(times1[baselinestr][currentindex1]) + ", " + str(checktime)
                    currentindex1 += 1
                while currentindex2 < len(times2[baselinestr]) and times2[baselinestr][currentindex2] < checktime - 0.05/86400.0:
                    currentindex2 += 1
                if currentindex1 == 0:
                    #print "File 1 didn't start baseline " + baselinestr + " early enough, continuing..."
                    checktime += 5.0/1440.0
                    continue
                if currentindex2 == 0:
                    #print "File 2 didn't start baseline " + baselinestr + " early enough, continuing..."
                    checktime += 5.0/1440.0
                    continue
                if currentindex1 >= len(times1[baselinestr]):
                    #print "File 1 ran out of times too early for baseline " + baselinestr + ", continuing..."
                    checktime += 5.0/1440.0
                    continue
                if currentindex2 >= len(times2[baselinestr]):
                    #print "File 2 ran out of times too early for baseline " + baselinestr + ", continuing..."
                    checktime += 5.0/1440.0
                    continue
                print currentindex2, currentindex1, len(times1[baselinestr]), len(times2[baselinestr])
                if math.fabs(times1[baselinestr][currentindex1] - times2[baselinestr][currentindex2])*86400 > 0.5:
                    print "This time doesn't match well enough (%.2f, %.2f) - moving forward half a minute" % (times1[baselinestr][currentindex1]*86400, times2[baselinestr][currentindex2]*86400)
                    checktime += 0.5/1440.0
                    continue
                delaydiff = (refwavelength2*wvals2[baselinestr][currentindex2]-refwavelength1*wvals1[baselinestr][currentindex1])/speedOfLight
                print "%s: at checktime %.6f, delay diff is %.10f" % (baselinestr, checktime, delaydiff)
                if lasttime < -999:
                    lasttime = times1[baselinestr][currentindex1]
                    lastdelaydiff = delaydiff
                    checktime += 5.0/1440.0
                    continue #Now we have something to estimate rate from...
                if times1[baselinestr][currentindex1] == lasttime:
                    checktime += 1
                    continue
                measuredtime = (times1[baselinestr][currentindex1] + lasttime)/2.0
                measureddelaydiff = (delaydiff + lastdelaydiff)/2.0
                measuredratediff = (delaydiff - lastdelaydiff)/(86400*(times1[baselinestr][currentindex1] - lasttime))
                ntimeturns = math.fabs(bestinttime*measuredratediff*reffreq)
                nbandturns = math.fabs(measureddelaydiff*refbw)
                mintimedecorr = 1.0
                if ntimeturns > 1.0:
                    mintimedecorr = 1.0/(2*math.pi*(math.floor(math.fabs(ntimeturns)) + 0.5))
                minbanddecorr = 1.0
                if nbandturns > 1.0:
                    minbanddecorr = 1.0/(2*math.pi*(math.floor(math.fabs(nbandturns)) + 0.5))
                maxcorrelation = mintimedecorr * minbanddecorr
                totalturns = math.sqrt(nbandturns*nbandturns + ntimeturns*ntimeturns)
                doflag = False
                print "For baseline %s at time %.6f, delay was %.3f microseconds, rate was %.3g s/s" % (baselinestr, checktime, measureddelaydiff*1e6, measuredratediff)
                print "  ntimeturns was %.3f and nbandturns was %.3f, so maximum remaining correlated signal fraction is %.3f" % (ntimeturns, nbandturns, maxcorrelation)
                if maxcorrelation > correlationlimit:
                    doflag = True
                    #print "For baseline %s at time %.6f, delay was %.3f microseconds, rate was %.3g s/s" % (baselinestr, checktime, measureddelaydiff*1e6, measuredratediff)
                    #print "  ntimeturns was %.3f and nbandturns was %.3f, so minimum expected total decoorelation is %.3f" % (ntimeturns, nbandturns, maxcorrelation)
                flags[baselinestr].append([measuredtime, doflag])
                lastdelaydiff = delaydiff
                lasttime = times1[baselinestr][currentindex1]
                checktime += 5.0/1440.0
            flags[baselinestr].append([maxtime+bestinttime/86400.0, False])
            lastflagtime = mintime-bestinttime/86400.0
            lastflagval = False
            for f in flags[baselinestr]:
                flagtime = f[0]
                flagval = f[1]
                if flagval or lastflagval:
                    print "FLAGGING " + baselinestr
                    uvflg.antenna[1] = a
                    uvflg.baseline[1] = b
                    uvflg.timerang[1:5] = fracdayToAIPSTime(lastflagtime)
                    uvflg.timerang[5:9] = fracdayToAIPSTime(flagtime+1./86400.)
                    uvflg()
                    flaggedtime += f[0]-lastflagtime
                else:
                    unflaggedtime += f[0]-lastflagtime
                lastflagtime = f[0]
                lastflagval = f[1]
    return flaggedtime/(unflaggedtime+flaggedtime)

####### CALCULATE SOURCE POSITION CHANGE BASED ON PMPAR FIT FILE ###############
def calcPosOffsetVsPmpar(vexfile, source, fitfile):
    if not os.path.exists(vexfile):
        print vexfile + " does not exist!"
        sys.exit()
    if not os.path.exists(fitfile):
        print fitfile + " does not exist!"
        sys.exit()
    mjd = -1
    orgdecrad = -99
    orgrarad = -99
    vexlines = open(vexfile).readlines()
    for i in range(len(vexlines)):
        line = vexlines[i]
        if mjd < 0 and "MJD" in line:
            mjd = float(line.split()[-1])
        if "nominal_start" in line:
            nominalstart = line.strip()[-10:]
        if "nominal_stop" in line:
            nominalstop = line.strip()[-10:]
        if source in line and "def" in line:
            srcline = vexlines[i+3]
            rastr = srcline.split()[2][:-2]
            decstr = srcline.split()[5][:-2]
            rastr = rastr[:2] + ':' + rastr[3:5] + ':' + rastr[6:]
            decstr = decstr[:2] + ':' + decstr[3:5] + ':' + decstr[6:]
            orgrarad = astro_utils.stringToRad(rastr, True)
            orgdecrad = astro_utils.stringToRad(decstr, False)
            break
    if mjd < 0 or orgrarad < 0:
        print "Vex file " + vexfile + " not parsed properly"
        sys.exit()
    startdayfrac = float(nominalstart[:2])/24.0 + float(nominalstart[3:5])/1440.0 + \
                   float(nominalstart[6:8])/86400.0
    stopdayfrac = float(nominalstop[:2])/24.0 + float(nominalstop[3:5])/1440.0 + \
                  float(nominalstop[6:8])/86400.0
    if stopdayfrac < startdayfrac:
        stopdayfrac += 1
    mjd += (startdayfrac+stopdayfrac)/2.0

    #sutable = uvdata.table('SU', 1)
    #found = False
    #for row in sutable:
    #    if row.source.strip() == srcname.strip():
    #        orgrarad = row.raepo*math.pi/180.0
    #        orgdecrad = row.decepo*math.pi/180.0
    #        found = True
    #        break
    #if not found:
    #    print "Couldn't find source " + srcname + " in SU table!"
    #    sys.exit()
    fitlines = open(fitfile).readlines()
    tmpout = open("tmpfit.in", "w")
    for line in fitlines:
        if "=" in line:
            tmpout.write(line)
        if "RA" in line:
            raline = line
        if "Dec" in line:
            decline = line
    tmpout.write("%.5f 0:0:0.0 0.1 0:0:0.0 0.1\n" % mjd)
    tmpout.close()
    os.system("pmpar tmpfit.in")
    elines = open("pmpar_e").readlines()
    if not len(elines) == 1:
        print "Some problem with pmpar - pmpar_e is too long!"
        sys.exit()
    rastr = raline.split('=')[-1].split()[0]
    decstr = decline.split('=')[-1].split()[0]
    newrarad = astro_utils.stringToRad(rastr, True)
    newdecrad = astro_utils.stringToRad(decstr, False)
    splite = elines[0].split()
    newrarad += float(splite[5])*math.pi/(180*60*60*1000*math.cos(newdecrad))
    newdecrad += float(splite[6])*math.pi/(180*60*60*1000)
    return newrarad, newdecrad, orgrarad, orgdecrad

####### ENSURE THAT NO HIGHER SEQNO FILES EXIST, THEN RETURN DESIRED OBJECT ####
def zapAndCreateUVData(srcname, klass, disk, seqno):
    for i in range(seqno, 100):
        outputdata = AIPSUVData(srcname,klass,disk,i)
        if outputdata.exists():
            outputdata.zap()
    outputdata = AIPSUVData(srcname,klass,disk,seqno)
    return outputdata

####### ATLOD FOR RPFITS FILES FROM THE LBA ####################################
def atlod(atfile, uvdata, sources):
    if uvdata.exists():
        if interaction.yesno("Delete existing UV dataset " + uvdata.name + \
                             "? (No will abort pipeline)"):
            uvdata.zap()
        else:
            sys.exit()
    atlod = AIPSTask('atlod', version = aipsver)
    atlod.datain      = atfile
    atlod.outdata     = uvdata
    atlod.sources[1:] = sources
    atlod.douvcomp    = -1
    atlod.aparm[1:]   = [0]
    atlod.aparm[3]    = 1
    atlod.aparm[7]    = 1
    atlod()

####### FITLD FOR NON-CORRELATOR UV FILES ######################################
def fitld_uvfits(uvfitsfile, aipsdata, sources):
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 0
    fitld.ncount = 1
    if sources != "":
        fitld.sources[1:len(sources)+1] = sources
    fitld.outdata = aipsdata
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits
    fitld()
    os.system("rm -f " + tempinfits)

####### FITLD FOR CORRELATOR UV FILES, NO WEIGHT DISCARDING ####################
def fitld_corr_noweightthresh(file, aipsdata, sources):
    fitld = AIPSTask('fitld', version = aipsver)
    splitfile = file.split(':')
    fitld.doconcat = 0
    if len(splitfile) > 1:
        fitld.doconcat = 1
    fitld.optype = ''
    fitld.ncount = 0
    fitld.dotable = 1
    fitld.douvcomp = 1
    fitld.qual = -1
    fitld.bchan = 0
    fitld.echan = 0
    fitld.bif = 0
    fitld.eif = 0
    fitld.clint = 0.166667
    fitld.digicor = 1
    fitld.wtthresh = 0.0
    fitld.sources[1:len(sources)+1] = sources
    fitld.outdata = aipsdata
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    for f in splitfile:
        fitld.datain = f
        fitld()

####### FITLD FOR CORRELATOR UV FILES, NO WEIGHT DISCARDING ####################
def fitld_vlba(file, aipsdata, sources, wtthreshhold=0.4, rdate='', cltablemin=0.166667):
    fitld = AIPSTask('fitld', version = aipsver)
    splitfile = file.split(':')
    fitld.doconcat = 0
    if len(splitfile) > 1:
        fitld.doconcat = 1
    print fitld.doconcat
    fitld.optype = ''
    fitld.ncount = 0
    fitld.dotable = 1
    fitld.douvcomp = 1
    fitld.qual = -1
    fitld.bchan = 0
    fitld.echan = 0
    fitld.bif = 0
    fitld.eif = 0
    fitld.clint = cltablemin
    fitld.refdate = rdate
    fitld.digicor = 1
    fitld.wtthresh = wtthreshhold
    fitld.sources[1:len(sources)+1] = sources
    fitld.outdata = aipsdata
    fitld.antname[1] = 'VLBA'
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    for f in splitfile:
        fitld.datain = f
        fitld()

####### FITLD FOR CORRELATOR UV FILES ##########################################
def fitld_corr(file, aipsdata, sources, antennalist='', wtthreshold=0.0, rdate='', cltablemin=0.166667):
    fitld = AIPSTask('fitld', version = aipsver)
    splitfile = file.split(':')
    splitantenna = antennalist.split(',')
    fitld.doconcat = 0
    if len(splitfile) > 1:
        fitld.doconcat = 1
    fitld.optype = ''
    fitld.ncount = 0
    fitld.dotable = 1
    fitld.douvcomp = -1
    fitld.qual = -1
    fitld.bchan = 0
    fitld.echan = 0
    fitld.bif = 0
    fitld.eif = 0
    fitld.clint = cltablemin
    fitld.digicor = 1
    fitld.wtthresh = wtthreshold
    fitld.refdate = rdate
    fitld.sources[1:len(sources)+1] = sources
    fitld.outdata = aipsdata
    i = 1
    for a in splitantenna:
        fitld.antname[i] = a
        i += 1
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    for f in splitfile:
        fitld.datain = f
        fitld()

####### FITLD FOR IMAGE FILES ##################################################
def fitld_image(file, aipsdata):
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.dotable = 1
    fitld.douvcomp = 1
    fitld.digicor = -1
    fitld.datain = file
    fitld.outdata = aipsdata
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    fitld()  

####### SORT INTO TB ORDER #####################################################
def uvsrt(inputuvdata, outputuvdata):
    uvsrt = AIPSTask('uvsrt', version = aipsver)
    uvsrt.indata  = inputuvdata
    uvsrt.outdata = outputuvdata
    uvsrt.sort    = 'TB'
    uvsrt()
    indxr = AIPSTask('indxr', version = aipsver)
    indxr.indata = outputuvdata
    indxr.cparm[1:] = [0]
    indxr.cparm[3] = 0.1666666666666666667
    indxr()

####### SORT INTO TB ORDER #####################################################
def indxr(inputuvdata):
    indxr = AIPSTask('indxr', version = aipsver)
    indxr.indata = inputuvdata
    indxr.cparm[1:] = [0]
    indxr.cparm[3] = 0.16666666666666667
    indxr()

####### RENUMBER ANTENNAS TO CONFORM WITH NORMAL VLBA NUMBERING ################
def renumber_vlba(uvdata):
    reseq = AIPSTask('reseq', version = aipsver)
    reseq.indata = uvdata
    reseq.outdata = uvdata
    reseq.antennas[1:] = [0]
    reseq.baseline[1:] = [0]
    antennas = ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC']
    allocated = []
    count = 0
    unknowncount = 0
    antable = uvdata.table('AN', 1)
    for row in antable:
        print row
        a = row.anname.strip()
        antnum = row.nosta
        reseq.baseline[count+1] = antnum
        found = False
        for i in range(len(antennas)):
            if antennas[i] == a:
                 reseq.antennas[count+1] = i+1
                 allocated.append(a)
                 found = True
        if not found:
            unknowncount += 1
            print "Didn't know how to renumber %s! Setting it to %d" % \
                  (a, unknowncount)
            reseq.antennas[count+1] = len(antennas) + unknowncount

        count += 1
    #for i in range(len(antennas) - len(uvdata.antennas)):
    #    reseq.baseline[len(uvdata.antennas)+i+1] = len(uvdata.antennas)+i+1
    #    count = 0
    #    for a in antennas:
    #        count += 1
    #        if not a in allocated:
    #            reseq.antennas[len(uvdata.antennas)+i+1] = count
    #            allocated.append(a)
    #print reseq.baseline
    #print reseq.antennas
    reseq()
    print "uvdata has been resequenced - now:"
    for row in antable:
        print row

####### DBCON N CATALOG ENTRIES ################################################
def dbcon(inputuvdatas, outputuvdata):
    junkos = []
    errormsg = "ANTENNA tables don't match! Can't DBCON with doarray=1.  " + \
               "This is set so that there are no subarrays, and difmap " + \
               "can still process the data.  Write a manual dbcon to get " + \
               "around this if subarrays are ok."
    ocount = 1
    junko = AIPSUVData('JUNKO','JUNKO',1,1)
    if junko.exists():
        junko.zap()
    dbcon = AIPSTask('dbcon', version = aipsver)
    dbcon.reweight[1:] = [0]
    dbcon.dopos[1:][1:] = [0]
    dbcon.doarray = 1
    dbcon.fqtol = -1
    missingan1 = []
    inputuvdatas.sort(key=lambda a: len(a.antennas), reverse=True)
    for i in range(len(inputuvdatas)-1):
        for j in range(i+1,len(inputuvdatas)):
            antable1 = WizAIPSUVData(inputuvdatas[i]).table('AN', 1)
            antable2 = WizAIPSUVData(inputuvdatas[j]).table('AN', 1)
            if len(inputuvdatas[i].antennas) < len(inputuvdatas[j].antennas):
                print errormsg
                sys.exit()
            for a in inputuvdatas[j].antennas:
                index1 = -1
                index2 = -1
                for row in antable1:
                    if row.anname.strip() == a.strip():
                        index1 = row.nosta
                        break
                for row in antable2:
                    if row.anname.strip() == a.strip():
                        index2 = row.nosta
                        break
                if index1 < 0:
                    missingan1.append(a)
                if index1 != index2 and index1 >= 0 and index2 >= 0:
                    print a
                    print index1, index2
                    print inputuvdatas[i].antennas
                    for row in antable1:
                        print row
                    print inputuvdatas[j].antennas
                    for row in antable2:
                        print row
                    print errormsg
                    sys.exit()
                if index1 >= 0 and index2 >= 0:
                    for row1 in antable1:
                        if row1.anname.strip() == a.strip():
                            staxof = row1.staxof
                            x = row1.stabxyz[0]
                            y = row1.stabxyz[1]
                            z = row1.stabxyz[2]
                            break
                    for row in antable2:
                        if row.anname.strip() == a.strip():
                            print "UPDATING " + a
                            staxdiff = row.staxof - staxof
                            xdiff = row.stabxyz[0] - x
                            ydiff = row.stabxyz[1] - y
                            zdiff = row.stabxyz[2] - z
                            if math.fabs(staxdiff) > 0.5 or math.fabs(xdiff) > 0.5 or \
                               math.fabs(ydiff) > 0.5 or math.fabs(zdiff) > 0.5:
                                print row
                                print row1
                                print errormsg
                                sys.exit()
                            row.staxof = staxof
                            row.stabxyz[0] = x
                            row.stabxyz[1] = y
                            row.stabxyz[2] = z
                            row.update()
            for a in missingan1:
                for row in antable2:
                    if row.anname.strip() == a.strip():
                        antable1.append(row)
            antable1.close()
            antable2.close()
    # Merge the two antenna tables manually
    dbcon.indata  = inputuvdatas[0]
    dbcon.in2data = inputuvdatas[1]
    print inputuvdatas[0].antennas
    print inputuvdatas[1].antennas
    dbcon.outdata = junko
    dbcon()
    for uvdata in inputuvdatas[2:]:
        ocount += 1
        dbcon.indata  = junko
        dbcon.in2data = uvdata
        print junko.antennas
        print uvdata.antennas
        junko = AIPSUVData('JUNKO','JUNKO',1,ocount)
        if junko.exists():
            junko.zap()
        dbcon.outdata = junko
        dbcon()
    print junko.antennas
    uvsrt = AIPSTask('UVSRT', version = aipsver)
    uvsrt.indata = junko
    uvsrt.outdata = outputuvdata
    uvsrt.sort = 'TB'
    uvsrt()
    for i in range(ocount):
        junko = AIPSUVData('JUNKO','JUNKO',1,i+1)
        junko.zap()

####### MERGE PULSAR BINS ######################################################
def mergebins(inputdatafiles, outputdatafile, noaverage):
    junk1 = AIPSUVData('JUNK1','JUNK1',1,1)
    junk2 = AIPSUVData('JUNK2','JUNK2',1,1)
    junko = AIPSUVData('JUNKO','JUNKO',1,1)
    if junk1.exists():
        junk1.zap()
    if junk2.exists():
        junk2.zap()
    if junko.exists():
        junko.zap()
    dbcon = AIPSTask('dbcon', version = aipsver)
    dbcon.reweight[1:] = [0]
    dbcon.dopos[1:][1:] = [0]
    dbcon.doarray = 1
    dbcon.fqtol = -1
    fitld_corr(inputdatafiles[0], junk1, [''])
    for infile in inputdatafiles[1:]:
        fitld_corr(inputdatafiles[0], junk2, [''])
        dbcon.indata  = junk1
        dbcon.in2data = junk2
        dbcon.outdata = junko
        dbcon()
        junk1.zap()
        junk2.zap()
        junko.rename('JUNK1','JUNK1',1) #Don't know why this only wants 3...
        junko = AIPSUVData('JUNKO','JUNKO',1,1)
    uvsrt = AIPSTask('UVSRT', version = aipsver)
    uvsrt.indata = junk1
    uvsrt.outdata = junk2
    uvsrt.sort = 'TB'
    uvsrt()
    if noaverage:
        writedata(junk2, outputdatafile, False)
    else:
        uvavg = AIPSTask('UVAVG', version = aipsver)
        uvavg.indata = junk2
        uvavg.outdata = junko
        uvavg()
        writedata(junko, outputdatafile, False)
        junko.zap()
    junk1.zap()
    junk2.zap()

####### LOAD AN ANTAB FORMAT FILE ##############################################
def antab(uvdata, antabfile, tyver, gcver):
    antab = AIPSTask('antab', version = aipsver)
    antab.indata = uvdata
    antab.calin = antabfile
    antab.tyver = tyver
    antab.gcver = gcver
    antab()

####### DIFFERENCE TWO DATASETS ################################################
def diffuv(avgdata1, avgdata2, diffdata):
    difuv = AIPSTask('difuv', version = aipsver)
    difuv.indata  = avgdata1
    difuv.in2data = avgdata2
    difuv.outdata = diffdata
    difuv.solint  = 0
    difuv.optype  = ''
    difuv()

####### TABLE/FREQUENCY/POLARISATION MERGING AND SORTING #######################
def vlbafix(uvdataset):
    vlbamcal = AIPSTask('vlbamcal', version = aipsver)
    vlbamcal.indata = uvdataset
    vlbamcal()

    vlbafix = AIPSTask('vlbafix', version = aipsver)
    vlbafix.indata   = uvdataset
    vlbafix.clint    = 0.25
    vlbafix.subarray = 0
    vlbafix()

####### PARALLACTIC ANGLE CORRECTION ###########################################
def clcor_pang(uvdataset, clversion):
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdataset
    clcor.gainver = clversion
    clcor.gainuse = clversion+1
    clcor.opcode = 'PANG'
    clcor.clcorprm[1] = 1
    clcor.clcorprm[2:] = [0]
    clcor()

####### SOURCE POSITION CORRECTION #############################################
def shift_source(uvdataset, source, rashift, decshift, clversion):
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdataset
    clcor.sources[1] = source
    clcor.gainver = clversion
    clcor.gainuse = clversion
    clcor.opcode = 'ANTC'
    clcor.clcorprm[1:] = [0]
    clcor.clcorprm[5] = rashift/1000.0
    clcor.clcorprm[6] = decshift/1000.0
    clcor.clcorprm[8] = 0.0
    clcor.clcorprm[9] = 0.0
    clcor.infile = ""
    clcor()

####### ANTENNA POSITION CORRECTION ############################################
def shift_antenna(uvdataset, antenna, xshift, yshift, zshift, clversion):
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdataset
    clcor.antenna[1] = antenna
    clcor.gainver = clversion
    clcor.gainuse = clversion
    clcor.opcode = 'ANTC'
    clcor.clcorprm[1:] = [0]
    clcor.clcorprm[1] = xshift
    clcor.clcorprm[2] = yshift
    clcor.clcorprm[3] = zshift
    clcor.clcorprm[5] = 0.0
    clcor.clcorprm[6] = 0.0
    clcor.clcorprm[8] = 0.0
    clcor.clcorprm[9] = 0.0
    clcor.infile = ""
    clcor()

####### UVFIX (SINGLE SOURCE) ##################################################
def uvfix(uvdataset, outdataset, rashift, decshift):
    uvfix = AIPSTask('uvfix', version = aipsver)
    uvfix.indata = uvdataset
    uvfix.outdata = outdataset
    uvfix.shift[1] = rashift/1000.0
    uvfix.shift[2] = decshift/1000.0
    uvfix.uvfixprm[1] = 1
    uvfix()

####### SETUP TEMP LOG DIRECTORY FOR A SPECIFIC IONEX MODEL ####################
def setup_iono(logdir, imod):
    templogdir = logdir
    if templogdir[-1] == '/':
        templogdir = templogdir[:-1]
    templogdir += '.saved/'
    os.system("mkdir " + templogdir)
    for f in os.listdir(logdir):
        if f.find(imod) >= 0:
            os.system("cp " + logdir + '/' + f + " " + templogdir)
    os.system("gunzip " + templogdir + "*.gz")
    os.system("gunzip " + templogdir + "*.Z") 
    return templogdir
    
####### IONOSPHERIC CORRECTIONS USING TECOR ####################################
def correct_iono(uvdataset, tecordirectory, clversion, follow=0.2):
    tecor = AIPSTask('tecor', version = aipsver)
    for t in uvdataset.tables:
        if 'CL' in t[1] and t[0] == clversion+1:
            if interaction.yesno("Delete existing CL table versiuon " + \
                                 str(clversion+1) + "? "):
                uvdataset.zap_table('CL', clversion+1)
            else:
                sys.exit()
            break
    tecor.indata = uvdataset
    files = sorted(os.listdir(tecordirectory))
    print files
    ionextypes = ["jplg","esag","codg","upcg","igsg","c1pg","c2pg","u2pg","e1pg"]
    selectedionextype = ""
    for ionextype in ionextypes:
        for filename in files:
            if ionextype in filename and filename[-1] == 'i':
                selectedionextype = ionextype
                break
        if selectedionextype != "":
            break
    if selectedionextype == "":
        print "Couldn't find a suitable unzipped ionex file - aborting!"
        sys.exit()
    numfiles = 0
    for filename in files:
        if selectedionextype in filename and filename[-1:] == 'i':
            if numfiles == 0:
                tecor.infile = tecordirectory + filename
            elif (filename.split('.')[-1] < tecor.infile.split('.')[-1]):
                tecor.infile = tecordirectory + filename
            elif (filename.split('.')[-1] == tecor.infile.split('.')[-1]) and \
                 ((tecordirectory + filename) < tecor.infile):
                tecor.infile = tecordirectory + filename
            numfiles = numfiles + 1
    tecor.nfiles = numfiles
    tecor.subarr = 0
    tecor.antennas[1:] = [0]
    tecor.gainver = clversion
    tecor.gainuse = (clversion+1)
    tecor.aparm[2:] = [0]
    tecor.aparm[1] = 1
    tecor.aparm[2] = follow
    if numfiles > 0:
        print 'Running TECOR with ' + str(numfiles) + ' files - first one is ' + tecor.infile
        print "Follow is " + str(follow)
        tecor() #Do it!
    else:
        print 'Could not find any suitable ionex files - aborting!!!'
        sys.exit()

####### GENERATE DIFFERENTIAL TECOR CORRECTIONS FROM CL TABLES AND COPY IN  ####
def insert_differential_iono(calibsuvdata, calinver, targetuvdata, targetinver, 
                             targetoutver):
    wizcaldata = WizAIPSUVData(calibsuvdata)
    wiztargetdata = WizAIPSUVData(targetuvdata)
    cltable1 = calibsuvdata.table('CL', calinver)
    cltable2 = targetuvdata.table('CL', targetinver)
    cdisp1 = []
    cdisp2 = []
    tdisp1 = []
    tdisp2 = []
    print "Storing iono corrections for calibrator"
    for row in cltable1:
        cdisp1.append(row.disp_1)
        cdisp2.append(row.disp_2)
    print "Storing iono corrections for target"
    for row in cltable2:
        tdisp1.append(row.disp_1)
        tdisp2.append(row.disp_2)
    wizcltable = wiztargetdata.table('CL', targetoutver)
    increment = len(cltable1)/10
    count = 0
    for row in wizcltable:
        if count % increment == 0:
            print str((10*count)/increment) + "% done"
        row.disp_1 = tdisp1[count] - cdisp1[count]
        row.disp_2 = tdisp2[count] - cdisp2[count]
        row.update()
        count += 1
    wizcltable.close()
'''
####### DIFFERENCE TWO IONOSPHERIC SOLUTIONS AND WRITE TO A FILE ###############
def diff_iono(calibuvdata, targetuvdata, tabledir, clversion, ionexmodel, 
              follow):
    wizcaldata = WizAIPSUVData(calibuvdata)
    wiztargetdata = WizAIPSUVData(targetuvdata)
    cltable1 = calibuvdata.table('CL', clversion)
    cltable2 = targetuvdata.table('CL', clversion)
    #newcltable = wiztargetdata.attach_table('CL', clversion+1, 
    #                                        no_term=cltable2.keywords['NO_TERM'])
    #newcltable.keywords['NO_ANT'] = cltable2.keywords['NO_ANT']
    #newcltable.keywords['NO_POL'] = cltable2.keywords['NO_POL']
    #newcltable.keywords['NO_IF'] =  cltable2.keywords['NO_IF']
    cdisp1 = []
    cdisp2 = []
    tdisp1 = []
    tdisp2 = []
    cantno = []
    tantno = []
    ctime = []
    ttime = []
    print "Storing iono corrections for calibrator"
    for row in cltable1:
        cdisp1.append(row.disp_1)
        cdisp2.append(row.disp_2)
        cantno.append(row.antenna_no)
        ctime.append(row.time)
    print "Storing iono corrections for target"
    for row in cltable2:
        tdisp1.append(row.disp_1)
        tdisp2.append(row.disp_2)
        tantno.append(row.antenna_no)
        ttime.append(row.time)
    print "Updating cl table, version " + str(clversion)
    wizcltable = wiztargetdata.table('CL', clversion)
    toffmul = 0
    coffmul = 0
    increment = len(cltable1)/10
    cl1len = len(cltable1)
    cl2len = len(cltable2)
    if cl2len > cl1len:
        tablelendiff = cl2len - cl1len
        print "Warning - CL table lengths differ! Will skip over target rows as necessary"
        toffmul = 1
    if cl1len > cl2len:
        tablelendiff = cl1len - cl2len
        print "Warning - CL table lengths differ! Will skip over cal rows as necessary"
        coffmul = 1
    count = 0
    offset = 0
    for row in wizcltable:
        if count % increment == 0:
            print str((10*count)/increment) + "% done"
        tcount = toffmul*offset + count
        ccount = coffmul*offset + count
        while cantno[ccount] != tantno[tcount] or \
              ctime[ccount] != ttime[tcount]:
            print "skipping because cal ant is " + str(cantno) + \
                  "target ant is " + str(tantno) + ", cal time is " + \
                  str(ctime[ccount]) + ", target time is " + str(ttime[tcount])
            offset += 1
            print "offset is now " + str(offset)
            tcount = toffmul*offset + count
            ccount = coffmul*offset + count
        row.disp_1 = tdisp1[tcount] - cdisp1[ccount]
        row.disp_2 = tdisp2[tcount] - cdisp2[ccount]
        row.update()
#        newcltable.append(row2)
        count += 1
#    newcltable.close()
    wizcltable.close()
    print "Eventual offset was " + str (offset) + ", while initial difference " + \
           "in table lengths was " + str(tablelendiff)
    outfilename = tabledir + ('differentialtecor.%s.%d.CL' % \
                              (ionexmodel, int(follow*100)))
    writetable(targetuvdata, 'CL', clversion, outfilename)
    calibuvdata.table('CL', clversion).zap()
    targetuvdata.table('CL', clversion).zap()

####### COPY IN DIFFERENTIAL TECOR CORRECTIONS TO A CL TABLE ###################
def insert_differential_iono(targetuvdata, dtecorfile, clversion):
    loadtable(targetuvdata, dtecorfile, clversion+1)
    wizuvdata = WizAIPSUVData(targetuvdata)
    cltable = wizuvdata.table('CL', clversion)
    difftable = targetuvdata.table('CL', clversion+1)
    for row, drow in zip(cltable, difftable):
        row.disp_1 = row.disp_1 + drow.disp_1
        row.disp_2 = row.disp_2 + drow.disp_2
        row.update()
    cltable.close()
    targetuvdata.table('CL', clversion+1).zap()
'''
####### EOP CORRECTIONS USING CLCOR ############################################
def correct_eops(uvdataset, eopsdirectory, clversion):
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdataset
    clcor.gainver = clversion
    clcor.gainuse = clversion + 1
    clcor.infile = eopsdirectory + '/usno_finals.erp'
    clcor.opcode = 'EOPS'
    clcor.clcorprm[1] = 1
    clcor.clcorprm[2] = 5
    print "About to run clcor on gainver %d, gainuse %d" % (clcor.gainver, clcor.gainuse)
    clcor()

####### FINDING THE MAXIMUM TABLE NUMBER #######################################
def get_max_tableno(uvdataset, tabletype):
    tables = uvdataset.tables
    maxtableversion = 0
    for table in tables:
	if table[1][-2:] == tabletype:
            if table[0] > maxtableversion:
                maxtableversion = table[0]
    return maxtableversion

####### LOADING A TABLE ########################################################
def loadtable(uvdataset, tablename, tableversion):
    tables = uvdataset.tables
    tabletype = tablename.split('.')[-1].upper()
    for table in tables:
        if table[1][-2:] == tabletype and table[0] >= tableversion:
            if yesno("Can i delete " + table[1][-2:] + " table version " + \
                     str(table[0]) + "? (No will abort pipeline)"):
                uvdataset.table(tabletype, table[0]).zap()
            else:
                sys.exit()
    tbin = AIPSTask('tbin', version = aipsver)
    tbin.outdata = uvdataset
    tbin.intext = tablename
    tbin()

####### WRITING A TABLE ########################################################
def writetable(uvdataset, tabletype, tableversion, filename):
    tbout = AIPSTask('tbout', version = aipsver)
    tbout.indata  = uvdataset
    tbout.inext   = tabletype
    tbout.invers  = tableversion
    tbout.outtext = filename
    tbout()

####### ZAPPING A TABLE ########################################################
def deletetable(uvdataset, tabletype, tableversion):
    try:
        uvdataset.table(tabletype, tableversion).zap()
        return True
    except IOError:
        print "Apparently " + tabletype + " version " + str(tableversion) + \
              " didn't exist..."
        return False

####### SPLIT OFF ALL TARGETS, DBCON TOGETHER AS NEEDED ########################
def splitanddbcon(uvdatasets, numtargets, numfields, numfieldsources, 
                  fieldsourcenames, clversion, aipsclass, limitedtargetlist, 
                  doconcat, skipping=False, doband=True):
    targetnames    = []
    targetcount    = []
    orgtargetindex = []
    aipsnames      = []
    alreadycleared = []
    for i in range(numtargets):
        if not (len(limitedtargetlist) == 1 and limitedtargetlist[0] == '') \
           and not fieldsourcenames[i][0][0] in limitedtargetlist:
            for j in range(numfields[i]):
                for k in range(numfieldsources[i][j]):
                    if not fieldsourcenames[i][j][k] in targetnames:
                        targetnames.append(fieldsourcenames[i][j][k])
                        targetcount.append(0)
                        orgtargetindex.append(i)
            continue
        if not skipping:
            for c in range(50): #Clear any old split catalog entries
                for j in range(numfields[i]):
                    for target in fieldsourcenames[i][j]:
                        aipsname = target[:12]
                        if not aipsname in alreadycleared:
                            if c == 1:
                                print "About to delete " + aipsname + " version " + str(c)
                            splitdata = AIPSUVData(aipsname, aipsclass, 1, c)
                            if splitdata.exists():
                                splitdata.zap()
                            elif c == 1:
                                print "Doesn't exist"
                            if len(target) >= 16: # Crappy old style name
                                aipsname = target[4:16]
                            splitdata = AIPSUVData(aipsname, aipsclass, 1, c)
                            if splitdata.exists():
                                splitdata.zap()
                            if c == 49:
                                alreadycleared.append(aipsname)
        for j in range(numfields[i]):
            for k in range(numfieldsources[i][j]):
                aipsname = fieldsourcenames[i][j][k][:12]
                if len(fieldsourcenames[i][j][k]) >= 16: # Crappy old style name
                    aipsname = fieldsourcenames[i][j][k][4:16]
                if not fieldsourcenames[i][j][k] in targetnames:
                    targetnames.append(fieldsourcenames[i][j][k])
                    targetcount.append(0)
                    orgtargetindex.append(i)
                    targetindex = len(targetnames)-1
                else:
                    targetindex = targetnames.index(fieldsourcenames[i][j][k])
                targetcount[targetindex] += 1
                if not skipping:
                    try:
                        splittoseq(uvdatasets[k], clversion, aipsclass, fieldsourcenames[i][j][k], 
                                   targetcount[targetindex], True, doband)
                        if len(fieldsourcenames[i][j][k]) >= 16: # Crappy old style name
                            splitdata = AIPSUVData(fieldsourcenames[i][j][k][:12], aipsclass, 1,
                                                   targetcount[targetindex])
                            splitdata.rename(aipsname, aipsclass, targetcount[targetindex])
                    except RuntimeError:
                        print "No valid visibilities for " + fieldsourcenames[i][j][k]
                        targetcount[targetindex] -= 1
    alreadyconcatted = []
    for i in range(len(targetnames)):
        print "Doing " + targetnames[i] + " which has a count of " + str(targetcount[i])
        if targetcount[i] > 1 and doconcat and not skipping:
            aipsname = targetnames[i][:12]
            if len(targetnames[i]) >= 16: # Crappy old style name
                aipsname = targetnames[i][4:16]
            inputuvdatas = []
            if not aipsname in alreadyconcatted:
                print "Concatenating " + aipsname
                for j in range(targetcount[i]):
                    inputuvdatas.append(AIPSUVData(aipsname, aipsclass, 1, j+1))
                outputuvdata = AIPSUVData(aipsname, aipsclass, 1, targetcount[i]+1)
                dbcon(inputuvdatas, outputuvdata)
                for j in range(targetcount[i]):
                    inputuvdatas[j].zap()
                outputuvdata.rename(aipsname, aipsclass, 1)
                alreadyconcatted.append(aipsname)
        else:
            aipsname = targetnames[i][:12]
            if len(targetnames[i]) >= 16: # Crappy old style name
                aipsname = targetnames[i][4:16]
        aipsnames.append(aipsname)
    return targetnames, aipsnames, targetcount, orgtargetindex


####### SPLIT OFF A TIMERANGE, AVERAGE IN FREQUENCY (DEFAULT NO CORRECTIONS) ###
def splat(uvdataset, chanavg, timerange, outdataset, clversion=-1, solint=0):
    splat = AIPSTask('splat', version = aipsver)
    splat.indata = uvdataset
    splat.timerang[1:] = timerange
    splat.aparm[1] = 3
    splat.aparm[5] = 1
    splat.channel = chanavg
    splat.flagver = 1
    splat.chinc = chanavg
    splat.outdata = outdataset
    if clversion > 0:
        splat.docal = 2
        splat.gainuse = clversion
    splat.solint = solint
    splat()

####### AVERAGE IN TIME ########################################################
def uvavg(uvdataset, avgsecs, outdataset):
    uvavg = AIPSTask('UVAVG', version = aipsver)
    uvavg.indata = uvdataset
    uvavg.outdata = outdataset
    uvavg.yinc = avgsecs - 0.1
    uvavg()

####### AUTOCORRELATION CORRECTION #############################################
def accor(uvdataset, solint=-2):
    accor = AIPSTask('accor', version = aipsver)
    accor.indata = uvdataset
    accor.solint = solint
    accor()

####### AMPLITUDE CALIBRATION ##################################################
def apcal(uvdataset, snver):
    apcal = AIPSTask('apcal', version = aipsver)
    apcal.indata = uvdataset
    apcal.snver = snver
    apcal.opcode = ''
    apcal.aparm[1] = 1
    apcal.dofit[1:] = [-1]
    apcal.solint = 0
    apcal.invers = 1
    apcal()

###### SN TABLE MEDIAN WINDOW FILTERING ########################################
def snmwfclip(uvdataset, timemins, ampdev, phasedev, delaydev, ratedev, snver, 
              refant, coherentratesmoothing=False, ratesmoothingmins=3.0):
    snsmo = AIPSTask('snsmo', version = aipsver)
    snsmo.indata = uvdataset
    snsmo.samptype = 'MWF'
    snsmo.bparm[1:] = [0]
    snsmo.doblank = -1
    snsmo.dobtween = 1
    snsmo.cparm[1:] = [0]
    if ampdev > 0.0:
        snsmo.cparm[1] = timemins/60.0
        snsmo.cparm[6] = ampdev
    if phasedev > 0.0:
        snsmo.cparm[2] = timemins/60.0
        snsmo.cparm[7] = phasedev
    if delaydev > 0.0:
        snsmo.cparm[4] = timemins/60.0
        snsmo.cparm[9] = delaydev
    if ratedev > 0.0:
        snsmo.cparm[3] = timemins/60.0
        snsmo.cparm[8] = ratedev
    if coherentratesmoothing:
        snsmo.samptype='BOX'
        snsmo.smotype='VLRI'
        snsmo.bparm[3] = ratesmoothingmins/60.0
    snsmo.inver = snver
    snsmo.outver = snver+1
    snsmo.refant = refant
    snsmo()

###### SN TABLE SMOOTHING ######################################################
def snsmo(uvdataset, smotype, timemins, ampdev, phasedev, delaydev, snver,
          refant, passthrough=False, ratedev=0.0):
    snsmo = AIPSTask('snsmo', version = aipsver)
    snsmo.indata = uvdataset
    snsmo.samptype = 'MWF'
    snsmo.bparm[1:] = [0]
    snsmo.doblank = -1
    if passthrough:
        snsmo.doblank = 1
        snsmo.dobtween = 1
        snsmo.bparm[1] = timemins/60.0
        snsmo.bparm[2] = timemins/60.0
        snsmo.bparm[3] = timemins/60.0
        snsmo.bparm[4] = timemins/60.0
        snsmo.bparm[5] = timemins/60.0
    snsmo.smotype = smotype
    snsmo.cparm[1:] = [0]
    if ampdev > 0.0:
        snsmo.cparm[1] = timemins/60.0
        snsmo.cparm[6] = ampdev
    if phasedev > 0.0:
        snsmo.cparm[2] = timemins/60.0
        snsmo.cparm[7] = phasedev
    if delaydev > 0.0:
        snsmo.cparm[4] = timemins/60.0
        snsmo.cparm[9] = delaydev
    if ratedev > 0.0:
        snsmo.cparm[3] = timemins/60.0
        snsmo.cparm[8] = ratedev
    snsmo.inver = snver
    snsmo.outver = snver+1
    snsmo.refant = refant
    snsmo()

##### SN TABLE EDITING #########################################################
def snedt(uvdataset, snver):
    snedt = AIPSTask('snedt', version = aipsver)
    snedt.indata = uvdataset
    snedt.inext = 'sn'
    snedt.invers = snver
    snedt.dodelay = 1
    snedt.freqid = -1
    snedt.dotwo = 1
    snedt()

##### COPY A TABLE #############################################################
def tacop(indataset, tabletype, inver, outdataset, outver):
    tacop = AIPSTask('tacop', version = aipsver)
    tacop.indata = indataset
    tacop.inext = tabletype
    tacop.ncount = 1
    tacop.invers = inver
    tacop.outdata = outdataset
    tacop.outvers = outver
    tacop()

##### Apply a SN table and generate a CL table #################################
def applysntable(uvdataset, snversion, interpol, clversion, refant, 
                 sourcelist=[], opcode = 'CALI', cutoff=0):
    clcal = AIPSTask('clcal', version = aipsver)
    clcal.indata = uvdataset
    clcal.sources[1:] = clcal.soucode = clcal.calsour[1:] = ''
    for i in range(len(sourcelist)):
        clcal.sources[i+1] = sourcelist[i]
    clcal.opcode = opcode
    clcal.interpol = interpol
    clcal.cutoff = cutoff
    clcal.samptype = ''
    clcal.bparm[1:] = [0]
    clcal.doblank = 1
    clcal.dobtween = -1
    clcal.smotype = ''
    clcal.snver = snversion
    clcal.inver = snversion
    clcal.gainver = clversion
    clcal.gainuse = clversion+1
    clcal.refant = refant
    clcal()

##### Merge SN tables ##########################################################
def mergesntables(uvdataset, startsnver, sncount, refant):
    clcal = AIPSTask('clcal', version = aipsver)
    clcal.indata = uvdataset
    clcal.opcode = 'MERG'
    clcal.snver = startsnver
    clcal.inver = startsnver + sncount - 1
    clcal.gainuse = startsnver + sncount
    clcal.refant = refant
    clcal()

def multipageps2tiledpng(directory, psfile, ncolumns):
    orgwd = os.getcwd()
    os.chdir(directory)
    pnmfiles = glob.glob("*.pnm")
    if len(pnmfiles) > 0:
        print "There exist pnm files in this directory (%s); can't make tiled png!" % directory
    os.system("pstopnm %s" % psfile)
    pnmfiles = glob.glob("*.ppm")
    pngfiles = []
    for pnmfile in pnmfiles:
        os.system("pnmtopng %s > %s.png" % (pnmfile, pnmfile[:-4]))
        pngfiles.append(pnmfile[:-4] + ".png")
    os.system("rm -f *.ppm")
    outputfile = psfile[:-(len(psfile.split('.')[-1])+1)] + '.png'
    if len(pngfiles) > 1:
        runline = "montage -tile %dx " % (ncolumns)
        for pngfile in pngfiles:
            runline += pngfile + " "
        runline += outputfile
        print runline
        os.system(runline)
        for pngfile in pngfiles:
            os.system("rm -f " + pngfile)
    else:
        os.system("mv %s %s" % (pngfiles[0], outputfile))
    os.chdir(orgwd)

##### Find the two timeranges corresponding to the two pulsar scan groups ######
def get_dataset_mjd_midpoint(uvdataset):
    listr = AIPSTask('listr', version = aipsver)
    listr.indata = uvdataset
    listr.optype = 'SCAN'
    listr.docrt = 0
    listr.outprint = 'junk'
    listr()
    listrjunk = open('junk')
    listrmsgs = listrjunk.readlines()
    listrjunk.close()
    os.remove('junk')
    validrange = False
    timerangelines = []
    for line in listrmsgs:
        if "Scan" in line and "Source" in line and "Qual" in line:
            validrange = True
        elif "Source summary" in line:
            validrange = False
        elif validrange and len(line.split()) > 7:
            # This is a line with a timerange
            timerangelines.append(line)
    if len(timerangelines) == 0:
        print "Didn't find any timerange lines!"
        sys.exit()
    starttimer = timerangelines[0].strip().split()[5]
    endtimer   = timerangelines[-1].strip().split()[7]

    startdd = int(starttimer.split('/')[0])
    enddd   = int(endtimer.split('/')[0])
    starthh = int(starttimer.split('/')[1].split(':')[0])
    endhh   = int(endtimer.split('/')[1].split(':')[0])
    startmm = int(starttimer.split('/')[1].split(':')[1])
    endmm   = int(endtimer.split('/')[1].split(':')[1])
    startss = int(starttimer.split('/')[1].split(':')[2])
    endss   = int(endtimer.split('/')[1].split(':')[2])
 
    startseconds = startdd*86400 + starthh*3600 + startmm*60 + startss
    endseconds   = enddd*86400   + endhh*3600   + endmm*60   + endss

    obsdate = uvdataset.header.date_obs
    year = int(obsdate[0:4])
    month = int(obsdate[5:7])
    day = int(obsdate[8:10])
    mjd = year*367 - int(7*(year + int((month + 9)/12))/4) + \
          int(275*month/9) + day - 678987

    return mjd + float(startseconds + endseconds)/(2.0*86400.0)

##### Find the two timeranges corresponding to the two pulsar scan groups ######
def get_pulsar_scangroup_timer(uvdataset, srcname, returnscanlist=False):
    listr = AIPSTask('listr', version = aipsver)
    listr.indata = uvdataset
    listr.optype = 'SCAN'
    listr.docrt = 0
    listr.outprint = 'junk'
    listr()
    listrjunk = open('junk')
    listrmsgs = listrjunk.readlines()
    listrjunk.close()
    os.remove('junk')
    validrange = False
    timerangelines = []
    for line in listrmsgs:
        if srcname in line: print line
        if "Scan" in line and "Source" in line and "Qual" in line:
            validrange = True
        elif "Source summary" in line:
            validrange = False
        elif validrange and srcname[:8] in line:
            # This is a line with a timerange
            timerangelines.append(line)
    if len(timerangelines) == 0:
        print "Didn't find any timerange lines!"
        sys.exit()

    if returnscanlist:
        toreturn = []
        for line in timerangelines:
            begin = line.strip().split()[5]
            end   = line.strip().split()[7]
            blist = [int(begin[0])]
            for b in begin.split('/')[-1].split(':'):
                blist.append(int(b))
            elist = [int(end[0])]
            for e in end.split('/')[-1].split(':'):
                elist.append(int(e))
            toreturn.append([blist, elist])
        return toreturn

    numscans = len(timerangelines)
    q1timer     = timerangelines[(numscans-1)/4].strip().split()[7]
    middletimer = timerangelines[(numscans-1)/2].strip().split()[7]
    q3timer     = timerangelines[(3*(numscans-1))/4].strip().split()[7]

    splittimes1 = q1timer.split('/')[-1].split(':')
    splittimes2 = middletimer.split('/')[-1].split(':')
    splittimes3 = q3timer.split('/')[-1].split(':')

    toreturn1 = [int(q1timer[0])]
    for s in splittimes1:
        toreturn1.append(int(s))
    toreturn2 = [int(middletimer[0])]
    for s in splittimes2:
        toreturn2.append(int(s))
    toreturn3 = [int(q3timer[0])]
    for s in splittimes3:
        toreturn3.append(int(s))


    return toreturn1, toreturn2, toreturn3

##### Find a timerange when the main calibrator is up ##########################
def get_ampcal_timer(uvdataset, srcname, scanno):
    timerang = []
    listr = AIPSTask('listr', version = aipsver)
    listr.indata = uvdataset
    listr.optype = 'SCAN'
    listr.docrt = 0
    listr.outprint = 'junk'
    listr()
    listrjunk = open('junk')
    listrmsgs = listrjunk.readlines()
    listrjunk.close()
    os.remove('junk')
    
    atline = 0
    while len(listrmsgs[atline].split()) < 5 or \
              not listrmsgs[atline].split()[0] == "Scan":
        atline = atline + 1
    for i in range(scanno):
        atline = atline + 1
        while atline < len(listrmsgs):
            if len(listrmsgs[atline].split()) <= 1:
                atline += 1
            else:
                if listrmsgs[atline].split()[1] == srcname:
                     print "Found " + srcname + " for scan " + \
                           str(i+1) + "/" + str(scanno)
                     break
                #print "Source was " + listrmsgs[atline].split()[1]
                #print listrmsgs[atline+1]
                atline += 1
    if atline >= len(listrmsgs):
        print "Couldn't find scan " + str(scanno) + " on source " + srcname
        print "Sorry, but I'll have to abort the pipeline"
        sys.exit()
    print listrmsgs[atline]
    starttime = listrmsgs[atline][41:51]
    endtime = listrmsgs[atline][56:66]
    print "Starttime of ampcal scan is " + starttime
    print "Endtime of ampcal scan is " + endtime

    startday  = int(starttime[0])
    starthour = int(starttime[2:4])
    startmin  = int(starttime[5:7])
    startsec  = int(starttime[8:10])
    endday    = int(endtime[0])
    endhour   = int(endtime[2:4])
    endmin    = int(endtime[5:7])
    endsec    = int(endtime[8:10])
    #Clip some, if we think we can get away with it
    timediff  = (endday-startday)*86400 + (endhour-starthour)*3600 + (endmin-startmin)*60 + endsec - startsec
    if timediff >= 70:
        endsec   -= 10
        timediff -= 10
        if endsec < 0:
            endsec += 60
            endmin -= 1
            if endmin < 0:
                endmin += 60
                endhour -= 1
                if endhour < 0:
                    endhour += 24
                    endday -= 1
        tosubtract = 0
        if timediff > 600:
            tosubtract = 240
        elif timediff > 360:
            tosubtract = 120
        elif timediff > 180:
            tosubtract = 60
        elif timediff > 120:
            tosubtract = 30
        elif timediff > 90:
            tosubtract = 20
        elif timediff > 70:
            tosubtract = 10
        startday  = endday
        starthour = endhour
        startmin  = endmin
        startsec  = endsec - (timediff - tosubtract)
        while startsec < 0:
            startmin -= 1
            startsec += 60
        while startmin < 0:
            startmin += 60
            starthour -= 1
        if starthour < 0:
            starthour += 24
            startday -= 1
    timerang.append(startday)
    timerang.append(starthour)
    timerang.append(startmin)
    timerang.append(startsec)
    timerang.append(endday)
    timerang.append(endhour)
    timerang.append(endmin)
    timerang.append(endsec)
    return timerang

##### Do the pulse cal corrections #############################################
def pccor(uvdataset, srcname, snversion, ampcalscanno, refant):
    pccor = AIPSTask('pccor', version = aipsver)
    pccor.indata = uvdataset
    pccor.timerang[1:] = get_ampcal_timer(uvdataset, srcname, ampcalscanno)
    print "Running PCCOR over timerange:"
    print pccor.timerang[1:]
    pccor.snver = snversion
    pccor.inver = 0
    pccor.refant = refant
    pccor.cutoff = 0
    pccor.delcorr = 0
    pccor()

##### Unflag a specified antenna ###############################################
def unflagantenna(uvdataset, inflagver, outflagver, antenna):
    taflg = AIPSTask('taflg', version = aipsver)
    taflg.indata = uvdataset
    taflg.inext = 'FG'
    taflg.invers = inflagver
    taflg.outvers = outflagver
    taflg.optype = '='
    taflg.aparm[1] = 4
    taflg.aparm[2] = 1
    taflg.bparm[1] = 1
    taflg.bparm[2] = antenna
    taflg()

##### User-specified additional flags ##########################################
def unflag(uvdataset, outflagver, antenna, timerange=[-1,-1,-1,-1,-1,-1,-1,-1],
           reason=' ',source='-',freqid=-1,bif=-1,eif=-1,bchan=-1,
           echan=-1,stokes='----'):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.outfgver = outflagver
    uvflg.antenna[1] = -antenna
    uvflg.timerang[1:9] = timerange
    uvflg.freqid = freqid
    uvflg.aparm[1:] = [0]
    uvflg.reason = reason
    uvflg.bif = bif
    uvflg.eif = eif
    uvflg.bchan = bchan
    uvflg.echan = echan
    uvflg.stokes = stokes
    uvflg.sources[1] = source    
    uvflg.opcode='UFLG'
    uvflg()

##### User-specified additional flags ##########################################
def userflag(uvdataset, outflagver, filename):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = filename
    uvflg.outfgver = outflagver
    uvflg()

##### User-specified additional flags ##########################################
def flagedgechannels(uvdataset, outflagver):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.outfgver = outflagver
    uvflg.bchan = 1
    uvflg.echan = 1
    uvflg()
    try:
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth[0]/\
                uvdataset.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth/\
                uvdataset.table('FQ', 1)[0].ch_width
    uvflg.bchan = nchan
    uvflg.echan = nchan
    uvflg()

##### Elevation-based flagging #################################################
def elevationflag(uvdataset, outflagver, elevationlimit):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = ''
    uvflg.outfgver = outflagver
    uvflg.aparm[1] = -90
    uvflg.aparm[2] = elevationlimit
    uvflg.aparm[3:] = [0]
    uvflg.opcode = 'FLAG'
    uvflg.reason = 'ELEVATION'
    uvflg()

##### Shadowing flagging #################################################
def shadowflag(uvdataset, outflagver, shadowdiameter, xtalkbl):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = ''
    uvflg.outfgver = outflagver
    uvflg.aparm[1:] = [0]
    uvflg.aparm[5] = shadowdiameter  # > 0 flag for shadowing; shadow diameter in m
    uvflg.aparm[6] = xtalkbl         # flag for cross-talk; baseline (BL) in m
    uvflg.opcode = 'FLAG'
    uvflg.reason = 'SHADOWING'
    uvflg()

##### Fringe rate flagging #####################################################
def fringerateflag(uvdataset, outflagver, suppressionfactor):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = ''
    uvflg.outfgver = outflagver
    uvflg.aparm[1:] = [0]
    uvflg.aparm[3] = suppressionfactor
    uvflg.opcode = 'FLAG'
    uvflg.reason = 'FRINGERATE'
    uvflg()

##### Read a list of times to quack from a file and do it ######################
def quackfromfile(uvdata, quackfile):
    if not os.path.exists(quackfile):
        print "Quack file listing " + quackfile + " does not exist - aborting!"
        sys.exit()
    input = open(quackfile)
    lines = input.readlines()
    input.close()
    for line in lines:
        splitline = line.split()
        if not len(splitline) == 3:
            print "Bad quack line " + line + " - will ignore!"
        count = 1
        found = False
        for a in uvdata.antennas:
            if a.strip() == splitline[0]:
                found = True
                break
            count += 1
        if found:
            quack(uvdata, [count], splitline[1], float(splitline[2]))
        else:
            print "Bad quack antenna " + splitline[0] + " - will ignore!"

##### Run quack on the specified antennas for the specified time range #########
def quack(uvdata, antennas, type, seconds):
    quack = AIPSTask('quack', version = aipsver)
    quack.indata = uvdata
    quack.sources[1:] = ''
    quack.timerang[1:] = [0,0,0,0,0,0,0,0]
    quack.antennas[1:] = [0]
    for i in range(len(antennas)):
        quack.antennas[i+1] = antennas[i]
    quack.flagver = 1
    quack.opcode = type #E.g. BEG OR END OR ENDB
    quack.reason = 'QUACK'
    quack.aparm[1] = 0
    quack.aparm[2] = seconds/60.0
    quack.aparm[3] = 0.5
    quack()

##### Fringe fit without calibrating weights ###################################
def fring_noweights(uvdataset, snversion, clversion, solintmins, inttimesecs, 
                    srcname,refant,doband=False,snr=7.5):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = uvdataset
    fring.docal = 100
    fring.gainuse = clversion
    fring.calsour[1] = srcname
    fring.flagver = 1
    fring.doband = -1
    if doband:
        fring.doband = 1
        fring.bpver = 1
    fring.refant = refant
    fring.aparm[1:] = [3,0,0,0,0,2,snr,0,1,0]
    fring.dparm[1:] = [0,600,40,inttimesecs,0,0,0,0]
    print "Inttime in seconds is " + str(inttimesecs)
    fring.solint = solintmins
    fring.snver = snversion
    print "Just set snver to " + str(fring.snver)
    fring()

##### Leakage calculation #################################################
def leakagecalc(uvdataset, lpcalsource, modelimage, lpcaltablefilename, \
                refant, ampcalmins, phasecalmins, leakagecalscan, clversion, \
                doband = False, outputfilename = "", uvrange=[0,0], weightit=0,
                ifmultiplier=0):
    # First split the leakage cal source
    splitdata = AIPSUVData("LPCAL", "JUNK", 1, 1)
    if splitdata.exists():
        splitdata.zap()
    for i in range(10):
        splitdata = AIPSUVData(lpcalsource, "LPCAL", 1, i)
        if splitdata.exists():
            splitdata.zap()
        splitdata = AIPSUVData(lpcalsource, "LPJNK", 1, i)
        if splitdata.exists():
            splitdata.zap()
        splitdata = AIPSUVData("LPCAL", "CALIB", 1, i)
        if splitdata.exists():
            splitdata.zap()
    splitdata = AIPSUVData(lpcalsource, "LPCAL", 1, 1)
    split2data = AIPSUVData(lpcalsource, "LPJNK", 1, 1)
    multidata = AIPSUVData("LPCAL", "MULTI", 1, 1)
    #morifdata = AIPSUVData("LPCAL", "MORIF", 1, 1)
    if split2data.exists():
        split2data.zap()
    if multidata.exists():
        multidata.zap()
    #if morifdata.exists():
    #    morifdata.zap()
    domulti = False
    beginif = -1
    endif = -1
    combineifs = False
    dopol = 0
    dolocalband = True
    if doband:
        dolocalband = False

    # Then run a bandpass cal on the leakage source, if none has been applied yet
    if dolocalband:
        if tableexists(uvdataset, "BP", -1):
            print "Not allowed to do local bandpasscal if a BP table already exists!"
            sys.exit()
        bpass(uvdataset, lpcalsource, clversion, leakagecalscan, modelimage, 0) # 0 = No dopol yet

    splittoseq(uvdataset, clversion, "LPCAL", lpcalsource, 1, domulti, True, # We're doing a bandpass one way or another
               beginif, endif, combineifs, dopol)
    if dolocalband: 
       uvdataset.zap_table("BP", -1) 
    splitdata.rename("LPCAL", "JUNK", 1) 

    ## Run MORIF if requested
    #if ifmultiplier > 1:
    #    morif = AIPSTask("MORIF", version = aipsver)
    #    morif.indata = splitdata
    #    morif.outdata = morifdata
    #    morif.npiece = ifmultiplier
    #    morif()
    #    splitdata.zap()
    #    splitdata = morifdata

    # Convert to multi
    multi = AIPSTask("MULTI", version = aipsver)
    multi.indata = splitdata
    multi.outdata = multidata
    multi.aparm[1] = 0.16667
    multi()

    # Then selfcal the leakage cal source
    localsn = 1
    localcl = 1
    doamp = False
    dostokesi = False
    soltype = "L1R"
    calibsnr = 5
    averageifs = False
    singlesource_calib(splitdata, modelimage, localsn, refant, doamp, phasecalmins,
                       dostokesi, soltype, calibsnr, averageifs, uvrange, weightit)
    plottops(splitdata, "SN", localsn, "PHAS", 0, 2, 4, outputfilename + ".phase1.ps")
    tacop(splitdata, "SN", localsn, multidata, localsn)
    applysntable(multidata, localsn, '2PT', localcl, refant)
    localcl += 1
    splittoseq(multidata, localcl, "LPJNK", lpcalsource, 1, domulti, False, # Bandpass was already applied
               beginif, endif, combineifs, dopol)
    doamp = True
    localsn += 1
    sourcesn = 1
    singlesource_calib(split2data, modelimage, sourcesn, refant, doamp, ampcalmins,
                       dostokesi, soltype, calibsnr, averageifs, uvrange, weightit)
    plottops(split2data, "SN", sourcesn, "PHAS", 0, 2, 4, outputfilename + ".phase2.ps")
    plottops(split2data, "SN", sourcesn, "AMP", 0, 2, 4, outputfilename + ".amp.ps")
    tacop(split2data, "SN", sourcesn, multidata, localsn)
    applysntable(multidata, localsn, '2PT', localcl, refant)
    localsn += 1
    localcl += 1

    # Then run LPCAL
    lpcal = AIPSTask("lpcal", version = aipsver)
    lpcal.indata = multidata
    lpcal.calsour[1] = lpcalsource
    lpcal.timerang[1:] = get_ampcal_timer(multidata, lpcalsource, leakagecalscan)
    lpcal.docalib = 1
    lpcal.gainuse = localcl
    lpcal.doband = -1 # Bandpass was already applied at the split stage
    lpcal.in2data = modelimage
    lpcal.ncomp[1] = 1 #FIXME: Need to generalise this
    lpcal.nmaps = 1
    lpcal()

    # Split and apply the calibration, and write out the calibrated data, if requested
    calibrateddata = AIPSUVData(lpcalsource, "CALD", 1, 1)
    if calibrateddata.exists():
        calibratedata.zap()
    if not outputfilename == "":
        splitseqno = 1
        splitmulti = False
        splitband = False
        splitbeginif = 1
        splitendif = 0
        combineifs = False
        leakagedopol = 2
        splittoseq(multidata, localcl, 'CALD', lpcalsource, splitseqno, splitmulti, 
                   splitband, splitbeginif,  splitendif, combineifs, leakagedopol)
        writedata(calibrateddata, outputfilename, True)
        calibrateddata.zap()

    # Then write out the AN table and delete the temporary datasets
    writetable(multidata, "AN", 1, lpcaltablefilename)
    splitdata.zap()
    split2data.zap()
    multidata.zap()
    splitdata = AIPSUVData("LPCAL", "CALIB", 1, i)
    if splitdata.exists():
        splitdata.zap()

##### Cross-polar delay calibration ############################################
def xpoldelaycal(uvdataset, clversion, refant, xpolsource, xpolscanno, xpolmodel,
                 solintmins, inttimesecs, sntablefilename, delaywin=400, ratewin=40):
    # Create temp datasets
    tempdata = AIPSUVData("XPOLUNSW", "JUNK", 1, 1)
    if tempdata.exists():
        tempdata.zap()
    tempswapdata = AIPSUVData("XPOLSWAP", "JUNK", 1, 1)
    if tempswapdata.exists():
        tempswapdata.zap()

    # Create the temp dataset
    uvcop = AIPSTask("uvcop", version = aipsver)
    uvcop.indata = uvdataset
    uvcop.outdata = tempdata
    uvcop.sources[1] = xpolsource
    uvcop.timerang[1:] = get_ampcal_timer(uvdataset, xpolsource, xpolscanno)
    uvcop.antennas[1] = refant
    uvcop.flagver = 1
    uvcop()

    # Clear the SN tables in this dataset
    tempdata.zap_table("SN", -1)

    print tempdata.tables
    # Run fring on this dataset
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = tempdata
    fring.aparm[1] = 2
    fring.aparm[2:] = [0]
    fring.aparm[6] = 2
    fring.docalib = 2
    fring.gainuse = clversion
    fring.snver = 1
    fring.in2data = xpolmodel
    fring.solint = solintmins
    fring.dparm[1:] = [0,delaywin,ratewin,inttimesecs,0,0,0,1]
    fring()

    # Apply this fring to temp dataset
    applysntable(tempdata, 1, '2PT', clversion, refant)

    # Run SWPOL to swap the antenna polarisations
    swpol = AIPSTask("SWPOL", version = aipsver)
    swpol.indata = tempdata
    swpol.outdata = tempswapdata
    swpol.antennas[1] = refant
    swpol.docalib = 1
    swpol()

    # Copy NX table to swapdata
    tacop(tempdata, "NX", 1, tempswapdata, 1)

    # Run FRING on the Xpol data
    print tempswapdata.tables
    fring.indata = tempswapdata
    fring.gainuse = clversion + 1
    fring.snver = 2
    fring()

    # Run POLSN on the output
    polsn = AIPSTask("POLSN", version = aipsver)
    polsn.indata = tempswapdata
    polsn.invers = 2
    polsn.outvers = 3
    polsn()

    # Write out the SN table
    writetable(tempswapdata, 'SN', 3, sntablefilename)

    # Clean up
    tempdata.zap()
    tempswapdata.zap()

##### Split IFs into more IFs ##################################################
def morif(inputdata, outputdata, ifmultiplier):
    morif = AIPSTask("MORIF", version = aipsver)
    morif.indata = inputdata
    morif.outdata = outputdata
    morif.npiece = ifmultiplier
    morif()

##### Fringe fit ###############################################################
def fring(uvdataset, snversion, clversion, solintmins, inttimesecs, srcname,
          refant,doband=False,snr=7.5,sumifs=False,modeldata=None,sumrrll=False,
          uvrange=[0,0],zerorates=False, delaywin=400, ratewin=40, doexhaustive=False,
          halfbandwidth=False, dispersivefit=False, dopol=0):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = uvdataset
    fring.docalib = 2
    fring.gainuse = clversion
    fring.calsour[1] = srcname
    fring.flagver = 1
    fring.doband = -1
    fring.dopol = dopol
    fring.smooth[1] = 5 # Hanning smoothing, width 4 channels
    
    if fring.gainuse < 0:
        fring.docalib = -1

    if not modeldata == None:
        fring.in2data = modeldata
    if doband:
        fring.doband = 1
        fring.bpver = 1
    try:
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth[0]/\
                uvdataset.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth/\
                uvdataset.table('FQ', 1)[0].ch_width
    if nchan > 8:
        fring.bchan = 2
        fring.echan = nchan - 1
    else:
        fring.bchan = 1
        fring.echan = nchan
    fring.refant = refant
    fring.uvrange[1] = uvrange[0]
    fring.uvrange[2] = uvrange[1]
    fring.aparm[1:] = [3,0,0,0,0,2,snr,0,1,0]
    #fring.aparm[1:] = [3,0,0,0,0,2,snr,0,0,0]
    if sumrrll:
        fring.aparm[3] = 1
    if sumifs:
        if halfbandwidth:
            fring.aparm[5] = 3
        else:
            fring.aparm[5] = 1
    else:
        if halfbandwidth:
            print "Can't have halfbandwidth without also having sumIFs"
            sys.exit()
    if doexhaustive:
        fring.aparm[9] = 1
    if dispersivefit:
        fring.aparm[10] = 1
    fring.dparm[1:] = [0,delaywin,ratewin,inttimesecs,0,0,0,0]
    if zerorates:
        fring.dparm[8] = 1
    print "Delay window is " + str(delaywin)
    print "Rate window is " + str(ratewin)
    print fring.dparm[3]
    print "Inttime in seconds is " + str(inttimesecs)
    fring.solint = solintmins
    fring.snver = snversion
    print "Just set snver to " + str(fring.snver)
    fring()

##### Fringe fit on an in-beam calibrator #######################################
def inbeam_fring(uvdataset, calmodel, snversion, clversion, solintmins,
                 inttimesecs, srcname, refant, dostokesi):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = uvdataset
    fring.docal = 2
    fring.gainuse = clversion
    fring.calsour[1] = srcname
    fring.flagver = 1
    fring.doband = 1
    fring.bpver = 1
    fring.refant = refant
    if dostokesi:
        fring.aparm[1:] = [3,0,1,0,0,2,5,0,1,0]
    else:
        fring.aparm[1:] = [3,0,0,0,0,2,5,0,1,0]
    fring.dparm[1:] = [0,-1,4,inttimesecs,0,0,0,0]
    fring.solint = solintmins
    fring.snver = snversion
    fring()

##### Bandpass corrections #####################################################
def bpass(uvdataset, srcname, clversion, ampcalscanno, ampcalmodeldata=None, \
          dopol = 0):
    bpass = AIPSTask('bpass', version = aipsver)
    bpass.indata = uvdataset
    bpass.calsour[1] = srcname
    bpass.timerang[1:] = get_ampcal_timer(uvdataset, srcname, ampcalscanno)
    bpass.docalib = 2
    bpass.gainuse = clversion
    bpass.flagver = 1
    bpass.outver = 1
    bpass.solint = 0
    bpass.dopol = dopol
    if ampcalmodeldata == None:
        bpass.in2name = ''
        bpass.in2class = ''
    else:
        bpass.in2data = ampcalmodeldata
    bpass.invers = 0
    bpass.ncomp[1:] = [0]
    bpass.flux = 0
    bpass.smodel[1:] = [0]
    bpass.bpassprm[1:] = [0]
    bpass.bpassprm[5] = 1
    bpass.bpassprm[9] = 1
    bpass.bpassprm[10] = 1
    print "CL version is " + str(clversion) + ", calsour is " + \
          str(bpass.calsour[1]) + ", timer is " + str(bpass.timerang)
    bpass()

##### Polynomial-based bandpass correction #####################################
def cpass(uvdataset, srcname, clversion, ampcalscanno, ampcalmodeldata=None, npoly=10):
    cpass = AIPSTask('cpass', version = aipsver)
    cpass.indata = uvdataset
    cpass.calsour[1] = srcname
    cpass.timerang[1:] = get_ampcal_timer(uvdataset, srcname, ampcalscanno)
    cpass.docalib = 2
    cpass.gainuse = clversion
    cpass.flagver = 1
    cpass.outver = 1
    cpass.solint = 0
    if ampcalmodeldata == None:
        cpass.in2name = ''
        cpass.in2class = ''
    else:
        cpass.in2data = ampcalmodeldata
    cpass.invers = 0
    cpass.ncomp[1:] = [0]
    cpass.flux = 0
    cpass.smodel[1:] = [0]
    cpass.bpassprm[1:] = [0]
    cpass.bpassprm[10] = 1
    cpass.bpassprm[11] = 1
    cpass.cparm[1:] = [0]
    cpass.cparm[1] = npoly
    cpass.cparm[2] = 80
    cpass.cparm[3] = 0.005
    cpass.cparm[5] = 2
    cpass.cparm[8] = 1
    cpass()

##### Get the number of channels in a uv file ##################################
def getNumChannels(uvdata):
    numchannels = -1
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            thisnumchannels = int((row.total_bandwidth[0] / row.ch_width[0]) + 0.5)
        except (AttributeError, TypeError):
            thisnumchannels = int((row.total_bandwidth / row.ch_width) + 0.5)
        if numchannels < 0:
            numchannels = thisnumchannels
        elif thisnumchannels != numchannels:
            print "Different num channels for different IFs! Returning -1"
            return -1
    return numchannels

##### Get the channel bandwidth for a uv file ##################################
def getchannelbandwidthhz(uvdata):
    chanbw = -1
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            thischanbw = row.ch_width[0]
        except (AttributeError, TypeError):
            thischanbw = row.ch_width
        if chanbw < 0:
            chanbw = thischanbw
        elif thischanbw != chanbw:
            print "Different channel widths for different IFs! Returning -1"
            return -1
    return chanbw

##### Get the IF bandwidth for a uv file #######################################
def getNumIFs(uvdata):
    numifs = 0
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            for c in row.ch_width:
                numifs += 1
        except (AttributeError, TypeError):
            numifs += 1
    return numifs

##### Get the IF bandwidth for a uv file #######################################
def getIFbandwidthhz(uvdata):
    totalbw = -1
    fqtable = uvdata.table('FQ', 1)
    for row in fqtable:
        try:
            thistotalbw = row.total_bandwidth[0]
        except (AttributeError, TypeError):
            thistotalbw = row.total_bandwidth
        if totalbw < 0:
            totalbw = thistotalbw
        elif thistotalbw != totalbw:
            print "Different IFs have different widths! Returning -1"
            return -1
    return totalbw

##### Get the integration time for a uv file ###################################
def getintegrationtime(uvdata):
    dtsum = AIPSTask("dtsum", version = aipsver)
    dtsum.indata = uvdata
    dtsum.docrt = 132
    dtsum.outprint = ''
    #dtsum.go()
    #dtsummessage = dtsum.message()
    dtsum.docrt = 0
    dtsum.outprint = 'junk.txt'
    dtsum()
    dtsummessage = open('junk.txt').readlines()
    inttime = -1.0
    for line in dtsummessage:
        if "Data integration time" in line:
            splitline = line.split('=')
            inttime = float(splitline[1].strip())
    if inttime < 0:
        print "Couldn't find the integration time!"
    return inttime

##### Split some data to multiple channels, applying calibration ###############
def splitmulti(uvdataset, clversion, outklass, srcname, seqno, doband=True):
    splittoseq(uvdataset, clversion, outklass, srcname, seqno, True, doband)

##### Split some data to single channel, applying calibration ##################
def split(uvdataset, clversion, outklass, srcname):
    splittoseq(uvdataset, clversion, outklass, srcname, 1, False, True)

##### Split some data to single channel, applying calibration, choose seqno ####
def splittoseq(uvdataset, clversion, outklass, srcname, seqno, domulti=False,
    doband=True, beginif=-1, endif=-1, combineifs=False, dopol=0):
    split = AIPSTask('split', version = aipsver)
    split.indata = uvdataset
    split.sources[1] = srcname
    split.bchan = 1
    try:
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth[0]/\
                uvdataset.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth/\
                uvdataset.table('FQ', 1)[0].ch_width
    split.echan = nchan
    if not beginif < 0:
        split.bif = beginif
    if not endif < 0:
        split.eif = endif
    if domulti:
        split.aparm[1] = 1
        split.nchav = 1
    else:
        split.aparm[1] = 2
        split.echan = nchan-1
        if combineifs:
            split.aparm[1] = 3
        if nchan > 8:
            split.bchan = 2
            split.echan = nchan-2
    split.dopol = dopol
    split.docalib = 2
    split.gainuse = clversion
    if doband:
        split.doband = 1
        split.bpver = 1
    else:
        split.doband = -1
    split.flagver = 1
    split.outclass = outklass
    split.outseq = seqno
    split.outdisk = 1
    split.douvcomp = -1
    split()

##### Split some data to single channel, applying calibration with no bandpass #
def splittoseqnoband(uvdataset, clversion, outklass, srcname, seqno):
    splittoseq(uvdataset, clversion, outklass, srcname, seqno, False, False)

##### Split some data to single channel, applying calibration with no bandpass #
def splittoseqnobandnoweight(uvdataset, clversion, outklass, srcname, seqno):
    split = AIPSTask('split', version = aipsver)
    split.indata = uvdataset
    split.sources[1] = srcname
    split.bchan = 1
    try:
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth[0]/\
                uvdataset.table('FQ', 1)[0].ch_width[0]
    except (AttributeError, TypeError):
        nchan = uvdataset.table('FQ', 1)[0].total_bandwidth/\
                uvdataset.table('FQ', 1)[0].ch_width
    split.bchan = 1
    split.echan = nchan-1
    if nchan > 8:
        split.bchan = 2
        split.echan = nchan-2
    split.aparm[1] = 2
    split.docalib = 100
    split.flagver = 1
    split.gainuse = clversion
    split.doband = -1
    split.outclass = outklass
    split.outseq = seqno
    split.outdisk = 1
    split.douvcomp = -1
    split()

##### Write out uvdata or image ################################################
def writedata(dataset, outputfile, autostomp):
    if os.path.exists(outputfile):
        if autostomp:
            os.remove(outputfile)
        else:
            if yesno("Remove disk file " + outputfile + " to proceed " + \
                     "(no will abort pipeline)? "):
                os.remove(outputfile)
            else:
                sys.exit()

    fittp = AIPSTask('fittp', version = aipsver)
    fittp.indata = dataset
    fittp.dataout = outputfile
    fittp()

def multisourcefit(image, predictedrms, maxsources, outfile):
    sad = AIPSTask('sad', version = aipsver)
    sad.indata = image
    sad.doresid = -1
    sad.ngauss = maxsources
    sad.cparm[1] = 10.0*predictedrms
    sad.cparm[2] = 4.0*predictedrms
    sad.sort = 'S'
    sad.docrt = 132
    sad.fitout = outfile
    sad.outver = -1
    sad.doall = 1
    sad.dowidth[1][1] = 1.0
    sad.dparm[1:] = [0]
    sad.dparm[4] = 200
    sad.dparm[6] = 10
    sad()

##### Image a single-source file, optionally cleaning using autobox ############
def widefieldimage(uvdataset, srcname, numcells, cellmas, doclean, stopflux,
                   taperml, rashiftmas, decshiftmas, numcc, cleanradius):
    imagr = AIPSTask('imagr', version = aipsver)
    imagr.indata = uvdataset
    imagr.sources[1] = srcname
    imagr.stokes = 'I'
    imagr.outname = srcname
    imagr.flux = stopflux
    imagr.outseq = 1
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
            imagr.im2parm[1:6] = [1, 6.5, 9, 0.005, 0]
    else:
        imagr.niter = 0
    imagr.dotv = -1
    imagr()

##### Use imean to get stats from an image and return ##########################
def getimagestats(image, cutofffrac): # [peak, rms, peakx, peaky]
    imean = AIPSTask('imean', version = aipsver)
    imean.indata = image
    xguardpix = int(image.header.naxis[0]*cutofffrac)
    yguardpix = int(image.header.naxis[1]*cutofffrac)
    print "Avoiding a guard band of %dx%d pixels around edge of image..." % \
          (xguardpix, yguardpix)
    imean.blc[1:] = [xguardpix, yguardpix]
    imean.trc[1:] = [image.header.naxis[0]-xguardpix,image.header.naxis[1]-yguardpix]
    imean.dohist = -1
    imean.doinvers = -1
    imean.outtext = ''
    imean.nboxes = 0
    #imean.go()
    try:
        imean()
    except RuntimeError:
        print "IMEAN failed - this image must be bogus!"
        return [0.0, 1.0, image.header.naxis[0]/2, image.header.naxis[1]/2, "", ""]
    imeanlines = imean.message()
    numlines = len(imeanlines)
    atline = 0
    while atline < numlines and imeanlines[atline].split()[-1] != "histogram":
        atline += 1
    if atline >= numlines:
        atline = 0
        noise1 = -1
    else:
        val = imeanlines[atline].split()[4]
        if val == "****":
            val = imeanlines[atline].split()[3]
        noise1 = float(val)
    while atline < numlines and imeanlines[atline].split()[-1] != "pixels":
        atline += 1
    val = imeanlines[atline].split()[4]
    if val == "JY/BEAM":
        val = imeanlines[atline].split()[3]
    noise2 = float(val)
    while atline < numlines and imeanlines[atline].split()[1] != "Skypos:":
        atline += 1
    while atline < numlines and not "Maximum" in imeanlines[atline]:
        atline += 1
    val1 = imeanlines[atline].split()[2]
    val2 = imeanlines[atline].split()[4]
    val3 = imeanlines[atline].split()[5]
    if val1 == "at":
        val1 = (imeanlines[atline].split()[1]).split('=')[-1]
        val2 = imeanlines[atline].split()[3]
        val3 = imeanlines[atline].split()[4]
    peak = float(val1)
    peakx = int(val2)
    peaky = int(val3)
    while atline < numlines and not "Skypos" in imeanlines[atline]:
        atline += 1
    splitline = imeanlines[atline].split()
    rastr = splitline[3] + ":" + splitline[4] + ":" + splitline[5]
    decstr = splitline[7] + ":" + splitline[8] + ":" + splitline[9]
    if noise1 < 0:
        return [peak, noise2, peakx, peaky, rastr, decstr]
    return [peak, noise1, peakx, peaky, rastr, decstr]

##### Image a single-source file and make a ps plot ############################
def image(uvdataset, cellmas, numcells, numcc, flux, srcname, plotfile,
          doallbands, nointeraction, stokesi):
    imagr = AIPSTask('imagr', version = aipsver)
    imagr.indata = uvdataset
    imagr.sources[1] = srcname
    imagr.stokes = 'I'
    if len(srcname) > 12:
        imagr.outname = srcname[:12]
    else:
        imagr.outname = srcname
    imagr.outseq = 1
    imagedata = AIPSImage(imagr.outname, 'ICL001', 1, 1)
    if imagedata.exists():
        imagedata.zap()
    imagr.cellsize[1] = cellmas/1000.0
    imagr.cellsize[2] = cellmas/1000.0
    imagr.imsize[1:] = [numcells, numcells]
    imagr.nfield = 1
    imagr.do3dimag = -1
    imagr.uvwtfn = 'CS'
    imagr.robust = -2
    imagr.uvbxfn = 1
    imagr.xtype = 5
    imagr.ytype = 5
    imagr.niter = numcc
    imagr.nboxes = 1
    imagr.clbox[1:][1:] = [0]
    imagr.gain = 0.1
    imagr.flux = flux/1000.0
    imagr.minpatch = 51
    if nointeraction:
        imagr.dotv = -1
    else:
        imagr.dotv = 1
    imagr.grch = 1
    imagr()
    cntr = AIPSTask('cntr', version = aipsver)
    cntr.indata = imagedata
    cntr.ltype = -3
    cntr.dotv = 0
    cntr()
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = imagedata
    lwpla.plver = 1
    lwpla.invers = 1
    lwpla.outfile = plotfile
    lwpla()
    if doallbands:
        imagr.dotv = 0
        for i in range(4):
            imagr.outseq = i+1
            if stokesi:
                klass = 'icl001'
                imagedata = AIPSImage(srcname, klass, 1, 2*(i+1))
                if imagedata.exists():
                    imagedata.zap()
                imagr.bif=i+1
                imagr.eif=i+1
                imagr.stokes='I'
                imagr()
            else:
                klass = 'RCL001'
                imagedata = AIPSImage(srcname, klass, 1, 2*(i+1))
                if imagedata.exists():
                    imagedata.zap()
                imagr.bif=i+1
                imagr.eif=i+1
                imagr.stokes='RR'
                imagr()
                klass = 'LCL001'
                imagr.stokes='LL'
                imagedata = AIPSImage(srcname, klass, 1, 2*(i+1)+1)
                if imagedata.exists():
                    imagedata.zap()
                imagr()

##### Correct manually provided delay errors using CLCOR   #####################
def clcordelaysfromfile(uvdata, filename, clversion):
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdata
    clcor.gainuse = clversion
    clcor.gainver = clversion
    clcor.opcode = 'CLOC'
    clcor.clcorprm[1:] = [0]
    clcor.clcorprm[7] = 1
    anrow_noss = {}
    antable = uvdata.table('AN', 1)
    for row in antable:
        anrow_noss[row.anname.strip()] = row.nosta
    for line in open(filename):
        if len(line) == 0 or line[0] == "#":
            continue
        splitline = line.split()
        if not len(splitline) == 6:
            print "Badly formed line " + line + " in " + filename
            sys.exit()
        try:
            antindex = anrow_noss[splitline[0]]
        except KeyError:
            print "Couldn't find " + splitline[0] + " in the AN table!"
            sys.exit()
        clcor.antennas[1] = antindex
        print "Correcting " + splitline[0] + ", which is index " + str(antindex)
        clcor.bif = int(splitline[1])
        clcor.eif = int(splitline[1])
        clcor.stokes = splitline[2]
        reftime = splitline[4].split(',')
        if len(reftime) != 4:
            print "Bad reftime " + splitline[3]
        clcor.clcorprm[1] = float(splitline[5])*24.0
        clcor.clcorprm[2] = float(splitline[3])
        clcor.clcorprm[3] = int(reftime[0])
        clcor.clcorprm[4] = int(reftime[1])
        clcor.clcorprm[5] = int(reftime[2])
        clcor.clcorprm[6] = int(reftime[3])
        clcor()

##### Correct delay errors from an amp cal scan using CLCOR#####################
def clcordelays(uvdataset, snversion, clversion, refant, numifs, numstokes):
    clcor = AIPSTask('clcor', version = aipsver)
    stokes = ['RR','LL']
    clcor.indata = uvdataset
    clcor.gainuse = clversion
    clcor.gainver = clversion
    clcor.opcode = 'CLOC'
    clcor.clcorprm[1:] = [0]
    clcor.clcorprm[7] = 1
    wizdataset = WizAIPSUVData(uvdataset)
    sntable = wizdataset.table('SN', snversion)
    for i in range(len(uvdataset.antennas)):
        clcor.antennas[1] = i+1
        found = False
        for row in sntable:
            if row.antenna_no == i+1:
                found = True
                for j in range(numifs):
                    clcor.bif = j+1
                    clcor.eif = j+1
                    clcor.stokes = stokes[0]
                    clcor.clcorprm[2] = row.delay_1[j]*1000000000.0
                    clcor()
                    if numstokes == 2:
                        clcor.stokes = stokes[1]
                        clcor.clcorprm[2] = row.delay_2[j]*1000000000.0
                        clcor()
                break
        if not found:
            print "Couldn't find solution for antenna " + str(i+1)

##### Create a flag table for a dataset for when a model will be too resolved ##
def createmodelflagtable(uvdataset, modelimage, minflux, flagver):
    tempvis = AIPSUVData("JNKVIS","JNKVIS",1,1)
    if tempvis.exists():
        tempvis.zap()
    divvis = AIPSUVData("DIVVIS","DIVVIS",1,1)
    if divvis.exists():
        divvis.zap()
    uvsub = AIPSTask("uvsub", version = aipsver)
    uvsub.indata = uvdataset
    uvsub.outdata = tempvis
    uvsub.smodel[1] = 1.0
    uvsub.opcode = "MODL"
    uvsub()
    uvsub.smodel[1] = 0
    uvsub.in2data = modelimage
    uvsub.indata = tempvis
    uvsub.outdata = divvis
    uvsub.opcode = "DIV"
    uvsub()
    clip = AIPSTask("clip", version = aipsver)
    clip.indata = divvis
    clip.aparm[1] = 1.0/minflux
    clip.aparm[9] = 1
    clip.outfgver = 1
    clip()
    tacop = AIPSTask("tacop", version = aipsver)
    tacop.indata = divvis
    tacop.outdata = uvdataset
    tacop.inext = "FG"
    tacop.inver = 1
    tacop.outver = flagver
    tacop()
    divvis.zap()
    tempvis.zap()

##### Selfcal a single-source file using CALIB #################################
def singlesource_calib(uvdataset, modelimage, snver, refant, doamp, solmins,
                       dostokesi, soltype, calibsnr, averageifs, uvrange=[0,0],
                       weightit=0,flagwheremodelbelow=-1.0,normalise=False):
    calib = AIPSTask('calib', version = aipsver)
    calib.indata = uvdataset
    if modelimage == None:
        calib.smodel[1] = 1.0
        calib.smodel[2:] = [0]
    elif isinstance(modelimage, float):
        print "Calibrating with a delta function of amp " + str(modelimage)
        calib.smodel[1] = modelimage
        calib.smodel[2:] = [0]
    else:
        calib.in2data = modelimage
    if flagwheremodelbelow > 0.0: 
        if modelimage == None or isinstance(modelimage, float):
            print "Can't flag with no image model!"
            sys.exit()
        createmodelflagtable(uvdataset, modelimage, flagwheremodelbelow, 1)
        calib.flagver = 1
    calib.uvrange[1] = uvrange[0]
    calib.uvrange[2] = uvrange[1]
    print "UV range is " + str(calib.uvrange)
    calib.cmethod = 'DFT'
#    calib.cmethod = 'FFT'
#    calib.cmodel = 'IMAG'
    calib.refant = refant
    calib.solint = solmins
    calib.weightit = weightit
    if solmins < 0:
        print "I got a solution interval of " + str(solmins) + " - aborting!"
        sys.exit()
    calib.aparm[1:] = [0]
    calib.aparm[1] = 3
    calib.aparm[7] = calibsnr
    calib.soltype = soltype
    if dostokesi:
        calib.aparm[3] = 1
    if averageifs:
        calib.aparm[5] = 1
    if doamp:
        calib.solmode = 'A&P'
        if normalise:
            calib.normaliz = 1
    else:
        calib.solmode = 'P'
        if normalise:
            print "You can't normalise phase-only solutions"
    calib.snver = snver
    calib()

##### Average all baselines together to a specified averaging time #############
def tbavg(uvdatain, uvdataout, solmins, rashiftmas=0.0, decshiftmas=0.0):
    tbavg = AIPSTask("tbavg", version = aipsver)
    if not uvdatain.exists():
        print "TBAVG: Input uv data doesn't exist!"
        sys.exist()
    if uvdataout.exists():
        print "TBAVG: Output uv data already exists!"
        sys.exit()
    tbavg.indata = uvdatain
    tbavg.outdata = uvdataout
    tbavg.solint = solmins*60
    tbavg.shift[1] = rashiftmas/1000.0
    tbavg.shift[2] = decshiftmas/1000.0
    tbavg()


##### Selfcal using CALIB ######################################################
def calib(uvdataset, modelimage, srcname, snver, clversion, refant, doamp, 
          solmins, dostokesi, soltype, calibsnr=5.0):
    calib = AIPSTask('calib', version = aipsver)
    calib.indata = uvdataset
    calib.calsour[1] = srcname
    calib.docal = 2
    calib.gainuse = clversion
    calib.flagver = 1
    calib.doband = 1
    calib.bpver = 1
    calib.in2data = modelimage
    calib.cmethod = 'DFT'
    calib.refant = refant
    calib.solint = solmins
    calib.aparm[1:] = [0]
    calib.aparm[1] = 3
    calib.aparm[7] = calibsnr
    calib.soltype = soltype
    if dostokesi:
        calib.aparm[3] = 1
    if doamp:
        calib.solmode = 'A&P'
    else:
        calib.aparm[1] = 2
        calib.solmode = 'P'
    calib.snver = snver
    calib()

##### Fix clock offsets and rates ##############################################
def fixclocks(clockfile, uvdata, clversion):
    antdict = {}
    count = 1
    for ant in uvdata.antennas:
        antdict[ant] = count
        count += 1
    clcor = AIPSTask('clcor', version = aipsver)
    clcor.indata = uvdata
    clcor.sources[1:] = ''
    clcor.freqid = clcor.subarray = 0
    clcor.antennas[1:] = [0]
    clcor.gainver = clversion
    clcor.gainuse = clversion
    clcor.opcode = 'CLOC'
    clcor.baddisk[1:] = [0]
    clcor.infile = ''
    clockin = open(clockfile)
    clocklines = clockin.readlines()
    clockin.close()
    for line in clocklines:
        if line == "" or line[0] == '#':
            continue
        splitline = line.split()
        clcor.timerang[1:] = [0]
        clcor.clcorprm[1:] = [0]
        clcor.stokes = '  '
        clcor.antennas[1] = antdict[splitline[0]]
        clcor.bif = int(splitline[3])
        clcor.eif = int(splitline[4])
        clcor.clcorprm[7] = 1
        if splitline[1] == "RATE":
            clcor.clcorprm[1] = float(splitline[2]) #ns/day
        else:
            clcor.clcorprm[2] = float(splitline[2]) #ns
        if splitline[5] == "RR" or splitline[5] == "LL":
            clcor.stokes = splitline[5]
        if len(splitline) > 6:
            if len(splitline) == 14:
                for i in range(8):
                    clcor.timeran[i+1] = int(splitline[i+6])
            else:
                print "Dodgy clock line:"
                print "  " + line
                print "Aborting!"
                sys.exit()
        clcor()

##### Make a plot from the amplitudes and phases of two datasets ###############
def plotvistime(amparray1, phsarray1, amparray2, phsarray2, timearray, numstokes,
                numchan, numif, baselinestr, pol, ifnum, channel, prefix, plottype, 
                label1, label2, notitle, iserr, dogray):
    if len(amparray1) != len(amparray2) or len(amparray1) != len(timearray):
        print "Len(amparray1) is %d, len(amparray2) is %d, len(timearray) is %d" % \
              (len(amparray1), len(amparray2), len(timearray))
        return False
   
    if dogray:
        line1type = "k-+"
        line2type = "k--x"
        pylab.gray()
    else:
        line1type = "g-+"
        line2type = "r-x"
    if not plottype == "ps" and not plottype == "png":
        print "Unknown plottype " + plottype + " - changing to png!"
        plottype = "png"
    toff = []
    for time in timearray:
        toff.append(24.0*(time - int(time)))
    if iserr:
        pylab.xlabel("Time (Hours since MJD %d)" % (timearray[0]))
        pylab.ylabel("Correlation coefficient")
        if not notitle:
            pylab.title(pol + " for IF " + str(ifnum) + ", channel " + str(channel) + \
                        " of baseline " + baselinestr)
        pylab.plot(toff,amparray1,line1type,label=label1) #Actually the real diff array
        pylab.plot(toff,amparray2,line2type,label=label2) #Actually the imag diff array
        pylab.legend(loc='best')
        pylab.savefig(prefix + ".IF" + str(ifnum) + '.chan' + str(channel) + "-" + pol + "." + baselinestr + ".errortime." + plottype)
    else:
        pylab.subplot(211)
        pylab.xlabel("Time (Hours since MJD %d)" % (timearray[0]))
        pylab.ylabel("Amplitude")
        if not notitle:
            pylab.title(pol + " for IF " + str(ifnum) + ", channel " + str(channel) + \
                        " of baseline " + baselinestr)
        pylab.plot(toff,amparray1,line1type,label=label1)
        pylab.plot(toff,amparray2,line2type,label=label2)
        pylab.legend(loc='best')
        pylab.subplot(212)
        pylab.xlabel("Time (Hours since MJD %d)" % (timearray[0]))
        pylab.ylabel("Phase (degrees)")
        pylab.plot(toff,phsarray1,line1type)
        pylab.plot(toff,phsarray2,line2type)
        pylab.savefig(prefix + ".IF" + str(ifnum) + '.chan' + str(channel) + "-" + pol + "." + baselinestr + ".vistime." + plottype)
        pylab.clf()
        pylab.subplot(111)
        pylab.xlabel("Time (Hours since MJD %d)" % (timearray[0]))
        pylab.ylabel("Phase difference (degrees)")
        phsdiff = []
        for i in range(len(phsarray1)):
            phsdiff.append(phsarray1[i] - phsarray2[i])
            if phsdiff[-1] > 180.0:
                phsdiff[-1] -= 180.0
        pylab.plot(toff,phsdiff,line1type)
        maxphsdiff = max(phsdiff)
        minphsdiff = min(phsdiff)
        if -minphsdiff > maxphsdiff:
             maxphsdiff = -minphsdiff
        print "Maximum phase difference for IF " + str(ifnum) + ' chan' + str(channel) + " " + pol + " on baseline " + baselinestr + " was " + str(maxphsdiff)
        pylab.savefig(prefix + ".IF" + str(ifnum) + '.chan' + str(channel) + "-" + pol + "." + baselinestr + ".phsdifference." + plottype)
    pylab.clf()
    return True

##### Make a postscript of a bandpass ##########################################
def plotbandpass(uvdata, bpver, plotbptable, plotsperpage, outputfile, clversion=1, ifs=[0,0], smooth=0, chans=[0,0], stokes=""):
    possm = AIPSTask('possm', version = aipsver)
    possm.indata = uvdata
    possm.stokes = stokes
    if isinstance(ifs, (list,)) and len(ifs) == 2:
        possm.bif = ifs[0]
        possm.eif = ifs[1]
    elif isinstance(ifs, (int,)):
        possm.bif = ifs
        poss.eif = ifs
    else:
        print "IFs parameter must be either a len(2) list for bif and eif, or else a single integer, is", ifs
        sys.exit()
    if isinstance(chans, (list,)) and len(chans) == 2:
        possm.bchan = chans[0]
        possm.echan = chans[1]
    else:
        print "Chans parameter must be a len(2) list for bchan and echan, was", chans
        sys.exit()
    if clversion > 0:
        possm.docalib = 1
        possm.gainuse = clversion
    else:
        possm.docalib = 0
    if plotbptable:
        if bpver <= 0:
            print "bpver must be >= 1 in order to plot BP table"
            sys.exit()
        possm.bpver = bpver
        possm.aparm[8] = 2
    else:
        if bpver > 0:
            possm.doband = 1
            possm.bpver = bpver
        possm.aparm[8] = 0
        possm.smooth[1] = 13 # Hanning
        possm.smooth[2] = smooth
    possm.nplots = plotsperpage
    possm.dotv = 0
    possm.aparm[9] = 1
    possm()
    if os.path.exists(outputfile):
        os.system("rm -f " + outputfile)
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.invers = 9999
    lwpla.outfile = outputfile
    lwpla.lpen = 5
    lwpla.dparm[5] = 0
    lwpla()
    for table in uvdata.tables:
        if table[1] == 'AIPS PL':
            deletetable(uvdata, 'PL', table[0])

##### Make a postscript of a given SN or CL table ##############################
def plottops(uvdata, tabletype, tablever, plotvariable, nifs, npols, plotsperpage, outputfile, doautoorientation=True):
    snplt = AIPSTask('snplt', version = aipsver)
    snplt.indata = uvdata
    snplt.inver  = tablever
    snplt.inext  = tabletype
    snplt.optype = plotvariable
    snplt.bif = 1
    snplt.eif = nifs
    snplt.stokes = ''
    if npols == 1:
        snplt.stokes = 'R'
    snplt.nplots = plotsperpage
    snplt.ltype = -3
    snplt.dotv = -1
    snplt.grch = 1
    snplt()
    if os.path.exists(outputfile):
        os.system("rm -f " + outputfile)
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = uvdata
    lwpla.plver = 1
    lwpla.invers = 9999
    lwpla.outfile = outputfile
    lwpla.lpen = 5
    lwpla.dparm[5] = 1
    if not doautoorientation:
        lwpla.dparm[5] = 0
    lwpla()
    for table in uvdata.tables:
        if table[1] == 'AIPS PL':
            deletetable(uvdata, 'PL', table[0])

##### Check if a given table exists in a UV dataset ############################
def tableexists(uvdata, tabletype, tablever):
    for table in uvdata.tables:
        if table[1] == 'AIPS ' + tabletype:
            if table[0] == tablever or tablever < 0:
                return True
    return False

##### Make a plot from two POSSM dump files ####################################
def plotpossmresults(lines1, lines2, numstokes, numchan, ifnum, baselinestr, prefix, plottype, label1, label2, notitle, dogray):
    if len(lines1) != len(lines2):
        return False
    if  dogray:
        pylab.gray()
        line1type = "k-+"
        line2type = "k--x"
    else:
        line1type = "g-+"
        line2type = "r-x"
    if not plottype == "ps" and not plottype == "png":
        print "Unknown plottype " + plottype + " - changing to png!"
        plottype = "png"
    at = 0
    for i in range(numstokes):
        while at < len(lines1) and not (len(lines1[at].split()) > 0 and \
                                        len(lines2[at].split()) > 0 and \
                                        lines1[at].split()[0] == "Channel" and \
                                        lines2[at].split()[0] == "Channel"):
            at += 1
        if at == len(lines1):
            return False
        at += 1
        freq = []
        amp1 = []
        amp2 = []
        erroramp = []
        errorphs = []
        phase1 = []
        phase2 = []
        for j in range(numchan):
            split1line = lines1[at].split()
            split2line = lines2[at].split()
            pol = split1line[2]
            freq.append(float(split1line[3]))
            if split1line[2] != split2line[2]:
                print "Warning - STOKES does not agree: " + pol + ", " + split2line[2]
            if split1line[3] != split2line[3]:
                print "Warning - frequency does not agree: " + split1line[3]+ ", " + split2line[3]
            re1 = float(split1line[5])*math.cos(float(split1line[6])*math.pi/180.0)
            im1 = float(split1line[5])*math.sin(float(split1line[6])*math.pi/180.0)
            re2 = float(split2line[5])*math.cos(float(split2line[6])*math.pi/180.0)
            im2 = float(split2line[5])*math.sin(float(split2line[6])*math.pi/180.0)
            amp1.append(float(split1line[5]))
            amp2.append(float(split2line[5]))
            phase1.append(float(split1line[6]))
            phase2.append(float(split2line[6]))
            if phase1[-1] > 180.0:
                phase1[-1] -= 360.0
            if phase2[-1] > 180.0:
                phase2[-1] -= 360.0
            erroramp.append(math.sqrt((re1-re2)*(re1-re2) + (im1-im2)*(im1-im2)))
            errorphs.append(math.atan2(im1-im2,re1-re2)*180.0/math.pi)
            at = at+1
        pylab.subplot(211)
        pylab.xlabel("Frequency")
        pylab.ylabel("Amplitude")
        if not notitle:
             pylab.title(pol + " for IF " + str(ifnum) + " of baseline " + baselinestr)
        handle1 = pylab.plot(freq,amp1,line1type,label=label1)
        handle2 = pylab.plot(freq,amp2,line2type,label=label2)
        pylab.legend(loc='best')
        pylab.subplot(212)
        pylab.xlabel("Frequency")
        pylab.ylabel("Phase (degrees)")
        pylab.plot(freq,phase1,line1type)
        pylab.plot(freq,phase2,line2type)
        pylab.savefig(prefix + "." + str(ifnum) + "-" + pol + "." + baselinestr + ".visfreq." + plottype)
        pylab.clf()
        pylab.subplot(211)
        pylab.xlabel("Frequency")
        pylab.ylabel("Amplitude")
        if not notitle:
             pylab.title(pol + " errors for IF " + str(ifnum) + " of baseline " + baselinestr)
        handle1 = pylab.plot(freq,erroramp,line1type,label=label1 + "-" + label2)
        pylab.legend(loc='best')
        pylab.subplot(212)
        pylab.xlabel("Frequency")
        pylab.ylabel("Phase (degrees)")
        pylab.plot(freq,errorphs,line1type)
        pylab.savefig(prefix + "." + str(ifnum) + "-" + pol + "." + baselinestr + ".errorfreq." + plottype)
        pylab.clf()
    return True

##### Use difmap to do a selfcal ###############################################
def difmapselfcal(inputfile, tabledir, modeldir, experiment, source):
    mastermodelfile = modeldir + '/' + source + '.mod'
    expmodelfile    = modeldir + '/' + source + '.' + experiment + '.mod'
    rrsnfile        = tabledir + '/' + source + '.selfcal.rr.sn'
    llsnfile        = tabledir + '/' + source + '.selfcal.ll.sn'
    isnfile         = tabledir + '/' + source + '.selfcal.sn'

    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE)
    difmap.stdin.write("obs " + inputfile + "\n")
    difmap.stdin.write("uvaver 10\n")
    difmap.stdin.write("mapsize 1024,1\n")
    difmap.stdin.write("uvweight 0\n")
    difmap.stdin.write("device /xs\n")
    difmap.stdin.write("select i\n")
    difmap.stdin.write("rmod " + mastermodelfile + "\n")
    difmap.stdin.write("vplot\n")
    difmap.stdin.write("radplot\n")
    loopresponse = ''
    for i in range(20):
        loopresponse = loopresponse + 'modelfit 20\nselfcal true,true,2.5\n'
    print "First make a model you are happy with - enter to dump that model." + \
          "  I suggest selfcal'ing at 2.5 min solution intervals"
    response = raw_input("Enter a difmap command - enter to go on to wmod")
    if response == 'loop': response = loopresponse
    while response != "":
        if response[:6] == 'uvaver':
            print "Averaging is not allowed!!!"
        else:
            difmap.stdin.write(response + "\n")
        response = raw_input("Enter a difmap command - enter to go on to wmod")
        if response == 'loop': response = loopresponse
    difmap.stdin.write("wmod " + expmodelfile + "\n")
    difmap.stdin.write("select rr\n")
    loopresponse = ''
    for i in range(20):
        loopresponse = loopresponse + 'selfcal true,true,2.5\n'
    print "Selfcal RR (I recommend solution intervals of 2.5 minutes) until " + \
          "chi-squared is ok. NO MODEL-FITTING!!!"
    response = raw_input("Enter a difmap command - enter to go on to cordump: ")
    if response == 'loop': response = loopresponse
    while response != "":
        if response[:6] == 'uvaver':
            print "Averaging is not allowed!!!"
        else:
            difmap.stdin.write(response + "\n")
        response = raw_input("Enter a difmap command - enter to dump solutions")
        if response == 'loop': response = loopresponse
    difmap.stdin.write("cordump " + rrsnfile + "\n")
    difmap.stdin.write("select ll\n")
    print "Now selfcal for LL (I recommend solution intervals of 2.5 minutes)" + \
          " until chi-squared is ok. NO MODEL-FITTING!!!"
    response = raw_input("Enter a difmap command - enter to go on to cordump: ")
    if response == 'loop': response = loopresponse
    while response != "":
        if response[:6] == 'uvaver':
            print "Averaging is not allowed!!!"
        else:
            difmap.stdin.write(response + "\n")
        response = raw_input("Enter a difmap command - enter to dump solutions")
        if response == 'loop': response = loopresponse
    difmap.stdin.write("cordump " + llsnfile + "\n")
    difmap.stdin.write("exit\n\n")
    difmap.wait()
    os.system("dasncon " + rrsnfile + " " + llsnfile + " " + isnfile + "\n")

def getimagerms(source, uvfitsfile, resultsfile):
    tempout = open("junkresults.txt", "w")
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=tempout)
    difmap.stdin.write('obs ' + uvfitsfile + '\n')
    difmap.stdin.write('select i\n')
    difmap.stdin.write('uvweight 0,-1\n')
    difmap.stdin.write('mapsize 2048,1\n')
    difmap.stdin.write('print "imagermsguess", imstat(rms)\n')
    difmap.stdin.write('exit\n\n')
    difmap.wait()
    tempin = open("junkresults.txt")
    for line in tempin:
        if "imagermsguess" in line:
            resultsout = open(resultsfile, "a")
            resultsout.write("%s %.2f\n" % (source, 1000.0*float(line.split()[1])))
            resultsout.close()
            break

##### Write a script to map one band of a pulsar ###############################
def write_difmappsrscript(imagename, bands, difmap, pixsize, finepix,npixels=1024):
    imagename = imagename + "." + bands
    difmap.stdin.write("clrmod true\n")
    difmap.stdin.write("unshift\n")
    difmap.stdin.write("shift -peakx,-peaky\n")
    difmap.stdin.write("mapsize " + str(npixels) + "," + str(finepix) + "\n")
    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("addcmp pkflux, true, finepeakx, finepeaky, true, 0, " + \
                       "false, 1, false, 0, false, 0, 0, 0\n")
    difmap.stdin.write("mapsize " + str(npixels) + "," + str(pixsize) + "\n")
    difmap.stdin.write("modelfit 50\n")
    difmap.stdin.write("rmsflux = imstat(rms)\n")
    difmap.stdin.write("restore\n")
    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("ilevs = pkflux/rmsflux\n")
    difmap.stdin.write("lowlev = 300.0/ilevs\n")
    difmap.stdin.write("loglevs lowlev\n")
    #difmap.stdin.write("device " + imagename + ".clean.ps/PS\n")
    psfile = imagename.split('/')[-1] +  ".clean.ps"
    os.system("rm -f " + psfile)
    difmap.stdin.write("device %s/PS\n" % psfile)
    difmap.stdin.write("mappl cln\n")
    if os.path.exists(imagename):
        os.system("rm -f " + imagename)
    difmap.stdin.write("wmap " + imagename + "\n")
    difmap.stdin.write("wmod " + imagename + "mod\n")
    difmap.stdin.write("unshift\n")

##### Use difmap to map a target ###############################################
def difmap_maptarget(uvfile, imagefile, nointeraction, stokesi, pixsize=1.0, 
                     mapsize=1024, uvweightstr="0,-1", uvaverstr='20,True', dogaussian=False, 
                     beginif=1, endif=4, ifrange="", finalmapsize=1024, finepix=0.2):
    """
    Note that for VLBI search (with no previous VLBI detection), uvaverstr should be set to <=30 to avoid smearing effect.
    """
    inputmsg = "Enter a difmap command for the LL data - enter to go to fitting"
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE)
    if pixsize/2.0 < finepix:
        finepix = pixsize/2.0
    difmap.stdin.write("float pkflux\n")
    difmap.stdin.write("float peakx\n")
    difmap.stdin.write("float peaky\n")
    difmap.stdin.write("float finepeakx\n")
    difmap.stdin.write("float finepeaky\n")
    difmap.stdin.write("float rmsflux\n")
    difmap.stdin.write("integer ilevs\n")
    difmap.stdin.write("float lowlev\n")
    difmap.stdin.write("obs " + uvfile + "\n")
    difmap.stdin.write("mapsize " + str(mapsize) + "," + str(pixsize) + "\n")
    difmap.stdin.write("uvweight " + uvweightstr + "\n")
    difmap.stdin.write("uvaver " + uvaverstr + "\n")
    difmap.stdin.write("mapcolor none\n")
    if nointeraction:
        difmap.stdin.write("device /null\n")
    else:
        difmap.stdin.write("device /xs\n")
    if not nointeraction:
        if stokesi:
            difmap.stdin.write("select i\n")
            response = raw_input(inputmsg)
            while response != "":
                difmap.stdin.write(response + "\n")
                response = raw_input(inputmsg)
        else:
            difmap.stdin.write("select ll\n")
            response = raw_input(inputmsg)
            while response != "":
                difmap.stdin.write(response + "\n")
                response = raw_input(inputmsg)
            difmap.stdin.write("select rr\n")
            inputmsg = "Enter a difmap command for the RR data - " + \
                       "enter to go to fitting"
            response = raw_input(inputmsg)
            while response != "":
                difmap.stdin.write(response + "\n")
                response = raw_input(inputmsg)
    difmap.stdin.write("select i\n")      
    difmap.stdin.write("save " + uvfile + ".junk\n")
    difmap.stdin.write("select i\n")
    difmap.stdin.write("peakx = peak(x,abs)\n")
    difmap.stdin.write("peaky = peak(y,abs)\n")
    difmap.stdin.write("shift -peakx,-peaky\n")
    difmap.stdin.write("mapsize 1024," + str(finepix) + "\n")
    difmap.stdin.write("finepeakx = peak(x,abs)\n")
    difmap.stdin.write("finepeaky = peak(y,abs)\n")
    if stokesi:
        pols = ['i']
    else:
        pols = ['rr','ll']
    for p in pols:
        for i in range(beginif, endif+1):
            print "Doing IF " + str(i) + ", pol " + p
            difmap.stdin.write("select " + p + "," + str(i) + "\n")
            if dogaussian:
                write_difmapmapscript(imagefile, p + "." + str(i), difmap,pixsize,finepix,finalmapsize)
            else:
                write_difmappsrscript(imagefile, p + "." + str(i), difmap,pixsize,finepix,finalmapsize)
        if not p == "i":
            print "Doing IF " + str(beginif) + "-" + str(endif) + ", pol " + p
            difmap.stdin.write("select " + p + "," + str(beginif) + "," + str(endif) + '\n')
            write_difmapmapscript(imagefile, p + ".a", difmap,pixsize,finepix)
    if endif > 0:
        selectstring = "select i," + str(beginif) + "," + str(endif)
    elif ifrange != "":
        selectstring = "select i," + ifrange
    else:
        selectstring = "select"
    difmap.stdin.write(selectstring + "\n")
    if dogaussian:
        write_difmapmapscript(imagefile, "ii.a", difmap,pixsize,finepix,finalmapsize)
        write_difmapradplotscript(imagefile, "ii.a", difmap)
    else:
        write_difmappsrscript(imagefile, "ii.a", difmap,pixsize,finepix,finalmapsize)
    #difmap.stdin.write("device %s.ps/PS\n" % imagefile)
    #difmap.stdin.write("mappl cln\n")
    difmap.stdin.write("exit\n\n")
    difmap.wait()
    imagefilesplit = imagefile.split('/')
    if len(imagefilesplit) > 1:
        localimagefile = imagefile.split('/')[-1]
        outputdir = imagefile[:-len(localimagefile)]
        os.system("mv -f %s* %s" % (localimagefile, outputdir))

##### Write a script to map one band of a target ###############################
def write_difmapmapscript(imagename, bands, difmap, pixsize, finepix,npixels=1024):
    imagename = imagename + "." + bands
    difmap.stdin.write("clrmod true\n")
    difmap.stdin.write("unshift\n")
    difmap.stdin.write("shift -peakx,-peaky\n")
    difmap.stdin.write("mapsize " + str(npixels) + "," + str(finepix) + "\n")
    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("addcmp pkflux, true, finepeakx, finepeaky, true, 0.5, " + \
                       "true, 1, true, 0, true\n")
    difmap.stdin.write("mapsize " + str(npixels) + "," + str(pixsize) + "\n")
    difmap.stdin.write("modelfit 50\n")
    difmap.stdin.write("rmsflux = imstat(rms)\n")
    difmap.stdin.write("restore\n")
    difmap.stdin.write("pkflux = peak(flux,abs)\n")
    difmap.stdin.write("ilevs = pkflux/rmsflux\n")
    difmap.stdin.write("lowlev = 300.0/ilevs\n")
    difmap.stdin.write("loglevs lowlev\n")
    #difmap.stdin.write("device " + imagename + ".clean.ps/PS\n")
    psfile = imagename.split('/')[-1] +  ".clean.ps"
    os.system("rm -f " + psfile)
    difmap.stdin.write("device %s/PS\n" % psfile)
    difmap.stdin.write("mappl cln\n")
    if os.path.exists(imagename):
        os.system("rm -f " + imagename)
    difmap.stdin.write("wmap " + imagename + "\n")
    difmap.stdin.write("unshift\n")

##### Write a script to make a radplot showing how good the target was #########
def write_difmapradplotscript(imagename, bands, difmap):
    difmap.stdin.write("clrmod true\n")
    difmap.stdin.write("unshift\n")
    difmap.stdin.write("shift -peakx,-peaky\n")
    # Add a clean window in roughly the right place


##### Use difmap to map a source ###############################################
def difmap_mapsource(uvfile, imagefile, pixsize, mapsize, coarsepixsize):
    inputmsg = "Enter a difmap command for the LL data - enter to go to fitting"
    finepix = 0.2
    if pixsize/2.0 < finepix:
        finepix = pixsize/2.0
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE)
    difmap.stdin.write("float pkflux\n")
    difmap.stdin.write("float peakx\n")
    difmap.stdin.write("float peaky\n")
    difmap.stdin.write("float finepeakx\n")
    difmap.stdin.write("float finepeaky\n")
    difmap.stdin.write("float rmsflux\n")
    difmap.stdin.write("integer ilevs\n")
    difmap.stdin.write("float lowlev\n")
    difmap.stdin.write("obs " + uvfile + "\n")
    difmap.stdin.write("select i\n")
    difmap.stdin.write("mapsize " + str(mapsize) + "," + str(coarsepixsize) + "\n")
    difmap.stdin.write("uvweight 0,-1\n")
    difmap.stdin.write("mapcolor none\n")
    difmap.stdin.write("device /none\n")
    difmap.stdin.write("peakx = peak(x,abs)\n")
    difmap.stdin.write("peaky = peak(y,abs)\n")
    difmap.stdin.write("shift -peakx,-peaky\n")
    difmap.stdin.write("mapsize 1024," + str(finepix) + "\n")
    difmap.stdin.write("finepeakx = peak(x,abs)\n")
    difmap.stdin.write("finepeaky = peak(y,abs)\n")
    write_difmapmapscript(imagefile, difmap, pixsize, finepix)
    difmap.stdin.write("exit\n\n")
    difmap.wait()

##### Use JMFIT to get position/error estimate for a loaded AIPS image #########
def imagrjmfit(imagedata, jmfitfile, xpos, ypos):
    targetstatout = open(jmfitfile + ".stats", "w")
    outdata = AIPSImage('JUNK', 'IMG', 1, 2)
    if outdata.exists():
        outdata.zap()
    obsdate = imagedata.header.date_obs
    year = int(obsdate[0:4])
    month = int(obsdate[5:7])
    day = int(obsdate[8:10])
    mjd = year*367 - int(7*(year + int((month + 9)/12))/4) + \
          int(275*month/9) + day - 678987
    jmfit = AIPSTask('jmfit', version = aipsver)
    jmfit.indata = imagedata
    jmfit.outdata = outdata
    jmfit.niter = 4000
    jmfit.ngauss = 1
    jmfit.fwidth[1][1] = 0
    jmfit.fwidth[1][2] = 0
    jmfit.fwidth[1][3] = 0
    jmfit.domax[1:] = [1]
    jmfit.dopos[1][1] = 1
    jmfit.dopos[1][2] = 1
    jmfit.dowidth[1][1] = 1
    jmfit.dowidth[1][2] = 1
    jmfit.dowidth[1][3] = 1
    jmfit.blc[1:] = [xpos-48, ypos-48]
    jmfit.trc[1:] = [xpos+48, ypos+48]
    jmfit.go()
    jmfitmessage = jmfit.message()

##### Make a contour plot (ps and jpg) of a 256x256 region of a cleaned image ##
def make_contour_plot(cleanimage, psfile, xcentrepix, ycentrepix, jpgfile=None,
                      snr=0.0):
    if os.path.exists(psfile):
        os.system("rm -f " + psfile)
    kntr = AIPSTask('kntr', version = aipsver)
    kntr.docont = 1
    kntr.dogrey = 0
    kntr.dovect = 0
    kntr.indata = cleanimage
    kntr.in2data = cleanimage
    ralo  = xcentrepix-128
    rahi  = xcentrepix+127
    declo = ycentrepix-128
    dechi = ycentrepix+127
    if ralo < 1:
        ralo = 1
        rahi = 256
    if rahi > cleanimage.header.naxis[0]:
        rahi = cleanimage.header.naxis[0]
        ralo = cleanimage.header.naxis[0] - 255
    if declo < 1:
        declo = 1
        dechi = 256
    if dechi > cleanimage.header.naxis[1]:
        dechi = cleanimage.header.naxis[1]
        declo = cleanimage.header.naxis[1] - 255
    kntr.blc[1:] = [ralo, declo]
    kntr.trc[1:] = [rahi, dechi]
    kntr.functype = 'LN'
    if snr > 0:
        if snr < 10:
            kntr.plev = 20
            kntr.levs[1:9] = [-4,-3,-2,-1,1,2,3,4]
        elif snr < 20:
            kntr.plev = 10
            kntr.levs[1:19] = [-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9]
        elif snr < 40:
            kntr.plev = 5
            kntr.levs[1:11] = [-8,-4,-2,-1,1,2,4,8,16,32]
        elif snr < 80:
            kntr.plev = 2
            kntr.levs[1:12] = [-8,-4,-2,-1,1,2,4,8,16,32,64]
        else:
            kntr.plev = 1
            kntr.levs[1:13] = [-8,-4,-2,-1,1,2,4,8,16,32,64,128]
    kntr.ltype = -3
    kntr.dotv = 0
    kntr()
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = cleanimage
    lwpla.plver = 1
    lwpla.invers = 1
    lwpla.outfile = psfile
    lwpla.dparm[5] = 1
    lwpla()
    cleanimage.table('PL', 1).zap()

    if not jpgfile == None:
        if os.path.exists(jpgfile):
            os.system("rm -f " + jpgfile)
        # Convert the postscript to a jpeg
        os.system("gs -sDEVICE=jpeg -dBATCH -dNOPAUSE -sOutputFile=%s %s" % (jpgfile, psfile))

##### Make a plot of amplitude vs uv distance for a dataset ####################
def aipsradplot(uvdata, outfilename):
    maxplver = 0
    for table in uvdata.tables:
        if table[1] == 'AIPS PL':
            if table[0] > maxplver:
                maxplver = table[0]
    uvplt = AIPSTask('uvplt', version = aipsver)
    uvplt.indata = uvdata
    uvplt.ltype = -3
    uvplt.dotv = -1
    uvplt()
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = uvdata
    lwpla.plver = maxplver+1
    lwpla.invers = maxplver+1
    lwpla.outfile = outfilename
    lwpla.dparm[5] = 1
    lwpla()
    uvdata.table('PL', maxplver+1).zap()

##### Make an plot (ps and jpg) of a 256x256 region of a cleaned image #########
def make_greyscale_plot(cleanimage, psfile, xcentrepix, ycentrepix, jpgfile=None):
    if os.path.exists(psfile):
        os.system("rm -f " + psfile)
    kntr = AIPSTask('kntr', version = aipsver)
    kntr.dogrey = 1
    kntr.dovect = 0
    kntr.indata = cleanimage
    kntr.in2data = cleanimage
    ralo  = xcentrepix-128
    rahi  = xcentrepix+127
    declo = ycentrepix-128
    dechi = ycentrepix+127
    if ralo < 1:
        ralo = 1
        rahi = 256
    if rahi > cleanimage.header.naxis[0]:
        rahi = cleanimage.header.naxis[0]
        ralo = cleanimage.header.naxis[0] - 255
    if declo < 1:
        declo = 1
        dechi = 256
    if dechi > cleanimage.header.naxis[1]:
        dechi = cleanimage.header.naxis[1]
        declo = cleanimage.header.naxis[1] - 255

    kntr.blc[1:] = [ralo, declo]
    kntr.trc[1:] = [rahi, dechi]
    kntr.functype = 'NG'
    kntr.docont = 0
    kntr.ltype = -3
    kntr.dotv = 0
    kntr()
    lwpla = AIPSTask('lwpla', version = aipsver)
    lwpla.indata = cleanimage
    lwpla.plver = 1
    lwpla.invers = 1
    lwpla.outfile = psfile
    lwpla()
    cleanimage.table('PL', 1).zap()

    if not jpgfile == None:
        if os.path.exists(jpgfile):
            os.system("rm -f " + jpgfile)
        # Convert the postscript to a jpeg
        os.system("gs -sDEVICE=jpeg -dBATCH -dNOPAUSE -sOutputFile=%s %s" % (jpgfile, psfile))

##### Use JMFIT to get a position and error estimate ###########################
def nonpulsarjmfit(imagefile, jmfitfile, target, centrerapixel=-1, 
                   centredecpixel=-1, fitwidth=True,doextended=False, 
                   loadedfile=None, pixwindow=96):
    imagedata = AIPSImage('JUNK', 'IMG', 1, 1)
    outdata = AIPSImage('JUNK', 'IMG', 1, 2)
    targetstatout = open(jmfitfile + ".stats", "w")
    if outdata.exists():
        outdata.zap()
    if loadedfile == None:
        if imagedata.exists():
            imagedata.zap()
        fitld = AIPSTask('fitld', version = aipsver)
        fitld.datain = imagefile
        fitld.outdata = imagedata
        fitld.dotable = 1
        fitld.douvcomp = -1
        fitld()
    else:
        imagedata = loadedfile

    ccflux = 0.0
    try:
        cctable = imagedata.table('CC', 1)
        for row in cctable:
            ccflux += row.flux
        print "CC table flux is " + str(ccflux)
    except ValueError:
        print "No CC table: leaving CC flux as 0 mJy"

    obsdate = imagedata.header.date_obs
    year = int(obsdate[0:4])
    month = int(obsdate[5:7])
    day = int(obsdate[8:10])
    mjd = year*367 - int(7*(year + int((month + 9)/12))/4) + \
          int(275*month/9) + day - 678987

    jmfit = AIPSTask('jmfit', version = aipsver)
    jmfit.indata = imagedata
    jmfit.outdata = outdata
    jmfit.niter = 4000
    jmfit.ngauss = 1
    jmfit.fwidth[1][1] = 0
    jmfit.fwidth[1][2] = 0
    jmfit.fwidth[1][3] = 0
    jmfit.dopos[1][1] = 1
    jmfit.dopos[1][2] = 1
    jmfit.domax[1:] = [0]
    jmfit.dowidth[1][1] = 0
    jmfit.dowidth[1][2] = 0
    jmfit.dowidth[1][3] = 0
    if fitwidth:
        jmfit.domax[1:] = [1]
        jmfit.dowidth[1][1] = 1
        jmfit.dowidth[1][2] = 1
        jmfit.dowidth[1][3] = 1
    if centrerapixel < 0:
        centrerapixel = int(imagedata.header.naxis[0])/2
    if imagedata.header.naxis[0] < 128:
        print "Image too small in RA axis (%d) - increase size!" % (imagedata.header.naxis[0])
        sys.exit()
    if centredecpixel < 0:
        centredecpixel = int(imagedata.header.naxis[1])/2
    if imagedata.header.naxis[1] < 128:
        print "Image too small in Dec axis (%d) - increase size!" % (imagedata.header.naxis[1])
        sys.exit()
    ralo  = centrerapixel-pixwindow/2
    rahi  = centrerapixel+pixwindow/2
    declo = centredecpixel-pixwindow/2
    dechi = centredecpixel+pixwindow/2
    if ralo < 1:
        ralo = 1
        rahi = pixwindow+1
    if rahi > imagedata.header.naxis[0]:
        rahi = imagedata.header.naxis[0]
        ralo = imagedata.header.naxis[0] - pixwindow
    if declo < 1:
        declo = 1
        dechi = pixwindow+1
    if dechi > imagedata.header.naxis[1]:
        dechi = imagedata.header.naxis[1]
        declo = imagedata.header.naxis[1] - pixwindow
    jmfit.blc[1:] = [ralo, declo]
    jmfit.trc[1:] = [rahi, dechi]
    jmfit.go()
    jmfitmessage = jmfit.message()
    msgindex = 0
    exciselinenos = []
    for i in range(len(jmfitmessage)-1):
        if jmfitmessage[len(jmfitmessage)-(i+1)][:5] != 'JMFIT':
            exciselinenos.append(len(jmfitmessage)-(i+1))
        if "channelbeam" in jmfitmessage[len(jmfitmessage)-(i+1)]:
            exciselinenos.append(len(jmfitmessage)-(i+1))
    for e in exciselinenos:
        jmfitmessage = jmfitmessage[:e] + jmfitmessage[e+1:]
    while jmfitmessage[msgindex].find('X-ref') == -1 and \
              msgindex < len(jmfitmessage):
        msgindex = msgindex + 1
    if msgindex >= len(jmfitmessage):
        print "Did not find JMFIT solution!"
        targetstatsout.close()
        os.remove(targetstatsout)
        return
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
    while not "Beam=" in jmfitmessage[msgindex] and msgindex < len(jmfitmessage):
        msgindex = msgindex + 1
    if msgindex >= len(jmfitmessage):
        print "Did not find intrinsic beam for " + imagefile
        targetstatsout.close()
        os.remove(targetstatsout)
        return
    beamsplit = jmfitmessage[msgindex].split()
    if not len(beamsplit) == 10:
        print "*****WARNING*****: Poorly formatted beam string!"
        if "=" in beamsplit[-8]:
            beammaj = float(beamsplit[-8].split('=')[-1])
        else:
            beammaj = float(beamsplit[-8])*pixsize
    else:
        beammaj = float(beamsplit[-8])*pixsize
    beammin = float(beamsplit[-6])*pixsize
    beampa  = float(beamsplit[-2])*pixsize
    msgindex += 1
    while jmfitmessage[msgindex].find('********* Solution from JMFIT') == -1 \
              and msgindex < len(jmfitmessage):
        msgindex = msgindex + 1
    if msgindex >= len(jmfitmessage):
        print "Did not find JMFIT solution!"
        targetstatsout.close()
        os.remove(targetstatsout)
        return
    if '*' in jmfitmessage[msgindex+5].split()[3]:
        xpospixels = -99999.9
    else:
        xpospixels = float(jmfitmessage[msgindex+5].split()[3])
    if '*' in jmfitmessage[msgindex+6].split()[3]:
        ypospixels = -99999.9
    else:
        ypospixels     = float(jmfitmessage[msgindex+6].split()[3])
    extraline = 0
    if "RASHIFT" in jmfitmessage[msgindex+15]:
        extraline = 1
    sourcerasplit  = jmfitmessage[msgindex+7].split()
    sourcedecsplit = jmfitmessage[msgindex+8].split()
    fluxsplit      = jmfitmessage[msgindex+3].split()
    intfluxsplit   = jmfitmessage[msgindex+4].split()
    fitmajsplit    = jmfitmessage[msgindex+12].split()
    fitminsplit    = jmfitmessage[msgindex+13].split()
    fitpasplit     = jmfitmessage[msgindex+14].split()
    dfitmajsplit   = jmfitmessage[msgindex+extraline+23].split()
    dfitminsplit   = jmfitmessage[msgindex+extraline+24].split()
    dfitpasplit    = jmfitmessage[msgindex+extraline+25].split()
    try:
        fitmaj = float(fitmajsplit[4])
        fitmin = float(fitminsplit[4])
    except ValueError:
        print "Bad value for one or more of major axis " +  fitmajsplit[4] + \
              ", minor axis " + fitminsplit[4] + ", both will be set to zero"
        fitmaj = 0.0
        fitmin = 0.0
    try:
        dfitmaj = float(dfitmajsplit[3])
        dfitmin = float(dfitminsplit[3])
        dfitpa  = float(dfitpasplit[3])
        dfitmajminmax = [float(dfitmajsplit[4]), float(dfitmajsplit[5])]
        dfitminminmax = [float(dfitminsplit[4]), float(dfitminsplit[5])]
        dfitpaminmax  = [float(dfitpasplit[4]), float(dfitpasplit[5])]
    except ValueError:
        print "Bad value for one or more of deconvolved axis/pa"
        print dfitmajsplit
        print dfitminsplit
        print dfitpasplit
        print "All will be set to 0.0"
        dfitmaj = 0.0
        dfitmin = 0.0
        dfitpa  = 0.0
        dfitmajminmax = [0.0,0.0]
        dfitminminmax = [0.0,0.0]
        dfitpaminmax  = [0.0,0.0]
    targetstatout.write("Source " + target + ": MJD " + str(mjd) + ":\n")
    if fluxsplit[4] == '+/-':
        flux = float(fluxsplit[3][1:])*1000.0
        rms  = float(fluxsplit[5])*1000.0
    else:
        flux = float(fluxsplit[4])*1000.0
        rms  = float(fluxsplit[6])*1000.0
    if intfluxsplit[3] == '+/-':
        intflux = float(intfluxsplit[2].split('=')[1])*1000.0
        intfluxerr = float(intfluxsplit[4])*1000.0
    else:
        intflux = float(intfluxsplit[3])*1000.0
        intfluxerr = float(intfluxsplit[5])*1000.0
    if xpospixels > int(imagedata.header.naxis[0])/2 + 48 \
       or xpospixels < int(imagedata.header.naxis[0])/2 - 48 \
       or ypospixels > int(imagedata.header.naxis[1])/2 + 48 \
       or ypospixels < int(imagedata.header.naxis[1])/2 - 48:
        print "Fitted position was outside boundary - setting flux to zero"
        flux = 0.0
    if rms < 0.00000000000000000001:
        rms = 999999999999999999999
    targetstatout.write("Peak flux (mJy/beam): " + str(flux) + "\n")
    targetstatout.write("Integrated flux (mJy):" + str(intflux) + " +/- " + str(intfluxerr) + "\n")
    targetstatout.write("CC flux (mJy):        " + str(1000*ccflux) + "\n")
    targetstatout.write("S/N:                  " + str(flux/rms) + "\n")
    rasub = 0
    if len(centrerasplit[2]) > 4: #Large ref pix, two normally separate tokens are together
        rasub = 1
    decsub = 0
    if len(centredecsplit[2]) > 4: #Large ref pix, two normally separate tokens are together
        decsub = 1
    centrerahours = float(centrerasplit[5-rasub]) + \
                    float(centrerasplit[6-rasub])/60.0 + \
                    float(centrerasplit[7-rasub])/3600.0
    centredecdegs = float(centredecsplit[5-decsub]) + \
                    float(centredecsplit[6-decsub])/60.0 + \
                    float(centredecsplit[7-decsub])/3600.0
    srcrahours = float(sourcerasplit[2]) +  float(sourcerasplit[3])/60.0 + \
                 float(sourcerasplit[4])/3600.0
    srcdecdegs = float(sourcedecsplit[2]) +  float(sourcedecsplit[3])/60.0 + \
                 float(sourcedecsplit[4])/3600.0
    raoff = (srcrahours-centrerahours)*15*60*60*1000*\
            math.cos(centredecdegs*math.pi/180.0)
    decoff = (srcdecdegs-centredecdegs)*60*60*1000
    targetstatout.write("Actual RA:            " + sourcerasplit[2] + ":" + \
                        sourcerasplit[3] + ":" + sourcerasplit[4] + "\n")
    targetstatout.write("Actual Dec:           " + sourcedecsplit[2] + ":" + \
                        sourcedecsplit[3] + ":" + sourcedecsplit[4] + "\n")
    targetstatout.write("Fit:                  " + str(fitmin*1000) + "x" + \
                        str(fitmaj*1000) + " at " + fitpasplit[4] + \
                        " degrees; beam " + str(beammin) + "x" + \
                        str(beammaj) + " at " + str(beampa) + " degrees\n")
    formatline = "Deconvolved fit:      %.2fx%.2f at %.1f degrees ; min/max" + \
                 "(%.2f/%.2f)_(%.2f/%.2f)_(%.1f/%.1f)\n"
    writeline = formatline % (dfitmin*1000, dfitmaj*1000, dfitpa, dfitminminmax[0],
                              dfitminminmax[1]*1000, dfitmajminmax[0]*1000, dfitmajminmax[1]*1000,
                              dfitpaminmax[0], dfitpaminmax[1])
    targetstatout.write(writeline)
    raerr = float(sourcerasplit[6])*1000
    decerr = float(sourcedecsplit[6])*1000
    targetstatout.write("Est. RA error (mas):  " + \
                        str(raerr*math.cos(centredecdegs*math.pi/180.0)*15.0) \
                        + "\n")
    targetstatout.write("Est. RA error (hms):  " + str(raerr) + "\n")
    targetstatout.write("Est. Dec error (mas): " + str(decerr) + "\n")
    targetstatout.close()

    return flux, flux/rms
    
##### Use JMFIT to get a position and error estimate ###########################
def jmfit(imagefile, jmfitfile, target, stokesi, nifs = 4, pixwindow=20, exactmjd=-1):
    directory = os.path.dirname(imagefile)

    todo = []
    beammaj = 0
    beammin = 0
    beampa = 0
    if stokesi:
        for i in range(nifs):
            todo.append('i.' + str(i+1))
        todo.append('ii.a')
    else:
        for i in range(nifs):
            todo.append('rr.' + str(i+1))
            todo.append('ll.' + str(i+1))
        todo.append('rr.a')
        todo.append('ll.a')
        todo.append('ii.a')
    imagedata = AIPSImage('JUNK', 'IMG', 1, 1)
    outdata = AIPSImage('JUNK', 'IMG', 1, 2)
    mjd = 99999
    centrerasplit = " : : : : :00:00:00.0000".split(':')
    centredecsplit = " : : : : :00:00:00.0000".split(':')
    if stokesi:
        targetstatout = open(jmfitfile + ".stokesi.stats", "w")
    else:
        targetstatout = open(jmfitfile + ".stats", "w")
    for ifpol in todo:
        if not os.path.exists(imagefile + '.' + ifpol):
            print "*******************************************"
            print "Could not find file " + imagefile + '.' + ifpol
            print "Setting S/N to 0"
            print "*******************************************"
            targetstatout.write("Pulsar " + target + ": MJD " + str(mjd) + \
                                ": Frequency " + ifpol[-1] + \
                                ", pol " + ifpol[:2].upper() + ":\n")
            targetstatout.write("Centre RA:            " + centrerasplit[5] + ":" + \
                                centrerasplit[6] + ":" + centrerasplit[7] + "\n")
            targetstatout.write("Centre Dec:           " + centredecsplit[5] + ":" + \
                                centredecsplit[6] + ":" + centredecsplit[7] + "\n")
            targetstatout.write("Flux (mJy):           0.0\n")
            targetstatout.write("S/N:                  0.0\n")
            targetstatout.write("RA offset (mas):      0.0\n")
            targetstatout.write("Dec offset (mas):     0.0\n")
            targetstatout.write("Actual RA:            " + centrerasplit[2] + ":" + \
                                centrerasplit[3] + ":" + centrerasplit[4] + "\n")
            targetstatout.write("Actual Dec:           " + centredecsplit[2] + ":" + \
                                centredecsplit[3] + ":" + centredecsplit[4] + "\n")
            targetstatout.write("Fit:                  9999.9x9999.9 at 0 " + \
                                "degrees; beam " + str(beammin) + "x" + \
                                str(beammaj) + " at " + str(beampa) + " degrees\n")
            targetstatout.write("Est. RA error (mas):  999999999.9\n")
            targetstatout.write("Est. RA error (hms):  999999999.9\n")
            targetstatout.write("Est. Dec error (mas): 999999999.9\n")
            continue
        if imagedata.exists():
            imagedata.zap()
        if outdata.exists():
            outdata.zap()

        tempinfits = directory + "/templink.fits"
        os.system("rm -f " + tempinfits)
        os.system("ln -s " + imagefile + '.' + ifpol + " " + tempinfits)
        fitld = AIPSTask('fitld', version = aipsver)
        fitld.datain = tempinfits
        fitld.outdata = imagedata
        fitld.dotable = 1
        fitld.douvcomp = -1
        fitld()

        if exactmjd > 0:
            mjd = exactmjd
        else:        
            obsdate = imagedata.header.date_obs
            year = int(obsdate[0:4])
            month = int(obsdate[5:7])
            day = int(obsdate[8:10])
            mjd = year*367 - int(7*(year + int((month + 9)/12))/4) + \
                  int(275*month/9) + day - 678987

        imean = AIPSTask("imean", version = aipsver)
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
        print xpix, ypix
        blc = [None,0,0]
        trc = [None,0,0]
        blc[1:] = [xpix-pixwindow/2, ypix-pixwindow/2]
        trc[1:] = [xpix+pixwindow/2, ypix+pixwindow/2]
        print blc
        print trc
        if blc[1] < 0: blc[1] = 1
        if blc[2] < 0: blc[2] = 1
        if trc[1] > int(imagedata.header.naxis[0]): trc[1] = int(imagedata.header.naxis[0])-1
        if trc[2] > int(imagedata.header.naxis[1]): trc[2] = int(imagedata.header.naxis[1])-1

        
        jmfit = AIPSTask('jmfit', version = aipsver)
        jmfit.indata = imagedata
        jmfit.outdata = outdata
        jmfit.niter = 4000
        jmfit.ngauss = 1
        jmfit.fwidth[1][1] = 0
        jmfit.fwidth[1][2] = 0
        jmfit.fwidth[1][3] = 0
        jmfit.domax[1:] = [1]
        jmfit.dopos[1][1] = 1
        jmfit.dopos[1][2] = 1
        jmfit.dowidth[1][1] = 1
        jmfit.dowidth[1][2] = 1
        jmfit.dowidth[1][3] = 1
        #jmfit.blc[1:] = [int(imagedata.header.crpix[0])-48,
        #                 int(imagedata.header.crpix[1])-48]
        #jmfit.trc[1:] = [int(imagedata.header.crpix[0])+48,
        #                 int(imagedata.header.crpix[1])+48]
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
            print "Did not find solution for " + imagefile + ", band " + ifpol
            continue
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
        while not "Beam=" in jmfitmessage[msgindex] and msgindex < len(jmfitmessage):
            msgindex = msgindex + 1
        if msgindex >= len(jmfitmessage):
            print "Did not find solution for " + imagefile + ", band " + ifpol
            continue
        beamsplit = jmfitmessage[msgindex].split()
        if not len(beamsplit) == 10:
            print "*****WARNING*****: Poorly formatted beam string!"
        beammaj = float(beamsplit[2])*pixsize
        beammin = float(beamsplit[4])*pixsize
        beampa  = float(beamsplit[8])*pixsize
        msgindex += 1
        while jmfitmessage[msgindex].find('********* Solution from JMFIT') == -1 \
                  and msgindex < len(jmfitmessage):
            msgindex = msgindex + 1
        if msgindex >= len(jmfitmessage):
            print "Did not find solution for " + imagefile + ", band " + ifpol
            continue
        sourcerasplit = jmfitmessage[msgindex+7].split()
        sourcedecsplit = jmfitmessage[msgindex+8].split()
        fluxsplit = jmfitmessage[msgindex+3].split()
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
        targetstatout.write("Pulsar " + target + ": MJD " + str(mjd) + \
                            ": Frequency " + ifpol[-1] + \
                            ", pol " + ifpol[:2].upper() + ":\n")
        targetstatout.write("Centre RA:            " + centrerasplit[5] + ":" + \
                            centrerasplit[6] + ":" + centrerasplit[7] + "\n")
        targetstatout.write("Centre Dec:           " + centredecsplit[5] + ":" + \
                            centredecsplit[6] + ":" + centredecsplit[7] + "\n")
        if fluxsplit[4] == '+/-':
            flux = float(fluxsplit[3][1:])*1000.0
            rms = float(fluxsplit[5])*1000.0
        else:
            flux = float(fluxsplit[4])*1000.0
            rms = float(fluxsplit[6])*1000.0
        targetstatout.write("Flux (mJy):           " + str(flux) + "\n")
        targetstatout.write("S/N:                  " + str(flux/rms) + "\n")
        centrerahours = float(centrerasplit[5]) +  float(centrerasplit[6])/60.0 \
                        + float(centrerasplit[7])/3600.0
        centredecdegs = float(centredecsplit[5]) +  float(centredecsplit[6])/60.0 \
                        + float(centredecsplit[7])/3600.0
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

##### Change the "source" info in one datasets header to match another #########
def match_headersource(template, toupdate):
    wizdata = WizAIPSUVData(toupdate)
    wizdata.header.object = template.header.object
    wizdata.header.crval = template.header.crval
    wizdata.header.update()

##### Normalise a UV dataset (divide through by model, weight appropriately) ###
def normaliseUVData(uvdata, clnimage, normdata, beginif=-1, endif=-1):
    tempdata = AIPSUVData("JUNKNORM", "NORM", 1, 1)
    if tempdata.exists():
        tempdata.zap()

    uvsub = AIPSTask('uvsub', version = aipsver)
    uvsub.indata = uvdata
    uvsub.in2data = clnimage
    if beginif > 0 or endif > 0:
        uvsub.outdata = tempdata
    else:
        uvsub.outdata = normdata
    uvsub.nmaps = 1
    uvsub.channel = 0
    uvsub.bif = 1
    uvsub.eif = 0
    if beginif > 0:
        uvsub.bif = beginif
    if endif > 0:
        uvsub.eif = endif
    uvsub.inver = 0
    uvsub.flux = 0
    uvsub.cmethod = 'DFT'
    uvsub.cmodel = 'COMP'
    uvsub.opcode = 'DIV'
    uvsub()

    if beginif > 0 or endif > 0:
        if endif < 0:
            endif = 9999
        numifs = getNumIFs(tempdata)
        uvflg = AIPSTask('uvflg', version = aipsver)
        uvflg.indata = tempdata
        for f in range(1,numifs+1):
            if f < beginif or f > endif:
                uvflg.bif = f
                uvflg.eif = f
                uvflg()
        splat = AIPSTask('splat', version = aipsver)
        splat.indata = tempdata
        splat.flagver = 1
        splat.outdata = normdata
        splat()

    #wtmod = AIPSTask('wtmod')
    #wtmod.indata = normdata
    #wtmod.outdata = normdata
    #wtmod.aparm[1] = snr
    #wtmod()

##### Generate a 3 character hash string from an str (max val 10000) ###########
def gethash(hashstr):
    hashval = abs(float(hashstr))
    while hashval < 1000.0:
        hashval = hashval*10
    hashval = int(hashval)
    hasha = hashval/576
    hashb = (hashval - 576*hasha)/24
    hashc = hashval - 576*hasha - 24*hashb
    return "%c%c%c" % (chr(ord('a') + hasha), chr(ord('a') + hashb), 
                       chr(ord('a') + hashc))

class calibrate_target_phase_with_two_colinear_phscals:
    def __init__(s, inbeamuvdata):
        s.inbeamuvdata = inbeamuvdata
        #s.compile_into_table()
        #s.solve_phase_ambiguity_semi_automatically()
    def read_inbeamselfcalp1_solutions(s):
        wizuvdata = WizAIPSUVData(s.inbeamuvdata)
        sntable = wizuvdata.table('SN', 0)
        s.numantennas = len(s.inbeamuvdata.antennas)
        times = np.array([])
        antenna_nos = np.array([])
        phi_degs = np.array([])
        row_nos = np.array([])
        j = 0
        for row in sntable:
            mag = row.real1[0]**2 + row.imag1[0]**2
            if mag < 1.1: ## just a sanity check, supposed to 1, but allow for some round-up error
                #row.imag1 = [0,0]
                times = np.append(times, row.time)
                antenna_nos = np.append(antenna_nos, row.antenna_no)
                phi_rad = math.atan2(row.imag1[0], row.real1[0])
                phi_deg = phi_rad/math.pi*180
                phi_degs = np.append(phi_degs, phi_deg)
                row_nos = np.append(row_nos, int(j))
            j += 1
        #print phi_degs
        return row_nos, antenna_nos, times, phi_degs
    def read_inbeamselfcalpn_solutions(s):
        wizuvdata = WizAIPSUVData(s.inbeamuvdata)
        sntable = wizuvdata.table('SN', 0)
        num_if  = int(sntable.keywords['NO_IF'])
        phi_degs = []
        for row in sntable:
            imag1 = np.array(row.imag1)
            real1 = np.array(row.real1)
            phi_rad = np.arctan2(imag1, real1)
            phi_deg = phi_rad/math.pi*180
            phi_degs.append(phi_deg)
        print phi_degs, len(phi_degs)
        return phi_degs
    def compile_into_table(s):
        from astropy.table import Table
        [row_nos, antenna_nos, times, phi_degs] = s.read_inbeamselfcalp1_solutions()
        s.t = Table([row_nos, antenna_nos, times, phi_degs], names=['row_no', 'antenna_no', 'time', 'phi'])
        print s.t
    def edit_inbeamselfcalpn_in_AIPS_and_write_out(s, correction_factor, snver, outputsntable):
        phi_degs = s.read_inbeamselfcalpn_solutions()
        wizuvdata = WizAIPSUVData(s.inbeamuvdata)
        sntable = wizuvdata.table('SN', snver)
        num_if  = int(sntable.keywords['NO_IF'])
        j = 0
        for row in sntable:
            #if j in s.t2['row_no']:
            #    index_in_t2 = s.t2['row_no']==j
            for i in range(num_if):
                phi = phi_degs[j][i]
                if phi != 45: #INDE points
                    phi *= correction_factor
                    phi_rad = phi/180*math.pi
                    real = math.cos(phi_rad)
                    imag = math.sin(phi_rad)
                    for k in [1,2]:    
                        exec("row.imag%d[i] = imag" % k)
                        exec("row.real%d[i] = real" % k)
                row.update()
            j += 1
        sntable.close()
        
        if os.path.exists(outputsntable):
            os.remove(outputsntable)
        writetable(s.inbeamuvdata, 'SN', snver, outputsntable)
        
    def interactively_solve_phase_ambiguity(s, plotpath2save):
        from astropy.table import Table
        import pickle
        for parameter in ['row_no', 'antenna_no', 'time', 'phi']:
            exec('%ss = np.array([])' % parameter)
        saved_phase_edit = plotpath2save+'/.corrected_phases_inbeam_selfcal'
        t = s.t
        if os.path.exists(saved_phase_edit):
            choise = raw_input('Do you want to use saved phase edit and continue the edit? Press y if affirmative: ')
            if choise == 'y' or choise=='yes':
                print "Continue on the saved phase edit."
                readfile = open(saved_phase_edit, 'r')
                t = pickle.load(readfile)
                readfile.close()
            else:
                print "Start new phase edit."
        
        phase_shifts = []
        for i in range(1,s.numantennas+1):
            index = t['antenna_no']==i
            eachAnt = t[index]
            if len(eachAnt)==0:
                continue
            if phase_shifts == 's':
                for parameter in ['row_no', 'antenna_no', 'time', 'phi']:
                    exec("%ss = np.append(%ss, eachAnt['%s'])" % (parameter, parameter, parameter))
                continue
            if max(eachAnt['phi'])>=90 or min(eachAnt['phi'])<=-90:
                print eachAnt['row_no']
                print("Here is the diagnostic plot for phase edit on a new antenna.\n")
                phase_shifts = []
                while not phase_shifts in ['s', 'n']:
                    s.plot_diagnostic_phi_versus_time_for_each_antenna1(eachAnt, i)
                    phase_shifts = raw_input("\nYou can press n to proceed to the next antenna and s to immediately exit.\nEnter an array of two integers N1,N2, separated by comma (e.g. 2,3), meaning shifting phases of the N1th and onward points (start from 0) by N2*360 degree (minus N2 accepted).\nNote that if -10<N2<-4, then the point would be deleted; when N2<-9, then the phases of the N1th and onward points will be deleted: ")
                    try: 
                        N1 = int(phase_shifts.split(',')[0].strip())
                        N2 = int(phase_shifts.split(',')[1].strip())
                    except ValueError:
                        if not phase_shifts in ['s', 'n']:  
                            print "Input format does not match the requirement. Start over again."
                        continue
                    if N2>-5:
                        eachAnt['phi'][N1:len(eachAnt)] += N2*360
                    elif N2>-10 and N2<-4:
                        eachAnt.remove_row(N1)
                    else:
                        eachAnt.remove_rows(range(N1,len(eachAnt)))
            for parameter in ['row_no', 'antenna_no', 'time', 'phi']:
                exec("%ss = np.append(%ss, eachAnt['%s'])" % (parameter, parameter, parameter))
            print row_nos
        s.t1 = Table([row_nos, antenna_nos, times, phis], names=['row_no', 'antenna_no', 'time', 'phi'])
        print s.t1
        writefile = open(saved_phase_edit, 'w')
        pickle.dump(s.t1, writefile)
        writefile.close()
    def copy_inbeamselfcal_sntable(s, sntable):
        oldsntable = sntable.replace('.icalib.p1.sn', '_original.icalib.p1.sn')
        os.system('cp %s %s' % (sntable, oldsntable))
        return oldsntable
    def load_final_inbeamselfcal_phase_edit_and_prepare_for_edit_in_AIPS(s, final_phase_edit, correction_factor):
        from astropy.table import Table
        import pickle
        readfile = open(final_phase_edit, 'r')
        t = pickle.load(readfile)
        readfile.close()
        phis = t['phi'] * correction_factor
        reals = np.cos(phis*math.pi/180)
        imags = np.sin(phis*math.pi/180)
        s.t2 = Table([t['row_no'], t['antenna_no'], t['time'], reals, imags], 
                                names=['row_no', 'antenna_no', 'time', 'real', 'imag'])
        print s.t2
    def edit_AIPS_sntable_and_write_out(s, snver, outputsntable):
        wizuvdata = WizAIPSUVData(s.inbeamuvdata)
        sntable = wizuvdata.table('SN', snver)
        num_if  = int(sntable.keywords['NO_IF'])
        j = 0
        for row in sntable:
            print j in s.t2['row_no']
            if j in s.t2['row_no']:
                index_in_t2 = s.t2['row_no']==j
                for k in [1,2]:
                    exec("row.imag%d = list(s.t2['imag'][index_in_t2]*np.ones(num_if))" % k)
                    exec("row.real%d = list(s.t2['real'][index_in_t2]*np.ones(num_if))" % k)
                row.update()
            j += 1
        sntable.close()
        
        if os.path.exists(outputsntable):
            os.remove(outputsntable)
        writetable(s.inbeamuvdata, 'SN', snver, outputsntable)
    def plot_diagnostic_phi_versus_time_for_each_antenna(s, antenna_no):
        import matplotlib.pyplot as plt
        index = s.t['antenna_no']==antenna_no
        plt.scatter(s.t['time'][index], s.t['phi'][index], marker='.')
        plt.ylim(-180,180)
        plt.xlabel('time (day)')
        plt.ylabel('phase (degree)')
        plt.title('phase-time evolution for antenna%d' % antenna_no)
        #plt.draw()
        plt.show()
        plt.clf()
    def plot_diagnostic_phi_versus_time_for_each_antenna1(s, t, antenna_no):
        import matplotlib.pyplot as plt
        plt.scatter(t['time'], t['phi'], marker='.')
        plt.xlabel('time (day)')
        plt.ylabel('phase (degree)')
        plt.title('phase-time evolution for antenna%d' % antenna_no)
        plt.show()
        plt.clf()
    def plot_phi_versus_time_for_each_antenna(s, plotpath2save, projectcode, before_phase_correction=True):
        import matplotlib.pyplot as plt
        if before_phase_correction:
            t = s.t
            before_or_after = 'before'
        else:
            t = s.t1
            before_or_after = 'after'
        for i in range(1,s.numantennas+1):
            index = t['antenna_no']==i
            if len(t[index])==0:
                continue
            plt.scatter(t['time'][index], t['phi'][index], marker='.')
            if before_phase_correction:
                plt.ylim(-180,180)
            plt.xlabel('time (day)')
            plt.ylabel('phase (degree)')
            plt.title('phase-time evolution for antenna%d' % i)
            plt.savefig("%s/%s_inbeamselfcal_phi_time_for_antenna%d_%s_phase_correction.eps" % (plotpath2save, projectcode, i, before_or_after))
            plt.clf()

class calibrate_target_phase_with_three_phscals:
    def __init__(s, targetname, cal1, cal2, cal3):
        pass
