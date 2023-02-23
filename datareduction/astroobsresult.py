#!/usr/bin/python

import math, astro_utils, os, sys
import numpy as np

class PmparInput:
    def __init__(self, filename=""):
        self.defaultepoch = 56000
        self.headerendlinenum = 0
        self.haveepoch = False
        self.numkeywords = 0
        self.astrometric_points = []
        self.numpoints = 0
        self.numepochs = 0
        self.pmparlines = []
        lastmjd = 0
        self.mjds = []
        self.mjdstrs = []
        self.rarad = []
        self.decrad = []
        self.perturbedrarad = []
        self.perturbeddecrad = []
        self.ra_err = []
        self.dec_err = []
        self.sysra_err = []
        self.sysdec_err = []
        self.perturbedra_err = []
        self.perturbeddec_err = []
        if not filename == "":
            self.pmparlines = open(filename).readlines()
            while '=' in self.pmparlines[self.headerendlinenum] or \
                  self.pmparlines[self.headerendlinenum].rstrip() == "" or \
                  self.pmparlines[self.headerendlinenum].lstrip()[0] == '#':
                if ('epoch' in self.pmparlines[self.headerendlinenum] or \
                    'Epoch' in self.pmparlines[self.headerendlinenum]) and \
                   (not self.pmparlines[self.headerendlinenum].lstrip()[0] == '#') :
                    self.haveepoch = True
                if '=' in self.pmparlines[self.headerendlinenum]:
                    self.numkeywords = self.numkeywords + 1
                print("Skipping line " + self.pmparlines[self.headerendlinenum][0:-1])
                self.headerendlinenum = self.headerendlinenum + 1
            if not self.haveepoch:
                self.numkeywords = self.numkeywords+1       
            startline = self.numkeywords
            for line in self.pmparlines[self.headerendlinenum:]:
                if not line.lstrip() == "" and not line.lstrip()[0] == '#':
                    self.astrometric_points.append(line)
                    self.numpoints = self.numpoints + 1
                    splitline = line.split()
                    mjd = float(splitline[0])
                    if abs(mjd-lastmjd) > 1.0:
                        self.numepochs = self.numepochs + 1
                        lastmjd = mjd
                        self.mjds.append(mjd)
                        self.mjdstrs.append(splitline[0])
                        self.rarad.append([])
                        self.decrad.append([])
                        self.perturbedrarad.append([])
                        self.perturbeddecrad.append([])
                        self.ra_err.append([])
                        self.dec_err.append([])
                        self.sysra_err.append([])
                        self.sysdec_err.append([])
                        self.perturbedra_err.append([])
                        self.perturbeddec_err.append([])
                    try:
                        rarad = astro_utils.stringToRad(splitline[1],True)
                        decrad = astro_utils.stringToRad(splitline[3],False)
                    except ValueError:
                        print("Bad position entry", splitline[1], splitline[3])
                        rarad = 0
                        decrad = 0
                    self.rarad[self.numepochs-1].append(rarad)
                    self.decrad[self.numepochs-1].append(decrad)
                    self.perturbedrarad[self.numepochs-1].append(rarad)
                    self.perturbeddecrad[self.numepochs-1].append(decrad)
                    self.ra_err[self.numepochs-1].append(float(splitline[2]))
                    self.dec_err[self.numepochs-1].append(float(splitline[4]))
                    self.sysra_err[self.numepochs-1].append(0.0)
                    self.sysdec_err[self.numepochs-1].append(0.0)
                    self.perturbedra_err[self.numepochs-1].append(float(splitline[2]))
                    self.perturbeddec_err[self.numepochs-1].append(float(splitline[4]))

    def getAverageError(self):
        for i in range(len(self.decrad)):
            print(len(self.decrad[i]))
        print(np.array(self.decrad[0]))
        print(np.array(self.decrad))
        meanraerr = np.mean(np.array(self.ra_err)*15*np.cos(np.array(self.decrad)))
        meandecerr = np.mean(np.array(self.dec_err))
        return math.sqrt(meanraerr*meanraerr + meandecerr*meandecerr)

    def setSystematicError(self, sysraerr, sysdecerr, rchisqscaling, mjdmatch=0.5):
        for i in range(len(self.mjdstrs)):
            mjd = self.mjds[i]
            scale = [-1,-1]
            for m in list(rchisqscaling.keys()):
                if math.fabs(mjd - float(m)) < mjdmatch:
                    scale = rchisqscaling[m]
                    break
            if scale[0] < 0:
                print("Couldn't match " + str(mjd) + " in rchisqscaling dict!")
                sys.exit()
            for j in range(len(self.sysra_err[i])):
                self.sysra_err[i][j] = sysraerr*scale[0]
                self.sysdec_err[i][j] = sysdecerr*scale[1]

    def perturbPositions(self, posoffsetsmas):
        self.perturbedrarad = []
        self.perturbeddecrad = []
        count = 0
        keys = list(posoffsetsmas.keys())
        for mjd in self.mjdstrs:
            self.perturbedrarad.append([])
            self.perturbeddecrad.append([])
            try:
                usekey = mjd
                for k in keys:
                    if math.fabs(float(k) - float(mjd)) < 0.1:
                        usekey = k
                offsetmas = posoffsetsmas[usekey]
                racount = 0
                for ra in self.rarad[count]:
                    self.perturbedrarad[count].append(ra + offsetmas[0]*math.pi/(math.cos(self.decrad[0][0])*180*60*60*1000.0))
                    if len(offsetmas) > 2:
                        self.perturbedra_err[count][racount] = math.sqrt(self.ra_err[count][racount]*self.ra_err[count][racount] + offsetmas[2]*offsetmas[2])
                    racount += 1
                deccount = 0
                for dec in self.decrad[count]:
                    self.perturbeddecrad[count].append(dec + offsetmas[1]*math.pi/(180*60*60*1000.0))
                    if len(offsetmas) > 2:
                        self.perturbeddec_err[count][deccount] = math.sqrt(self.dec_err[count][deccount]*self.dec_err[count][deccount] + offsetmas[3]*offsetmas[3])
                    deccount += 1
#                    rastr1, decstr1 = astro_utils.posradians2string(ra, dec)
#                    rastr2, decstr2 = astro_utils.posradians2string(ra, self.perturbeddecrad[count][-1])
            except KeyError:
                print("Could not find " + mjd + " in the position offset dictionary...")
            count += 1

    def generateOffsets(self):
        flatrarad = [item for sublist in self.rarad for item in sublist]
        flatdecrad = [item for sublist in self.decrad for item in sublist]
        meanrarad = sum(flatrarad) / len(flatrarad)
        meandecrad = sum(flatdecrad) / len(flatdecrad)
        offsetmas = {}
        for r,d,r_e,d_e,mjd in zip(self.rarad, self.decrad, self.ra_err, self.dec_err,self.mjdstrs):
            if len(r) > 1:
                print("This only works for single-fit-per-epoch type files!")
                sys.exit()
            offsetmas[mjd] = [(meanrarad-r[0])*3.6e6*180/math.pi, (meandecrad-d[0])*3.6e6*180/math.pi, r_e[0], d_e[0]]
            print(offsetmas[mjd], r[0], meanrarad)
        print(offsetmas)
        return offsetmas

    def addEpoch(self, jmfitresultarray, snrcut=0.0):
        self.numepochs += 1
        self.rarad.append([])
        self.decrad.append([])
        self.ra_err.append([])
        self.dec_err.append([])
        self.mjds.append(jmfitresultarray[0].mjd)
        self.mjdstrs.append("%.4f" % jmfitresultarray[0].mjd)
        for jmfitresult in jmfitresultarray:
            if jmfitresult.snr < snrcut:
                print("Skipping result with S/N " + str(jmfitresult.snr))
                continue
            self.rarad[-1].append(jmfitresult.rarad)
            self.decrad[-1].append(jmfitresult.decrad)
            self.ra_err[-1].append(jmfitresult.raerrhhh/1000.0)
            self.dec_err[-1].append(jmfitresult.decerrmas/1000.0)
        if len(self.rarad[-1]) == 0:
            self.numepochs -= 1
            self.rarad = self.rarad[:-1]
            self.decrad = self.decrad[:-1]
            self.ra_err = self.ra_err[:-1]
            self.dec_err = self.dec_err[:-1]
            self.mjds = self.mjds[:-1]
            self.mjdstrs = self.mjdstrs[:-1]
       
    def writeSysErrorEstimates(self, filename):
        tempout = open(filename, 'w')
        for j in range(self.numepochs):
            for k in range(len(self.rarad[j])):
                tempout.write(self.mjdstrs[j] + " " + str(self.sysra_err[j][k]) + " " + str(self.sysdec_err[j][k]) + "\n")
        tempout.close()

    def writePmparFile(self, filename, doperturbed=False,addsys=False):
        tempout = open(filename, 'w')
        tempout.writelines(self.pmparlines[0:self.headerendlinenum])
        if not self.haveepoch:
            tempout.write("epoch = %.1f\n\n" % self.defaultepoch)
        for j in range(self.numepochs):
            for k in range(len(self.rarad[j])):
                if doperturbed:
                    ra, dec = astro_utils.posradians2string(self.perturbedrarad[j][k], self.perturbeddecrad[j][k])
                    raerr = self.perturbedra_err[j][k]
                    decerr = self.perturbeddec_err[j][k]
                    if addsys:
                        raerr = math.sqrt(self.perturbedra_err[j][k]**2 + self.sysra_err[j][k]**2)
                        decerr = math.sqrt(self.perturbeddec_err[j][k]**2 + self.sysdec_err[j][k]**2)
                else:
                    ra, dec = astro_utils.posradians2string(self.rarad[j][k], self.decrad[j][k])
                    raerr = self.ra_err[j][k]
                    decerr = self.dec_err[j][k]
                    if addsys:
                        raerr = math.sqrt(self.ra_err[j][k]**2 + self.sysra_err[j][k]**2)
                        decerr = math.sqrt(self.dec_err[j][k]**2 + self.sysdec_err[j][k]**2)
                tempout.write(self.mjdstrs[j] + " " + ra + " " + str(raerr) + " " + \
                              dec + " " + str(decerr) + "\n")
        tempout.close()

class AstroObsResult:
    def __init__(self, obj, exp, lines, isexp):
        self.object = obj
        self.experiment = exp
        self.isexploratory = isexp
        self.gategain = -1.0
        self.lowerlimit = True
        self.mjd = float(lines[0].split()[3][:-1])
        self.raerrmas = float(lines[10].split()[-1])
        self.raerrhhh = float(lines[11].split()[-1])
        self.decerrmas = float(lines[12].split()[-1])
        self.snr = float(lines[4].split()[-1])
        self.flux = float(lines[3].split()[-1])
        self.ra = lines[7].split()[-1]
        self.dec = lines[8].split()[-1]
        self.rarad = 0
        self.decrad = 0
        if self.snr != 0.0:
            self.rarad = astro_utils.stringToRad(self.ra, True)
            self.decrad = astro_utils.stringToRad(self.dec, False)
        tempsplit = lines[9].split()
        self.fitmin = float(tempsplit[1].split('x')[0])
        self.fitmaj = float(tempsplit[1].split('x')[1])
        self.fitpa  = float(tempsplit[3])
        self.hasbeam = False
        if len(tempsplit) > 6 and "x" in tempsplit[6]:
            self.hasbeam = True
            self.beammin = float(tempsplit[6].split('x')[0])
            self.beammaj = float(tempsplit[6].split('x')[1])
            self.beampa  = float(tempsplit[8])
            self.beamparad = math.pi*self.beampa/180.0

    def updateRA(self, newrarad, newraraderr=-1):
        self.rarad = newrarad
        newra, xxx = astro_utils.posradians2string(newrarad, 1)
        self.ra = newra
        if newraraderr > 0:
            self.raerrmas = newraraderr*180*60*60*1000/math.pi
            self.raerrhms = self.raerrmas/(15*math.cos(self.decrad))

    def updateDec(self, newdecrad, newdecraderr=-1):
        self.decrad = newdecrad
        xxx, newdec = astro_utils.posradians2string(1, newdecrad)
        self.dec = newdec
        if newdecraderr > 0:
            self.decerrmas = newdecraderr*180*60*60*1000/math.pi

    def write(self, outputfile, overwrite=False):
        if os.path.exists(outputfile) and not overwrite:
            print("%s exists and overwrite is False - aborting" % outputfile)
            sys.exit()
        output = open(outputfile, "w")
        output.write("Pulsar %s: MJD %.4f: Frequency X, pol ?.:\n" % (self.object, self.mjd))
        output.write("%s%s\n" % ("Centre RA:".ljust(22), self.ra))
        output.write("%s%s\n" % ("Centre Dec:".ljust(22), self.dec))
        output.write("%s%.4f\n" % ("Flux (mJy):".ljust(22), self.flux))
        output.write("%s%.4f\n" % ("S/N:".ljust(22), self.snr))
        output.write("%s%.4f\n" % ("RA offset (mas)".ljust(22), 0.0))
        output.write("%s%.4f\n" % ("Dec offset (mas)".ljust(22), 0.0))
        output.write("%s%s\n" % ("Actual RA:".ljust(22), self.ra))
        output.write("%s%s\n" % ("Actual Dec:".ljust(22), self.dec))
        output.write("%s%.4fx%.4f at %.2f degrees\n" % ("Fit:".ljust(22), self.fitmin, self.fitmaj, self.fitpa))
        output.write("%s%.4f\n" % ("Est. RA error (mas):".ljust(22), self.raerrmas))
        output.write("%s%.4f\n" % ("Est. RA error (hms):".ljust(22), self.raerrhms))
        output.write("%s%.4f\n" % ("Est. Dec error (mas):".ljust(22), self.decerrmas))
        output.close()

class BootstrapResult:
    def __init__(self, filename):
        lines = open(filename).readlines()
        self.mode = {}
        self.interval1sig = {}
        self.interval2sig = {}
        self.interval3sig = {}
        self.usedparams = ["ra", "dec", "pm_ra", "pm_dec", "pm_tot", "pi", "distance", "v_t", "omega", "inc", "incmod90"]
        for line in lines:
            splitline = line.split()
            if len(splitline) < 2:
                print("Garbage line " + line)
                continue
            param = splitline[1][:-1].lower()
            if not param in self.usedparams:
                print("Got unknown parameter " + param + " when initialising bootstrapresult")
                print(line)
                sys.exit()
            if "mode" in line:
                self.mode[param] = float(splitline[-1])
            elif "68.0%" in line:
                self.interval1sig[param] = (float(splitline[-5]), float(splitline[-3]), float(splitline[-1]))
            elif "95.0%" in line:
                self.interval2sig[param] = (float(splitline[-5]), float(splitline[-3]), float(splitline[-1]))
            elif "99.7%" in line:
                self.interval3sig[param] = (float(splitline[-5]), float(splitline[-3]), float(splitline[-1]))
            else:
                print("Unparseable line " + line)
            
    def get1sigmaconfidence(self, param):
        if not param.lower() in self.usedparams:
            print("Got unknown parameter " + param)
            sys.exit()
        return self.interval1sig[param.lower()]

    def get2sigmaconfidence(self, param):
        if not param.lower() in self.usedparams:
            print("Got unknown parameter " + param)
            sys.exit()
        return self.interval2sig[param.lower()]

    def get3sigmaconfidence(self, param):
        if not param.lower() in self.usedparams:
            print("Got unknown parameter " + param)
            sys.exit()
        return self.interval3sig[param.lower()]

class AstroFitResult:
    def __init__(self, filename=None, nepochs=0):
        self.nepochs = nepochs
        if filename == None:
            self.epochmjdstring = "-"
            self.rastring = "-"
            self.decstring = "-"
            self.glstring = "-"
            self.gbstring = "-"
            self.pmrastring = "-"
            self.pmdecstring = "-"
            self.pmglstring = "-"
            self.pmgbstring = "-"
            self.pxstring = "-"
            self.diststring = "-"
            self.vtstring = "-"
            self.scatterx = 9999999999.99
            self.scattery = 9999999999.99
            self.rchisqstring = "-"
            self.rchisq = 9999999999.99
            self.px = -999999999.99
            self.pxerr = -999999999.99
        elif os.path.exists(filename):
            lines = open(filename).readlines()
            self.name = "???"
            toadd = 0
            if "Name" in lines[0]:
                self.name = lines[0].split('=')[1].strip()
                toadd = 1
            self.epochmjdstring = lines[0+toadd].split('=')[1].strip()
            self.rastring = lines[1+toadd].split('=')[1].strip()
            self.decstring = lines[2+toadd].split('=')[1].strip()
            self.glstring = lines[3+toadd].split('=')[1].strip()
            self.gbstring = lines[4+toadd].split('=')[1].strip()
            self.pmrastring = lines[5+toadd].split('=')[1].strip().split('(')[0].strip()
            self.pmdecstring = lines[6+toadd].split('=')[1].strip()
            self.pmglstring = lines[7+toadd].split('=')[1].strip()
            self.pmgbstring = lines[8+toadd].split('=')[1].strip()
            self.pxstring = lines[9+toadd].split('=')[1].strip()
            self.px = float(self.pxstring.split()[0])
            self.pxerr = float(self.pxstring.split()[2])
            self.diststring = lines[10+toadd].split('=')[1].strip()
            self.vtstring = lines[11+toadd].split('=')[1].strip()
            self.scatterx = float(lines[14+toadd].split()[2])
            self.scattery = float(lines[15+toadd].split()[2])
            self.rchisqstring = lines[16+toadd].split()[3]
            self.rchisq = 9999999999.99
            if not self.rchisqstring == "inf":
                self.rchisq = float(self.rchisqstring)
            self.rchisqstring = self.rchisqstring
            self.rarad = astro_utils.stringToRad(self.rastring.split()[0], True)
            self.decrad = astro_utils.stringToRad(self.decstring.split()[0], False)
            self.pmramas = float(self.pmrastring.split()[0])
            self.pmdecmas = float(self.pmdecstring.split()[0])
            self.pmramaserr = float(self.pmrastring.split()[2])
            self.pmdecmaserr = float(self.pmdecstring.split()[2])
        else:
            print("Cannot find file " + filename)

    def getApproxPositionAtDate(self, mjd):
        yeardiff = (mjd - float(self.epochmjdstring))/365.25
        rarad = self.rarad + ((self.pmramas*yeardiff*math.pi/(180*60*60*1000))/math.cos(self.decrad))
        decrad = self.decrad + self.pmdecmas*yeardiff*math.pi/(180*60*60*1000)
        return rarad, decrad

class BasicJmfitResult:
    def __init__(self, obj, exp, filename, snrlimit=6.5, clearlyinflux=4.0, limitfile=None, limitscale=0.0):
        lines = open(filename).readlines()
        self.object = obj
        self.experiment = exp
        self.upperlimit = True
        self.mjd = float(lines[0].split()[3][:-1])
        self.ra = lines[3].split()[-1]
        self.dec = lines[4].split()[-1]
        self.rarad = astro_utils.stringToRad(self.ra, True)
        self.decrad = astro_utils.stringToRad(self.dec, False)
        self.raerrmas = float(lines[6].split()[-1])
        self.raerrhhh = float(lines[7].split()[-1])
        self.decerrmas = float(lines[8].split()[-1])
        self.snr = float(lines[2].split()[-1])
        self.peakflux = float(lines[1].split()[-1])
        tempsplit = lines[5].split()
        self.fitmin = float(tempsplit[1].split('x')[0])
        self.fitmaj = float(tempsplit[1].split('x')[1])
        self.fitpa  = float(tempsplit[3])
        self.hasbeam = False
        self.estimatedimagerms = -1.0
        self.peakpixel = -1.0
        self.iscomplex = False
        if len(tempsplit) > 6 and "x" in tempsplit[6]:
            self.hasbeam = True
            self.beammin = float(tempsplit[6].split('x')[0])
            self.beammaj = float(tempsplit[6].split('x')[1])
            self.beampa  = float(tempsplit[8])
        if self.snr > snrlimit or self.peakflux > clearlyinflux:
            self.upperlimit = False
        if not limitfile == None:
            limitlines = open(limitfile).readlines()
            self.estimatedimagerms = float(limitlines[-1].split(':')[-1])*limitscale
            self.measuredimagerms = float(limitlines[-2].split(':')[-1])*limitscale
            self.peakpixel = float(limitlines[-3].split(':')[-1])
            if self.peakpixel > 10.0*self.estimatedimagerms and self.upperlimit:
                print("Source %s was a non-detection, but peak pixel (%.3f) is %.1f x the expected image rms!" % \
                      (self.object, self.peakpixel, self.peakpixel/self.estimatedimagerms))
            if self.upperlimit and self.peakpixel > 12*self.estimatedimagerms:
                self.iscomplex = True

class ExtendedJmfitResult:
    def __init__(self, obj, exp, filename, snrlimit=6.5, clearlyinflux=4.0, limitfile=None, limitscale=0.0):
        lines = open(filename).readlines()
        self.object = obj
        self.experiment = exp
        self.upperlimit = True
        self.mjd = float(lines[0].split()[3][:-1])
        self.ra = lines[5].split()[-1]
        self.dec = lines[6].split()[-1]
        self.rarad = astro_utils.stringToRad(self.ra, True)
        self.decrad = astro_utils.stringToRad(self.dec, False)
        self.raerrmas = float(lines[9].split()[-1])
        self.raerrhhh = float(lines[10].split()[-1])
        self.decerrmas = float(lines[11].split()[-1])
        self.snr = float(lines[4].split()[-1])
        self.peakflux = float(lines[1].split()[-1])
        intflux = lines[2].split(':')[-1].strip()
        if "+/-" in intflux:
            self.intflux = float(intflux.split('+/-')[0].strip())
            self.intfluxerr = float(intflux.split('+/-')[1].strip())
        else:
            self.intflux = float(intflux)
            self.intfluxerr = -1.0
        self.ccflux = float(lines[3].split(':')[-1].strip())
        tempsplit = lines[7].split(':')[-1].strip().split()
        self.fitmin = float(tempsplit[0].split('x')[0])
        self.fitmaj = float(tempsplit[0].split('x')[1])
        self.fitpa  = float(tempsplit[2])
        self.hasbeam = True
        self.estimatedimagerms = -1.0
        self.peakpixel = -1.0
        self.iscomplex = False
        self.beammin = float(tempsplit[5].split('x')[0])
        self.beammaj = float(tempsplit[5].split('x')[1])
        self.beampa  = float(tempsplit[7])
        tempsplit = lines[8].split(':')[-1].strip().split()
        self.deconvmin = float(tempsplit[0].split('x')[0])
        self.deconvmaj = float(tempsplit[0].split('x')[1])
        self.deconvpa  = float(tempsplit[2])
        self.deconvminplusminus = [0,0]
        self.deconvmajplusminus = [0,0]
        self.deconvpaplusminus = [0,0]
        plusminsplit = tempsplit[-1][7:].split('_')
        self.deconvminplusminus[0] = float(plusminsplit[0].split('/')[0][1:])
        self.deconvminplusminus[1] = float(plusminsplit[0].split('/')[1][:-1])
        self.deconvmajplusminus[0] = float(plusminsplit[1].split('/')[0][1:])
        self.deconvmajplusminus[1] = float(plusminsplit[1].split('/')[1][:-1])
        self.deconvpaplusminus[0] = float(plusminsplit[2].split('/')[0][1:])
        self.deconvpaplusminus[1] = float(plusminsplit[2].split('/')[1][:-1])
        if self.snr > snrlimit or self.peakflux > clearlyinflux:
            self.upperlimit = False
        if not limitfile == None:
            limitlines = open(limitfile).readlines()
            if len(limitlines) < 2:
                print("Some problem with limit file " + limitfile)
            else:
                self.estimatedimagerms = float(limitlines[-1].split(':')[-1])*limitscale
                self.measuredimagerms = float(limitlines[-2].split(':')[-1])*limitscale
                self.peakpixel = float(limitlines[-3].split(':')[-1])
                if self.peakpixel > 10.0*self.estimatedimagerms and self.upperlimit:
                    print("Source %s was a non-detection, but peak pixel (%.3f) is %.1f x the expected image rms!" % \
                          (self.object, self.peakpixel, self.peakpixel/self.estimatedimagerms))
                if self.upperlimit and self.peakpixel > 12*self.estimatedimagerms:
                    self.iscomplex = True

class PmparOutput:
    def __init__(self, pmparefile="pmpar_e", pmpartfile="pmpar_t"):
        pmparlines = open(pmparefile).readlines()
        self.meas_mjd = []
        self.meas_ra = []
        self.meas_dec = []
        self.meas_ra_err = []
        self.meas_dec_err = []
        for line in pmparlines:
            self.meas_mjd.append(float(line.split()[0]))
            self.meas_ra.append(float(line.split()[1]))
            self.meas_ra_err.append(float(line.split()[2]))
            self.meas_dec.append(float(line.split()[3]))
            self.meas_dec_err.append(float(line.split()[4]))
        lines = open("pmpar_t").readlines()
        self.pmparMJDs = []
        self.pmparRAs = []
        self.pmparDECs = []
        for line in lines:
            self.pmparMJDs.append(float(line.split()[0]))
            self.pmparRAs.append(float(line.split()[1]))
            self.pmparDECs.append(float(line.split()[2]))


