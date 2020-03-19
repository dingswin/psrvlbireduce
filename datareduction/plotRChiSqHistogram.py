#!/usr/bin/env python
import os,sys,glob,math,yaml
import astro_utils
import numpy
import astroobsresult
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from scipy.stats import chi2

######################################################################################################
# SUBROUTINES
######################################################################################################

# Utility routine to check if a file exists and optionally bail out if it doesn't
def checkFileExists(f, abortifnot=True):
    if not os.path.exists(f):
        print f + " doesn't exist!"
        if abortifnot:
            print "Aborting!"
            sys.exit()

# Routine to plot a histogram of reduced chi squared values, on a log scale from 0.1 to 100
def plotRChiSqHistogram(rchisqlist, outfilename):
    logchisq = []
    for r in rchisqlist:
        logchisq.append(math.log10(r))
    nbins = 16
    n, bins, patches = pyplot.hist(logchisq, nbins, normed=1)
    # Also plot the reduced chi squared for an 8 epoch obs (11 d.o.f.)
    dof = 11
    x = numpy.linspace(-0.6, 0.8, 100)
    c = dof*numpy.power(10, x)
    y = 16*chi2.pdf(c, dof)
    pyplot.plot(x, y)
    pyplot.xlabel('Log10(rchisquared)')
    pyplot.ylabel('Probability')
    #pyplot.axis([0.1, 10, 0, 0.1])
    pyplot.savefig(outfilename)
    pyplot.clf()

#############24 (options, junk) = parser.parse_args()#########################################################################################
# MAIN CODE
######################################################################################################
basedir = os.environ('PSRVLBAUXDIR')
##  THE following were for the microarcsecond-based approach
##separationsearchvals = numpy.linspace(0.0, 30.0, 11) # microarcseconds per deprojected arcmin
##snrsearchvals = numpy.linspace(0.0, 40.0, 21) # microarcseconds per 100/snr
#separationsearchvals = numpy.linspace(8.0, 8.0, 1) # microarcseconds per deprojected arcmin
#snrsearchvals = numpy.linspace(8.0, 8.0, 1) # microarcseconds per 100/snr

## Now for fractions of a beam, rather than microarcseconds
#separationsearchvals = numpy.linspace(0.0008, 0.0011, 4) # fractions of a beam per deprojected arcmin
#snrsearchvals = numpy.linspace(0.004, 0.008, 5) # fractions of a beam per 100/snr
separationsearchvals = numpy.linspace(0.001, 0.001, 1)  # fractions of a beam per deprojected arcmin
snrsearchvals = numpy.linspace(0.006, 0.006, 1) # fractions of a beam per 100/snr

# Run over the basic results (no systematic contribution)
pulsars = sorted(glob.glob("J*"))
basicfits = glob.glob("J*/pulsar/full/*.simple.lsqfit")
numfits = len(basicfits)
print "I have %d basic fits" % (numfits)
errorsum = 0
rchisqlist = []
for bff in basicfits:
    fitresult = astroobsresult.AstroFitResult(bff)
    rchisq = fitresult.rchisq
    rchisqlist.append(rchisq)
    errorsum += math.exp(2*math.fabs(math.log(rchisq)))
print errorsum/numfits
besterrorsum = errorsum/numfits
plotRChiSqHistogram(rchisqlist, "systematicshistograms/rchisqhistogram.nosystematic.png")

deprojout = open("deprojected-separations.txt", "w")
orgrchisqlist = list(rchisqlist)
# Loop over a bunch of values for sys correction based on median(posrefsep/sin(elevation)) and 100/inbeamsnr, testing rchisqsum
for i in separationsearchvals:
    for j in snrsearchvals:
        rchisqlist = []
        errorsum = 0.0
        for pulsar in pulsars:
            # Set up some needed variables for this pulsar
            pulsardir = basedir + pulsar + "/pulsar/full/"
            os.chdir(pulsardir)
            configfile = "../../" + pulsar + ".yaml"
            simplefitfile = pulsar + ".simple.lsqfit"
            if not os.path.exists(simplefitfile):
                continue ########################TEMPORARY #####################
            pmparinputfile = "local.jmfit.pmpar.in"
            checkFileExists(configfile)
            checkFileExists(pmparinputfile)
            checkFileExists(simplefitfile)
            config = yaml.load(open(configfile))
            primaryinbeam = config["primaryinbeam"]
            try:
                posreference = [config["positionreference"]]
            except KeyError:
                posreference = primaryinbeam.split(',')
            pulsarfit = astroobsresult.AstroFitResult(simplefitfile)
            posreferenceseparation = 0.0
            for p in posreference:
                posreferencepmparinfile = "../../" + p +"-divided/full/local.jmfit.pmpar.in"
                if not os.path.exists(posreferencepmparinfile):
                    print "Missing " + posreferencepmparinfile
                    continue
                checkFileExists(posreferencepmparinfile)
                posreferencefit = astroobsresult.PmparInput(posreferencepmparinfile)
                posreferenceseparation += (astro_utils.posdiff(pulsarfit.rarad, pulsarfit.decrad, posreferencefit.rarad[0][0], posreferencefit.decrad[0][0])*180*60/math.pi)**2
            posreferenceseparation = math.sqrt(posreferenceseparation)

            # Parse the standard input file
            pmparinput = astroobsresult.PmparInput(pmparinputfile)

            # For each observation, calculate median deprojected separation and inbeam snr, 
            # use that to get the sys error estimate for this obs, and add it
            statsfiles = sorted(glob.glob("*" + pulsar + "*.stats"))
            epochscales = {}
            allerrors = []
            for statsfile in statsfiles:
                finalobscode = statsfile.split('_')[0]
                orgobscode = finalobscode[:-1]
                if "bd152c" in orgobscode or "bd152d" in orgobscode: 
                    sumfile = "/home/deller/svn_codebase/psrpi/schedules/exploratory/observed/" + orgobscode + ".sum"
                else:
                    sumfile = "/home/deller/svn_codebase/psrpi/schedules/astrometric/observed/" + orgobscode + ".sum"
                checkFileExists(sumfile)
                inrange = False
                deprojected = []
                for line in open(sumfile).readlines():
                    if "SCAN  DAY START UT  SOURCE     TYPE  STATIONS" in line:
                        inrange = True # Beginning of the scan list
                    if "TIME RANGE OF RECORDINGS and TOTAL BYTES" in line:
                        break # Past the end of the scan list
                    if inrange and pulsar in line:
                        splitline = line.split()
                        for e in splitline[5:]:
                            deprojected.append(posreferenceseparation/math.sin(float(e)*math.pi/180.0))
                deprojectedseparation = sorted(deprojected)[len(deprojected)/2]
                deprojout.write("%s: %.2f %.1f\n" % (pulsar, deprojectedseparation, posreferenceseparation))
                inbeamresultfiles = []
                inbeamsnr = 0
                if primaryinbeam == pulsar:
                    inbeamresultfiles.append(statsfile)
                else:
                    for inbeam in primaryinbeam.split(','):
                        inbeam = inbeam.strip()
                        inbeamresultfiles.append("../../" + inbeam + "-divided/full/" + finalobscode + "_" + inbeam + "_full.clean.jmfit.stokesi.stats")
                for inbeamresultfile in inbeamresultfiles:
                    checkFileExists(inbeamresultfile, False)
                    inbeamresult = astroobsresult.AstroObsResult(primaryinbeam, orgobscode, open(inbeamresultfile).readlines(), False)
                    inbeamsnr += inbeamresult.snr**2
                inbeamsnr = math.sqrt(inbeamsnr)
                pulsarresult = astroobsresult.AstroObsResult(pulsar, orgobscode,  open(statsfile).readlines(), False)
                #syserror = (i*deprojectedseparation + j*100.0/inbeamsnr)/1e6 # arcseconds
                ## syserror is at this stage along the major axis of the beam.  Scale by the projected beam onto RA and dec axes
                #epochscales[inbeamresult.mjd] = [syserror, syserror] #RA, dec; in arcseconds still (conversion to seconds for RA takes place below)
                #allerrors.append(syserror)
                # Now using "fractions of a beam" rather than "arcseconds"
                raprojection = 1.0/math.sqrt((math.sin(pulsarresult.beamparad)*math.sin(pulsarresult.beamparad) / (pulsarresult.beammaj*pulsarresult.beammaj)) + (math.cos(pulsarresult.beamparad)*math.cos(pulsarresult.beamparad) / (pulsarresult.beammin*pulsarresult.beammin))) # in mas
                decprojection =  1.0/math.sqrt((math.sin(pulsarresult.beamparad)*math.sin(pulsarresult.beamparad) / (pulsarresult.beammin*pulsarresult.beammin)) + (math.cos(pulsarresult.beamparad)*math.cos(pulsarresult.beamparad) / (pulsarresult.beammaj*pulsarresult.beammaj))) # in mas
                #print raprojection, decprojection
                beamfraction = i*deprojectedseparation + j*100.0/inbeamsnr
                epochscales[inbeamresult.mjd] = [beamfraction * raprojection / 1000.0, beamfraction * decprojection / 1000.0]  #RA, dec; in arcseconds still (conversion to seconds for RA takes place below)
                allerrors.append(epochscales[inbeamresult.mjd])

            # Add the systematic error and write out the revised pmpar input file
            #medianerror = 1000.0*sorted(allerrors)[len(allerrors)/2] # A bit approximate, but near enough since there is little variation epoch to epoch
            #os.system("rm -f withsystematic/erroradded.txt")
            #os.system('echo "%.3f" >> withsystematic/erroradded.txt' % medianerror)
            if not os.path.exists("withsystematic"):
                os.mkdir("withsystematic")
            pmparinput.setSystematicError(1.0/(15.0*math.cos(pulsarfit.decrad)), 1.0, epochscales)
            pmparinput.writePmparFile("withsystematic/local.jmfit.pmpar.in", False, True)
            os.system("rm -f withsystematic/newerroradded.txt")
            for a in allerrors:
                os.system('echo "%.6f %.6f\n" >> withsystematic/newerroradded.txt' % (a[0], a[1]))

            # Run pmpar and get the reduced chi squared
            os.chdir("withsystematic")
            os.system("rm -f junk.pmparout")
            os.system("pmpar local.jmfit.pmpar.in > junk.pmparout")
            os.system("rm -f scaling.txt")
            os.system('echo "%.1f %.1f" > scaling.txt' % (i, j))
            fitresult = astroobsresult.AstroFitResult("junk.pmparout")
            rchisqlist.append(fitresult.rchisq)
            errorsum += math.exp(2*math.fabs(math.log(fitresult.rchisq)))

        # Print the result and plot the histogram
        os.chdir("../../../")
        print i, j, ': ', errorsum/numfits
        if errorsum/numfits < besterrorsum:
            bestrchisqlist = list(rchisqlist)
            besterrorsum = errorsum/numfits
            besti = i
            bestj = j
        plotRChiSqHistogram(rchisqlist, "/home/deller/data/vlbi/psrpi/final/systematicshistograms/rchisqhistogram.%d-%d.png" % (int(i*10000), int(j*10000)))

print "The overall best result was obtained with i=%.5f, j=%.5f" % (besti,bestj)
print "The value achieved was " + str(besterrorsum)
print "Originally, the rchisq 25, 50, 75%% values were %.2f, %.2f, %.2f"  % (sorted(orgrchisqlist)[(25*len(orgrchisqlist))/100], sorted(orgrchisqlist)[(50*len(orgrchisqlist))/100], sorted(orgrchisqlist)[(75*len(orgrchisqlist))/100])
print "Now, the rchisq 25, 50, 75%% values are %.2f, %.2f, %.2f"  % (sorted(bestrchisqlist)[(25*len(bestrchisqlist))/100], sorted(bestrchisqlist)[(50*len(bestrchisqlist))/100], sorted(bestrchisqlist)[(75*len(bestrchisqlist))/100])
deprojout.close()
