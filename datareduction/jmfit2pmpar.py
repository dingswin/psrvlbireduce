#!/usr/bin/env python2

import sys, glob, re, math, os
from optparse import OptionParser

# function which parses one stats file
def processFile(solout, statsfile, vexfile,dosinglefit, dosinglerrll, experiment, expcode=""):
    statsin = open(statsfile, 'r')
    statslines = statsin.readlines()
    statsin.close()
    runto = len(statslines)/statsinc - 1
    count = 0
    expsnr = 0.0
    # Load up the vex file and figure out a centroid
    # DO IT!
    date = float((statslines[0].split()[3])[:-1])
    try:
        date += dayoffsets[expcode]
    except KeyError:
        if not expcode == "":
            print "No dayoffset for expcode " + expcode + " of experiment " + experiment + " - aborting"
            sys.exit()
#    if (dosinglefit or (experiment == 'bd141' and (expcode == 'b' or expcode == 'g'))) and not dosinglerrll:
    if dosinglefit or (experiment == 'bd141' and (expcode == 'g' or (expcode == 'b' and not dosinglerrll))):
        print "Taking the single summed value for experiment bd141" + expcode + "!!"
        count = runto*statsinc
        runto = 1
    elif dosinglerrll:
        print "Taking the single RR and single LL vals"
        count = (runto - 2)*statsinc
        runto = 2
    elif runto == 6 or runto == 10:
        print "Looks like a new style file with bonus allRR and allLL - skipping these"
        runto -= 2
    for i in range(runto):
        snr    = float(statslines[count + 4][21:-1])
        ra     = statslines[count + 7][21:-1]
        dec    = statslines[count + 8][21:-1]
        raerr  = float(statslines[count + 11][21:-1])/1000.0
        decerr = float(statslines[count + 12][21:-1])/1000.0
        if snr > 0.0:
            expsnr = expsnr + snr/runto
        if (snr > snrcut):
            solout.write("%10.5f %s %010.8f %s %010.8f\n" % (date, ra, raerr, dec, decerr))
        else:
            print("Skipping solution for " + statslines[count][:-1] + " since snr was " + \
                  str(snr) + " and cutoff was " + str(snrcut) + "\n")
            solout.write("#SNR of the following was only " + str(snr) + "\n")
            solout.write("#%10.5f %s %010.8f %s %010.8f\n" % (date, ra, raerr, dec, decerr))
        count = count + statsinc
    print "Average SNR for " + experiment + expcode + " was " + str(expsnr)

# GLOBAL VARS
statsinc = 13 #Number of lines per psrstats entry
bd141dayoffsets = {'a': 0.5, 'b': 0.333, 'c': 0.083, 'd': 1.01, 'e': 0.57, 'f': 0.52, 'g': 0.09, 'h': 1.055, 'i': 0.5, 'j': 0.54}
v190dayoffsets = {'a': 0.916667, 'e': 0.375, 'g': 0.041667, 'k': 0.375, 'm': 0.2916667, 'n': 0.08333, 'o': 0.5417}
defaultdayoffsets = {'a': 0.5, 'b': 0.5, 'c': 0.5, 'd': 0.5, 'e': 0.5, 'f': 0.5, 'g': 0.5, 'h': 0.5, 'i': 0.5, 'j': 0.5, 'k': 0.5, 'l': 0.5, 'm': 0.5}

# Set up parser
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-p", "--pulsardetails", dest="pulsardetails", 
                  default="1023+0038,1023+0024,14.329362",
                  help="Pulsar name[,ref name,DM]")
parser.add_option("-b", "--basedir", dest="basedir", 
                  default="/export/home/marathon2/data/1023+0038/",
                  help="Base directory")
parser.add_option("-s","--suffix", dest="suffix", default="",
                  help="Suffix to append to all filenames")
parser.add_option("-e","--experiment", dest="experiment", default="bd141",
                  help="Root experiment code (ie no a,b,c etc)")
parser.add_option("-c","--codes", dest="codes", default="a,b,c,d,e,f,g,h,i,j",
                  help="Comma-separated list of exp. subcodes eg a,b,c")
parser.add_option("--epoch", dest="epoch", default="55000",
                  help="Epoch for proper motion")
parser.add_option("--snrcut", dest="snrcut", default="4",
                  help="SNR cutoff for transferring solutions")
parser.add_option("--stokesi", dest="stokesi", default=False, 
                  action="store_true",
                  help="Use fits to stokes I (c.f. RR, LL) for position")
parser.add_option("--dosinglefit", dest="dosinglefit", default=False,
                  action="store_true", help="Take the single i,1,2,3,4 fit")
parser.add_option("--dosinglerrll", dest="dosinglerrll", default=False,
                  action="store_true", help="Take the single rr1234,ll1234 fit")
parser.add_option("--posttecormode", dest="posttecormode", default="",
                  help="For results using post-cal tecor specify eg jplg.40")
parser.add_option("--dodir", dest="dodir", default=False, action="store_true",
                  help="Just take all stats files in the cwd as the list")

# Parse command line options
(options, junk) = parser.parse_args()
details         = options.pulsardetails.split(',')
pulsar          = details[0]
ref             = ""
dm              = ""
snrcut          = float(options.snrcut)
epoch           = options.epoch
basedir         = options.basedir
suffix          = options.suffix
experiment      = options.experiment
expcodes        = options.codes.split(',')
stokesi         = options.stokesi
dosinglefit     = options.dosinglefit
dosinglerrll    = options.dosinglerrll
posttecormode   = options.posttecormode
dodir           = options.dodir
doposttecor = False
if not posttecormode == "":
    doposttecor = True
    #suffix = posttecormode + suffix
if dodir:
    pulsarfile = os.getcwd() + '/local.jmfit.pmpar.in'
else:
    pulsarfile = basedir + '/solutions/' + pulsar + ".jmfit.pmpar.in" + suffix
if len(details) > 1:
    ref = details[1]
if len(details) > 2:
    dm = details [2]
if experiment == 'bd141':
    dayoffsets = bd141dayoffsets
elif experiment == 'v190':
    dayoffsets = v190dayoffsets
else:
    print "Setting day offsets to default!"
    dayoffsets = defaultdayoffsets

# Check command line options
if len(junk) > 0:
    print "You can't supply any command line arguments except the - " + \
          "and -- ones! Aborting"
    sys.exit()

# Write the header of the pmpar file
solout = open(pulsarfile, 'w')
solout.write("name = J" + pulsar + "\n\n")
if not ref == "":
    solout.write("ref = J" + ref + "\n\n")
if not dm == "":
    solout.write("dm = " + str(dm) + "\n\n")
solout.write("epoch = " + epoch + "\n\n")
solout.write("# Positions\n\n")

statslist = []
if dodir:
    filelist = os.listdir(os.getcwd())
    for f in filelist:
        if len(f) > 5 and f[-6:] == ".stats":
            statslist.append(f)
    for s in statslist:
        print "Processing " + s
        print "dosinglerrll is " + str(dosinglerrll)
        processFile(solout, s, dosinglefit, dosinglerrll, dayoffsets, experiment)
    solout.close()
    sollines = open(pulsarfile).readlines()
    solout = open(pulsarfile, 'w')
    linelist = []
    for line in sollines:
        splitline = line.split()
        if len(splitline) == 5 and not line[0] == "#" and float(splitline[0]) > 40000 and float(splitline[0]) < 70000:
            linelist.append(line)
        else:
            solout.write(line)
    #linelist.sort()
    sortedlinelist = sorted(linelist, key=lambda line: line.split()[0])
    for line in sortedlinelist:
        solout.write(line)
    solout.close()
    sys.exit()

for expcode in expcodes:
    oldstokesi = stokesi
    if experiment + expcode == 'v190k':
        stokesi = False
    if experiment + expcode == 'v190m' or \
       experiment + expcode == 'v190n' or \
       experiment + expcode == 'v190o':
       pulsar = 'J0437-4715'
    expdir    = basedir + '/' + experiment + expcode + '/'
    if doposttecor:
        if not experiment == "v190" and not "bb269" in experiment:
            if stokesi:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            "_pulsar.difmap.jmfit." + posttecormode + \
                            ".0.stokesi.stats" + suffix
            else:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            "_pulsar.difmap.jmfit." + posttecormode + \
                            ".0.stats" + suffix
        else:
            if stokesi:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            ("_%s.difmap.jmfit."%(pulsar)) + posttecormode + \
                            ".0.stokesi.stats" + suffix
            else:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            ("_%s.difmap.jmfit."%(pulsar)) + posttecormode + \
                            ".0.stats" + suffix
    else:
        if not experiment == "v190" and not "bb269" in experiment:
            if stokesi:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            "_pulsar.difmap.jmfit.stokesi.stats" + suffix
            else:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            "_pulsar.difmap.jmfit.stats" + suffix
        else:
            if stokesi:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            ("_%s.gated.difmap.jmfit.stokesi.stats"%(pulsar)) + suffix
            else:
                statsfile = expdir + experiment.upper() + expcode.upper() + \
                            ("_%s.gated.difmap.jmfit.stats"%(pulsar)) + suffix
    processFile(solout, statsfile, dosinglefit, dosinglerrll, dayoffsets, experiment, expcode)
    stokesi = oldstokesi
solout.close()
