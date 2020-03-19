#!/usr/bin/env ParselTongue
################################################################################
# General python imports
################################################################################
import sys, os, string, math, yaml, glob
import vlbatasks, astroobsresult, astro_utils, wmom
import matplotlib.pyplot as pyplot
from time import gmtime, strftime
from optparse import OptionParser
from matplotlib.patches import Ellipse

try:
    svndir = os.environ['PSRPISVNROOT']
except KeyError:
    print "PSRPISVNROOT is not defined - aborting!"
    sys.exit(1)

snrlimit = 10

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-p", "--pulsar", dest="pulsar", default="",
                  help="Experiment name")
(options, junk) = parser.parse_args()
pulsar          = options.pulsar
configdir       = svndir + '/configs/'

if pulsar == "":
    parser.error("You must supply a pulsar name")

targetconfig = yaml.load(open(configdir + pulsar + '.yaml'))
expconfigfiles = glob.glob(configdir + 'bd*yaml') + glob.glob(configdir + 'br*yaml')
usedexps = []
for e in expconfigfiles:
    clines = open(e).readlines()
    for line in clines:
        if pulsar in line:
            usedexps.append(e.split('/')[-1].split('.')[0])
            continue
if len(usedexps) == 0:
    print "Couldn't find any experiments matching pulsar " + pulsar
    sys.exit()
print usedexps

results = []
vexpositions = []
for exp in usedexps:
    expconfigfile = configdir + exp + '.yaml'
    expconfig     = yaml.load(open(expconfigfile))
    expdir        = expconfig['rootdir'] + '/' + exp + '/'
    print expdir
    sourcefile = expdir + exp + '.source'
    if not os.path.exists(sourcefile):
        print "Skipping ", exp, "as ", sourcefile, "doesn't exist"
        continue
    sourcefilelines = open(sourcefile).readlines()
    phsrefname = ""
    for j in range(len(sourcefilelines)):
        if pulsar in sourcefilelines[j]:
            phsrefname = sourcefilelines[j+1].split(':')[-1].strip()
    if phsrefname == "":
        print "Couldn't find the phase reference source for " + pulsar
        sys.exit()
    vexphsrefname = phsrefname
    if "IBC" in phsrefname or "NJ" in phsrefname or "FJ" in phsrefname:
        vexphsrefname = pulsar + "PT"
    statsfile = expdir + exp + '_' + phsrefname + '_ibshiftdiv.jmfit.stats'
    vexfile = expdir + exp + ".vex"
    if not os.path.exists(statsfile) or not os.path.exists(vexfile):
        print "Skipping ", exp, statsfile, vexfile
        continue
    result = astroobsresult.ExtendedJmfitResult(pulsar, exp, statsfile)
    if result.snr < snrlimit:
        print "Skipping ", exp, statsfile, "due to SNR limit"
        continue
    vexlines = open(vexfile).readlines()
    vexstring = ""
    for i in range(len(vexlines)):
        line = vexlines[i]
        if "source_name" in line and vexphsrefname in line:
            j = 1
            while vexlines[i+j][0] == "#" or not "ra = " in vexlines[i+j] and i+j < len(vexlines):
                j += 1
            vexstring = vexlines[i+j]
        i += 1
    if vexstring == "":
        print "Couldn't find " + phsrefname + ' in vex file'
        sys.exit()
    vexra = vexstring.split()[2][:-1]
    vexdec = vexstring.split()[5][:-1]
    results.append(result)
    vexpositions.append([astro_utils.stringToRad(vexra, True), astro_utils.stringToRad(vexdec, False)])

if len(results) == 0:
    print "Couldn't find any matching results!"
    sys.exit()

meanvexra = 0
meanvexdec = 0
for v in vexpositions:
    print v
    meanvexra += v[0]/len(vexpositions)
    meanvexdec += v[1]/len(vexpositions)

raoffsets = []
decoffsets = []
raerrs = []
decerrs = []
raweights = []
decweights = []
rad2mas = 180*60*60*1000/math.pi
for r in results:
    raoff = (r.rarad - meanvexra)*rad2mas
    decoff = (r.decrad - meanvexdec)*rad2mas
    raoffsets.append(raoff)
    decoffsets.append(decoff)
    raerrs.append(r.raerrmas)
    decerrs.append(r.decerrmas)
    raweights.append(1.0/r.raerrmas**2)
    decweights.append(1.0/r.decerrmas**2)

wmeanra, werrra, wstdra = wmom.wmom(raoffsets, raweights, None, False, True)
wmeandec, werrdec, wstddec = wmom.wmom(decoffsets, decweights, None, False, True)

print wmeanra, werrra, wstdra
print wmeandec, werrdec, wstddec

fig = pyplot.figure()
ax = fig.add_subplot(111, aspect='equal')
l1 = pyplot.errorbar(raoffsets, decoffsets, raerrs, decerrs, linestyle='none', elinewidth=2, color='b')
l2 = pyplot.scatter([0], [0], s=40, c='r', marker='o')
l3 = pyplot.scatter([wmeanra], [wmeandec], s=40, c='blue', marker = 'x')
errorellipse = Ellipse(xy=(wmeanra, wmeandec), width=2*wstdra, height=2*wstddec, angle=0)
ax.add_artist(errorellipse)
errorellipse.set_alpha(0.25)
errorellipse.set_facecolor('blue')
pyplot.xlabel("Offset in R.A. (mas)")
pyplot.ylabel("Offset in Dec. (mas)")
pyplot.savefig(phsrefname + ".reverseinbeamcorrections.png")
print "Mean error in RA is " + str(wmeanra) + " +/- " + str(wstdra) + " mas"
print "Mean error in dec is " + str(wmeandec) + " +/- " + str(wstddec) + " mas"
