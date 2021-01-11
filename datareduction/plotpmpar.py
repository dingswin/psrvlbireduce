#!/usr/bin/env python
################################################################################
## converted from pylabplot_pmpar.py by Adam Deller, 26 Aug 2011
## Hao Ding, Mar 2019
################################################################################
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import mspsrpifun, howfun
from optparse import OptionParser

# Run pmpar, optionally with proper motion subtracted, get results
# fit parameters with random extracted positions, while plot with consistent pm and refpos
def pmpar(pmparexec, filename, position_lines, line_index_range, nopmsubtract=False, nopi=False): 
    #first pmpar run to prepare
    nepoch = len(position_lines)
    tmpfile = '.plotpmpar.tmp'
    os.system("%s %s > %s" % (pmparexec, filename, tmpfile))
    [RA, Dec, epoch, pi, mu_a, mu_d, error_pi, error_mu_a, error_mu_d, l, b, rchsq] = mspsrpifun.readpmparout(tmpfile)
    RA = howfun.deg2dms(RA)
    Dec = howfun.deg2dms(Dec)
    writefile = open(tmpfile, 'w')
    writefile.write('epoch = %s\n' % epoch)
    writefile.write('RA = %s\n' % RA)
    writefile.write('Dec = %s\n' % Dec)
    writefile.write('mu_a = %s\n' % mu_a)
    writefile.write('mu_d = %s\n' % mu_d)
    if nopi:
        writefile.write('pi = 0\n')
    else:
        writefile.write('pi = %s\n' % pi)
    if line_index_range[0] < 0:
        line_index_range[0] += nepoch
    if line_index_range[1] < 0:
        line_index_range[1] += nepoch + 1
    #print type(line_index_range[0]), type(line_index_range[1])
    for i in range(line_index_range[0], line_index_range[1]):
        writefile.write(position_lines[i] + '\n')
    writefile.close()
    filename = tmpfile
    # second pmpar run to plot
    runline = "%s %s > /dev/null" % (pmparexec, filename)
    if not nopmsubtract:
        runline += " -om"
    os.system(runline)
    tin = open("pmpar_t")
    tlines = tin.readlines()
    tin.close()
    ein = open("pmpar_e")
    elines = ein.readlines()
    ein.close()
    ttimes = np.array([])
    tras   = np.array([])  
    tdecs  = np.array([])
    etimes = np.array([])
    eras   = np.array([])
    edecs  = np.array([])
    eraerr = np.array([])
    edecerr= np.array([])
    pras   = np.array([])
    pdecs  = np.array([])
    for line in tlines:
        splitline = line.split()
        ttimes = np.append(ttimes, float(splitline[0]))
        tras   = np.append(tras, float(splitline[1]))
        tdecs  = np.append(tdecs, float(splitline[2]))
    for line in elines:
        splitline = line.split()
        etimes = np.append(etimes, float(splitline[0]))
        eras   = np.append(eras, float(splitline[1]))
        edecs  = np.append(edecs, float(splitline[3]))
        eraerr = np.append(eraerr, float(splitline[2]))
        edecerr= np.append(edecerr, float(splitline[4]))
        pras   = np.append(pras, float(splitline[5]))
        pdecs  = np.append(pdecs, float(splitline[6]))
    return ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs

def pmparin2position_lines(pmparin):
    positions = np.array([])
    lines = open(pmparin).readlines()
    for line in lines:
        if 'epoch' in line:
            epoch = line.split('=')[1].strip()
        if line.count(':') == 4 and (not line.strip().startswith('#')):   
            positions = np.append(positions, line)
            positions.sort()
    return positions, epoch
################################################################################
## Main code
################################################################################
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--filename", dest="filename", default="",
                  help="The pmpar input file")
parser.add_option("--nopmsubtract", dest="nopmsubtract", default=False,
                  action="store_true", help="Don't subtract PM for " + \
                  "vs. time plots")
parser.add_option("--plottype", dest="plottype", default="pdf",
                  help="The file type for the plot (def. png)")
parser.add_option("--pmparexec", dest="pmparexec", default="pmpar",
                  help="Path to the pmpar executable")
parser.add_option("-b", "--bootstrapruns", dest="bootstrapruns",default=-1,
                  help="Plot the contents of this bootstrap file; if bootstrapruns==-1, no bootstrapplot is made")
parser.add_option("-r", "--rangeofepoch", dest="rangeofepoch", default=[0,-1],
                  help="set a epoch range to plot, whereas the fitting would still go through all the input epochs")
#parser.add_option("-t", "--addtext", dest="addtext", default=False,
#                  action="store_true", help="provide MJD notation")
(options, junk) = parser.parse_args()
pmparexec       = options.pmparexec
filename        = options.filename
target          = filename.split('.')[0].strip()
plottype        = options.plottype
nopmsubtract    = options.nopmsubtract
bootstrapruns   = int(options.bootstrapruns)
range2plot      = options.rangeofepoch
exec('range2plot = %s' % range2plot)
if len(junk) > 0:
    parser.error("You can only supply the allowed options")
if filename=="":
    parser.error("You must supply a filename to run on with -f or --filename")
## PARSE .PMPAR.IN FILE
[position_lines, epoch] = pmparin2position_lines(filename)
nepoch = len(position_lines)
### FIRST DO RA VS DEC
for pm in [True, False]:
    ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs = pmpar(pmparexec, filename, position_lines, range2plot, pm)
    plt.plot(tras,tdecs)
    plt.errorbar(eras, edecs, xerr=eraerr, yerr=edecerr, fmt='.', markersize=3)
    plt.xlabel('relative RA. (mas)')
    plt.ylabel('relative Decl. (mas)')
    plt.title('Sky position evolution')
    ax = plt.gca()
    ax.invert_xaxis()
    if pm: 
        plotfile = "radec" + target + "." + plottype
    else:
        plotfile = "radec_nopm" + target + "." + plottype
    plt.savefig(plotfile)
    plt.clf()
pmparout = target + '.pmpar.out'
os.system("%s %s > %s" % (pmparexec, filename, pmparout))
[RA0, Dec0, junk1, junk2, mu_a0, mu_d0, junk3, junk4, junk5, junk6, junk7, junk8] = mspsrpifun.readpmparout(pmparout)
parameters = [RA0, Dec0, mu_a0, mu_d0]
### NOW DO THE RA VS TIME WITH PM SUBTRACTED
if not nopmsubtract:
    ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs = pmpar(pmparexec, filename, position_lines, range2plot)
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

ax1.plot(ttimes, tras)
ax1.set_xlabel('time (MJD)')
ax1.set_ylabel('RA. offset (mas)')
ax1.set_title('RA-time (proper motion removed)')
ax1.errorbar(etimes, eras, yerr=eraerr, fmt='.', markersize=5, capsize=3)
#plt.savefig(plotfile_RAtime)
#plt.clf()

### NOW DO THE DEC VS TIME WITH PM SUBTRACTED
if nopmsubtract:
    plotfile_ra_dec_time = "ra_dec_time." + plottype
    if not target == "":
        plotfile_ra_dec_time = "ra_dec_time." + target + "." + plottype
else:
    plotfile_ra_dec_time = "ra_dec_time_nopm." + plottype
    if not target == "":
        plotfile_ra_dec_time = "ra_dec_time_nopm." + target + "." + plottype
ax2.plot(ttimes, tdecs)
ax2.set_xlabel('time (MJD)')
ax2.set_ylabel('Decl. offset (mas)')
ax2.set_title('Dec-time (proper motion removed)')
ax2.errorbar(etimes, edecs, yerr=edecerr, fmt='.', markersize=5, capsize=3)
plt.tight_layout()
plt.savefig(plotfile_ra_dec_time)
plt.clf()
#####################################################################################
## bootstrap begins here
####################################################################################
if bootstrapruns < 0:
    sys.exit()
## PLOT STACKED BOOTSTRAPPING RA/DEC-TIME (PM REMOVED) PLOTS
print "\nrunning bootstrap and making stacked RA/Dec-time (PM removed) plots...\n"
# get the measurements first (without pm and ref position)
count = 0
#bootstrapruns = 2000
pulsitions = target + '.pmpar.in.bootstrap'
if not target=='':
    plotfile = 'bootstrap_nopm_%s_%druns.%s' % (target, bootstrapruns, plottype)
else:
    plotfile = 'bootstrap_nopm_%druns.%s' % (bootstrapruns, plottype)
fig = plt.figure(1)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ttimes, tras2, tdecs2, junk1, junk2, junk3, junk4, junk5, pras, pdecs = pmpar(pmparexec, filename, position_lines, range2plot, True, True)
while count < bootstrapruns:
    indices = np.random.randint(0, nepoch, nepoch)
    if len(np.unique(indices)) < 5:
        continue
    writefile = open(pulsitions, 'w')
    writefile.write('epoch = %s\n' % epoch)
    for i in indices:
        writefile.write(position_lines[i])
    writefile.close()
    ## PLOT NOW
    ttimes, tras1, tdecs1, junk1, junk2, junk3, junk4, junk5, pras, pdecs = pmpar(pmparexec, pulsitions, position_lines, range2plot, True, False)
    tras = tras1 - tras2 # bootstrap sky evolution - best fit model, both including pm
    tdecs = tdecs1 - tdecs2
    ax1.plot(ttimes, tras, alpha=0.01, color='k', linewidth=0.1)
    ax2.plot(ttimes, tdecs, alpha=0.01, color='k', linewidth=0.1)
    print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / bootstrapruns)) + " \r")
    count += 1
ax1.errorbar(etimes, eras, yerr=eraerr, fmt='.', markersize=5, capsize=3) # use the inputs from pre-bootstrap
ax2.errorbar(etimes, edecs, yerr=edecerr, fmt='.', markersize=5, capsize=3) 
ax1.set_ylabel('RA. offset (mas)')
ax2.set_ylabel('Decl. offset (mas)')
ax2.set_xlabel('time (MJD)')
#ax1.set_ylim([-1,1])
#ax2.set_ylim([-1,1])
plt.savefig(plotfile)
