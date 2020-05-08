#!/usr/bin/env python
################################################################################
## converted from pylabplot_pmpar.py by Adam Deller, 26 Aug 2011
## Hao Ding, May 2020
## The plot in the XTEJ1810-197 paper is made with the following command:
## plotpmparXTEJ1810.py -f XTEJ1810-197.pmpar.in.use.A1 -s XTEJ1810-197.pmpar.in.use.A1.non_phscal -r [2,-1] -b 500
## 
################################################################################
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
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
    [start_index, end_index] = line_index_range
    if start_index < 0:
        start_index += nepoch
    if end_index < 0:
        end_index += nepoch + 1
    #print type(line_index_range[0]), type(line_index_range[1])
    for i in range(start_index, end_index):
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
    return positions, epoch
################################################################################
## Main code
################################################################################
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--filename", dest="filename", default="",
                  help="The pmpar input file")
parser.add_option("-s", "--secondfile", dest="filename1", default="",
                  help="The second pmpar input file for the data points of non-phscal setup")
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
(options, junk) = parser.parse_args()
pmparexec       = options.pmparexec
filename        = options.filename
filename1       = options.filename1
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
fig = plt.figure()
gs = gridspec.GridSpec(6,14)
fig1 = fig.add_subplot(gs[:6, :4])
ttimes, tras, tdecs, etimes, eras, edecs, eraerr, edecerr, pras, pdecs = pmpar(pmparexec, filename, position_lines, [0,-1], True)
fig1.plot(tras,tdecs, color='k', alpha=0.3)
fig1.errorbar(eras, edecs, xerr=eraerr, yerr=edecerr, fmt='.', markersize=3)
fig1.set_xlabel('relative RA. (mas)')
fig1.set_ylabel('relative Decl. (mas)')
fig1.tick_params(axis='x', direction='in', bottom=True)
#fig1.set_title('Sky position evolution')
#ax = fig1.gca()
fig1.invert_xaxis()
#ax1.axvline(x=most_probable_PI, c='black', linestyle='--', linewidth=0.5)
ttimes, tras, tdecs, etimes,  eras,  edecs,  eraerr,  edecerr,  pras,  pdecs = pmpar(pmparexec, filename, position_lines, range2plot, False)
if filename1 != '':
    [position_line1s, junk] = pmparin2position_lines(filename1)
    junk1, junk2, junk3, etimes1, eras1, edecs1, eraerr1, edecerr1, pras1, pdecs1 = pmpar(pmparexec, filename1, position_line1s, range2plot, False)

pmparout = target + '.pmpar.out'
os.system("%s %s > %s" % (pmparexec, filename, pmparout))
[RA0, Dec0, junk1, junk2, mu_a0, mu_d0, junk3, junk4, junk5, junk6, junk7, junk8] = mspsrpifun.readpmparout(pmparout)
parameters = [RA0, Dec0, mu_a0, mu_d0]
#####################################################################################
## bootstrap begins here
####################################################################################
## PLOT STACKED BOOTSTRAPPING RA/DEC-TIME (PM REMOVED) PLOTS
print "\nrunning bootstrap and making stacked RA/Dec-time (PM removed) plots...\n"
# get the measurements first (without pm and ref position)
count = 0
#bootstrapruns = 2000
pulsitions = target + '.pmpar.in.bootstrap'
if not target=='':
    plotfile = 'bootstrap_nopm_%s_%druns_plus_radec.%s' % (target, bootstrapruns, plottype)
else:
    plotfile = 'bootstrap_nopm_%druns_plus_radec.%s' % (bootstrapruns, plottype)
#fig = plt.figure(1)
fig2 = fig.add_subplot(gs[:3, 4:14])
fig3 = fig.add_subplot(gs[3:6, 4:14])
#fig3 = fig.add_subplot(gs[3:6, 4:14], sharex=fig2)
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
    fig2.plot(ttimes, tras, alpha=0.01, color='k', linewidth=0.1)
    fig3.plot(ttimes, tdecs, alpha=0.01, color='k', linewidth=0.1)
    print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / bootstrapruns)) + " \r")
    count += 1
fig2.errorbar(etimes, eras, yerr=eraerr, fmt='.', markersize=5, capsize=3) # use the inputs from pre-bootstrap
fig2.errorbar(etimes1, eras1, yerr=eraerr1, fmt='.', markersize=5, capsize=4, color='y', alpha=0.4) # use the inputs from pre-bootstrap
errorbar0 = fig3.errorbar(etimes, edecs, yerr=edecerr, fmt='.', markersize=5, capsize=3) 
errorbar1 = fig3.errorbar(etimes1, edecs1, yerr=edecerr1, fmt='.', markersize=5, capsize=4, color='y', alpha=0.4) 
fig3.legend([errorbar0, errorbar1], ['virtual-cal frame', 'J1819 frame'], loc='lower right')
fig2.set_ylabel('RA. offset (mas)')
fig2.yaxis.set_label_position('right')
fig2.yaxis.tick_right()
plt.setp(fig2.get_xticklabels(), visible=False)
#fig2.set_xticklabels([])
fig2.tick_params(axis='x', direction='in', bottom=True)
fig3.tick_params(axis='x', direction='in', bottom=True, top=True)
fig3.set_ylabel('Decl. offset (mas)')
fig3.yaxis.set_label_position('right')
fig3.yaxis.tick_right()
fig3.set_xlabel('time (MJD)')
fig2.set_ylim([-1.2,1.2])
#fig3.set_ylim([-1,1])
nbins = len(fig2.get_xticklabels())
fig2.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='lower'))
nbins = len(fig3.get_xticklabels())
fig3.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
plt.subplots_adjust(hspace=0, wspace=3)
#gs.tight_layout(fig)
#fig.tight_layout()
plt.savefig(plotfile)
