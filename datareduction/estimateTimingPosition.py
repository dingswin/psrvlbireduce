#!/usr/bin/env python
#####################################################################################
### functions:
### read timing fit, gaussianly generate astrometric parameters, and pmparfit the
### position at targeted time
### usage: estimateTimingPosition.py -f timingfit_inputfile -t epoch
### all timingfit_inputfile under /fred/oz002/hding/mspsrpi/processing/J1012+5307/pmparesults/timing_results/
#####################################################################################
import os,glob,sys,yaml,math,pickle
import howfun,mspsrpifun
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.time import Time
from optparse import OptionParser
np.set_printoptions(suppress=True)

usage = "usage: %prog []\n-f --inputfile\n-t --targetepoch\n\n-n --noplotbootstrap\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-t", "--targetepoch", dest="targetepoch", default="",
                   help="epoch in MJD of VLBI epoch to transform to")
parser.add_option("-n", "--noplotbootstrap", dest="noplotbootstrap", default=False,
                   action="store_true",help="don't plot and save bootstrap samples")
parser.add_option("-f", "--inputfile", dest="inputfile", default="", 
                   help="timingfit file as input")
(options, junk) = parser.parse_args()
targetepoch      = options.targetepoch
noplotbootstrap   = options.noplotbootstrap
inputfile      = options.inputfile

if targetepoch=='' or inputfile=='':
    print usage
    sys.exit()
if not os.path.exists(inputfile):
    print "input file not found; aborting"
    sys.exit()
## read inputfile #######################################################
estimates = np.array([])
errors = np.array([])
lines = open(inputfile).readlines()
for line in lines:
    if 'epoch' in line:
        epoch = line.split('=')[-1].strip()
    for estimate in ['RA', 'Dec']:
        if estimate in line:    
            line1 = line.split('=')[-1].strip()
            exec("%s = line1.split('+-')[0].strip()" % estimate)
            exec("err_%s = float(line1.split('+-')[1].strip())/3600" % estimate)
            exec("estimates = np.append(estimates, howfun.dms2deg(%s))" % estimate)
            exec("errors = np.append(errors, err_%s)" % estimate)
    for estimate in ['pi', 'mu_a', 'mu_d']:
        if estimate in line:
            line1 = line.split('=')[-1].strip()
            exec("%s = line1.split('+-')[0].strip()" % estimate)
            exec("err_%s = float(line1.split('+-')[1].strip().split(' ')[0])" % estimate)
            exec("estimates = np.append(estimates, float(%s))" % estimate)
            exec("errors = np.append(errors, err_%s)" % estimate)
print estimates,errors
## get predicted position ##################################################
tmpfile = '.pmpar.tmp'
predicted_position = '.timing_predicted_position.tmp'
writefile = open(tmpfile, 'w')
writefile.write('epoch = %s\n' % epoch)
writefile.write('RA = %s\n' % RA)
writefile.write('Dec = %s\n' % Dec)
writefile.write('pi = %s\n' % pi)
writefile.write('mu_a = %s\n' % mu_a)
writefile.write('mu_d = %s\n' % mu_d)
writefile.close()
os.system("pmpar %s -z %s > %s" % (tmpfile, targetepoch, predicted_position))
lines = open(predicted_position).readlines()
for line in lines:
    RA0 = line.split('  ')[1].strip()
    Dec0 = line.split('  ')[2].strip()
    RA0 = howfun.dms2deg(RA0)
    Dec0 = howfun.dms2deg(Dec0)
print RA0, Dec0, "\n"

## generate random astrometric parameters and predict position at targetepoch with pmpar #########
bootstrapruns = 10000
count = 0
bs_pmparinput = '.boostrap.pmpar.in'
bs_positions = 'bs_timing_positions'
if os.path.exists(bs_positions):
    os.remove(bs_positions)
while count < bootstrapruns:
    # generate random parameters
    [RA_bs, Dec_bs, pi_bs, mu_a_bs, mu_d_bs] = np.random.normal(estimates, errors)
    RA_bs = howfun.deg2dms(RA_bs)
    Dec_bs = howfun.deg2dms(Dec_bs)
    # write pmpar.in
    writefile = open(bs_pmparinput, 'w')
    writefile.write('epoch = %s\n' % epoch)
    writefile.write('RA = %s\n' % RA_bs)
    writefile.write('Dec = %s\n' % Dec_bs)
    writefile.write('pi = %f\n' % pi_bs)
    writefile.write('mu_a = %f\n' % mu_a_bs)
    writefile.write('mu_d = %f\n' % mu_d_bs)
    writefile.close()
    # predict positions with pmpar
    os.system('pmpar %s -z %s >> %s' % (bs_pmparinput, targetepoch, bs_positions))
    print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / bootstrapruns)) + " \r") 
    count += 1
## parse bs_positions ###########################################################################
RAs = np.array([])
Decs = np.array([])
lines = open(bs_positions).readlines()
for line in lines:
    line1 = line.split('  ')
    RAs = np.append(RAs, howfun.dms2deg(line1[1].strip()))
    Decs = np.append(Decs, howfun.dms2deg(line1[2].strip()))
RAs.sort()
Decs.sort()
samplevolume = len(RAs)
print samplevolume
# estimate uncertainty and write out
fileWrite = open(predicted_position, 'w')
for confidencelevel in [0.9973, 0.6827]: # 1sigma, 3sigma->0.9973
    [valueRA1, errorRA1] = howfun.sample2estimate(RAs, confidencelevel)
    [valueDec1, errorDec1] = howfun.sample2estimate(Decs, confidencelevel)
    print("confidencelevel = %f" % confidencelevel)
    print("RA = %s +- %f" % (howfun.deg2dms(valueRA1), errorRA1*3600))
    print("Dec = %s +- %f" % (howfun.deg2dms(valueDec1), errorDec1*3600))
    fileWrite.write("confidencelevel = %f:\n" % confidencelevel)
    fileWrite.write("epoch = %s\n" % targetepoch)
    fileWrite.write("RA = %s +- %f\n" % (howfun.deg2dms(valueRA1), errorRA1*3600))
    fileWrite.write("Dec = %s +- %f\n" % (howfun.deg2dms(valueDec1), errorDec1*3600))
fileWrite.close()
## plot and save probability density functions for RAs and Decs
if noplotbootstrap:
    sys.exit()
print "\nmaking bootstrap plots...\n"
for parameter in [RAs, Decs]:
    plt.hist(parameter, 1000, density=True, facecolor='g', alpha=0.75)
    minimum = min(parameter)
    maximum = max(parameter)
    x = np.arange(minimum,maximum,(maximum-minimum)/1000)
    if (parameter == RAs).any():
        y = norm.pdf(x, valueRA1, errorRA1)
    if (parameter == Decs).any():
        y = norm.pdf(x, valueDec1, errorDec1)
    plt.plot(x,y)
    plt.ylabel('probability density')
    if (parameter == RAs).any():
        plt.xlabel('Right Ascension (hour)')
        plt.savefig('timing_predicted_RA_pdf.eps')
    if (parameter == Decs).any():
        plt.xlabel('Declination (degree)')
        plt.savefig('timing_predicted_Dec_pdf.eps')
    plt.clf()
print "\nwe successfully fly over\n"
