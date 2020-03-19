#!/usr/bin/env python
#####################################################################################
### functions:
### 1) compile statsfiles for target to  pmpar.in.preliminary, 
### 2) estimate systematics using empirical functions, then generate pmpar.in file
### 3) run bootstraps from pmpar.in, make PDFs, estimate uncertainties and make plots
### tailored for MSPSRPI
### composed by Hao Ding
#####################################################################################
import os,glob,sys,yaml,math,pickle
import howfun,mspsrpifun
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#plt.rcParams.update({'font.size': 22})
#from matplotlib import rc
#rc('text', usetex=True)
from scipy.stats import norm
from scipy import constants
from astropy.time import Time
from optparse import OptionParser
np.set_printoptions(suppress=True)

usage = "usage: %prog []\n-t --target\n-n --nosys\n-b --bootstrap\n-p --plotbootstrap\n-i --binno\n-m --montecarlo\n-d --diagnosebootstrap\n-s --estnd\n-o --bootstrapprior\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-t", "--target", dest="target", default="",
                   help="target on which positions are gathered into pmpar.in file")
parser.add_option("-n", "--nosys", dest="nosys", default=False,
                   action="store_true",help="no systematic error will be estimated")
parser.add_option("-b", "--bootstrap", dest="bootstrap", default=False,
                   action="store_true",help="run bootstrapping after estimating systematic errors")
parser.add_option("-p", "--plotbootstrap", dest="plotbootstrap", default=False,
                   action="store_true",help="plot and save bootstrap samples")
parser.add_option("-d", "--diagnosebootstrap", dest="diagnosebootstrap", default="",
                   help="diagnose bootstrap by removing one epoch from fitting")
parser.add_option("-i", "--binno", dest="binno", default=-1, help="for binned gated data, there are multiple gated*stats, choose which bin to compile")
parser.add_option("-m", "--montecarlo", dest="montecarlo", default=False,
                   action="store_true", help="use montecarlo instead of bootstrap for parameter estimation, this function still needs verification")
parser.add_option("-s", "--estnd", dest="estnd", default=False,
                   action="store_true", help="use eliminative monte carlo to estimate n-dimensional function of known-pdfs variables")
parser.add_option("-o", "--bootstrapprior", dest="bootstrapprior", default='', help="adopt timing or other prioris in bootstrap")
(options, junk) = parser.parse_args()
targetname      = options.target
nosys           = options.nosys
bootstrap       = options.bootstrap
plotbootstrap   = options.plotbootstrap
diagnosebootstrap = options.diagnosebootstrap
binno = int(options.binno)
montecarlo = options.montecarlo
estnd = options.estnd
if bootstrap == True:
    bootstrapprior = options.bootstrapprior

auxdir    = os.environ['PSRVLBAUXDIR']
configdir = auxdir + '/configs/'
targetdir = auxdir + '/processing/' + targetname
if not os.path.exists(targetdir):
    print("%s doesn't exist; aborting\n" % targetdir)
    sys.exit()
statsfiles=glob.glob(r'%s/*/*.gated.difmap.jmfit.stokesi.stats' % targetdir) #find statsfile for each epoch
if binno != -1:
    statsfiles=glob.glob(r'%s/*/bin_astrometric_run_results/bin%d/*stats' % (targetdir, binno))
statsfiles.sort() #sort statsfile
## find pulsar name automatically and name the position file    
mjds       = np.array([])
decyears   = np.array([])
error0RAs  = np.array([]) #error0 -> formal errors; error1 ->systematic errors; error -> total errors
error0Decs = np.array([]) 
SNprIBCs   = np.array([]) #SNprIBC -> S/N of primary IBC
beamLAs    = np.array([]) #beamLA -> beam Long axis length; beamSA -> beam short axis length;
beamSAs    = np.array([]) 
beamPAs    = np.array([]) #beamPA -> beam position angle;
RAs        = np.array([])
Decs       = np.array([])
epoch      = 57700

pmparesultsdir = targetdir + '/pmparesults/'
if not os.path.exists(pmparesultsdir):
    os.system('mkdir %s' % pmparesultsdir)
pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in.preliminary'
fileWrite = open(pulsitions, 'w')
fileWrite.write("### pmpar format\n")
fileWrite.write("#name = " + targetname + "\n")
fileWrite.write("#ref = " + targetname + "\n")
fileWrite.write("epoch = %f\n" % epoch)
fileWrite.write("#pi = 0\n")
fileWrite.write("# decimalyear RA +/- Dec +/-\n")
#find position info from each statsfile and compile them into pulsitions
for statsfile in statsfiles:
    expno = statsfile.split('/')[-2]
    if binno != -1:
        expno = statsfile.split('/')[-4]
    print expno
    vexfile = targetdir + "/" + expno + "/" + expno + ".vex"
    if not os.path.exists(vexfile):
        print vexfile + " does not exist - aborting!"
        sys.exit()
    
    ## read obs_time from vex file
    vexlines = open(vexfile).readlines()
    for line in vexlines:
        if 'MJD' in line:
            MJD = int(line.split(':')[-1])
        searchContent = 'source=' + targetname            
        if searchContent in line:
            startime=line.split(';')[0]
            break
    
    # translate the obs_time to decimalyear format
    startime=startime.split('d')[-1]
    starthour=float(startime[0:2])    
    startmin=float(startime[3:5])
    print starthour,startmin
    dayfrac = (startmin/60+starthour)/24
    MJD = float(MJD) + dayfrac
    mjds = np.append(mjds,MJD)
    print mjds
    print "average MJD is: %f" % np.average(mjds)
    MJD = Time(MJD,format='mjd')
    decyear = format(MJD.decimalyear,'.4f')
    decyears = np.append(decyears, decyear)
    print decyears
    # dump decyears for plotscatter use
    tmpfile = targetdir + '/.variable_decyears.tmp'
    writefile = open(tmpfile, 'w')
    pickle.dump(decyears, writefile)
    writefile.close()

    ## read positions and errors from stats file
    statsline = open(statsfile).readlines()[-7:]
    for line in statsline:
        if 'Actual RA' in line:
            RA = line.split('RA:')[-1].strip()
            RAs = np.append(RAs, RA)
            print RA
        if 'Actual Dec' in line:
            Dec = line.split('Dec:')[-1].strip() 
            Decs = np.append(Decs, Dec)
            print Dec
        if 'RA error (hms)' in line:
            error0RA = float(line.split('(hms):')[-1].strip())*0.001 #unit: s
            error0RAs = np.append(error0RAs, error0RA)
            print error0RAs
        if 'Dec error (mas)' in line:
            error0Dec = float(line.split('(mas):')[-1].strip())*0.001 #unit: "
            error0Decs = np.append(error0Decs, error0Dec)
            print error0Decs
    fileWrite.write("%s %s %.7f %s %.6f\n" % (decyear, RA, error0RA, Dec, error0Dec))
fileWrite.close()
nepoch = len(RAs)
##### estimate systematic errors #####################################################################
#####################################################################################################
if nosys:
    print "\nsystematics NOT estimated\n"
    sys.exit()

print("\nestimating systematics...\n")
## look for primary IBC
expconfigfile = configdir + targetname + '.yaml'
if not os.path.exists(expconfigfile):
    parser.error("Experiment config file %s does not exist!" % expconfigfile)
expconfig = yaml.load(open(expconfigfile))
prIBCname = expconfig['primaryinbeam']
print prIBCname
## get primary IBC position, S/N and synthesized beam size from its stats file
prIBCstatsfiles = glob.glob(r'%s/*/*%s.difmap.jmfit.stokesi.stats' % (targetdir, prIBCname))
prIBCstatsfiles.sort()
print prIBCstatsfiles
for prIBCstatsfile in prIBCstatsfiles:
    lines = open(prIBCstatsfile).readlines()[-10:]
    for line in lines:
        if 'S/N' in line:
            SNprIBC = line.split(':')[-1].strip()
            SNprIBC = float(SNprIBC)
            SNprIBCs = np.append(SNprIBCs, SNprIBC)
            print SNprIBCs
        if 'Actual RA' in line:
            prIBCRA = line.split('RA:')[-1].strip() #just leave the prIBCRA/Dec for the last epoch, that is enough
            print prIBCRA
        if 'Actual Dec' in line:
            prIBCDec = line.split('Dec:')[-1].strip()
            print prIBCDec
        if 'beam' in line:
            line = line.split('beam')[-1].strip().split(' ')
            beamPA   = float(line[-2])
            beamPAs  = np.append(beamPAs, beamPA)
            beamsize = line[0]
            beamLA   = float(beamsize.split('x')[-1])
            beamLAs  = np.append(beamLAs, beamLA)
            beamSA   = float(beamsize.split('x')[0])
            beamSAs  = np.append(beamSAs, beamSA)
            print beamPAs,beamSAs,beamLAs
## calculate separation between target and prIBC
sep = howfun.separation(RA,Dec,prIBCRA,prIBCDec) #unit: arcmin
print sep           
## parse sum file for elevation, and calculate averaged csc(el)
sumfiles = glob.glob(r'%s/*/*.sum' % targetdir)
sumfiles.sort()
print sumfiles
Av_cscEls = np.array([])
for sumfile in sumfiles:
    TelAv_csc_Els = np.array([])
    lines = open(sumfile).readlines()
    for line in lines:
        if targetname in line:
            if ':' in line.split(targetname)[0]:
                elevations = line.split('-')[-1].strip().split('    ')
                try:
                    elevations = map(float, elevations)
                except ValueError:
                    continue
                elevations = np.asarray(elevations)
                elevations_rad = elevations*math.pi/180
                csc_elevations = (np.sin(elevations_rad))**(-1)
                TelAv_csc_El = np.average(csc_elevations) #telescope averaged csc(el) for each scan
                TelAv_csc_Els = np.append(TelAv_csc_Els, TelAv_csc_El)
    Av_cscEl = np.average(TelAv_csc_Els) #averaged csc(el) over telescopes and scans
    Av_cscEls = np.append(Av_cscEls,Av_cscEl)
    print Av_cscEls
## calculate delta_sys, or systematic errors in unit of beam, for each epochs
parameterA = 0.001
parameterB = 0.6
delta_sys = parameterA*sep*Av_cscEls+parameterB/SNprIBCs
print delta_sys
## calculate RA and Dec extension of synthesized beams
beamRAs = np.array([])
beamDecs = np.array([])
for i in range(len(beamPAs)):
    [beamRA, beamDec] = howfun.deprojectbeam2xy(beamLAs[i],beamSAs[i],beamPAs[i]) #full-width deprojection on RA/Dec from beam
    beamRAs  = np.append(beamRAs, beamRA)
    beamDecs = np.append(beamDecs, beamDec)
print beamRAs,beamDecs
## calculate systematic errors with beamRAs_rect, beamDecs_rect and delta_sys
sysErrRAs  = beamRAs*delta_sys #in mas
sysErrDecs = beamDecs*delta_sys #in mas
print sysErrRAs,sysErrDecs
## generate pmpar.in including the systematic errors
sysErrRAs = howfun.mas2ms(sysErrRAs,Dec) #in ms
print sysErrRAs
sysErrRAs  = sysErrRAs/1000 #ms -> s
sysErrDecs = sysErrDecs/1000 #mas ->"
errorRAs   = (error0RAs**2 + sysErrRAs**2)**0.5
errorDecs  = (error0Decs**2 + sysErrDecs**2)**0.5
print errorRAs,errorDecs

pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in'
fileWrite = open(pulsitions, 'w')
fileWrite.write("### pmpar format\n")
fileWrite.write("#name = " + targetname + "\n")
fileWrite.write("#ref = " + targetname + "\n")
fileWrite.write("epoch = %f\n" % epoch)
fileWrite.write("#pi = 0\n")
fileWrite.write("# decimalyear RA +/- Dec +/-\n")
for i in range(nepoch):
#for i in [0, 1, 2, 5, 6, 7, 9]: # for J1537+1155
    fileWrite.write("%s %s %.7f %s %.6f\n" % (decyears[i], RAs[i], errorRAs[i], Decs[i], errorDecs[i]))
fileWrite.close()
## get parallax, mu_a and mu_d
pmparout = pmparesultsdir + '/' + targetname + '.pmpar.out'
os.system("pmpar %s > %s" % (pulsitions, pmparout))
[PI0, mu_a0, mu_d0, RA0, Dec0, rchisq] = howfun.pmparout2paras(pmparout)
D0 = 1/PI0
print PI0,mu_a0,mu_d0,RA0,Dec0
##### bootstrapping begins here #####################################################################
####################################################################################################
if not bootstrap:
    sys.exit()
print "\nrunning bootstrapping...\n"
errorRAs_h  = errorRAs/3600 # s -> h
errorDecs_d = errorDecs/3600 # " -> d
print("%.10f %.9f" % (errorRAs[0], errorDecs[0]))

pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in.bootstrap'
if diagnosebootstrap != "":
    pmparesults_bootstrap = pmparesultsdir + '/diagnose_' + targetname + '.pmparesults.bootstrap'
    if os.path.exists(pmparesults_bootstrap):
        os.remove(pmparesults_bootstrap)
else:
    pmparesults_bootstrap = pmparesultsdir + '/' + targetname + '.pmparesults.bootstrap'
if montecarlo:
    pmparesults_bootstrap = pmparesultsdir + '/' + targetname + '.pmparesults.montecarlo'
    pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in.montecarlo'

count = 0
bootstrapruns = 100 #including monte-carlo
if diagnosebootstrap != '':
    bootstrapruns = 6000

if not montecarlo:
    print("\nrunning %d bootstrap runs...\n" % bootstrapruns)
    while count < bootstrapruns:
        random_indices = np.random.randint(0, nepoch, nepoch)
        if len(np.unique(random_indices)) < 4: #use at least 4 different positions for astrometric fit
            continue
        if diagnosebootstrap != '':
            if int(diagnosebootstrap) in random_indices: #locate the cause of double-peak structure 
                continue
        fileWrite = open(pulsitions, 'w')
        fileWrite.write("epoch = %f\n" % epoch)
        # use timing prioris
        if bootstrapprior != '':
            fileWrite.write("%s\n" % bootstrapprior)
        for i in random_indices:
            fileWrite.write("%s %s %.7f %s %.6f\n" % (decyears[i], RAs[i], errorRAs[i], Decs[i], errorDecs[i]))
        fileWrite.close()
        ## run pmpar and get pmparesults.bootstrap
        os.system('pmpar %s >> %s' % (pulsitions, pmparesults_bootstrap))
        print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / bootstrapruns)) + " \r")
        count += 1
else:
    print("\nrunning %d monte-carlo runs...\n" % bootstrapruns)
    while count < bootstrapruns:
        RAs_MC_h = np.random.normal(howfun.dms2deg(RAs), errorRAs_h)
        Decs_MC_d = np.random.normal(howfun.dms2deg(Decs), errorDecs_d)
        RAs_MC = howfun.deg2dms(RAs_MC_h)
        Decs_MC = howfun.deg2dms(Decs_MC_d)

        fileWrite = open(pulsitions, 'w')
        fileWrite.write("epoch = %f\n" % epoch)
        for i in range(nepoch):
            fileWrite.write("%s %s %.7f %s %.6f\n" % (decyears[i], RAs_MC[i], errorRAs[i], Decs_MC[i], errorDecs[i]))
        fileWrite.close()
        ## run pmpar and get pmparesults.bootstrap
        os.system('pmpar %s >> %s' % (pulsitions, pmparesults_bootstrap))
        print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / bootstrapruns)) + " \r")
        count += 1
    
## read pmparesults.bootstrap and get PIs, mu_as, mu_ds
randomparas = pmparesultsdir + '/.' + targetname + '.pmparesults.bootstrap.swp' #buffer file
if bootstrapprior != '':
    randomparas_suffix = re.sub('[=.]', '', bootstrapprior) 
    randomparas += '.for.plot.' + randomparas_suffix
print randomparas
print("shorten %s...\n" % pmparesults_bootstrap)

#if bootstrapprior != '':
#    for parameter in ['pi', 'mu_a', 'mu_d']:
#        if parameter in bootstrapprior:


mu_as = np.array([])
mu_ds = np.array([])
PIs   = np.array([])
RA1s = np.array([])
Dec1s = np.array([])
fileWrite = open(randomparas, 'w')
lines = open(pmparesults_bootstrap).readlines()
for line in lines:
    if ('mu_a' in line) or ('mu_d' in line) or ('pi' in line) or ('RA' in line) or ('Dec  ' in line):
        line = line.split('+-')[0].strip()
        fileWrite.write(line + '\n')
        line1 = line.split('=')[-1].strip()
        if 'RA' in line:
            RA1s = np.append(RA1s, howfun.dms2deg(line1))
        if 'Dec' in line:
            Dec1s = np.append(Dec1s, howfun.dms2deg(line1))
        if 'mu_a' in line:
            mu_as = np.append(mu_as, float(line1))
        if 'mu_d' in line:
            mu_ds = np.append(mu_ds, float(line1))
        if 'pi' in line:
            PIs = np.append(PIs, float(line1))
fileWrite.close()
os.system('cat %s > %s' % (randomparas, pmparesults_bootstrap))
Ds = 1./PIs
#VTs = (mu_as**2+mu_ds**2)**0.5*Ds*mspsrpifun.estimate_uncertainty.A
VT0 = (mu_a0**2+mu_d0**2)**0.5*D0*mspsrpifun.estimate_uncertainty.A

## estimate PI, mu_a and mu_d from the bootstrap samples
PIs.sort()
mu_as.sort()
mu_ds.sort()
RA1s.sort()
Dec1s.sort()
#VTs.sort()

if estnd:
    err_volume_3D, ratios_3D, estimates_3D, errors_3D = howfun.nD_sample2estimates(1, PIs, mu_as, mu_ds)
    print err_volume_3D, ratios_3D, estimates_3D, errors_3D
    estimates_3D_min = abs(estimates_3D) - errors_3D
    estimates_3D_max = abs(estimates_3D) + errors_3D
    VT_min = ((estimates_3D_min[1])**2+(estimates_3D_min[2])**2)**0.5/estimates_3D_max[0]*mspsrpifun.estimate_uncertainty.A
    VT_max = ((estimates_3D_max[1])**2+(estimates_3D_max[2])**2)**0.5/estimates_3D_min[0]*mspsrpifun.estimate_uncertainty.A
    print "VT = %f + %f - %f" % (VT0, VT_max-VT0, VT0-VT_min)

bootstrap_estimates_output = pmparesultsdir + '/' + targetname + '.bootstrap.estimates.out'
fileWrite = open(bootstrap_estimates_output, 'w')
fileWrite.write("estimates obtained with %d bootstrap runs:\n" % len(PIs))
for confidencelevel in [0.9973, 0.9545, 0.6827]: # 3sigma, 2sigma, 1sigma
    print len(PIs)
    most_probable_PI = howfun.sample2most_probable_value(PIs, 500)
    most_probable_mu_a = howfun.sample2most_probable_value(mu_as, 300)
    most_probable_mu_d = howfun.sample2most_probable_value(mu_ds, 150)
    for estimate in ['PI', 'mu_a', 'mu_d']:
        exec("[value_%s, error_%s] = howfun.sample2estimate(%ss, confidencelevel)" % (estimate, 
        estimate, estimate))
        exec("error_%s_symm = howfun.sample2uncertainty(%ss, %s0, confidencelevel)" % (
        estimate, estimate, estimate))
        exec("error_%s_symm1 = howfun.sample2uncertainty(%ss, most_probable_%s, confidencelevel)" % (estimate, estimate, estimate))
        exec("median_%s = howfun.sample2median(%ss)" % (estimate, estimate))
        if confidencelevel==0.9545:
            exec('%s_upper_limit_pdf_plot = value_%s + error_%s + 0.2*abs(error_%s)' % (estimate, estimate, estimate, estimate))
            exec('%s_lower_limit_pdf_plot = value_%s - error_%s - 0.2*abs(error_%s)' % (estimate, estimate, estimate, estimate))
    for estimate in ['RA', 'Dec']:
        exec("[value_%s1, error_%s1] = howfun.sample2estimate(%s1s, confidencelevel)" % (estimate, 
        estimate, estimate))
        exec("error_%s1_symm = howfun.sample2uncertainty(%s1s, %s0, confidencelevel)" % (
        estimate, estimate, estimate))
        if confidencelevel==0.9545:
            exec('%s_upper_limit_pdf_plot = value_%s1 + error_%s1 + 0.2*abs(error_%s1)' % (estimate, estimate, estimate, estimate))
            exec('%s_lower_limit_pdf_plot = value_%s1 - error_%s1 - 0.2*abs(error_%s1)' % (estimate, estimate, estimate, estimate))
    errorDec1_symm_mas = error_Dec1_symm*3600*1000
    errorRA1_symm_ms = error_RA1_symm*3600*1000
    RA0_dms = howfun.deg2dms(RA0)
    Dec0_dms = howfun.deg2dms(Dec0)
    errorRA1_symm_mas = howfun.ms2mas(errorRA1_symm_ms, Dec0_dms)
    print("confidencelevel = %f" % confidencelevel)
    print("PI = %f +- %f (mas)" % (value_PI, error_PI))
    print("PI_symmetric = %f +- %f (mas)" % (PI0, error_PI_symm))
    print("mu_a_symm = %f +- %f (mas/yr)" % (mu_a0, error_mu_a_symm))
    print("mu_d_symm = %f +- %f (mas/yr)" % (mu_d0, error_mu_d_symm))
    print("std_D = %f" % (np.std(Ds)))
    print("mu_a = %f +- %f (mas/yr)" % (value_mu_a, error_mu_a))
    print("mu_d = %f +- %f (mas/yr)" % (value_mu_d, error_mu_d))
    print("RA = %s +- %f (mas)" % (RA0_dms, errorRA1_symm_mas))
    print("Dec = %s +- %f (mas)" % (Dec0_dms, errorDec1_symm_mas))
    fileWrite.write(70*"="+"\n")
    fileWrite.write("confidencelevel = %f:\n" % confidencelevel)
    fileWrite.write("epoch = %f\n" % epoch)
    fileWrite.write("RA = %s +- %f (mas)\n" % (RA0_dms, errorRA1_symm_mas))
    fileWrite.write("Dec = %s +- %f (mas)\n" % (Dec0_dms, errorDec1_symm_mas))
    fileWrite.write("pi = %f + %f - %f (mas)\n" % (PI0, value_PI + error_PI - PI0, PI0 - value_PI + error_PI))
    fileWrite.write("D_best_fit = %f + %f -%f (kpc)\n" % (D0, 1/(value_PI - error_PI) - D0, D0 - 1/(value_PI + error_PI)))
    fileWrite.write("mu_a = %f + %f - %f (mas/yr) #mu_a=mu_ra*cos(dec)\n" % (mu_a0, value_mu_a + error_mu_a - mu_a0, mu_a0 - value_mu_a + error_mu_a))
    fileWrite.write("mu_d = %f + %f - %f (mas/yr)\n" % (mu_d0, value_mu_d + error_mu_d - mu_d0, mu_d0 - value_mu_d +error_mu_d))
    #fileWrite.write("v_t = %f + %f - %f (km/s)\n" % (VT0, value_VT +error_VT -VT0, VT0 -value_VT +error_VT))
    fileWrite.write("PI_symm_from_best_fit = %f +- %f (mas)\n" % (PI0, error_PI_symm))
    fileWrite.write("mu_a_symm_from_best_fit = %f +- %f (mas/yr)\n" % (mu_a0, error_mu_a_symm))
    fileWrite.write("mu_d_symm_from_best_fit = %f +- %f (mas/yr)\n" % (mu_d0, error_mu_d_symm))
    #fileWrite.write("D_symm = %f +- %f -%f (kpc)\n" % (D0, error_D_symm))
    #fileWrite.write("v_t_symm = %f +-%f (km/s)\n" % (VT0, error_VT_symm))
    fileWrite.write("PI_sample = %f +- %f (mas)\n" % (value_PI, error_PI))
    fileWrite.write("mu_a_sample = %f +- %f (mas/yr)\n" % (value_mu_a, error_mu_a))
    fileWrite.write("mu_d_sample = %f +- %f (mas/yr)\n" % (value_mu_d, error_mu_d))
    fileWrite.write("median_PI = %f\n" % median_PI)
    fileWrite.write('median_mu_a = %f\n' % median_mu_a)
    fileWrite.write("median_mu_d = %f\n" % median_mu_d)
    fileWrite.write("most_probable_PI = %f + %f - %f (mas)\n" % (most_probable_PI, value_PI+error_PI-most_probable_PI, most_probable_PI - (value_PI-error_PI)))
    fileWrite.write("most_probable_mu_a = %f + %f - %f (mas/yr)\n" % (most_probable_mu_a, value_mu_a+error_mu_a-most_probable_mu_a, most_probable_mu_a - (value_mu_a-error_mu_a)))
    fileWrite.write("most_probable_mu_d = %f + %f - %f (mas/yr)\n" % (most_probable_mu_d, value_mu_d+error_mu_d-most_probable_mu_d, most_probable_mu_d - (value_mu_d-error_mu_d)))
    fileWrite.write("D = %f + %f -%f (kpc)\n" % (1/most_probable_PI, 1/(value_PI - error_PI) - 1/most_probable_PI, 1/most_probable_PI - 1/(value_PI + error_PI)))
    VT1 = (most_probable_mu_a**2+most_probable_mu_d**2)**0.5/most_probable_PI*mspsrpifun.estimate_uncertainty.A
    VT1_upper = ((abs(value_mu_a)+error_mu_a)**2+(abs(value_mu_d)+error_mu_d)**2)**0.5/(value_PI-error_PI)*mspsrpifun.estimate_uncertainty.A
    VT1_lower = ((abs(value_mu_a)-error_mu_a)**2+(abs(value_mu_d)-error_mu_d)**2)**0.5/(value_PI+error_PI)*mspsrpifun.estimate_uncertainty.A
    fileWrite.write("VT = %f + %f - %f (km/s)\n" % (VT1, VT1_upper-VT1, VT1-VT1_lower))
    fileWrite.write("PI_symm = %f +- %f (mas)\n" % (most_probable_PI, error_PI_symm1))
    fileWrite.write("mu_a_symm = %f +- %f (mas)\n" % (most_probable_mu_a, error_mu_a_symm1))
    fileWrite.write("mu_d_symm = %f +- %f (mas)\n" % (most_probable_mu_d, error_mu_d_symm1))
    fileWrite.write("D_symm = %f + %f -%f (kpc)\n" % (1/most_probable_PI, 1/(most_probable_PI - error_PI_symm1) - 1/most_probable_PI, 1/most_probable_PI - 1/(most_probable_PI + error_PI_symm1)))
    VT1_symm_upper = ((abs(most_probable_mu_a)+error_mu_a_symm1)**2+(abs(most_probable_mu_d)+error_mu_d_symm1)**2)**0.5/(most_probable_PI-error_PI_symm1)*mspsrpifun.estimate_uncertainty.A
    VT1_symm_lower = ((abs(most_probable_mu_a)-error_mu_a_symm1)**2+(abs(most_probable_mu_d)-error_mu_d_symm1)**2)**0.5/(most_probable_PI+error_PI_symm1)*mspsrpifun.estimate_uncertainty.A
    fileWrite.write("VT_symm = %f + %f - %f (km/s)\n" % (VT1, VT1_symm_upper-VT1, VT1-VT1_symm_lower))
fileWrite.close()

## plot and save probability density functions for PIs, mu_as and mu_ds
if not plotbootstrap:
    sys.exit()
print "\nbaking bootstrap plots...\n"
bins = 400
diagnoseBSstring = '_noepoch' + diagnosebootstrap
if diagnosebootstrap == '':
    diagnoseBSstring = ''
for parameter in [PIs, mu_as, mu_ds, RA1s, Dec1s]:
    plt.rc('text', usetex=True)
    plt.hist(parameter, bins, density=True, facecolor='g', alpha=0.75)
    #plt.rcParams.update({'font.size': 22})
    minimum = min(parameter)
    maximum = max(parameter)
    if maximum == minimum:
        continue
    x = np.arange(minimum,maximum,(maximum-minimum)/1000)
    if (parameter == PIs).any():
        y = norm.pdf(x, value_PI, error_PI)
    if (parameter == mu_as).any():
        y = norm.pdf(x, value_mu_a, error_mu_a)
    if (parameter == mu_ds).any():
        y = norm.pdf(x, value_mu_d, error_mu_d)
    if (parameter == RA1s).any():
        y = norm.pdf(x, value_RA1, error_RA1)
    if (parameter == Dec1s).any():
        y = norm.pdf(x, value_Dec1, error_Dec1)
    plt.plot(x,y)
    plt.ylabel('probability density')
    if (parameter == PIs).any():
        plt.xlabel('parallax (mas)')
        plt.xlim(PI_lower_limit_pdf_plot, PI_upper_limit_pdf_plot)
        plt.savefig('%s/parallax_%sbins_pdf%s.eps' % (pmparesultsdir, bins, diagnoseBSstring))
    if (parameter == mu_as).any():
        plt.xlabel(r'$\rm \mu_{\alpha}~(mas~yr^{-1})$')
        plt.xlim(mu_a_lower_limit_pdf_plot, mu_a_upper_limit_pdf_plot)
        plt.savefig('%s/mu_a_%sbins_pdf%s.eps' % (pmparesultsdir, bins, diagnoseBSstring))
    if (parameter == mu_ds).any():
        plt.xlabel(r'$\rm \mu_{\delta}cos(\delta)~(mas~yr^{-1})$')
        plt.xlim(mu_d_lower_limit_pdf_plot, mu_d_upper_limit_pdf_plot)
        plt.savefig('%s/mu_d_%sbins_pdf%s.eps' % (pmparesultsdir, bins, diagnoseBSstring))
    if (parameter == RA1s).any():
        plt.xlabel(r'Right Ascension')
        plt.savefig('%s/RA_%sbins_pdf%s.eps' % (pmparesultsdir, bins, diagnoseBSstring))
    if (parameter == Dec1s).any():
        plt.xlabel(r'Declination')
        plt.savefig('%s/Dec_%sbins_pdf%s.eps' % (pmparesultsdir, bins, diagnoseBSstring))
    plt.clf()

## plot 3+2 grid of the 5 PDFs
fig = plt.figure(figsize=[8,3])
gs = gridspec.GridSpec(2, 6)


ax1 = fig.add_subplot(gs[:2, :2])
#ax1.plot(range(0,10), range(0,10))
ax1.hist(PIs, 500, density=True, facecolor='g', alpha=0.75)
#[bin_counts, bin_values, junk] = ax1.hist(PIs, 500, density=True, facecolor='g', alpha=0.75)
#index_max_value = np.argmax(bin_counts)
#most_probable_value = bin_values[index_max_value]
#minimum = min(PIs)
#maximum = max(PIs)
#x = np.arange(minimum,maximum,(maximum-minimum)/1000)
#y = norm.pdf(x, value_PI, error_PI)
#ax1.plot(x,y)
ax1.set_xlabel('parallax (mas)')
ax1.set_xlim(PI_lower_limit_pdf_plot, PI_upper_limit_pdf_plot)
#ax1.set_xlim(0.9, 1.5)
ax1.set_ylabel('probability density')
ax1.axvline(x=most_probable_PI, c='black', linestyle='--', linewidth=0.5)
ax1.axvline(x=value_PI+error_PI, c='black', linestyle='-.', linewidth=0.5)
ax1.axvline(x=value_PI-error_PI, c='black', linestyle='-.', linewidth=0.5)

ax2 = fig.add_subplot(gs[:2, 2:4])
ax2.hist(mu_as, 300, density=True, facecolor='g', alpha=0.75)
#minimum = min(mu_as)
#maximum = max(mu_as)
#x = np.arange(minimum,maximum,(maximum-minimum)/1000)
#y = norm.pdf(x, value_mu_a, error_mu_a)
#ax2.plot(x,y)
#ax2.set_xlim(mu_a_lower_limit_pdf_plot, mu_a_upper_limit_pdf_plot)
ax2.set_xlim(2.45, 2.9)
ax2.set_xlabel(r'$\rm \mu_{\alpha}~(mas~yr^{-1})$')
ax2.axvline(x=most_probable_mu_a, c='black', linestyle='dashed', linewidth=0.5)
ax2.axvline(x=value_mu_a+error_mu_a, c='black', linestyle='-.', linewidth=0.5)
ax2.axvline(x=value_mu_a-error_mu_a, c='black', linestyle='-.', linewidth=0.5)

ax3 = fig.add_subplot(gs[:2, 4:6])
ax3.hist(mu_ds, 150, density=True, facecolor='g', alpha=0.75)
#minimum = min(mu_ds)
#maximum = max(mu_ds)
#x = np.arange(minimum,maximum,(maximum-minimum)/1000)
#y = norm.pdf(x, value_mu_d, error_mu_d)
#ax3.plot(x,y)
ax3.set_xlim(mu_d_lower_limit_pdf_plot, mu_d_upper_limit_pdf_plot)
ax3.set_xlabel(r'$\rm \mu_{\delta}cos\delta~(mas~yr^{-1})$')
ax3.axvline(x=most_probable_mu_d, c='black', linestyle='dashed', linewidth=0.5)
ax3.axvline(x=value_mu_d+error_mu_d, c='black', linestyle='-.', linewidth=0.5)
ax3.axvline(x=value_mu_d-error_mu_d, c='black', linestyle='-.', linewidth=0.5)

gs.tight_layout(fig) #rect=[0, 0.1, 1, 1])
plt.savefig('%s/PDF_grid.eps' % pmparesultsdir)
print "\nwe successfully fly over\n"
