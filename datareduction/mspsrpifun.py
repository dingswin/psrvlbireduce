import os,glob,sys,yaml,math,pickle,datetime
import howfun,mspsrpifun
import numpy as np
from optparse import OptionParser
np.set_printoptions(suppress=True)

## group 1a ################################################################################
def prepare_path_source(targetname):
    auxdir    = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    expconfigfile = configdir + '/' + targetname + '.yaml'
    targetdir = auxdir + '/processing/' + targetname
    if not os.path.exists(targetdir):
        print("%s doesn't exist; aborting\n" % targetdir)
        sys.exit()
    expconfig = yaml.load(open(expconfigfile))
    prIBCname = expconfig['primaryinbeam']
    phscalname = target2cals(targetname)[0]
    print phscalname, prIBCname
    return auxdir, configdir, targetdir, phscalname, prIBCname

def target2cals(targetname): #get phscal and bandpass cal given targetname
    auxdir = os.environ['PSRVLBAUXDIR']
    targetdir = auxdir + '/processing/' + targetname
    sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
    if sourcefiles == []:
        print("source files not found; abort")
        sys.exit()
    sourcefiles.sort()
    sourcefile = sourcefiles[0]
    lines = open(sourcefile).readlines()
    for line in lines:
        if 'BANDPASS' in line:
            bpcal = line.split(':')[-1].strip()
        if 'PHSREF' in line:
            phscal = line.split(':')[-1].strip()
    cals = [phscal, bpcal]
    return cals
    
def nonpulsar_statsfiles2positions(statsfiles): # tailored for nonpulsarjmfit format
    statsfiles.sort()
    RAs = np.array([])
    Decs = np.array([])
    for statsfile in statsfiles:
        lines = open(statsfile).readlines()
        for line in lines:
            if 'Actual RA' in line:
                temp = line.split('RA:')[-1].strip()
                RAs = np.append(RAs, temp)
            if 'Actual Dec' in line:
                temp = line.split('Dec:')[-1].strip()
                Decs = np.append(Decs, temp)
    return RAs, Decs

def statsfiles2positions(statsfiles):
    pass

def dms_positions2stat(RAs, Decs):
    [RA_average, RA_std_mas] = dms_array2stat(RAs)
    [Dec_average, Dec_std_mas] = dms_array2stat(Decs)
    Dec_av_rad = Dec_average*math.pi/180
    RA_std_mas = RA_std_mas*15*math.cos(Dec_av_rad)
    RA_average_dms = howfun.deg2dms(RA_average)
    Dec_average_dms = howfun.deg2dms(Dec_average)
    return [RA_average_dms, RA_std_mas], [Dec_average_dms, Dec_std_mas]
def dms_array2stat(array):
    if len(array) < 2:
        print("len(list)<2, thus unable to carry out statistics; abort")
        sys.exit()
    degrees = howfun.dms2deg(array)
    average = np.average(degrees)
    std_degree = np.std(degrees)
    std_mas = std_degree*3600*1000
    return average, std_mas

def target2positionscatter(targetname):
    #####################################################################################
    ## compile stats files (for phscal and primary IBC) to get systematics
    ## tailored for MSPSRPI datasets
    ## by Hao Ding on 5 May 2019
    #####################################################################################
    
    ## get paths and sourcenames ########################################################
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    ## find and parse statsfiles, reach RAs/Decs #########################
    phscalstatsfiles = glob.glob(r'%s/*/*_ibshiftdiv_difmap.jmfit.stats' % targetdir) #find statsfile for each epoch
    [phscalRAs, phscalDecs] = nonpulsar_statsfiles2positions(phscalstatsfiles)
    print phscalRAs, phscalDecs
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.divided.difmap.jmfit.stats' % targetdir)
    [prIBC_RAs, prIBC_Decs] = nonpulsar_statsfiles2positions(prIBCstatsfiles)
    print prIBC_RAs, prIBC_Decs
    ## statistics about RAs/Decs ##########################################
    return dms_positions2stat(phscalRAs, phscalDecs), dms_positions2stat(prIBC_RAs, prIBC_Decs)

def targetname2decyears(targetname):
    os.system('generatepmparin.py -t %s' % targetname)
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    tmpfile = targetdir + '/.variable_decyears.tmp'
    readfile = open(tmpfile, 'r')
    decyears = pickle.load(readfile)
    readfile.close()
    return decyears
def diffposition(RAs, refRA, Decs, refDec):
    diffDecs = howfun.dms2deg(Decs) - howfun.dms2deg(refDec)
    diffRAs = howfun.dms2deg(RAs) - howfun.dms2deg(refRA)
    diffDecs_mas = diffDecs*3600*1000
    diffRAs_ms = diffRAs*3600*1000
    refDec_rad = howfun.dms2deg(refDec)*math.pi/180
    diffRAs_mas = diffRAs_ms*15*math.cos(refDec_rad)
    return diffRAs_mas, diffDecs_mas
def prepareplotscatter(targetname): #generate input for MATLAB plotting
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    phscalstatsfiles = glob.glob(r'%s/*/*_ibshiftdiv_difmap.jmfit.stats' % targetdir) #find statsfile for each epoch
    [phscalRAs, phscalDecs] = nonpulsar_statsfiles2positions(phscalstatsfiles)
    [refRAphscal, refDECphscal] = srcposition(targetname, phscalname)
    [diffRAphscal, diffDECphscal] = diffposition(phscalRAs, refRAphscal, phscalDecs, refDECphscal)
    print diffRAphscal, diffDECphscal
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.divided.difmap.jmfit.stats' % targetdir)
    [prIBC_RAs, prIBC_Decs] = nonpulsar_statsfiles2positions(prIBCstatsfiles)
    [refRAprIBC, refDECprIBC] = srcposition(targetname, prIBCname)
    [diffRAprIBC, diffDECprIBC] = diffposition(prIBC_RAs, refRAprIBC, prIBC_Decs, refDECprIBC) 
    decyears = targetname2decyears(targetname)
    print decyears
    # now generate input file for MATLAB plotting
    nepoch = len(decyears)
    outfile = targetdir + '/pmparesults/' + phscalname + '_positionscatter.txt'
    writefile = open(outfile, 'w')
    writefile.write('!time(yr) RA(mas) Dec(mas) (relative to %s %s)\n' % (refRAphscal, refDECphscal))
    for i in range(nepoch):
        writefile.write('%s %f %f\n' % (decyears[i], diffRAphscal[i], diffDECphscal[i]))
    writefile.close()
    
    outfile = targetdir + '/pmparesults/' + prIBCname + '_positionscatter.txt'
    writefile = open(outfile, 'w')
    writefile.write('!time(yr) RA(mas) Dec(mas) (relative to %s %s)\n' % (refRAprIBC, refDECprIBC))
    for i in range(nepoch):
        writefile.write('%s %f %f\n' % (decyears[i], diffRAprIBC[i], diffDECprIBC[i]))
    writefile.close()
    print("two txt files are prepared for MATLAB plotting.")
## group 1b #############################################################################
def readpulsition(targetname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)    
    pmparout = targetdir + '/pmparesults/' + targetname + '.pmpar.out'
    return readpmparout(pmparout)
def readpmparout(pmparout):
    lines = open(pmparout).readlines()
    for line in lines:
        if 'epoch' in line:
            epoch = line.split('=')[1].strip()
        for estimate in ['RA', 'Dec  ', 'mu_a', 'mu_d', 'pi']:
            if estimate in line:
                print line.split('=')[-1].split('+')[0].strip()
                exec("%s = line.split('=')[-1].split('+')[0].strip()" % estimate.strip())
    return RA, Dec, epoch, pi, mu_a, mu_d
def bootstrapRADECerr(targetname, HowManySigma):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    BSresults = targetdir + '/pmparesults/' + targetname + '.bootstrap.estimates.out'
    errRAs = np.array([])
    errDecs = np.array([])
    lines = open(BSresults).readlines()
    for line in lines:
        if 'RA' in line:
            errRA = line.split('+-')[-1].split('(mas)')[0].strip()
            errRAs = np.append(errRAs, float(errRA))
        if 'Dec' in line:
            errDec = line.split('+-')[-1].split('(mas)')[0].strip()
            errDecs = np.append(errDecs, float(errDec))
    if HowManySigma == 1:
        errRA = errRAs[0]
        errDec = errDecs[0]
    elif HowManySigma == 2:
        errRA = errRAs[1]
        errDec = errDecs[1]
    else:
        print('the second parameter should be either 1 or 2 (sigma); aborting')
        sys.exit()
    err = np.array([errRA, errDec])
    return err
def abspsrposition(targetname, HowManySigma):
    [phscalstats, prIBCstats] = target2positionscatter(targetname)
    phscalRADEC = [phscalstats[0][0], phscalstats[1][0]]
    phscalRADECstd = np.array([phscalstats[0][1], phscalstats[1][1]])
    prIBC_RADEC = [prIBCstats[0][0], prIBCstats[1][0]]
    prIBC_RADECstd = np.array([prIBCstats[0][1], prIBCstats[1][1]])
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    [psrRA, psrDec, epoch, junk1, junk2, junk3] = readpulsition(targetname)
    psrRADEC0 = [psrRA, psrDec]
    rfcinfos = rfcposition(phscalname)
    phscal_absRADEC = [rfcinfos[0], rfcinfos[2]]
    errRAphscal = float(rfcinfos[1])
    errDECphscal = float(rfcinfos[3])
    ## 1) using prIBC shift ####################
    prIBCrefRADEC = srcposition(targetname, prIBCname)
    psrRADEC1 = howfun.dms2deg(psrRADEC0) + howfun.dms2deg(prIBC_RADEC) - howfun.dms2deg(prIBCrefRADEC)
    ## 2) using phscal shift ###################
    phscal_refRADEC = srcposition(targetname, phscalname)
    psrRADEC2 = howfun.dms2deg(psrRADEC0) + howfun.dms2deg(phscal_refRADEC) - howfun.dms2deg(phscalRADEC)
    ## 3) update absolute position for phscal #####
    psrRADEC2 = psrRADEC2 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
    psrRADEC2 = howfun.deg2dms(psrRADEC2)
    psrRADEC1 = psrRADEC1 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
    psrRADEC1 = howfun.deg2dms(psrRADEC1)
    ## error estimation #####################################################
    err_psr2prIBC = bootstrapRADECerr(targetname, HowManySigma)
    err_abs_phscal = np.array([errRAphscal, errDECphscal])
    err1 = (err_psr2prIBC**2 + err_abs_phscal**2 + prIBC_RADECstd**2)**0.5
    err2 = (err_psr2prIBC**2 + err_abs_phscal**2 + phscalRADECstd**2)**0.5
    return psrRADEC1, err1, psrRADEC2, err2, epoch
def exportabspsrposition(targetname, HowManySigma):
    [psrRADEC1, err1, psrRADEC2, err2, epoch] = abspsrposition(targetname, HowManySigma)
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    err1[0] = howfun.mas2ms(err1[0], str(psrRADEC1[1]))
    err2[0] = howfun.mas2ms(err2[0], str(psrRADEC2[1]))
    err1 = err1/1000
    err2 = err2/1000
    outfile = targetdir + '/' + targetname + '.abs.position'
    writefile = open(outfile, 'w')
    writefile.write("epoch = %s\n" % epoch)
    writefile.write("using primary IBC shift, we get:\n")
    writefile.write("RA = %s +- %f\n" % (psrRADEC1[0], err1[0]))
    writefile.write("Dec = %s +- %f\n" % (psrRADEC1[1], err1[1]))
    writefile.write("using phscal shift, we get:\n")
    writefile.write("RA = %s +- %f\n" % (psrRADEC2[0], err2[0]))
    writefile.write("Dec = %s +- %f\n" % (psrRADEC2[1], err2[1]))
    writefile.close()

## group 2 #############################################################################################
def targetbeams(targetname): #also get SNprIBCs, using prIBC statsfile,
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    beamLAs    = np.array([]) #beamLA -> beam Long axis length; beamSA -> beam short axis length;
    beamSAs    = np.array([]) 
    beamPAs    = np.array([]) #beamPA -> beam position angle;
    SNprIBCs    = np.array([])
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
            if 'beam' in line:
                line = line.split('beam')[-1].strip().split(' ')
                beamPA   = float(line[-2])
                beamPAs  = np.append(beamPAs, beamPA)
                beamsize = line[0]
                beamLA   = float(beamsize.split('x')[-1])
                beamLAs  = np.append(beamLAs, beamLA)
                beamSA   = float(beamsize.split('x')[0])
                beamSAs  = np.append(beamSAs, beamSA)
    return beamPAs, beamSAs, beamLAs, SNprIBCs

def srcposition(targetname, srcname): #src should be non-target calibrator
    os.system('srcposition.py %s %s' % (targetname, srcname))
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    tmpfile = targetdir + '/.variables.tmp'
    readfile = open(tmpfile, 'r')
    position = pickle.load(readfile)
    readfile.close()
    return position
def srcseparation(targetname, src1, src2):
    [RA1, Dec1] = srcposition(targetname, src1)
    [RA2, Dec2] = srcposition(targetname, src2)
    sep = howfun.separation(RA1,Dec1,RA2,Dec2) #unit: arcmin
    return sep

def Av_cscEls(targetname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
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
                    elevations = map(float, elevations)
                    elevations = np.asarray(elevations)
                    elevations_rad = elevations*math.pi/180
                    csc_elevations = (np.sin(elevations_rad))**(-1)
                    TelAv_csc_El = np.average(csc_elevations) #telescope averaged csc(el) for each scan
                    TelAv_csc_Els = np.append(TelAv_csc_Els, TelAv_csc_El)
        Av_cscEl = np.average(TelAv_csc_Els) #averaged csc(el) over telescopes and scans
        Av_cscEls = np.append(Av_cscEls,Av_cscEl)
    return Av_cscEls
def delta_sys(targetname, src1, src2): #empirical delta_syss for sysErrs between src1 and src2 in the astrometric experiment for targetname
    parameterA = 0.001
    parameterB = 0.6
    sep = srcseparation(targetname, src1, src2)
    Av_cscEls1 = Av_cscEls(targetname)
    [beamPAs, beamSAs, beamLAs, SNprIBCs] = targetbeams(targetname)
    delta_sys = parameterA*sep*Av_cscEls1+parameterB/SNprIBCs
    return delta_sys
def beam_in_RA_Dec(targetname):
    [beamPAs, beamSAs, beamLAs, SNprIBCs] = targetbeams(targetname)
    beamRAs = np.array([])
    beamDecs = np.array([])
    for i in range(len(beamPAs)):
        [beamRA, beamDec] = howfun.deprojectbeam2xy(beamLAs[i],beamSAs[i],beamPAs[i]) #full-width deprojection on RA/Dec from beam
        beamRAs  = np.append(beamRAs, beamRA)
        beamDecs = np.append(beamDecs, beamDec)
    return beamRAs, beamDecs 
def sysErr(targetname, src1, src2):
    [beamRAs, beamDecs] = beam_in_RA_Dec(targetname)
    delta_sys1 = delta_sys(targetname, src1, src2)    
    sysErrRAs  = beamRAs*delta_sys1 #in mas
    sysErrDecs = beamDecs*delta_sys1 #in mas
    return sysErrRAs, sysErrDecs

## group 3: draw online absolute position and errors ###############################################
def parse_rfctxt(phscal, rfctxt):
    lines = open(rfctxt).readlines() 
    for line in lines:
        if phscal in line:
            words = line.split('  ')
    word1 = []
    for word in words:
        if word not in ['', ' ', '  ', '   ']:
            word1.append(word.strip())
    for sign in ['+', '-']:
        if sign in word1[2]:
            position = word1[2].split(sign)
            RA = position[0].strip()
            Dec = sign + position[1].strip()
            RA = RA.replace(' ', ':')
            Dec = Dec.replace(' ', ':')
    return RA, word1[3], Dec, word1[4]
def rfcposition(phscal):
    auxdir = os.environ['PSRVLBAUXDIR']
    downloaddir = auxdir + '/downloads/'
    year = datetime.datetime.today().year
    year = str(year)
    versions = [year+'d', year+'c', year+'b', year+'a']
    for version in versions:
        filename = 'rfc_' + version + '_cat.txt'
        ftpfile = 'http://astrogeo.org/vlbi/solutions/rfc_' + version + '/' + filename
        downloadfile = downloaddir + filename
        if not os.path.exists(downloadfile):
            os.system('wget -P %s %s' % (downloaddir, ftpfile))
        if os.path.exists(downloadfile):
            print("%s already exists" % filename)
            break
    return parse_rfctxt(phscal, downloadfile)
