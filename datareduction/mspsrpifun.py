##########################################################################################
## Codes written for astrometry data analysis
## Hao Ding
##########################################################################################


import os,glob,sys,yaml,math,pickle,datetime
from scipy import constants
from astropy import constants as AC
import howfun
import numpy as np
from optparse import OptionParser
from astropy.table import Table
#import matplotlib.pyplot as plt
from scipy.stats import norm
np.set_printoptions(suppress=True)

## group 1a ################################################################################
def expno2sources(expno, inverse_referencing=False):
    """
    Note
    ----
    It might be little bit confusing, when inverse_referencing=True, the inverse-referenced source (normally an AGN) is set as the prIBC.

    Return parameters
    -----------------
    targetname, foldername, phscalname, prIBCname
    """
    auxdir    = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    expconfigfile = configdir + '/' + expno + '.yaml'
    expconfig = yaml.load(open(expconfigfile))
    targetname = expconfig['targets'][0]
    configfile = configdir + '/' + targetname + '.yaml'
    config = yaml.load(open(configfile))
    if not inverse_referencing:
        prIBCname = config['primaryinbeam']
    else:
        prIBCname = config['primarytarget']    
    targetdir = expconfig['rootdir']
    foldername = targetdir.split('/')[-1]
    if foldername == '':
        foldername = targetdir.split('/')[-2]
    try:
        othertargetname = config['othertargetname']
    except KeyError:
        othertargetname = ''
    if othertargetname != '':
        phscalname = target2cals(othertargetname, expno)[0]
    else:
        if not expconfig['dodefaultnames']:
            phscalname = target2cals(targetname, expno)[0]
        else:
            phscalname = target2cals(foldername, expno)[0]
    return targetname, foldername, phscalname, prIBCname

def prepare_path_source(targetname, inverse_referencing=False):
    """
    Note
    ----
    It might be little bit confusing, when inverse_referencing=True, the inverse-referenced source (normally an AGN) is set as the prIBC.

    Input parameters
    ----------------

    Return parameters
    -----------------
    auxdir, configdir, targetdir, phscalname, prIBCname
    """
    auxdir    = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    configfile = configdir + '/' + targetname + '.yaml'
    if not os.path.exists(configfile):
        configfile = configdir + '/PULSAR.yaml'
        if not os.path.exists(configfile):
            print('target yaml file does not exist; aborting')
            sys.exit()
    targetdir = auxdir + '/processing/' + targetname
    if not os.path.exists(targetdir):
        print("%s doesn't exist; return False\n" % targetdir)
        return False
        sys.exit()
    config = yaml.load(open(configfile))
    if not inverse_referencing:
        prIBCname = config['primaryinbeam']
    else:
        prIBCname = config['primarytarget']    
    phscalname = target2cals(targetname)[0]
    print phscalname, prIBCname
    return auxdir, configdir, targetdir, phscalname, prIBCname

def target2cals(targetname, expno=''): #get phscal and bandpass cal given targetname
    """
    Given that the source file feature will be deprecated in 2021, this function will search in exp.yaml first.
    If relevant keys are not found, then the function will proceed to search in source files.
    """
    auxdir = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    targetdir = auxdir + '/processing/' + targetname
    if expno == '':
        vexfiles = glob.glob(r'%s/*/*.vex' % targetdir)
        vexfiles.sort()
        vexfile = vexfiles[0]
        expno = vexfile.split('/')[-2].strip()
    expconfigfile = configdir + expno + '.yaml'
    if not os.path.exists(expconfigfile):
        print('%s does not exist; aborting' % expconfigfile)
    expconfig = yaml.load(open(expconfigfile))
    cals = []
    try: 
        phscalnames = expconfig['phscalnames']
        if type(phscalnames) == str:
            cals.append(phscalnames)
        else:
            for phscalname in phscalnames:
                cals.append(phscalname)
    except KeyError:
        pass
    try:
        bpcals = expconfig['ampcalsrc']
        if type(bpcals) == str:
            cals.append(bpcals)
        else:
            for bpcal in bpcals:
                cals.append(bpcal)
    except KeyError:
        ## if 'inbeamnames' is unfound, neither should be phscalnames, then proceeding to sourcefiles ###
        sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
        if sourcefiles == []:
            print("source files not found; abort")
            sys.exit()
        sourcefiles.sort()
        sourcefile = sourcefiles[0]
        if expno != '':
            for sourcefile1 in sourcefiles:
                if expno in sourcefile1:
                    sourcefile = sourcefile1        
        lines = open(sourcefile).readlines()
        for line in lines:
            if 'BANDPASS' in line:
                bpcal = line.split(':')[-1].strip()
            if 'PHSREF' in line:
                phscal = line.split(':')[-1].strip()
        cals = [phscal, bpcal]
        if expno != '':
            for sourcefile1 in sourcefiles:
                if expno in sourcefile1:
                    sourcefile = sourcefile1        
    return cals

def target2phscals_and_inbeamcals(targetname, expno=''):
    """
    Given that the source file feature will be deprecated in 2021, this function will search in exp.yaml first.
    If relevant keys are not found, then the function will proceed to search in source files.
    """
    auxdir = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    targetdir = auxdir + '/processing/' + targetname
    if expno == '':
        vexfiles = glob.glob(r'%s/*/*.vex' % targetdir)
        vexfiles.sort()
        vexfile = vexfiles[0]
        expno = vexfile.split('/')[-2].strip()
    expconfigfile = configdir + expno + '.yaml'
    if not os.path.exists(expconfigfile):
        print('%s does not exist; aborting' % expconfigfile)
    expconfig = yaml.load(open(expconfigfile))
    cals = []
    try: 
        phscalnames = expconfig['phscalnames']
        if type(phscalnames) == str:
            cals.append(phscalnames)
        else:
            for phscalname in phscalnames:
                cals.append(phscalname)
    except KeyError:
        pass
    try:
        inbeamnames = expconfig['inbeamnames']
        if type(inbeamnames) == str:
            cals.append(inbeamnames)
        else:
            for inbeamname in inbeamnames:
                cals.append(inbeamname)
    except KeyError:
        ## if 'inbeamnames' is unfound, neither should be phscalnames, then proceeding to sourcefiles ###
        sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
        if sourcefiles == []:
            print("source files not found; abort")
            sys.exit()
        sourcefiles.sort()
        sourcefile = sourcefiles[0]
        if expno != '':
            for sourcefile1 in sourcefiles:
                if expno in sourcefile1:
                    sourcefile = sourcefile1        
        lines = open(sourcefile).readlines()
        cals = []
        for line in lines:
            if 'PHSREF' in line:
                phscal = line.split(':')[-1].strip()
                cals.append(phscal)
            if 'INBEAM' in line and 'NAME' in line:
                inbeamcal = line.split(':')[-1].strip()
                cals.append(inbeamcal)
    return cals

def print_out_calibrator_plan(targetname, expno='', savefig=False):
    [RA_t, Dec_t] = srcposition(targetname, targetname)
    cals = target2phscals_and_inbeamcals(targetname, expno)
    RAs_cals = np.array([])
    Decs_cals = np.array([])
    for cal in cals:
        [RA, Dec] = srcposition(targetname, cal)
        RAs_cals = np.append(RAs_cals, RA)
        Decs_cals = np.append(Decs_cals, Dec)
    plot_calibrator_plan(targetname, cals, RA_t, Dec_t, RAs_cals, Decs_cals, savefig)
def plot_calibrator_plan(targetname, cals, RA_t, Dec_t, RAs, Decs, savefig=False):
    import matplotlib.pyplot as plt
    RAs_h = howfun.dms2deg(RAs)
    Decs_deg = howfun.dms2deg(Decs)
    RA_t_h = howfun.dms2deg(RA_t)
    Dec_t_deg = howfun.dms2deg(Dec_t)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(RAs_h, Decs_deg)
    ax.set_ylabel('Declination (deg)')
    ax.set_xlabel('Right Ascension (h)')
    plt.gca().invert_xaxis()
    aspect_ratio = 1./15/math.cos(Dec_t_deg/180*math.pi)
    ax.set_aspect(aspect_ratio)
    ax.plot(RA_t_h, Dec_t_deg, 'rs')
    ax.annotate(targetname, (RA_t_h, Dec_t_deg))
    for i, cal in enumerate(cals):
        ax.annotate(cal, (RAs_h[i], Decs_deg[i]))
    if savefig:
        [junk1, junk2, targetdir, junk3, junk4] = prepare_path_source(targetname)
        directory2save = targetdir + '/pmparesults/'
        if not os.path.exists(directory2save):
            os.makedirs(directory2save)
        figname2save = directory2save + targetname + '_calibrator_plan.pdf'
        plt.savefig(figname2save)
    else:
        plt.show()
    plt.clf()

def nonpulsar_statsfile2position(statsfile):
    lines = open(statsfile).readlines()
    for line in lines:
        if 'Actual RA' in line:
            RA = line.split('RA:')[-1].strip() # no matter how many 'Actual RA's, we use the last one
        if 'Actual Dec' in line:
            Dec = line.split('Dec:')[-1].strip()
    return RA, Dec 
def nonpulsar_statsfiles2positions(statsfiles): # tailored for nonpulsarjmfit format
    statsfiles.sort()
    RAs = np.array([])
    Decs = np.array([])
    for statsfile in statsfiles:
        [RA, Dec] = nonpulsar_statsfile2position(statsfile)
        RAs = np.append(RAs, RA)
        Decs = np.append(Decs, Dec)
    return RAs, Decs

def statsfiles2positions(statsfiles):
    pass

def statsfiles2expnos(statsfiles):
    statsfiles.sort()
    expnos = []
    for statsfile in statsfiles:
        expno = statsfile.split('/')[-2]
        expno = expnos.append(expno)
    return expnos

def dms_positions2stat(RAs, Decs):
    """
    Function
    --------
    Consume calibrator positions (relative to another calibrator) and calculate the average position and the position scatter.
    
    Outlier exclusion
    -----------------
    3-sigma theshold is used to exclude outliers in an iterative way.

    Output parameters
    -----------------
    return [RA_average_dms, RA_std_mas], [Dec_average_dms, Dec_std_mas]
    """
    outlier_count = 0
    while outlier_count < 20: ## when estimating scatter, exclude 2sigma outliers in an iterative way
        [RA_average, RA_std_ms] = dms_array2stat(RAs)
        [Dec_average, Dec_std_mas] = dms_array2stat(Decs)
        Dec_av_rad = Dec_average*math.pi/180
        RA_std_mas = RA_std_ms*15*math.cos(Dec_av_rad)
        RA_average_dms = howfun.deg2dms(RA_average)
        Dec_average_dms = howfun.deg2dms(Dec_average)
        sigma = (RA_std_mas**2+Dec_std_mas**2)**0.5
        RA1s = np.array([])
        Dec1s = np.array([])
        outliers = []
        for i in range(len(RAs)):
            print RAs[i], Decs[i], RA_average_dms
            sep_mas = howfun.separation(str(RAs[i]), str(Decs[i]), str(RA_average_dms), str(Dec_average_dms)) #in arcmin
            sep_mas *= 60*1000
            if sep_mas <= 3*sigma:
                RA1s = np.append(RA1s, RAs[i])
                Dec1s = np.append(Dec1s, Decs[i])
            else:
                outliers.append(i)
                outlier_count += 1
        RAs = RA1s
        Decs = Dec1s
        print outlier_count
        if len(outliers) == 0:
            break
    [RA_average, RA_std_ms] = dms_array2stat(RAs)
    [Dec_average, Dec_std_mas] = dms_array2stat(Decs)
    Dec_av_rad = Dec_average*math.pi/180
    RA_std_mas = RA_std_ms*15*math.cos(Dec_av_rad)
    RA_average_dms = howfun.deg2dms(RA_average)
    Dec_average_dms = howfun.deg2dms(Dec_average)
    print("There are %d outliers" % outlier_count)
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

def targetdir2phscalstatsfiles(targetdir, exceptions=''):
    phscalstatsfiles = glob.glob(r'%s/*/*_ibshiftdiv_difmap.jmfit.stats' % targetdir) #find statsfile for each epoch
    if exceptions != '':
        for exception in exceptions:
            while any(exception in phscalstatsfile for phscalstatsfile in phscalstatsfiles):
                for phscalstatsfile in phscalstatsfiles:
                    if exception in phscalstatsfile:
                        phscalstatsfiles.remove(phscalstatsfile)
    phscalstatsfiles.sort()
    return phscalstatsfiles
def targetdir2prIBCstatsfiles(targetdir, exceptions=''):
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.difmap.jmfit*.stats' % targetdir)
    if exceptions != '':
        for exception in exceptions:
            while any(exception in prIBCstatsfile for prIBCstatsfile in prIBCstatsfiles):
                for prIBCstatsfile in prIBCstatsfiles:
                    if exception in prIBCstatsfile:
                        prIBCstatsfiles.remove(prIBCstatsfile)
    prIBCstatsfiles.sort()
    return prIBCstatsfiles
def target2positionscatter(targetname, exceptions=''):
    """
    #####################################################################################
    ## compile stats files (for phscal and primary IBC) to get systematics
    ## tailored for MSPSRPI datasets
    #####################################################################################
    """
    ## get paths and sourcenames ########################################################
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    ## find and parse statsfiles, reach RAs/Decs #########################
    phscalstatsfiles = targetdir2phscalstatsfiles(targetdir, exceptions)
    [phscalRAs, phscalDecs] = nonpulsar_statsfiles2positions(phscalstatsfiles)
    print(phscalRAs, phscalDecs)
    prIBCstatsfiles = targetdir2prIBCstatsfiles(targetdir, exceptions)
    [prIBC_RAs, prIBC_Decs] = nonpulsar_statsfiles2positions(prIBCstatsfiles)
    print(prIBC_RAs, prIBC_Decs)
    ## statistics about RAs/Decs ##########################################
    return dms_positions2stat(phscalRAs, phscalDecs), dms_positions2stat(prIBC_RAs, prIBC_Decs)

class plot_position_scatter:
    def __init__(s, targetname, exceptions=''):
        s.targetname = targetname
        s.exceptions = exceptions
        [auxdir, configdir, s.targetdir, s.phscalname, s.prIBCname] = prepare_path_source(s.targetname)
        #s.plot_prIBC_preselfcal_scatter()
        #s.plot_phscal_scatter_interpolated_from_prIBC_and_prIBC_preselfcal_scatter()
    def plot_phscal_scatter_interpolated_from_prIBC_and_prIBC_preselfcal_scatter(s, colorbarstep=0.1):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        [RAs_phscal, Decs_phscal] = s.get_phscal_positions_interpolated_from_prIBC()
        [refRA_phscal, refDec_phscal] = srcposition(s.targetname, s.phscalname)
        [diffRAs_phscal, diffDecs_phscal] = diffposition(RAs_phscal, refRA_phscal, Decs_phscal, refDec_phscal)
        [RAs_prIBC, Decs_prIBC] = s.get_prIBC_preselfcal_positions()
        [refRA_prIBC, refDec_prIBC] = srcposition(s.targetname, s.prIBCname)
        [diffRAs_prIBC, diffDecs_prIBC] = diffposition(RAs_prIBC, refRA_prIBC, Decs_prIBC, refDec_prIBC)
        decyears = s.statsfiles2decyears(s.prIBCstatsfiles)
        decyears = decyears.astype(np.float)
        #print decyears
        #print '\n\n\n'
        
        fig = plt.figure()
        gs = gridspec.GridSpec(4, 10)
        cm = plt.cm.get_cmap('viridis')
        a1 = fig.add_subplot(gs[:4, :4])
        a1.scatter(diffRAs_prIBC, diffDecs_prIBC, c=decyears, marker='.', s=150, cmap=cm, zorder=1)
        a1.set_ylabel('relative Decl. (mas)')
        a1.set_xlabel('relative RA. (mas)')
        a1.set_title('%s' % s.prIBCname)
        ax1 = plt.gca()
        ax1.invert_xaxis()
        
        [RA_stat_prIBC, Dec_stat_prIBC] = target2positionscatter(s.targetname, s.exceptions)[1]
        for parameter in ['RA', 'Dec']:
            exec('[avg_%s_dms_prIBC, std_%s_mas_prIBC] = %s_stat_prIBC' % (parameter, parameter, parameter))
        [diff_avgRA_prIBC, diff_avgDec_prIBC] = diffposition(avg_RA_dms_prIBC, refRA_prIBC, avg_Dec_dms_prIBC, refDec_prIBC)
        #plt.plot(diff_avgRA_prIBC, diff_avgDec_prIBC, 'rx', zorder=2)
        rect = plt.Rectangle((diff_avgRA_prIBC-std_RA_mas_prIBC, diff_avgDec_prIBC-std_Dec_mas_prIBC), 
                       std_RA_mas_prIBC*2, std_Dec_mas_prIBC*2, facecolor="0.91", alpha=0.1, zorder=0)
        ax1.add_patch(rect)
        
        a2 = fig.add_subplot(gs[:4, 5:10])
        a2_scatter = a2.scatter(diffRAs_phscal, diffDecs_phscal, c=decyears, marker='.', s=150, cmap=cm, zorder=1)
        a2.set_xlabel('relative RA. (mas)')
        ax2 = plt.gca()
        ax2.invert_xaxis()
        a2.set_title('%s' % s.phscalname)
        cbar = plt.colorbar(a2_scatter)
        #cbar.ax.set_yticklabels(['{:.1f}'.format(x) for x in np.arange(math.floor(min(decyears)*10)/10., math.ceil(max(decyears)*10)/10.+colorbarstep, colorbarstep)])
        cbar.set_label('time (yr)', rotation=90)
        
        [RA_stat_phscal, Dec_stat_phscal] = target2positionscatter(s.targetname, s.exceptions)[0]
        for parameter in ['RA', 'Dec']:
            exec('[avg_%s_dms_phscal, std_%s_mas_phscal] = %s_stat_phscal' % (parameter, parameter, parameter))
        [diff_avgRA_phscal, diff_avgDec_phscal] = diffposition(avg_RA_dms_phscal, refRA_phscal, avg_Dec_dms_phscal, refDec_phscal)
        #plt.plot(diff_avgRA_phscal, diff_avgDec_phscal, 'rx', zorder=2)
        rect = plt.Rectangle((diff_avgRA_phscal-std_RA_mas_phscal, diff_avgDec_phscal-std_Dec_mas_phscal), 
                       std_RA_mas_phscal*2, std_Dec_mas_phscal*2, facecolor="0.91", alpha=0.1, zorder=0)
        ax2.add_patch(rect)
        #gs.tight_layout(fig)
        plt.savefig("%s/pmparesults/combined_scatter.eps" % s.targetdir)
        plt.clf()

    def get_phscal_positions_interpolated_from_prIBC(s):
        s.phscalstatsfiles = targetdir2phscalstatsfiles(s.targetdir, s.exceptions)
        [phscal_RAs, phscal_Decs] = nonpulsar_statsfiles2positions(s.phscalstatsfiles)
        return phscal_RAs, phscal_Decs
    def plot_prIBC_preselfcal_scatter(s, annotatesize=10, colorbarstep=0.1):
        import matplotlib.pyplot as plt
        [RAs, Decs] = s.get_prIBC_preselfcal_positions()
        expnos = statsfiles2expnos(s.prIBCstatsfiles)
        [refRA_prIBC, refDec_prIBC] = srcposition(s.targetname, s.prIBCname)
        [diffRAs, diffDecs] = diffposition(RAs, refRA_prIBC, Decs, refDec_prIBC)
        decyears = s.statsfiles2decyears(s.prIBCstatsfiles)
        decyears = decyears.astype(np.float)
        cm = plt.cm.get_cmap('viridis')
        a = plt.scatter(diffRAs, diffDecs, c=decyears, marker='.', cmap=cm, zorder=1)
        for i, expno in enumerate(expnos):
            plt.annotate(expno, (diffRAs[i], diffDecs[i]), size=annotatesize)
        cbar = plt.colorbar(a)
        cbar.ax.set_yticklabels(['{:.1f}'.format(x) for x in np.arange(math.floor(min(decyears)*10)/10., math.ceil(max(decyears)*10)/10.+colorbarstep, colorbarstep)])
        cbar.set_label('time (yr)', rotation=90)
        plt.xlabel('relative Right Ascension (mas)')
        ax1 = plt.gca()
        ax1.invert_xaxis()
        plt.ylabel('relative Declination (mas)')
        plt.title('%s positions before self-calibration' % s.prIBCname)

        [RA_stat, Dec_stat] = target2positionscatter(s.targetname, s.exceptions)[1]
        for parameter in ['RA', 'Dec']:
            exec('[avg_%s_dms, std_%s_mas] = %s_stat' % (parameter, parameter, parameter))
        [diff_avgRA, diff_avgDec] = diffposition(avg_RA_dms, refRA_prIBC, avg_Dec_dms, refDec_prIBC)
        plt.plot(diff_avgRA, diff_avgDec, 'rx', zorder=2)
        
        rect = plt.Rectangle((diff_avgRA-std_RA_mas, diff_avgDec-std_Dec_mas), 
                       std_RA_mas*2, std_Dec_mas*2, facecolor="0.91", alpha=0.5, zorder=0)
        ax1.add_patch(rect)
        
        ## rfc2018 position of J1819-2036 for comparison
        [RA_rfc2018, Dec_rfc2018] = ['18:19:36.895534', '-20:36:31.57089']
        [diffRA_rfc2018, diffDec_rfc2018] = diffposition(RA_rfc2018, refRA_prIBC, Dec_rfc2018, refDec_prIBC)
        #plt.plot(diffRA_rfc2018, diffDec_rfc2018, 'r*', zorder=1)
        plt.savefig("%s/pmparesults/prIBC_preselfcal_scatter.pdf" % s.targetdir)
        plt.clf()
        print("Plot made, where the mean position is %s %s." % (avg_RA_dms, avg_Dec_dms))

    def get_prIBC_preselfcal_positions(s): 
        s.prIBCstatsfiles = targetdir2prIBCstatsfiles(s.targetdir, s.exceptions)
        [prIBC_RAs, prIBC_Decs] = nonpulsar_statsfiles2positions(s.prIBCstatsfiles)
        return prIBC_RAs, prIBC_Decs
    def statsfiles2decyears(s, statsfiles):
        decyears = np.array([])
        for statsfile in statsfiles:
            expno = statsfile.split('/')[-2]
            if expno == '':
                expno = statsfile.split('/')[-3]
            decyear = s.expno2decyear(expno)
            decyears = np.append(decyears, decyear)
        return decyears
    def statsfile2expno_and_decyear(s, statsfile):
        expno = statsfile.split('/')[-2]
        if expno == '':
            expno = statsfile.split('/')[-3]
        decyear = s.expno2decyear(expno)
        return expno, decyear 
    def expno2decyear(s, expno):
        from astropy.time import Time
        vexfile = s.targetdir + "/" + expno + "/" + expno + ".vex"
        if not os.path.exists(vexfile):
            print vexfile + " does not exist - aborting!"
            sys.exit()
        targetname = expno2sources(expno)[0]
        startime_search_key1 = 'source=' + targetname
        startime_search_key2 = 'source=TARGETPT'
        ## read obs_time from vex file
        vexlines = open(vexfile).readlines()
        for line in vexlines:
            if 'MJD' in line:
                MJD = int(line.split(':')[-1])
            if startime_search_key1 in line or startime_search_key2 in line:
                startime=line.split(';')[0]
                break
        # translate the obs_time to decimalyear format
        startime=startime.split('d')[-1]
        starthour=float(startime[0:2])    
        startmin=float(startime[3:5])
        print starthour,startmin
        dayfrac = (startmin/60+starthour)/24
        MJD = float(MJD) + dayfrac
        #mjds = np.append(mjds,MJD)
        #print mjds
        #print "average MJD is: %f" % np.average(mjds)
        MJD = Time(MJD,format='mjd')
        decyear = format(MJD.decimalyear,'.4f')
        return decyear 

def targetname2decyears(targetname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    tmpfile = targetdir + '/.variable_decyears.tmp'
    if not os.path.exists(tmpfile):
        os.system('generatepmparin.py -t %s' % targetname)
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
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.difmap.jmfit*.stats' % targetdir)
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
    rchsq = 0
    lines = open(pmparout).readlines()
    for line in lines:
        if 'epoch' in line:
            epoch = line.split('=')[1].strip()
        if 'Reduced' in line:
            rchsq = float(line.split('=')[1].strip())
        for estimate in ['l', 'b']:
            if (estimate in line) and 'degrees' in line:
                exec("%s = %s" % (estimate, line.split('=')[-1].strip().split(' ')[0]))
        for estimate in ['RA', 'Dec  ', 'mu_a', 'mu_d', 'pi']:
            if estimate in line:
                print line.split('=')[-1].split('+')[0].strip()
                exec("%s = line.split('=')[-1].split('+')[0].strip()" % estimate.strip())
                if estimate in ['RA', 'Dec  ']:
                    exec("%s = howfun.dms2deg(%s)" % (estimate.strip(), estimate.strip()))
    for line in lines:
        for estimate in ['mu_a', 'mu_d', 'pi']:
            if estimate in line:
                error = line.split('+-')[1].strip().split(' ')[0]
                exec("error_%s = %s" % (estimate, error))
                exec("%s = float(%s)" % (estimate, estimate))
    return RA, Dec, epoch, pi, mu_a, mu_d, error_pi, error_mu_a, error_mu_d, l, b, rchsq
def read_bayesian_position_and_error(targetname):
    """
    read from the output of the astrometryfit code written by Adam
    """
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    bayesian_results = targetdir + '/bayesian_results/' + targetname + '.bayesian.results'
    if not os.path.exists(bayesian_results):
        print('%s does not exist; aborting' % bayesian_results)
        sys.exit()
    lines = open(bayesian_results).readlines()
    for line in lines:
        if ' RA:' in line:
            contents = line.split(' RA:')[-1].strip().split('(+')
            RA = float(contents[0].strip())
            RA_err_two_sides = contents[1].split(')')[0].split('/ -')
            RA_err_up = float(RA_err_two_sides[0].strip())
            RA_err_down = float(RA_err_two_sides[1].strip())
            RA_err = max(RA_err_up, RA_err_down)
        if ' DEC:' in line:
            contents = line.split(' DEC:')[-1].strip().split('(+')
            DEC = float(contents[0].strip())
            DEC_err_two_sides = contents[1].split(')')[0].split('/ -')
            DEC_err_up = float(DEC_err_two_sides[0].strip())
            DEC_err_down = float(DEC_err_two_sides[1].strip())
            DEC_err = max(DEC_err_up, DEC_err_down)
    err = np.array([RA_err*np.cos(DEC), DEC_err]) #in rad
    err *= 180./np.pi*3600*1000 #in mas
    position = np.array([RA/15., DEC])
    position *= 180./np.pi #in hr and deg
    return position, err

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
            #print(errRAs)
        if 'Dec' in line:
            errDec = line.split('+-')[-1].split('(mas)')[0].strip()
            errDecs = np.append(errDecs, float(errDec))
            #print(errDecs)
    if not HowManySigma in [1,2]:
        print('the second parameter should be either 1 or 2 (sigma); aborting')
        sys.exit()
    errRA = errRAs[3-HowManySigma]
    errDec = errDecs[3-HowManySigma]
    err = np.array([errRA, errDec])
    return err
def abspsrposition(targetname, HowManySigma):
    """
    Notice for use:
    The function is designed to estimate the absolute position for the target anchored to a primary in-beam calibrator, 
    and indirectly to a main phase calibrator. This works for all PSRPI and MSRPI targets. However, it doesn't apply to other
    observing setup.
    Measured value: 
    This function adopts the position estimate that is kept by the pmpar.out file (not from bootstrap, but almost the same).
    Then it will be shifted when the target is tied to the main phscal, and sebsequently aligned to the catalog phscal position.
    Uncertainty: 
    This function reports two sets of results (scatter of positions for prIBC or phscal). The uncertainties 
    include the bootstrap uncertainty, prIBC position uncertainty and catalog uncertainty for the main phscal.
    """
    [phscalstats, prIBCstats] = target2positionscatter(targetname)
    phscalRADEC = [phscalstats[0][0], phscalstats[1][0]]
    phscalRADECstd = np.array([phscalstats[0][1], phscalstats[1][1]])
    prIBC_RADEC = [prIBCstats[0][0], prIBCstats[1][0]]
    prIBC_RADECstd = np.array([prIBCstats[0][1], prIBCstats[1][1]])
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    [psrRA, psrDec, epoch, junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8, junk9] = readpulsition(targetname)
    psrRADEC0 = [psrRA, psrDec]
    bayesian_results = targetdir + '/bayesian_results/' + targetname + '.bayesian.results'
    if os.path.exists(bayesian_results):
        [psrRADEC0, err_psr2prIBC] = read_bayesian_position_and_error(targetname)
    rfcinfos = rfcposition(phscalname)
    phscal_absRADEC = [rfcinfos[0], rfcinfos[2]]
    errRAphscal = float(rfcinfos[1]) #in mas
    errDECphscal = float(rfcinfos[3]) #in mas
    ## 1) using prIBC shift ####################
    prIBCrefRADEC = srcposition(targetname, prIBCname)
    print(psrRADEC0, prIBC_RADEC, prIBCrefRADEC)
    psrRADEC1 = psrRADEC0 + howfun.dms2deg(prIBC_RADEC) - howfun.dms2deg(prIBCrefRADEC)
    ## 2) using phscal shift ###################
    phscal_refRADEC = srcposition(targetname, phscalname)
    psrRADEC2 = psrRADEC0 + howfun.dms2deg(phscal_refRADEC) - howfun.dms2deg(phscalRADEC)
    ## 3) update absolute position for phscal #####
    psrRADEC2 = psrRADEC2 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
    psrRADEC2 = howfun.deg2dms(psrRADEC2)
    psrRADEC1 = psrRADEC1 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
    psrRADEC1 = howfun.deg2dms(psrRADEC1)
    ## error estimation #####################################################
    if not os.path.exists(bayesian_results):
        err_psr2prIBC = bootstrapRADECerr(targetname, HowManySigma)
    err_abs_phscal = np.array([errRAphscal, errDECphscal])
    err1 = (err_psr2prIBC**2 + err_abs_phscal**2 + prIBC_RADECstd**2)**0.5 #in mas
    print(err_psr2prIBC, err_abs_phscal, prIBC_RADECstd)
    err2 = (err_psr2prIBC**2 + err_abs_phscal**2 + phscalRADECstd**2)**0.5 #in mas
    return psrRADEC1, err1, psrRADEC2, err2, epoch
def exportabspsrposition(targetname, HowManySigma):
    [psrRADEC1, err1, psrRADEC2, err2, epoch] = abspsrposition(targetname, HowManySigma)
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    err1[0] = howfun.mas2ms(err1[0], str(psrRADEC1[1]))
    err2[0] = howfun.mas2ms(err2[0], str(psrRADEC2[1]))
    err1 = err1/1000 #in s
    err2 = err2/1000 #in arcsecond
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

## group 2a #############################################################################################
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

def modeltype(targetname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    vexfiles = glob.glob(r'%s/*/*.vex' % targetdir)
    vexfiles.sort()
    vexfile = vexfiles[0]
    expno = vexfile.split('/')[-1].split('.')[0].strip().lower()
    expconfigfile = configdir + '/' +  expno + '.yaml'
    if not os.path.exists(expconfigfile):
        print("%s does not exist; aborting\n" % expconfigfile)
        sys.exit()
    expconfig = yaml.load(open(expconfigfile))
    userno = expconfig['userno']
    useprelimmodels = expconfig['useprelimmodels']
    if useprelimmodels:
        modeltype = 'preliminary'
    else:
        modeltype = 'final'
    return modeltype, userno
def jmfitfromfile(targetname, srcname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    outputstatsfile = targetdir + '/.' + srcname + '.stats'
    if not os.path.exists(outputstatsfile):
        [modeltype1, junk1] = modeltype(targetname)
        srcmodel = auxdir + '/sourcemodels/' + modeltype1 + '/' + srcname + '.clean.fits'
        if not os.path.exists(srcmodel):
            print("%s does not exist; aborting\n" % srcmodel)
            sys.exit()
        os.system("jmfitfromfile.py %s %s" % (srcmodel, outputstatsfile))
    return outputstatsfile
def srcposition(targetname, srcname):
    if targetname != srcname:
        statsfile = jmfitfromfile(targetname, srcname)
        print("The position corresponding to the brightest spot on the map is:\n")
        return nonpulsar_statsfile2position(statsfile)
    else:
        print("The map center for the target is:\n")
        return srcposition_map_center(targetname, srcname)
def srcposition_map_center(targetname, srcname):
    [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
    tmpfile = targetdir + '/.' + srcname + '.map.center.position.tmp'
    if not os.path.exists(tmpfile):
        os.system('srcposition_map_center.py %s %s' % (targetname, srcname))
    readfile = open(tmpfile, 'r')
    position = pickle.load(readfile)
    readfile.close()
    print position
    return position
def srcseparation(targetname, src1, src2):
    [RA1, Dec1] = srcposition(targetname, src1)
    [RA2, Dec2] = srcposition(targetname, src2)
    sep = howfun.separation(RA1,Dec1,RA2,Dec2) #unit: arcmin
    print sep
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

## group 2b: epoch-based sysErr estimator #########################################################
class expno_sysErr:
    paraA = 0.001
    paraB = 0.6
    def __init__(s, expno, dualphscal='', paraA_rchsq=1e-4, inverse_referencing=False):
        s.inverse_referencing = inverse_referencing
        [s.targetname, s.foldername, s.phsrefname, s.prIBCname] = expno2sources(expno, s.inverse_referencing)
        s.expno = expno
        auxdir    = os.environ['PSRVLBAUXDIR']
        targetdir = auxdir + '/processing/' + s.foldername
        s.expdir = targetdir + '/' + expno
        s.paraA_rchsq = paraA_rchsq
        s.dualphscal = dualphscal
        s.dodefaultnames = s.read_dodefaultnames_given_expno(s.expno)
    def sysErr(s):
        [beamRA, beamDec] = s.beam_in_RA_Dec()
        delta_sys1 = s.delta_sys()    
        sysErrRA  = beamRA*delta_sys1 #in mas
        sysErrDec = beamDec*delta_sys1 #in mas
        sysErrRA_ms = howfun.mas2ms(sysErrRA, s.Dec_target)
        return sysErrRA, sysErrDec, sysErrRA_ms
    def delta_sys(s):
        sep = s.separation_between_target_and_virtual_cal_in_dualphscal_mode()
        #sep = s.target_prIBC_separation()
        Av_cscEl1 = s.Av_cscEl()
        [junk1, junk2, junk3, SNprIBC] = s.targetbeam()
        if s.dualphscal:
            delta_sys = s.paraA_rchsq*sep*Av_cscEl1 ## since not using IBC and both phscals are relatively bright
        else:
            delta_sys = s.paraA*sep*Av_cscEl1 + s.paraB/SNprIBC
        return delta_sys
    def target_prIBC_separation(s):
        [RA1, Dec1] = srcposition(s.foldername, s.foldername)
        s.Dec_target = Dec1
        [RA2, Dec2] = srcposition(s.foldername, s.prIBCname)
        sep = howfun.separation(RA1,Dec1,RA2,Dec2) #unit: arcmin
        return sep
    def beam_in_RA_Dec(s):
        [beamPA, beamSA, beamLA, junk1] = s.targetbeam()
        [beamRA, beamDec] = howfun.deprojectbeam2xy(beamLA,beamSA,beamPA) #full-width deprojection on RA/Dec from beam
        return beamRA, beamDec
    def read_dodefaultnames_given_expno(s, expno):
        auxdir    = os.environ['PSRVLBAUXDIR']
        configdir = auxdir + '/configs/'
        expconfigfile = configdir + '/' + expno + '.yaml'
        expconfig = yaml.load(open(expconfigfile))
        return expconfig['dodefaultnames']
    def targetbeam(s): #also get SNprIBCs, using prIBC statsfile,
        """
        Note
        ----
        1. statsfile for divided IBC fitsfile is not used, as the image S/N normally increases after dividing the model.
        """
        if not s.dodefaultnames:
            prIBCstatsfiles = glob.glob(r'%s/*%s.difmap.jmfit.stokesi.stats' % (s.expdir, s.prIBCname))
        else:
            prIBCstatsfiles = glob.glob(r'%s/*_inbeam0_0_divided.difmap.jmfit.stokesi.stats' % s.expdir)
        print prIBCstatsfiles
        if len(prIBCstatsfiles) != 1:
            if s.prIBCname == s.phsrefname:
                prIBCstatsfiles = glob.glob(r'%s/*%s_ibshiftdiv_difmap.jmfit.stats' % (s.expdir, s.prIBCname))
                if len(prIBCstatsfiles) != 1:
                    print "No statsfile or more than one statsfiles for the primary in-beam calibrator is found."
                    sys.exit()
        prIBCstatsfile = prIBCstatsfiles[0]
        lines = open(prIBCstatsfile).readlines()[-10:]
        for line in lines:
            if 'S/N' in line:
                SNprIBC = line.split(':')[-1].strip()
                SNprIBC = float(SNprIBC)
            if 'beam' in line:
                line = line.split('beam')[-1].strip().split(' ')
                beamPA   = float(line[-2])
                beamsize = line[0]
                beamLA   = float(beamsize.split('x')[-1])
                beamSA   = float(beamsize.split('x')[0])
        return beamPA, beamSA, beamLA, SNprIBC
    def Av_cscEl(s):
        sumfiles = glob.glob(r'%s/*.sum' % s.expdir)
        if len(sumfiles) != 1:
            print "There is no or more than one sumfiles found; abort"
            sys.exit()
        sumfile = sumfiles[0]
        TelAv_csc_Els = np.array([])
        targetname = s.targetname
        if not s.targetname[-1].isdigit():
            targetname = s.targetname[:-1]
        lines = open(sumfile).readlines()
        if s.dodefaultnames:
            keyword = 'TARGETPT'
        else:
            keyword = targetname
        for line in lines:
            if keyword in line:
                #if ':' in line.split(s.targetname)[0] and howfun.no_alphabet(line.split(s.targetname)[-1]):
                if ':' in line.split(keyword)[0]:
                    try:
                        elevations = line.split('-')[-1].strip().split('    ')
                        elevations = map(float, elevations)
                    except ValueError:
                        try:
                            elevations = line.split('-')[-1].strip().split('   ')
                            elevations = map(float, elevations)
                        except ValueError:
                            continue
                    elevations = np.asarray(elevations)
                    elevations_rad = elevations*math.pi/180
                    csc_elevations = (np.sin(elevations_rad))**(-1)
                    TelAv_csc_El = np.average(csc_elevations) #telescope averaged csc(el) for each scan
                    TelAv_csc_Els = np.append(TelAv_csc_Els, TelAv_csc_El)
        Av_cscEl = np.average(TelAv_csc_Els) #averaged csc(el) over telescopes and scans
        return Av_cscEl
    def separation_between_target_and_virtual_cal_in_dualphscal_mode(s):
        [RA1, Dec1] = srcposition(s.foldername, s.foldername)
        s.Dec_target = Dec1
        if s.dualphscal and s.phsrefname != s.prIBCname:
            a = find_virtual_calibrator_position_with_colinear_calibrators(s.targetname, s.phsrefname, s.prIBCname) 
            [sep, inbeamsn_ratio] = a.separation_between_virtual_phscal_and_sources()
        else:
            print("Not a dual-phscal configuration; go with normal separation")
            sep = s.target_prIBC_separation() #correct it back after test
        return sep

# group 3: draw online absolute position and errors ###############################################
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
    return RA, word1[3], Dec, word1[4] #uncertainties in mas
def rfcposition(phscal):
    auxdir = os.environ['PSRVLBAUXDIR']
    downloaddir = auxdir + '/downloads/'
    year = datetime.datetime.today().year
    lastyear = str(int(year) - 1)
    year = str(year)
    versions = [year+'d', year+'c', year+'b', year+'a', lastyear+'d', lastyear+'c', lastyear+'b']
    for version in versions:
        filename = 'rfc_' + version + '_cat.txt'
        #ftpfile = 'http://astrogeo.org/vlbi/solutions/rfc_' + version + '/' + filename
        ftpfile = 'http://astrogeo.org/sol/rfc/rfc_' + version + '/' + filename #updated in May 2021
        downloadfile = downloaddir + filename
        if not os.path.exists(downloadfile):
            os.system('wget -P %s %s' % (downloaddir, ftpfile))
        if os.path.exists(downloadfile):
            print("%s already exists" % filename)
            break
    separation_line = 80*'|'+'\n'
    print("%srfc%s position is used for %s. Double check the catalog position!\n%s" % (separation_line, version, phscal, separation_line))
    return parse_rfctxt(phscal, downloadfile)

## align position for target to the same reference position in ad-hoc experiment
def align_position_in_adhoc_experiment(RA, errRA, Dec, errDec, targetname, exceptions, dualphscal=False, dualphscalratio=1): #str, float in s, str, float in "...
    """
    1) This is a combination of the align_position_in_adhoc_experiment0/1. 
    2) It chooses either prIBC relative to phscal or phscal to prIBC with smaller scatter in a semi-automatic way.
    3) 'refRA/Dec' is the average position of the calibrator positions.
    """
    ## ref position (and its scatter) of the phsref Cal relative to prIBC for new observations
    [refRA0, refDec0] = target2positionscatter(targetname, exceptions)[0]
    [refRA0, errRefRA0_mas] = refRA0
    [refDec0, errRefDec0_mas] = refDec0
    ## ref position (and its scatter) of the prIBC relative to phscal for new observations
    [refRA1, refDec1] = target2positionscatter(targetname, exceptions)[1]
    [refRA1, errRefRA1_mas] = refRA1
    [refDec1, errRefDec1_mas] = refDec1
    ## find the set of positions with smaller scatter
    if errRefRA0_mas < errRefRA1_mas and errRefDec0_mas < errRefDec1_mas:
        which_scatter = 0
        refRA = refRA0
        refDec = refDec0
        errRefRA_mas = errRefRA0_mas
        errRefDec_mas = errRefDec0_mas
    elif errRefRA0_mas > errRefRA1_mas and errRefDec0_mas > errRefDec1_mas:
        which_scatter = 1
        refRA = refRA1
        refDec = refDec1
        errRefRA_mas = errRefRA1_mas
        errRefDec_mas = errRefDec1_mas
    else:
        print('The RA/Dec scatter of phscal relative to prIBC is %f/%f' % (errRefRA0_mas, errRefDec0_mas))
        print('The RA/Dec scatter of prIBC relative to phscal is %f/%f' % (errRefRA1_mas, errRefDec1_mas))
        which_scatter = raw_input("Which set of scatter do you want to use (0 for phscal, 1 for prIBC, q to quit): ")
        while which_scatter not in ['0', '1', 'q']:
            which_scatter = raw_input("Please choose a set of scatter you want to use -- 0 for phscal and 1 for prIBC (q to quit): ")
        if which_scatter == 'q':
            print('you will quit for now')
            sys.exit()
        which_scatter = int(which_scatter)
        exec('refRA = refRA%d' % which_scatter)
        exec('refDec = refDec%d' % which_scatter)
        exec('errRefRA_mas = errRefRA%d_mas' % which_scatter)
        exec('errRefDec_mas = errRefDec%d_mas' % which_scatter)
    ## here begins to add the error and move the positions
    errRefRA_ms = howfun.mas2ms(errRefRA_mas, refDec)
    errRefDec_as = errRefDec_mas/1000
    errRefRA_s = errRefRA_ms/1000
    if dualphscal:
        errRefRA_s *= dualphscalratio
        errRefDec_as *= dualphscalratio
    errRA2 = (errRA**2+errRefRA_s**2)**0.5 # in s
    errDec2 = (errDec**2+errRefDec_as**2)**0.5 # in "
    # now get the reference position of phsref Cal for old data 
    #expno = exceptions[-1]
    #[junk1, junk2, phscalname, prIBCname] = expno2sources(expno)
    #if phscalname == prIBCname:
    #    pass
    calibratornames = prepare_path_source(targetname)
    calibratorname = calibratornames[which_scatter+3]
    [RA1calibrator, Dec1calibrator] = srcposition(targetname, calibratorname) # this unfortunately only works when old/new observations share the same name of phscalname in data
    #[diffRA_mas, diffDec_mas] = diffposition(refRA, RA1calibrator, refDec, Dec1calibrator) 
    if which_scatter == 1: #prIBC relative to phscal
        [diffRA_mas, diffDec_mas] = diffposition(RA1calibrator, refRA, Dec1calibrator, refDec)
        if dualphscal:
            [diffRA_mas, diffDec_mas] = np.array(diffposition(RA1calibrator, refRA, Dec1calibrator, refDec))*dualphscalratio
    else: # which_scatter == 0 -- phscal relative to prIBC
        [diffRA_mas, diffDec_mas] = diffposition(refRA, RA1calibrator, refDec, Dec1calibrator)
        if dualphscal:
            [diffRA_mas, diffDec_mas] = np.array(diffposition(refRA, RA1calibrator, refDec, Dec1calibrator))*dualphscalratio
    diffRA_ms = howfun.mas2ms(diffRA_mas, refDec) 
    RA2 = howfun.dms2deg(RA) + diffRA_ms/1000./3600.
    RA2 = howfun.deg2dms(RA2)
    Dec2 = howfun.dms2deg(Dec) + diffDec_mas/1000./3600.
    Dec2 = howfun.deg2dms(Dec2)
    return RA2, errRA2, Dec2, errDec2
def align_position_in_adhoc_experiment1(RA, errRA, Dec, errDec, targetname, exceptions, dualphscal=False, dualphscalratio=1): #str, float in s, str, float in "...
    # ref position (and its scatter) of the prIBC relative to phscal for new observations
    [refRA1, refDec1] = target2positionscatter(targetname, exceptions)[1]
    [refRA, errRefRA_mas] = refRA1
    [refDec, errRefDec_mas] = refDec1
    #refRA=RA
    #errRefRA_mas=0.1
    #refDec=Dec
    #errRefDec_mas=0.1
    errRefRA_ms = howfun.mas2ms(errRefRA_mas, refDec)
    errRefDec_as = errRefDec_mas/1000
    errRefRA_s = errRefRA_ms/1000
    if dualphscal:
        errRefRA_s *= dualphscalratio
        errRefDec_as *= dualphscalratio
    errRA2 = (errRA**2+errRefRA_s**2)**0.5 # in s
    errDec2 = (errDec**2+errRefDec_as**2)**0.5 # in "
    # now get the reference position of phsref Cal for old data 
    #expno = exceptions[-1]
    #[junk1, junk2, phscalname, prIBCname] = expno2sources(expno)
    #if phscalname == prIBCname:
    #    pass
    [junk1, junk2, junk3, phscalname, prIBCname] = prepare_path_source(targetname)
    [RA1prIBC, Dec1prIBC] = srcposition(targetname, prIBCname) # this unfortunately only works when old/new observations share the same name of phscalname in data
    #[diffRA_mas, diffDec_mas] = diffposition(refRA, RA1prIBC, refDec, Dec1prIBC) 
    [diffRA_mas, diffDec_mas] = diffposition(RA1prIBC, refRA, Dec1prIBC, refDec)
    if dualphscal:
        [diffRA_mas, diffDec_mas] = np.array(diffposition(RA1prIBC, refRA, Dec1prIBC, refDec))*dualphscalratio
    diffRA_ms = howfun.mas2ms(diffRA_mas, refDec) 
    RA2 = howfun.dms2deg(RA) + diffRA_ms/1000./3600.
    RA2 = howfun.deg2dms(RA2)
    Dec2 = howfun.dms2deg(Dec) + diffDec_mas/1000./3600.
    Dec2 = howfun.deg2dms(Dec2)
    return RA2, errRA2, Dec2, errDec2
def align_position_in_adhoc_experiment0(RA, errRA, Dec, errDec, targetname, exceptions): #str, float in s, str, float in "...
    # ref position (and its scatter) of the phsref relative to prIBC for new observations
    [refRA1, refDec1] = target2positionscatter(targetname, exceptions)[0]
    [refRA, errRefRA_mas] = refRA1
    [refDec, errRefDec_mas] = refDec1
    #refRA=RA
    #errRefRA_mas=0.1
    #refDec=Dec
    #errRefDec_mas=0.1
    errRefRA_ms = howfun.mas2ms(errRefRA_mas, refDec)
    errRefDec_as = errRefDec_mas/1000
    errRefRA_s = errRefRA_ms/1000
    errRA2 = (errRA**2+errRefRA_s**2)**0.5 # in s
    errDec2 = (errDec**2+errRefDec_as**2)**0.5 # in "
    # now get the reference position of phsref Cal for old data 
    expno = exceptions[0]
    [junk1, junk2, phscalname, junk3] = expno2sources(expno)
    [RA1phscal, Dec1phscal] = srcposition(targetname, phscalname) # this unfortunately only works when old/new observations share the same name of phscalname in data
    [diffRA_mas, diffDec_mas] = diffposition(refRA, RA1phscal, refDec, Dec1phscal) 
    diffRA_ms = howfun.mas2ms(diffRA_mas, refDec) 
    RA2 = howfun.dms2deg(RA) + diffRA_ms/1000/3600
    RA2 = howfun.deg2dms(RA2)
    Dec2 = howfun.dms2deg(Dec) + diffDec_mas/1000/3600
    Dec2 = howfun.deg2dms(Dec2)
    return RA2, errRA2, Dec2, errDec2

class estimate_uncertainty:
    """
    estimate uncertainty with direct-fitting uncertainties
    """
    yr2d = 365.242199
    A = 1./1000/3600/180*math.pi*(1./yr2d/24/3600) #mas/yr to rad/s
    A = A*3.0857*10**16 # rad/s*kpc -> km/s, altogether mas/yr*kpc ->km/s
    pc2ly = 3.26163344
    pc2m = 3.08567758e16
    
    GdotG = 7.1e-14 #Gdot/G by Hofmann et al. 2010
    errGdotG = 7.6e-14 #1sigma error for a
    q = 10.44 #m_p/m_c by Daniel et al. in prep.
    errQ = 0.11
    m_c = 0.174 #m_c/m_sun by Antoniadis et al. 2016
    errM_c = 0.011
    T = constants.G*AC.M_sun.value/constants.c**3

    def __init__(s, targetname='', timingfityear=''): #s -> self
        s.targetname, s.timingfityear = targetname, timingfityear
    def workflow():
        """
        to initiate the code
        """
        #print readpulsition(targetname)
        [s.RA, s.Dec, s.epoch, s.pi, s.mu_a, s.mu_d, s.error0_pi, s.error0_mu_a, 
        s.error0_mu_d, s.l, s.b, funk1] = readpulsition(s.targetname)
        if s.timingfityear != '':
            if type(s.timingfityear) != str:
                s.timingfityear = str(s.timingfityear)
            print '\n' + s.timingfityear
            [s.epochTm, s.DM, s.Pb, s.estimatesTm, s.errorsTm] = s.readtimingfit(s.timingfityear)
            [s.RATm, s.DecTm, s.piTm, s.mu_aTm, s.mu_dTm] = s.estimatesTm
            [s.error_RATm, s.error_DecTm, s.error_piTm, s.error_mu_aTm, s.error_mu_dTm] = s.errorsTm
            print s.readtimingfit(s.timingfityear)
    def readtimingfit(s, timingfityear):
        [junk1,junk2,targetdir,junk3,junk4] = prepare_path_source(s.targetname)
        timingfitfile = targetdir + '/pmparesults/timing_results/' + str(timingfityear) + 'timingfit'
        estimates = np.array([])
        errors = np.array([])
        lines = open(timingfitfile).readlines()
        for line in lines:
            for estimate in ['epoch', 'DM', 'Pb']:
                if estimate in line:
                    line1 = line.split('=')[-1].strip().split(' ')[0].strip()
                    exec("%s = %s" % (estimate, line1))
            for estimate in ['RA', 'Dec']:
                if estimate in line:
                    line1 = line.split('=')[-1].strip()
                    exec("%s = line1.split('+-')[0].strip()" % estimate)
                    exec("err_%s = float(line1.split('+-')[1].strip())/3600" % estimate)
                    exec("estimates = np.append(estimates, howfun.dms2deg(%s))" % estimate)
                    exec("errors = np.append(errors, err_%s)" % estimate)
            for estimate in ['pi', 'mu_a', 'mu_d', 'P_bdot']:
                if estimate in line:
                    line1 = line.split('=')[-1].strip()
                    if not 'P_bdot' in line:
                        exec("%s = %s" % (estimate, line1.split('+-')[0].strip()))
                        exec("err_%s = float(line1.split('+-')[1].strip().split(' ')[0].strip())" % estimate)
                        exec("estimates = np.append(estimates, float(%s))" % estimate)
                        exec("errors = np.append(errors, err_%s)" % estimate)
                    else:
                        PbObs = line.split('=')[-1].strip()
                        s.PbObs = float(PbObs.split('+-')[0].strip())
                        s.errPbObs = float(PbObs.split('+-')[-1].strip())
        return epoch, DM, Pb, estimates, errors
    def transverse_velocity1(s, use_other_result=True): #use one-side difference
        if not use_other_result:
            return s.transverse_velocity1a(s.mu_a, s.mu_d, s.pi, s.error0_mu_a, s.error0_mu_d, s.error0_pi)
        if use_other_result:
            return s.transverse_velocity1a(s.mu_aTm, s.mu_dTm, s.piTm, s.error_mu_aTm, s.error_mu_dTm, s.error_piTm)
    def transverse_velocity1a(s, mu_a, mu_d, pi, error_mu_a, error_mu_d, error_pi):
        v0_t = (mu_a**2 + mu_d**2)**0.5/pi*s.A
        err0_vt_pi = (mu_a**2 + mu_d**2)**0.5/(pi - error_pi)*s.A - v0_t
        err0_vt_mu_a = ((mu_a+error_mu_a)**2 + mu_d**2)**0.5/pi*s.A - v0_t
        err0_vt_mu_d = (mu_a**2 + (mu_d+error_mu_d)**2)**0.5/pi*s.A - v0_t
        err0_vt = (err0_vt_pi**2 + err0_vt_mu_a**2 + err0_vt_mu_d**2)**0.5
        return v0_t, err0_vt
    def transverse_velocity2(s, use_other_result=True):
        if not use_other_result:
            return s.transverse_velocity2a(s.mu_a, s.mu_d, s.pi, s.error0_mu_a, s.error0_mu_d, s.error0_pi)
        if use_other_result:
            return s.transverse_velocity2a(s.mu_aTm, s.mu_dTm, s.piTm, s.error_mu_aTm, s.error_mu_dTm, s.error_piTm)
    def transverse_velocity2a(s, mu_a, mu_d, pi, error_mu_a, error_mu_d, error_pi):
        v0_t = (mu_a**2 + mu_d**2)**0.5/pi*s.A
        ratio_square = (error_pi/pi)**2
        ratio_square += (error_mu_a/mu_a)**2/(1+(mu_d/mu_a)**2)**2
        ratio_square += (error_mu_d/mu_d)**2/(1+(mu_a/mu_d)**2)**2
        ratio = ratio_square**0.5
        err0_vt =v0_t * ratio
        return v0_t, err0_vt
    def uncertainty_PbShk(s, howmanysigma=1, usetimingfit=True):
        CL = math.erf(howmanysigma/2**0.5) #CL: 3->0.9973, 2->0.9545, 1->0.6827
        if not usetimingfit: #use timing Pb and VLBI mu_a, mu_d and pi
            pass
        if usetimingfit:
            errorsTm = s.errorsTm*howmanysigma #Tm --> timing
            [error_RATm, error_DecTm, error_piTm, error_mu_aTm, error_mu_dTm] = errorsTm
            DTm = 1/s.piTm*1000*s.pc2ly #kpc ->ly
            error_piTm = error_piTm/1000/s.pc2ly
            muTm = np.array([s.mu_aTm, error_mu_aTm, s.mu_dTm, error_mu_dTm])
            muTm = muTm/1000./3600./180*math.pi #mas/yr -> rad/yr
            [mu_aTm, error_mu_aTm, mu_dTm, error_mu_dTm] = muTm
            Pb = s.Pb/s.yr2d 
            PbShk = (mu_aTm**2+mu_dTm**2)*DTm*Pb
            errPbShk = (2*mu_aTm*Pb*DTm*error_mu_aTm)**2
            errPbShk += (2*mu_dTm*Pb*DTm*error_mu_dTm)**2
            errPbShk += ((mu_aTm**2+mu_dTm**2)*Pb*error_piTm*DTm**2)**2
            errPbShk = errPbShk**0.5
            print 'PbShk = %f +- %f' % (PbShk, errPbShk)
            return PbShk, errPbShk
    def uncertainty_PbGal(s, howmanysigma=1, usetimingfit=True):
        if not usetimingfit: #use timing Pb and VLBI mu_a, mu_d and pi
            pass
        if usetimingfit:
            errorsTm = s.errorsTm*howmanysigma
            [error_RATm, error_DecTm, error_piTm, error_mu_aTm, error_mu_dTm] = errorsTm
            b = s.b/180*math.pi
            l = s.l/180*math.pi
            DTm = 1/s.piTm #kpc
            error_DTm = DTm*error_piTm/s.piTm
            PbGal1 = s.PbGal_vertical(DTm, b, s.Pb)
            PbGal1_plus = s.PbGal_vertical(DTm+error_DTm, b, s.Pb) - PbGal1
            PbGal1_minus = PbGal1 - s.PbGal_vertical(DTm-error_DTm, b, s.Pb)
            print PbGal1, PbGal1_minus, PbGal1_plus
            errPbGal1 = max(PbGal1_minus, PbGal1_plus)
            
            #now the second term, should be A2, simplified to A
            Theta = 233.34 #km/s; McGaugh, 2018
            error_Theta = 1.40 #km/s
            R = 8.122 #kpc; GRAVITY Collaboration, 2018
            error_R = 0.031 #kpc
            beta = DTm/R*math.cos(b) - math.cos(l)
            div1 = (math.sin(l))**2+beta**2
            A_G = -math.cos(b)*Theta**2/R*(math.cos(l)+beta/div1) #km^2/s^2/kpc
            A_G = A_G*10e6/10e3/s.pc2m #m/s^2
            #print A_G
            PbGal2 = A_G/constants.c*s.Pb #d/s
            PbGal2 = PbGal2*24*3600
            print PbGal2
            dA_dTheta = 2/Theta*A_G #m/s^2/(km/s)
            dBeta_dR = -DTm/R**2*math.cos(b) #1/kpc
            div2 = (math.sin(l))**2-beta**2
            dA_dR = -math.cos(b)*Theta**2/R*dBeta_dR*div2/div1**2 #km^2/s^2/kpc/kpc
            dA_dR = dA_dR*10e6/10e3/s.pc2m #m/s^2/kpc
            dA_dR = dA_dR - A_G/R #m/s^2/kpc
            dBeta_dD = math.cos(b)/R #1/kpc
            dA_dD = -math.cos(b)*Theta**2/R*dBeta_dD*div2/div1**2 #km^2/s^2/kpc/kpc
            dA_dD = dA_dD*10e6/10e3/s.pc2m #m/s^2/kpc
            errA_sq = (dA_dTheta*error_Theta)**2
            errA_sq += (dA_dR*error_R)**2
            errA_sq += (dA_dD*error_DTm)**2
            errA = errA_sq**0.5 #m/s^2
            errPbGal2 = errA/constants.c*s.Pb*24*3600 #s/s
            print "Pbdot_Gal1 = %f +- %f (fs/s)\nPbdot_Gal2 = %f +- %f (fs/s)" % (PbGal1*1e15, errPbGal1*1e15, PbGal2*1e15, errPbGal2*1e15)
            PbGal = PbGal1 + PbGal2
            errPbGal = (errPbGal1**2+errPbGal2**2)**0.5
            return PbGal, errPbGal

    def PbGal_vertical(s, D, b, Pb): # b in rad, D in kpc
        z = abs(math.sin(b))*D #kpc
        Kz = 2.27*z + 3.68*(1-math.exp(-4.31*z)) #10pm/s^2
        PbGal1_r = -Kz*abs(math.sin(b))/10**11/constants.c #the vertical term of Pbdot^{Gal}/Pb, in /s
        PbGal1 = PbGal1_r*Pb*24*3600
        return PbGal1
    
    def PbGW(s):
        q = s.q
        errQ = s.errQ
        m_c = s.m_c
        errM_c = s.errM_c
        T = s.T
        Pb = s.Pb
        C = -192*math.pi/5*T**(5./3)
        C *= (2*math.pi/Pb)**(5./3) # in (s/d)**(5/3)
        C *= (1./24/3600)**(5./3) # in 1
        PbGW = C*m_c**(5./3)*q/(q+1)**(1./3)
        dPdM = 5./3*PbGW/m_c
        dPdQ = 1./q-1./3/(q+1)
        dPdQ *= PbGW
        errPbGW_sq = (dPdM*errM_c)**2 + (dPdQ*errQ)**2
        errPbGW = errPbGW_sq**0.5
        return PbGW, errPbGW

    def PbEx(s, howmanysigma=1, usetimingfit=True):
        PbObs = s.PbObs #Desvignes et al. 2016
        errPbObs = s.errPbObs
        [PbGW, errPbGW] = s.PbGW()
        if usetimingfit:
            PbGal, errPbGal = s.uncertainty_PbGal(howmanysigma) 
            PbShk, errPbShk = s.uncertainty_PbShk(howmanysigma)
        PbEx = PbObs - PbGW - PbGal - PbShk
        errors = np.array([errPbObs, errPbGW, errPbGal, errPbShk])
        errPbEx = (sum(errors**2))**0.5
        return PbEx, errPbEx
    def Gdot2PbGdot(s): #Gdot/G from LLR to Pbdot_Gdot
        a = s.GdotG #Gdot/G by Hofmann et al. 2010
        errA = s.errGdotG #1sigma error for a
        q = s.q #m_p/m_c by Callanan et al. 1998
        errQ = s.errQ
        m = s.m_c #m_c/m_sun
        errM = s.errM_c
        
        fq = 2*q+1-1/(q+1)
        PbG = -2*a*s.Pb*(1-0.05*fq*m) #Pbdot_Gdot in d/yr
        PbG = PbG/s.yr2d
        print PbG
        dPdA = PbG/a
        dPdM = 2*0.05*a*s.Pb*fq/s.yr2d
        dPdQ = 2*0.05*a*s.Pb/s.yr2d*m
        dPdQ *= 2-1/(q+1)**2
        errPbG_sq = (dPdA*errA)**2
        errPbG_sq += (dPdM*errM)**2
        errPbG_sq += (dPdQ*errQ)**2
        errPbG = errPbG_sq**0.5
        return PbG, errPbG
    def Gdot2PbDp(s): #Gdot/G from LLR to Pbdot^dipole
        PbEx, errPbEx = s.PbEx()
        PbG, errPbG = s.Gdot2PbGdot()
        PbDp = PbEx-PbG
        errPbDp = (errPbEx**2+errPbG**2)**0.5
        return PbDp, errPbDp
    def Gdot2kD(s): #Gdot/G from LLR to Pbdot^dipole to kappa_D
        PbDp, errPbDp = s.Gdot2PbDp() 
        q = s.q #m_p/m_c by Callanan et al. 1998
        errQ = s.errQ
        m = s.m_c #m_c/m_sun
        errM = s.errM_c
        T = s.T
        print T
        B = 25/T/math.pi**2
        kD = -B/m**3*(q+1)/q**3*s.Pb*PbDp #d/s
        kD *= 24*3600
        dKdM = -3/m*kD
        dKdPbDp = kD/PbDp
        dKdQ = (1/(q+1)-3/q)*kD
        errK_sq = (dKdM*errM)**2
        errK_sq += (dKdPbDp*errPbDp)**2
        errK_sq += (dKdQ*errQ)**2
        errK = errK_sq**0.5
        return kD, errK

class estimate_stellar_radius_with_distance_magnitudes_and_optionally_T_eff:
    zero_mag_flux_dict = {'U':1780, 'B':4050, 'V':3635, 'R':3080, 
    'I':2420, 'J':1585, 'H':1020, 'K':640, 'L':236} #f_nu in jansky
    lamda_eff_dict = {'U':0.366, 'B': 0.436, 'V': 0.545, 'R':0.641, 'I':0.789, 'J':1.22, 'K':2.19, 'L':3.80} # in micrometer
    lamda_FWHM_dict = {'U':0.15, 'B':0.22, 'V':0.16, 'R':0.23, 'I':0.19, 'J':0.16, 'K':0.23, 'L':0.14}
    def __init__(s, targetname, filename): #s -> self
        s.targetname = targetname
        print s.read_magnitudes(targetname, filename)
    def read_magnitudes(s, targetname, filename):
        [junk1,junk2,targetdir,junk3,junk4] = prepare_path_source(targetname)
        magnifile = targetdir + '/optical/' + filename
        lines = open(magnifile).readlines()
        s.bands = []
        s.mags = np.array([])
        s.err_mags = np.array([])
        s.mag0_fluxes = np.array([])
        s.lamda_effs = np.array([])
        s.dlamdas = np.array([])
        for line in lines:
            for estimate in ['U ', 'B ', 'V ', 'R ', 'I ', 'J ', 'H ', 'K ', 'L ']:
                if estimate in line:
                    estimate = estimate.strip()
                    s.bands.append(estimate) 
                    s.mag0_fluxes = np.append(s.mag0_fluxes, s.zero_mag_flux_dict[estimate])
                    s.lamda_effs = np.append(s.lamda_effs, s.lamda_eff_dict[estimate])
                    s.dlamdas = np.append(s.dlamdas, s.lamda_FWHM_dict[estimate]/2) # HWHM dlamda/lamda_eff
                    line1 = line.split('=')[-1].strip()
                    exec("mag = %s" % line1.split('+-')[0].strip())
                    exec("err_mag = %s" % line1.split('+-')[-1].strip())
                    s.mags = np.append(s.mags, mag)
                    s.err_mags = np.append(s.err_mags, err_mag)
            for estimate in ['T_eff', 'dist']:
                if estimate in line:
                    line1 = line.split('=')[-1].strip().split(' ')[0].strip()
                    exec("s.%s = %s" % (estimate, line1.split('+-')[0]))
                    exec("s.err_%s = %s" % (estimate, line1.split('+-')[-1]))
        return s.dist, s.err_dist, s.T_eff, s.err_T_eff, s.bands, s.mags, s.err_mags, s.mag0_fluxes, s.lamda_effs, s.dlamdas
    def mags2fluxes(s, mags=0, err_mags=0):
        if mags == 0:
            mags = s.mags
        if err_mags == 0:
            err_mags = s.err_mags
        fluxes = s.mag0_fluxes*100**(-mags/5) #f_nu in jansky or 1e-26 W/m^2/Hz
        return fluxes
    def B_nu(s): #B_nu(nu, T)
        nu_effs = constants.c/s.lamda_effs # m/s/1e-6m
        nu_effs *= 1e6 # Hz
        print s.lamda_effs
        print nu_effs
        B_nus = 2*constants.h*nu_effs**3/constants.c**2
        B_nus /= np.exp(constants.h*nu_effs/constants.k/s.T_eff) - 1 #W/m^2/Hz
        B_nus *= math.pi #integral across solid angle
        B_nus *= 1e26 #jansky
        return B_nus
    def estimate_radii(s):
        D = s.dist #pc
        D *= estimate_uncertainty.pc2m #m
        fluxes = s.mags2fluxes()
        B_nus = s.B_nu()
        radii_sq = D**2*fluxes/B_nus
        radii = radii_sq**0.5 #m
        radii /= 1000 #km
        return radii

class read_name_and_position_from_catalog_then_covert_to_topcat_friendly_file:
    """
    read catalog with specific format and update its position to the latest SIMBAD ones.
    """
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, catalog='catalog1', tableformat='ascii'):
        s.catalog = s.path + catalog
        s.srcnames = s.read_srcname_from_catalog(s.catalog)
        print s.srcnames
        s.convert_simbad_coordinates_to_topcat_friendly_file(tableformat)
        #[s.srcnames, s.RAs, s.Decs] = s.read_name_and_position_from_catalog()
        #s.convert_to_topcat_friendly_file(tableformat)
    def read_name_and_position_from_catalog(s):
        srcnames = np.array([])
        RAs = np.array([])
        Decs = np.array([])
        #kinds = np.array([])
        lines = open(s.catalog).readlines()
        for line in lines:
            if not line.isspace():
                contents = line.split('   ')
                [srcname, kind, RA, Dec] = contents[0:4]
                RA = howfun.colonizedms(RA.strip())
                Dec = howfun.colonizedms(Dec.strip())
                print Dec
                RA = 15*howfun.dms2deg(RA)
                Dec = howfun.dms2deg(Dec)
                srcnames = np.append(srcnames, srcname.strip())
                #kinds = np.append(kinds, kind.strip())
                RAs = np.append(RAs, RA)
                Decs = np.append(Decs, Dec)
        return srcnames, RAs, Decs
    def read_srcname_from_catalog(s, catalog):
        srcnames = np.array([])
        lines = open(catalog).readlines()
        for line in lines:
            srcname = line.strip().split('  ')[0]
            srcnames = np.append(srcnames, srcname.strip())
        return srcnames
    def convert_to_topcat_friendly_file(s, tableformat):
        t = Table([s.srcnames, s.RAs, s.Decs], names=['srcname', 'RA', 'Dec']) 
        print t
        output = s.catalog + '.' + tableformat
        t.write(output, format=tableformat, overwrite=True)
    def simbad_coordinates_with_srcnames(s):
        from astroquery.simbad import Simbad
        s.queryresults = Simbad.query_objects(s.srcnames)
        print s.queryresults
        output = s.catalog + '_simbad_coordinates'
        s.queryresults.write(output, format='ascii', overwrite=True)
    def simbad_coordinate2deg(s, Decs): #or RA 
        if len(Decs) == 1:
            return s.simbad_coordinate2deg1(Decs)
        else:
            Decs_deg = np.array([])
            for Dec in Decs:
                Dec2 = s.simbad_coordinate2deg1(Dec)
                Decs_deg = np.append(Decs_deg, Dec2)
            return Decs_deg
    def simbad_coordinate2deg1(s, Dec):
        Dec_str = howfun.table_str(Dec)
        Dec1 = howfun.colonizedms(Dec_str)
        print Dec1
        return howfun.dms2deg(Dec1) 

    def convert_simbad_coordinates_to_topcat_friendly_file(s, tableformat='fits'):
        s.simbad_coordinates_with_srcnames()
        RAs_deg = 15*s.simbad_coordinate2deg(s.queryresults['RA'])
        Decs_deg = s.simbad_coordinate2deg(s.queryresults['DEC'])
        print RAs_deg, Decs_deg
        t = Table([s.srcnames, RAs_deg, Decs_deg], names=['srcname', 'RA', 'Dec']) 
        print t
        output = s.catalog + '_simbad_coordinates_only.' + tableformat
        t.write(output, format=tableformat, overwrite=True)
    def table_deg2dms(s, tablename='1arcsecond_Gaia_cross_identification.8sources', tableformat='ascii'):
        tablename = s.path + '/' + tablename
        t = Table.read(tablename, format=tableformat)
        t1 = t
        t1['RA'] = howfun.deg2dms(t['RA']/15)
        t1['Dec'] = howfun.deg2dms(t['Dec'])
        t1['ra_x'] = howfun.deg2dms(t['ra_x']/15)
        t1['dec_x'] = howfun.deg2dms(t['dec_x'])
        output = s.path + '/sources_dms'
        t1.write(output, format='ascii', overwrite=True)

class filter_the_list_of_Gaia_candidates_for_PRE_bursters_generated_from_topcat_with_cross_match_criterion:
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, catalog, inputformat='ascii', outputformat='ascii'):
        s.tablename = s.path + catalog
        s.inputformat, s.outputformat = inputformat, outputformat
        s.outputfile = s.tablename + '_filtered'
        s.filter_list_of_Gaia_candidates_for_PRE_bursters_generated_from_topcat_with_cross_match_criterion()
    def filter_list_of_Gaia_candidates_for_PRE_bursters_generated_from_topcat_with_cross_match_criterion(s):
        t=Table.read(s.tablename, format=s.inputformat)
        t.sort(['angDist'])
        index1 = abs(t['pmra'])/t['pmra_error']>3 
        index2 = abs(t['pmdec'])/t['pmdec_error']>3
        index = np.logical_or(np.array(index1), np.array(index2))
        s.t1 = t[index]
        s.t1.write(s.outputfile, format=s.outputformat, overwrite=True)

class Simbad_source_to_Gaia_count_in_1deg_radius_to_AGNs_crossmatch_radius:
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, srcname, AGN_number=1, count_radius=1, howmanysigma=3):
        """
        row_limit=-1 means no row limit for Gaia cone search.
        """
        s.srcname = srcname
        s.AGN_number = AGN_number
        howmanysigma1 = 1
        CL1 = 100*math.erf(howmanysigma1/2**0.5)
        [RA_deg, Dec_deg] = s.Simbad_srcname_to_position(s.srcname)
        print RA_deg, Dec_deg
        [Gaia_count, CpSqAs] = s.position_to_Gaia_count_in_1deg_radius(RA_deg, Dec_deg, count_radius)
        print "%d Gaia sources counted within %f-deg radius around %s, equivalent to %f count/arcsecond^2" % (Gaia_count, count_radius, srcname, CpSqAs)
        cross_match_radius1 = s.CpSqAs_and_AGN_number_to_cross_match_radius(CpSqAs, AGN_number, howmanysigma1)
        print "Cross-match radius for each background AGN to %s is %f arcsecond, in order to make sure no random field star cross-matched for the whole sample of %d AGNs at %f percent confidence level." % (srcname, cross_match_radius1, AGN_number, CL1)
        CL = 100*math.erf(howmanysigma/2**0.5)
        cross_match_radius = s.CpSqAs_to_cross_match_radius(CpSqAs, howmanysigma)
        print "Cross-match radius for each background AGN to %s is %f arcsecond, in order to make sure no random field star cross-matched for every background AGN at %f percent confidence level." % (srcname, cross_match_radius, CL)
        
        outfile = s.path + '/cross-match_radius_summary.txt'
        writefile = open(outfile, 'a')
        writefile.write("%d Gaia sources counted within %f deg radius around %s, equivalent to %f count/arcsecond^2\n" % (Gaia_count, count_radius, srcname, CpSqAs))
        writefile.write("Cross-match radius for each background AGN to %s is %f arcsecond, in order to make sure no random field star cross-matched for the whole sample of %d AGNs at %f percent confidence level.\n" % (srcname, cross_match_radius1, AGN_number, CL1))
        writefile.write("Cross-match radius for each background AGN to %s is %f arcsecond, in order to make sure no random field star cross-matched for every background AGN at %f percent confidence level.\n\n" % (srcname, cross_match_radius, CL))
        writefile.close()
    def Simbad_srcname_to_position(s, srcname):
        from astroquery.simbad import Simbad
        queryresult = Simbad.query_object(srcname)
        print queryresult
        #print type(queryresult['RA'])
        a = read_name_and_position_from_catalog_then_covert_to_topcat_friendly_file()
        RA_deg = 15*a.simbad_coordinate2deg(queryresult['RA'])
        Dec_deg = a.simbad_coordinate2deg(queryresult['DEC'])
        return RA_deg, Dec_deg
        
    def position_to_Gaia_count_in_1deg_radius(s, RA_deg, Dec_deg, radius=1.0):
        """
        row_limit=-1 means no row limit for Gaia cone search.
        """
        import astropy.units as u
        import timeit
        from astropy.coordinates import SkyCoord
        from astroquery.gaia import Gaia
        #Gaia.ROW_LIMIT = row_limit
        coord = SkyCoord(ra=RA_deg, dec=Dec_deg, unit=(u.degree, u.degree), frame='icrs')
        radius_unit = u.Quantity(radius, u.deg)
        start_time_cone_search = timeit.default_timer()
        j = Gaia.cone_search_async(coord, radius_unit)
        stop_time_cone_search = timeit.default_timer()
        s.time_spent_on_cone_search = stop_time_cone_search - start_time_cone_search
        print("%f min spent on the cone search within %f-deg radius of the target" % (s.time_spent_on_cone_search/60, float(radius)))
        s.results = j.get_results()
        Gaia_count = len(s.results)
        output = s.path + '/.source_list_' + str(radius) + 'deg_around_' + s.srcname.replace(' ','') + '.ascii'
        s.results.write(output, format='ascii', overwrite=True)
        count_per_sq_as = Gaia_count*1.0/(radius**2)/(3600**2)
        print "%d Gaia sources within %fdeg-radius circle of the target" % (Gaia_count,float(radius))
        #print AGN_count, radius
        return Gaia_count, count_per_sq_as
    def CpSqAs_and_AGN_number_to_cross_match_radius(s, CpSqAs, AGN_number, howmanysigma=1): 
        CL = math.erf(howmanysigma/2**0.5) #CL: 3->0.9973, 2->0.9545, 1->0.6827
        howmanySqAs = math.log(CL)/math.log(1-CpSqAs)
        r = math.sqrt(howmanySqAs/AGN_number) #unit: arcsecond
        return r
    def CpSqAs_to_cross_match_radius(s, CpSqAs, howmanysigma=5):
        CL = math.erf(howmanysigma/2**0.5)
        howmanySqAs = math.log(CL)/math.log(1-CpSqAs)
        r = math.sqrt(howmanySqAs) #unit: arcsecond
        return r

class single_out_quasars_with__cone_searched_Gaia_sources__and__AgnCrossId:
    """
    There are two ways to cross-match Gaia counterparts for backgound AGNs: (a) cross-match cone-searched Gaia sources with AgnCrossId*csv,
        (b) make my own ID_ra_dec catalog for AGNs using the AgnCrossId*csv, then do cone-search with my own codes.
    necessities for method (a): 1) decompressed AgnCrossId*.csv files; 2) cone searched Gaia sources in ascii format by my default
    """
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, cone_search_Gaia_sources):
        """
        cone_search_Gaia_sources needs to be just the file name, e.g. .source_list_1deg_around_CenX-4.ascii
        """
        s.cone_search_Gaia_sources = s.path + cone_search_Gaia_sources
        s.AgnCrossIDs = s.path + '/.numpy_list_of_AgnCrossIDs.pickle'
    def read_cone_searched_Gaia_sources(s):
        #import pandas as pd
        t = Table.read(s.cone_search_Gaia_sources, format='ascii')
        #print(t)
        #t1 = pd.DataFrame(np.array(t))
        #print(t)
        return t
    def filter_source_ID(s, t1, t_filter):
        source_ids = t_filter['source_id'] 
        common_ids = set(t1['source_id']).intersection(source_ids)
        ids = list(common_ids)
        #print(ids)
        return ids
    def merge_AgnCrossIDs(s):
        AgnCrossID_files = glob.glob(r'%s/edr3_agn_cross_id/AgnCrossId*.csv' % s.path)
        AgnCrossIDs = np.array([], dtype='int64')
        for AgnCrossID_file in AgnCrossID_files:
            t = Table.read(AgnCrossID_file, format='ascii')
            new_AgnCrossIDs = np.array(t['source_id'], dtype='int64')
            print(new_AgnCrossIDs)
            AgnCrossIDs = np.concatenate((AgnCrossIDs, new_AgnCrossIDs))
            print(AgnCrossIDs)
        writefile = open(s.AgnCrossIDs, 'w')
        pickle.dump(AgnCrossIDs, writefile)
        writefile.close()
    def make_AgnCrossID_RA_Dec_table_step1(s, start=0):
        if not os.path.exists(s.AgnCrossIDs):
            print('running merge_AgnCrossIDs...')
            s.merge_AgnCrossIDs()
        readfile = open(s.AgnCrossIDs, 'r')
        AgnCrossIDs = pickle.load(readfile)
        readfile.close()
        #print(AgnCrossIDs)
        #s.ID_pos = Table(names=['source_id', 'ra', 'dec'], dtype=['int64','float64','float64'])
        s.GaiaSourceFiles = glob.glob(r'%s/edr3_gaia_source/GaiaSource*.csv' % s.path)
        s.GaiaSourceFiles.sort()
        number_files = len(s.GaiaSourceFiles)
        count = start
        for GaiaSourceFile in s.GaiaSourceFiles[start:]:
            count += 1
            output_table = s.path + '/prepare_edr3_agn_astrometric/ID_astrometric_' + str(count).zfill(4)
            if os.path.exists(output_table):
                continue
            t = Table.read(GaiaSourceFile, format='ascii')
            GaiaIDs = np.array(t['source_id'], dtype='int64')
            #print(GaiaIDs)
            mask = np.in1d(GaiaIDs, AgnCrossIDs)
            #print(common_ids)
            #mask = [x in common_ids for x in GaiaIDs]
            AGNs = t[mask]
            AGNs_astrometric = Table([AGNs['source_id'], AGNs['ra'], AGNs['dec'], AGNs['parallax'], AGNs['parallax_error'], AGNs['pmra'], AGNs['pmra_error'], 
                AGNs['pmdec'], AGNs['pmdec_error'], AGNs['phot_g_mean_mag'], AGNs['bp_rp']], names=['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'phot_g_mean_mag', 'bp_rp'])
            print(AGNs_astrometric)
            AGNs_astrometric.write(output_table, format='ascii', overwrite=True)
            print('%d/%d has been finished.' % (count, number_files))
    def make_AgnCrossID_RA_Dec_table_step2(s):
        from astropy.table import vstack
        s.AGNs_astrometric_outputs = glob.glob(r'%s/prepare_edr3_agn_astrometric/ID_astrometric*' % s.path)
        s.AGNs_astrometric_outputs.sort()
        print(s.AGNs_astrometric_outputs)
        s.ID_pos = Table(names=['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'phot_g_mean_mag', 'bp_rp'],
            dtype=['int64', 'float64', 'float64', 'float64', 'float64', 'float64','float64', 'float64', 'float64', 'float64', 'float64'])
        count = 0
        number_files = len(s.AGNs_astrometric_outputs)
        for output in s.AGNs_astrometric_outputs:
            count += 1
            t1 = Table.read(output, format='ascii')
            s.ID_pos = vstack([s.ID_pos, t1])
            print('%d/%d has been finished.' % (count, number_files))
        final_output_table = s.path + '/AGNs_ID_astrometric.ascii'
        s.ID_pos.write(final_output_table, format='ascii', overwrite=True)
    
    def Simbad_srcname_to_position(s, srcname):
        from astroquery.simbad import Simbad
        queryresult = Simbad.query_object(srcname)
        print queryresult
        #print type(queryresult['RA'])
        a = read_name_and_position_from_catalog_then_covert_to_topcat_friendly_file()
        RA_deg = 15*a.simbad_coordinate2deg(queryresult['RA'])
        Dec_deg = a.simbad_coordinate2deg(queryresult['DEC'])
        return RA_deg, Dec_deg

    def cone_search_given__AGNs_ID_astrometric__catalog(s, targetname, cone_radius_in_deg):
        r = cone_radius_in_deg
        final_output_table = s.path + '/AGNs_ID_astrometric.ascii'
        if not os.path.exists(final_output_table):
            print('Please make the AGNs_ID_astrometric.ascii first; aborting')
            sys.exit()
        [RA_deg, Dec_deg] = s.Simbad_srcname_to_position(targetname)
        t = Table.read(final_output_table, format='ascii')
        print(t)
        Ds_in_deg = howfun.separations_deg(RA_deg, Dec_deg, t['ra'], t['dec'])
        print(Ds_in_deg)
        mask = Ds_in_deg < r
        print(mask)
        s.t1 = t[mask]
        output_table = s.path + '/background_Gaia_AGNs_for_' + targetname.replace(' ', '_') + '_within_' + str(r) + 'deg.ascii'
        s.t1.write(output_table, format='ascii', overwrite=True)
        return s.t1

    def single_out_quasars(s):
        if not os.path.exists(s.AgnCrossIDs):
            print('running merge_AgnCrossIDs...')
            s.merge_AgnCrossIDs()
        readfile = open(s.AgnCrossIDs, 'r')
        AgnCrossIDs = pickle.load(readfile)
        readfile.close()
        print(AgnCrossIDs)
        s.y1 = AgnCrossIDs
        s.t = s.read_cone_searched_Gaia_sources()
        print('Reading cone-searched Gaia results...')
        coneGaiaIDs = np.array(s.t['source_id'], dtype='int64')
        print(coneGaiaIDs)
        s.y2 = coneGaiaIDs
        s.common_ids = np.intersect1d(AgnCrossIDs, coneGaiaIDs)
        print('There are %d Gaia counterparts for background AGNs found.' % len(s.common_ids))
        mask = [x in s.common_ids for x in s.t['source_id']]
        s.background_AGNs = s.t[mask]
        #print('There are %d Gaia counterparts for background AGNs found.' % len(s.background_AGNs))
        output_table = s.path + '/.Gaia_sources_for_background_quasars'
        s.background_AGNs.write(output_table, format='ascii', overwrite=True)
        return s.background_AGNs
        


class catalog_of_Gaia_counterparts_for_AGNs_to_zero_parallax_point(object):
    """
    workflow is the default one-stop function
    """
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    global_mean_parallax_zero_point = -0.03 #Gaia DR2
    def __init__(s, srcname, magnitude_match_offset=0.02, use_high_precision_subset=False, color_match_offset=1):
        """
        use_high_precision_subset satisfies "table['parallax_error'<1]" (mas)
        """
        s.srcname,   s.MMO,                  s.use_high_precision_subset, s.CMO = \
            srcname, magnitude_match_offset, use_high_precision_subset,   color_match_offset
    def workflow(s):        
        import copy
        s.T = s.read_catalog_of_Gaia_counterparts_for_AGNs(s.srcname)
        s.T0 = copy.deepcopy(s.T)
        s.T = s.delete_rows_where_the_required_parameter_is_absent(s.T, 'parallax')
        s.T1 = copy.deepcopy(s.T)
        s.T = s.filter_out_sources_with_detected_proper_motions(s.T, 3, s.use_high_precision_subset)
        s.T2 = copy.deepcopy(s.T)
        s.mag_g = s.read_phot_g_mean_mag_of_target(s.srcname)
        s.T = s.get_like_magnitude_AGNs(s.T, s.mag_g, s.MMO)
        #s.T3 = copy.deepcopy(s.T)
        s.bp_rp = s.read_bp_rp_of_target(s.srcname)
        s.T = s.get_like_color_AGNs(s.T, s.bp_rp, s.CMO)
        s.T3 = copy.deepcopy(s.T)
        if len(s.T3)>1:
            [s.zero_parallax_point, s.err_ZPP, s.std_parallax] = s.calculate_zero_parallax_and_its_sigma(s.T3)
            print s.zero_parallax_point, s.err_ZPP, s.std_parallax
    def plot__S_MMO_relation(s, lg_MMO_step=0.1):
        """
        an extra function that plots S-MMO relation, where
        S stands for weighted standard deviation of parallaxes;
        MMO stands for relative half width of G-band magnitude used to define the magnitude filter.
        This function is not a part of the workflow. It calls workflow many times to make the plot.
        """
        import matplotlib.pyplot as plt
        lg_MMOs = Ss = Ns = np.array([]) #N stands for volume of s.T3
        #for i in np.arange(-3, lg_MMO_step, lg_MMO_step):
        for i in np.arange(-2.0, -1+lg_MMO_step, lg_MMO_step):
            s.MMO = 10**i
            s.workflow()
            N = len(s.T3)
            lg_MMOs = np.append(lg_MMOs, i)
            Ss = np.append(Ss, s.std_parallax)
            Ns = np.append(Ns, N)
            lg_Ns = np.log10(Ns)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i+3)/3.*1000)/10) + " \r")
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(lg_MMOs, Ss, color='b')
        #ax1.axvline(x=-1.78, c='black', linestyle='--', linewidth=0.5)
        #ax1.axvline(x=-1.89, c='black', linestyle='--', linewidth=0.5)
        ax1.axvline(x=-1.41, c='black', linestyle='--', linewidth=0.5)
        ax1.axvline(x=-1.45, c='black', linestyle='--', linewidth=0.5)
        #ax1.axvline(x=s.value_PI+s.error_PI, c='black', linestyle='-.', linewidth=0.5)
        ax1.set_ylabel(r'$s_{\pi_0}$ (mas)', color='b')
        ax1.set_xlabel(r'$\log_{10} \Delta m_G$')
        ax2 = ax1.twinx()
        ax2.set_ylabel(r'$\log_{10} N_\mathrm{quasar}$', color='r', alpha=0.4)
        ax2.plot(lg_MMOs, lg_Ns, color='r', alpha=0.4)
        #plt.ylabel(r'$s_{\pi_0}$ (mas)')
        #plt.xlabel(r'$\lg{\Delta m_G}$')
        #plt.show()
        plt.savefig('%s/S_MMO_relation_for_%s.pdf' % (s.path, s.srcname.replace(' ', '')))
        plt.clf()
    def plot__S_CMO_relation(s, lg_CMO_step=0.1):
        """
        an extra function that plots S-CMO relation, where
        S stands for weighted standard deviation of parallaxes;
        CMO stands for relative half width of bp_rp used to define the color filter.
        This function is not a part of the workflow. It calls workflow many times to make the plot.
        """
        import matplotlib.pyplot as plt
        lg_CMOs = Ss = Ns = np.array([]) #N stands for volume of s.T3
        for i in np.arange(-3, lg_CMO_step, lg_CMO_step):
        #for i in np.arange(-2.2, -1.8+lg_CMO_step, lg_CMO_step):
            s.CMO = 10**i
            s.workflow()
            N = len(s.T3)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i+3)/3.*1000)/10) + " \r")
            if N in [0, 1]:
                continue
            lg_CMOs = np.append(lg_CMOs, i)
            Ss = np.append(Ss, s.std_parallax)
            Ns = np.append(Ns, N)
            lg_Ns = np.log10(Ns)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(lg_CMOs, Ss, color='b')
        ax1.axvline(x=-1.98, c='black', linestyle='--', linewidth=0.5)
        ax1.set_ylabel(r'$s_{\pi_0}$ (mas)', color='b')
        ax1.set_xlabel(r'$\log_{10} \Delta m_\mathrm{bp-rp}$')
        ax2 = ax1.twinx()
        ax2.set_ylabel(r'$\log_{10} N_\mathrm{quasar}$', color='r', alpha=0.4)
        ax2.plot(lg_CMOs, lg_Ns, color='r', alpha=0.4)
        #plt.ylabel(r'$s_{\pi_0}$ (mas)')
        #plt.xlabel(r'$\lg{\Delta m_G}$')
        #plt.show()
        plt.savefig('%s/S_CMO_relation_for_%s.pdf' % (s.path, s.srcname.replace(' ', '')))
        plt.clf()
    
    def plot__S_MMO__and__S_CMO__relations_for_CygX2_and_CenX4(s, lg_MMO_step=0.1, lg_CMO_step=0.1):
        """
        only used to make the plot in the paper
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        plt.rc('text', usetex=True)
        lg_CMOs_CenX4 = lg_MMOs_CenX4 = lg_CMOs_CygX2 = SsM_CenX4 = SsC_CenX4 = SsC_CygX2 = \
            NsM_CenX4 = NsC_CenX4 = NsC_CygX2 = np.array([]) #N stands for volume of s.T3
        s.MMO = s.CMO = 1
        s.srcname = 'Cen X-4'
        s.workflow()
        Nmax_CenX4 = len(s.T3)
        for i in np.arange(-3, lg_MMO_step, lg_MMO_step): #S-MMO relation for Cen X-4
            s.MMO = 10**i
            s.workflow()
            N = len(s.T3)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i+3)/3.*1000)/10) + " \r")
            if N in [0, 1]:
                continue
            lg_MMOs_CenX4 = np.append(lg_MMOs_CenX4, i)
            SsM_CenX4 = np.append(SsM_CenX4, s.std_parallax)
            NsM_CenX4 = np.append(NsM_CenX4, N)
            lg_NsM_CenX4 = np.log10(NsM_CenX4)
            if N == Nmax_CenX4:
                break
        
        s.MMO = s.MMO = 1
        for i in np.arange(-3, 3, lg_CMO_step): #S-CMO relation for Cen X-4
            s.CMO = 10**i
            s.workflow()
            N = len(s.T3)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i+3)/3.*1000)/10) + " \r")
            if N in [0, 1]:
                continue
            lg_CMOs_CenX4 = np.append(lg_CMOs_CenX4, i)
            SsC_CenX4 = np.append(SsC_CenX4, s.std_parallax)
            NsC_CenX4 = np.append(NsC_CenX4, N)
            lg_NsC_CenX4 = np.log10(NsC_CenX4)
            if N == Nmax_CenX4:
                break
        
        s.MMO = s.CMO = 1
        s.srcname = 'Cyg X-2'
        s.workflow()
        Nmax_CygX2 = len(s.T3)
        for i in np.arange(-3, 3, lg_CMO_step): #S-MMO relation for Cen X-4
            s.CMO = 10**i
            s.workflow()
            N = len(s.T3)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i+3)/3.*1000)/10) + " \r")
            if N in [0, 1]:
                continue
            lg_CMOs_CygX2 = np.append(lg_CMOs_CygX2, i)
            SsC_CygX2 = np.append(SsC_CygX2, s.std_parallax)
            NsC_CygX2 = np.append(NsC_CygX2, N)
            lg_NsC_CygX2 = np.log10(NsC_CygX2)
            if N == Nmax_CygX2:
                break

        #fig = plt.figure(figsize=[8,6])
        fig = plt.figure()
        gs = gridspec.GridSpec(4, 4)
        #subplot1
        ax1 = fig.add_subplot(gs[:2, :2])
        ax1.plot(lg_MMOs_CenX4, SsM_CenX4, color='b')
        ax1.axvline(x=-1.45, c='black', linestyle='--', linewidth=0.5)
        ax1.set_ylabel(r'$s_{\pi_0}$ (mas)', color='b')
        ax1.set_xlabel(r'$\log_{10} \Delta m_G$')
        #ax1.set_title(r'$s_{\pi_0}-\Delta m_G$ relation for Cen~X$-$4')
        ax1.set_title(r'Cen~X$-$4')
        ax2 = ax1.twinx()
        ax2.set_ylabel(r'$\log_{10} N_\mathrm{quasar}$', color='r', alpha=0.4)
        ax2.plot(lg_MMOs_CenX4, lg_NsM_CenX4, color='r', alpha=0.4)
        #subplot2 
        ax3 = fig.add_subplot(gs[:2, 2:4])
        ax3.plot(lg_CMOs_CenX4, SsC_CenX4, color='b')
        ax3.set_ylabel(r'$s_{\pi_0}$ (mas)', color='b')
        ax3.set_xlabel(r'$\log_{10} \Delta m_\mathrm{B-R}$')
        #ax3.set_title(r'$s_{\pi_0}-\Delta m_\mathrm{bp-rp}$ relation for Cen~X$-$4')
        ax3.set_title(r'Cen~X$-$4')
        ax4 = ax3.twinx()
        ax4.set_ylabel(r'$\log_{10} N_\mathrm{quasar}$', color='r', alpha=0.4)
        ax4.plot(lg_CMOs_CenX4, lg_NsC_CenX4, color='r', alpha=0.4)
        #subplot3
        ax5 = fig.add_subplot(gs[2:4, 1:3])
        ax5.plot(lg_CMOs_CygX2, SsC_CygX2, color='b')
        ax5.set_ylabel(r'$s_{\pi_0}$ (mas)', color='b')
        ax5.set_xlabel(r'$\log_{10} \Delta m_\mathrm{B-R}$')
        #ax5.set_title(r'$s_{\pi_0}-\Delta m_\mathrm{bp-rp}$ relation for Cyg~X$-$2')
        ax5.set_title(r'Cyg~X$-$2')
        ax6 = ax5.twinx()
        ax6.set_ylabel(r'$\log_{10} N_\mathrm{quasar}$', color='r', alpha=0.4)
        ax6.plot(lg_CMOs_CygX2, lg_NsC_CygX2, color='r', alpha=0.4)
        
        gs.tight_layout(fig) #rect=[0, 0.1, 1, 1])
        plt.savefig('%s/S_CMO__and__S_MMO__relations_for_CenX4_and_CygX2.pdf' % s.path)
        plt.clf()

    def workflow_analysing_parallax2magnitude_or_parallax2color_relation(s, parameter_str):
        s.T2 = s.sort_table_according_to_magnitude_or_color(s.T2, parameter_str)
        t0 = s.get_parallax_parameter_relation(s.T2, parameter_str, 100)
        print t0
        s.plot_parallax_parameter_relation(t0, parameter_str)
    def read_catalog_of_Gaia_counterparts_for_AGNs(s, srcname):
        if srcname == 'Cyg X-2':
            catalogfile = s.path + '/CygX-2_Gaia_sources_for_background_quasars'
        elif srcname == 'Cen X-4':
            catalogfile = s.path + '/CenX-4_Gaia_sources_for_background_quasars'
        elif srcname == 'NP Ser':
            catalogfile = s.path + '/NpSer_Gaia_sources_for_background_quasars'
        else:
            print "The provided srcname is not recognized; aborting..."
            sys.exit()
        catalogtable = Table.read(catalogfile, format='ascii')
        return catalogtable
    def delete_rows_where_the_required_parameter_is_absent(s, catalogtable, parameter_str):
        values = catalogtable[parameter_str]
        i = 0
        for value in values:
            if not float(value)<1e10:
                catalogtable.remove_row(i)
                i -= 1
            i += 1
        return catalogtable
    def filter_out_sources_with_detected_proper_motions(s, input_table, threshold_SNR=3, use_high_precision_subset=False):
        """
        also can filter in high_precision_subset with condition "parallax_error<1 (mas)"
        """
        t = input_table
        index1 = abs(t['pmra'])/t['pmra_error'] < threshold_SNR
        index2 = abs(t['pmdec'])/t['pmdec_error'] < threshold_SNR
        index = np.logical_and(np.array(index1), np.array(index2))
        if use_high_precision_subset:
            #index3 = (t['parallax']-s.global_mean_parallax_zero_point)/t['parallax_error'] < 3
            index3 = t['parallax_error'] < 1
            index = np.logical_and(np.array(index), np.array(index3))
        return s.T[index]
    def sort_table_according_to_magnitude_or_color(s, catalogtable, parameter_str):
        T0 = s.delete_rows_where_the_required_parameter_is_absent(catalogtable, parameter_str)
        T0.sort(parameter_str)
        return T0
    def get_parallax_parameter_relation(s, sorted_catalogtable, parameter_str, number_of_sources_in_a_bin=10):
        ## first step: bin the table
        T0 = sorted_catalogtable
        number = number_of_sources_in_a_bin
        length = len(T0) 
        bins = length/number
        ## second step: get average_parallax and errors for each bin
        parallaxes = T0['parallax']
        errs_parallax = T0['parallax_error']
        valueXs = T0[parameter_str]
        Xs = np.array([]) #mid-point of color/magnitude in each bin
        Ps = np.array([]) #average parallaxes in each bin
        errs_P = np.array([]) # errors of average parallaxes in each bin
        j = 0
        for i in range(0,bins):
            if i == 0:
                num = number + length%number
            else:
                num = number
            start = j
            end = j + num - 1
            X = valueXs[start:end+1].mean()
            Xs = np.append(Xs, X)
            [P, err_P] = howfun.weightX(parallaxes[start:end+1], errs_parallax[start:end+1])
            Ps = np.append(Ps, P)
            errs_P = np.append(errs_P, err_P)
            j = end + 1
        t = Table([Xs, Ps, errs_P], names=[parameter_str, 'parallax', 'parallax_error']) 
        return t
    def plot_parallax_parameter_relation(s, X_P_err_P_table, parameter_str):
        import matplotlib.pyplot as plt
        Xs = X_P_err_P_table[parameter_str]
        Ps = X_P_err_P_table['parallax']
        errs_P = X_P_err_P_table['parallax_error']
        plt.plot(Xs, Ps)
        plt.ylabel('zero-parallax point')
        plt.xlabel(parameter_str)
        plt.savefig('%s/marginalized_%s_parallax_relation.eps' % (s.path, parameter_str))
        plt.clf()
    
    def read_phot_g_mean_mag_of_target(s, srcname):
        if srcname == 'NP Ser':
            mag_g = 17.032724
        elif srcname == 'Cen X-4':
            mag_g = 17.864574
        elif srcname == 'Cyg X-2':
            mag_g = 17.66 #14.726028
        return mag_g
    def read_bp_rp_of_target(s, srcname):
        if srcname == 'NP Ser':
            bp_rp = 1.56335 
        elif srcname == 'Cen X-4':
            bp_rp = 1.64234
        elif srcname == 'Cyg X-2':
            bp_rp = 0.70213
        return bp_rp
    def get_like_magnitude_AGNs(s, catalogtable, mag_g, magnitude_match_offset):
        MMO = magnitude_match_offset
        i = catalogtable['phot_g_mean_mag']<mag_g*(1+MMO)
        T1 = catalogtable[i]
        i = T1['phot_g_mean_mag']>mag_g*(1-MMO)
        T2 = T1[i]
        return T2
    def get_like_color_AGNs(s, catalogtable, bp_rp, color_match_offset):
        CMO = color_match_offset
        i = catalogtable['bp_rp'] < bp_rp + abs(CMO*bp_rp)
        T1 = catalogtable[i]
        i = T1['bp_rp'] > bp_rp - abs(CMO*bp_rp)
        T2 = T1[i]
        return T2
    def calculate_zero_parallax_and_its_sigma(s, catalogtable):
        parallaxes = catalogtable['parallax']
        errs_parallax = catalogtable['parallax_error']
        [av_parallax, std_parallax] = howfun.weighted_avg_and_std(parallaxes, errs_parallax)
        [average_parallax, integral_err] = howfun.weightX(parallaxes, errs_parallax)
        return average_parallax, integral_err, std_parallax

class catalog_of_Gaia_EDR3_counterparts_for_quasars_to_zero_parallax_point(catalog_of_Gaia_counterparts_for_AGNs_to_zero_parallax_point):
    def __init__(s, srcname, magnitude_match_offset=0.02, use_high_precision_subset=False, color_match_offset=1):
        super(catalog_of_Gaia_EDR3_counterparts_for_quasars_to_zero_parallax_point, s).__init__(srcname, magnitude_match_offset) #python2 way to use super
        s.srcname,   s.MMO,                  s.use_high_precision_subset, s.CMO = \
            srcname, magnitude_match_offset, use_high_precision_subset,   color_match_offset
    def workflow_EDR3(s):
        import copy
        s.T = s.read_catalog_of_Gaia_EDR3_counterparts_for_AGNs(s.srcname)
        s.T0 = copy.deepcopy(s.T)
        s.T = s.delete_rows_where_the_required_parameter_is_absent(s.T, 'parallax')
        s.T1 = copy.deepcopy(s.T)
        s.T = s.filter_out_sources_with_detected_proper_motions(s.T, 3, s.use_high_precision_subset)
        s.T2 = copy.deepcopy(s.T)

    def read_catalog_of_Gaia_EDR3_counterparts_for_AGNs(s, srcname, radius_in_deg=5):
        r = radius_in_deg
        catalogfile = s.path + '/background_Gaia_AGNs_for_' + srcname.replace(' ', '_') + '_within_' + str(r) + 'deg.ascii'
        if not os.path.exists(catalogfile):
            print("%s does not exist; aborting" % catalogfile)
            sys.exit()
        catalogtable = Table.read(catalogfile, format='ascii')
        return catalogtable

class plot_positions_of_like_magnitude_background_AGNs(catalog_of_Gaia_counterparts_for_AGNs_to_zero_parallax_point):
    def __init__(s, srcname, magnitude_match_offset=0.02):
        super(plot_positions_of_like_magnitude_background_AGNs, s).__init__(srcname, magnitude_match_offset) #python2 way to use super
        s.srcname = srcname
        [s.RAs, s.Decs] = s.read_RAs_Decs()
        [s.RA_t, s.Dec_t] = s.read_target_position()
        s.plot_RAs_Decs()
    def read_RAs_Decs(s):
        RAs = s.T3['RA']/15
        Decs = s.T3['DEC']
        return RAs, Decs
    def read_target_position(s, tablename='1arcsecond_Gaia_cross_identification.8sources', tableformat='ascii'):
        tablename = s.path + '/' + tablename
        t = Table.read(tablename, format=tableformat)
        srcname = s.srcname
        if srcname != 'NP Ser':
            index = t['srcname']==s.srcname
        else:
            index = t['srcname']=='GX 17+2'
        RA = t['RA'][index]/15
        Dec = t['Dec'][index]
        return RA, Dec
    def plot_RAs_Decs(s):
        import matplotlib.pyplot as plt
        av_RAs = s.RAs.mean()
        av_Decs = s.Decs.mean()
        plt.scatter(s.RAs, s.Decs, marker='.')
        plt.plot(av_RAs, av_Decs, 'rs')
        plt.plot(s.RA_t, s.Dec_t, 'g^')
        plt.ylabel('Declination (deg)')
        plt.xlabel('Right Ascension (hour)')
        plt.savefig('%s/positions_of_background_quasars_for_%s.eps' % (s.path, s.srcname.replace(' ', '')))
        plt.clf()
        
class plot_Gaia_sources_around_a_source:
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, srcname, radius=10, src_position=''):
        """
        radius in arcsecond;
        workflow is the default function that would make the plots for you.
        srcposition in the form of 'HH:MM:SS.SSSSS,dd:mm:ss.sss';
        if srcposition is not given, then position will be searched and obtained from SIMBAD
        """
        s.srcname, s.radius, s.src_position = srcname, radius, src_position
    def workflow(s):
        RA_deg, Dec_deg = s.Simbad_srcname_to_position(s.srcname)
        s.T = s.cone_search_Gaia_sources_within_given_radius(RA_deg, Dec_deg, s.radius)
        s.plot_RAs_Decs(s.srcname, s.radius, RA_deg, Dec_deg)
    def Simbad_srcname_to_position(s, srcname):
        from astroquery.simbad import Simbad
        if s.src_position != '':
            [RA, Dec] = s.src_position.split(',')
            RA_deg = 15*howfun.dms2deg(RA.strip())
            Dec_deg = howfun.dms2deg(Dec.strip())
        else:
            queryresult = Simbad.query_object(srcname)
            print queryresult
            #print type(queryresult['RA'])
            a = read_name_and_position_from_catalog_then_covert_to_topcat_friendly_file()
            RA_deg = 15*a.simbad_coordinate2deg(queryresult['RA'])
            Dec_deg = a.simbad_coordinate2deg(queryresult['DEC'])
        return RA_deg, Dec_deg
    def cone_search_Gaia_sources_within_given_radius(s, RA_deg, Dec_deg, radius):
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astroquery.gaia import Gaia
        coord = SkyCoord(ra=RA_deg, dec=Dec_deg, unit=(u.degree, u.degree), frame='icrs')
        radius_unit = u.Quantity(radius, u.arcsec)
        j = Gaia.cone_search_async(coord, radius_unit)
        return j.get_results()
    def plot_RAs_Decs(s, srcname, radius, RA_t_deg, Dec_t_deg):
        import matplotlib.pyplot as plt
        #av_RAs = s.RAs.mean()
        #av_Decs = s.Decs.mean()
        RA_t_dms, Dec_t_dms = howfun.deg2dms(RA_t_deg/15), howfun.deg2dms(Dec_t_deg)
        RAs_dms, Decs_dms = howfun.deg2dms(s.T['ra']/15), howfun.deg2dms(s.T['dec'])
        diffRAs, diffDecs = diffposition(RAs_dms, RA_t_dms, Decs_dms, Dec_t_dms) #in mas
        diffRAs /= 1000 #in arcsec
        diffDecs /= 1000 #in arcsec
        #comments = ['Ser X-1 (DSe)','DSw','DN','','','']
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(diffRAs, diffDecs, marker='.', s=200)
        #for i, comment in enumerate(comments):
        #    ax.annotate(comment, (diffRAs[i], diffDecs[i]), color='cyan', )
        #ax.annotate('Ser X-1 (DSe)', (diffRAs[0], diffDecs[0]), color='cyan', xytext=(3.3,-0.1))
        #ax.annotate('DSw', (diffRAs[1], diffDecs[1]), color='cyan', xytext=(-1.5,-0.5))
        #ax.annotate('DN', (diffRAs[2], diffDecs[2]), color='cyan', xytext=(0.4,2))
        ax.text(5.0,-0.1,'Ser X-1 (DSe)', size=15, color='cyan')
        ax.text(-1.5,-0.5,'DSw', size=15, color='cyan')
        ax.text(0.8,2.2,'DN', size=15, color='cyan')
        #ax.plot(0, 0, 'g^', markersize=2)
        ax.set_ylabel('relative Decl. (arcsec)', fontsize='xx-large')
        ax.set_ylim(-4,4)
        #ax.set_yticklabels(y_ticks, fontsize='x-large')
        plt.yticks([-2,0,2], fontsize='xx-large')
        ax.set_xlabel('relative RA. (arcsec)', fontsize='xx-large')
        ax.set_xlim(-9,8)
        #ax.set_xticklabels(x_ticks, fontsize='x-large')
        plt.xticks(np.arange(-8,8,2), fontsize='xx-large')
        plt.gca().invert_xaxis()
        ax.set_aspect('equal')
        plt.savefig('%s/Gaia_sources_within_%darcsec_radius_around_%s.eps' % (s.path, radius, srcname.replace(' ', '')), transparent=True)
        plt.clf()

class estimate_space_velocities_of_PRE_bursters:
    """
    read in proper motions and distances for PRE bursters and generate their apparent/peculiar transverse velocities;
    """
    path = '/fred/oz002/hding/AQLX-1/PREBursters_catalog/'
    R0 = 8.15 #(+-0.15kpc) Reid et al. 2019
    V_Sun2GC = 247 #(+-4km/s)
    #V0 = 236 #(+-7km/s)
    def __init__(s):
        s.load_rotation_curve_by_Reid2019()
        a = estimate_uncertainty()
        s.A = a.A
    def calculate_v_t_with_distance(s, D, mu_a, mu_d):
        a = estimate_uncertainty()
        v_t = (mu_a**2+mu_d**2)**0.5 * D * a.A
        return v_t
    def calculate_peculiar_transverse_velocity_in_the_Galaxy(s, d, l, b, mu_l, mu_b): #in kpc, deg, deg, mas/yr, mas/yr
        """
        It is originally designed for J1810 falling into the first quadrant of the Galactic coordinates.
        It might need extra inspection for the validity in other qudrants.
        This function will be improved later.
        Non-circular motion of the Sun is neglected.
        """
        l *= math.pi/180
        b *= math.pi/180
        v_b_obs = v_t * mu_b / (mu_b**2 + mu_l**2)**0.5
        v_l_obs = v_t * mu_l / (mu_b**2 + mu_l**2)**0.5
        v_b_gsr = v_b_obs - V_Sun2GC*math.sin(l)*math.sin(b)
        v_l_gsr = v_l_obs + V_Sun2GC*math.cos(l)
        alpha = math.atan((R0*math.cos(l)-d)/(R0*math.sin(l)))
        v_l_loc = v_l_gsr - V0*math.sin(alpha)
        v_b_loc = v_b_gsr + V0*math.cos(alpha)*math.sin(b)
        return v_b_loc, v_l_loc
    def distance_to_the_Galatic_center_based_on_D_and_l(s, D, l, b):
        """
        D in kpc, l in deg, D_GC in kpc.
        """
        l *= math.pi/180
        b *= math.pi/180
        D1 = D * math.cos(b)
        D_GC = s.R0**2 + D1**2 - 2*s.R0*D1*math.cos(l)
        return D_GC**0.5
    def distances_to_the_Galatic_center_based_on_D_and_l(s):
        """
        D_GC means the distance from the target to the Galactic center;
        given up upon in the end.
        """
        s.write_out_table_summaring_l_b_mu_l_mu_d()
        Dmaxs, Dmins, ls, bs = s.t['Dmax'], s.t['Dmin'], s.t1['l'], s.t1['b']
        Dmax_GCs = Dmin_GCs = np.array([])
        for i in range(len(ls)):
            D1_GC = s.distance_to_the_Galatic_center_based_on_D_and_l(Dmaxs[i], ls[i], bs[i])
            D2_GC = s.distance_to_the_Galatic_center_based_on_D_and_l(Dmins[i], ls[i], bs[i])
            Dmax_GC, Dmin_GC = max(D1_GC, D2_GC), min(D1_GC, D2_GC)
            Dmax_GCs, Dmin_GCs = np.append(Dmax_GCs, Dmax_GC), np.append(Dmin_GCs, Dmin_GC)
        s.t2 = Table([s.t['srcname'], Dmin_GCs, Dmax_GCs], names=['srcname', 'Dmin_GC', 'Dmax_GC'])
    def load_rotation_curve_by_Reid2019(s):
        rotation_curve = s.path + 'rotation_curve_Reid2019'
        s.RC = Table.read(rotation_curve, format='ascii')
    def interpolate_local_circular_velocity_with_the_rotation_curve_from_Reid2019(s, R):
        """
        R in kpc; R means the distance from the target to the Galactic center, same as D_GC;
        v0 represents the rotation velocity of a target around the Galactic center.
        """
        index = abs(s.RC['R'] - R) < 0.25
        RC1 = s.RC[index]
        if len(RC1) != 2:
            if len(RC1) == 1:
                v0 = RC1['v0'][0]
            else:
                print('%f is out of the range of the rotation curve by Reid et al. 2019.' % R)
                sys.exit()
        else:
            v0 = abs(RC1['R'][0]-R)/0.25*RC1['v0'][1] + abs(RC1['R'][1]-R)/0.25*RC1['v0'][0]
        return v0
    def v_l_pec__for_the_first_and_fourth_Galactic_quadrants(s, D, mu_l, l, b):
        """
        v_l_pec is the v_l relative to that of the neighbourhood of the target;
        D in kpc, mu_l in mas/yr, l in deg, b in deg;
        V_Sun2GC is regarded a constant (without accounting for its uncertainty).
        """
        R = s.distance_to_the_Galatic_center_based_on_D_and_l(D, l, b) 
        v0 = s.interpolate_local_circular_velocity_with_the_rotation_curve_from_Reid2019(R)
        l *= math.pi/180
        alpha = math.acos(s.R0*math.sin(l)/R)
        d1 = D - s.R0*math.cos(l)
        v_l_pec = mu_l * D * s.A + s.V_Sun2GC * math.cos(l) + d1/abs(d1)*v0*math.sin(alpha)
        return v_l_pec
    def v_l_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(s, l, b, mu_l, err_mu_l, Dmin, Dmax, D_step=0.1):
        v_l_pec_max = float('-inf')
        v_l_pec_min = float('inf')
        for mu_l1 in [mu_l-err_mu_l, mu_l+err_mu_l]:
            for D in np.arange(Dmin, Dmax+D_step, D_step):
                v_l_pec = s.v_l_pec__for_the_first_and_fourth_Galactic_quadrants(D, mu_l1, l, b)
                if v_l_pec > v_l_pec_max:
                    v_l_pec_max = v_l_pec
                if v_l_pec < v_l_pec_min:
                    v_l_pec_min = v_l_pec
        return v_l_pec_min, v_l_pec_max
    def v_b_pec__for_the_first_and_fourth_Galactic_quadrants(s, D, mu_b, l, b):
        """
        b in deg; mu_l in mas/yr; D in kpc.
        """
        R = s.distance_to_the_Galatic_center_based_on_D_and_l(D, l, b) 
        v0 = s.interpolate_local_circular_velocity_with_the_rotation_curve_from_Reid2019(R)
        l *= math.pi/180
        b *= math.pi/180
        cos_alpha = s.R0*math.sin(l)/R
        v_b_pec = mu_b*D*s.A - s.V_Sun2GC*math.sin(l)*math.sin(b) + v0*cos_alpha*math.sin(b)*math.sin(l)/abs(math.sin(l))
        return v_b_pec
    def v_b_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(s, l, b, mu_b, err_mu_b, Dmin, Dmax, D_step=0.1):
        v_b_pec_max = float('-inf')
        v_b_pec_min = float('inf')
        for mu_b1 in [mu_b-err_mu_b, mu_b+err_mu_b]:
            for D in np.arange(Dmin, Dmax+D_step, D_step):
                v_b_pec = s.v_b_pec__for_the_first_and_fourth_Galactic_quadrants(D, mu_b1, l, b)
                if v_b_pec > v_b_pec_max:
                    v_b_pec_max = v_b_pec
                if v_b_pec < v_b_pec_min:
                    v_b_pec_min = v_b_pec
        return v_b_pec_min, v_b_pec_max
    def v_t_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(s, l, b, mu_l, err_mu_l,\
            mu_b, err_mu_b, Dmin, Dmax, D_step=0.1):
        v_t_pec_max = float('-inf')
        v_t_pec_min = float('inf')
        for mu_l1 in [mu_l-err_mu_l, mu_l+err_mu_l]:
            for mu_b1 in [mu_b-err_mu_b, mu_b+err_mu_b]:
                for D in np.arange(Dmin, Dmax+D_step, D_step):
                    v_l_pec = s.v_l_pec__for_the_first_and_fourth_Galactic_quadrants(D, mu_l1, l, b)
                    v_b_pec = s.v_b_pec__for_the_first_and_fourth_Galactic_quadrants(D, mu_b1, l, b)
                    v_t_pec = (v_l_pec**2 + v_b_pec**2)**0.5
                    if v_t_pec > v_t_pec_max:
                        v_t_pec_max = v_t_pec
                    if v_t_pec < v_t_pec_min:
                        v_t_pec_min = v_t_pec
        return v_t_pec_min, v_t_pec_max
    def calculate_proper_motion_in_Galactic_coordinates_with_uncertainties(s, RA, Dec, mu_a, mu_d, err_mu_a, err_mu_d):
        """
        RA in h, Dec in deg, or both in 'xx:xx:xx.xxx';
        mu_a and mu_d in mas/yr
        """
        mu_a_l = howfun.upper_limit_or_lower_limit_with_larger_magnitude(mu_a, err_mu_a)
        mu_a_s = howfun.upper_limit_or_lower_limit_with_smaller_magnitude(mu_a, err_mu_a)
        mu_d_l = howfun.upper_limit_or_lower_limit_with_larger_magnitude(mu_d, err_mu_d)
        mu_d_s = howfun.upper_limit_or_lower_limit_with_smaller_magnitude(mu_d, err_mu_d)
        a = howfun.equatorial2galactic_coordinates(RA, Dec, mu_a, mu_d)
        mu_l, mu_b = a.equatorial_proper_motion2galactic_proper_motion()
        a_l = howfun.equatorial2galactic_coordinates(RA, Dec, mu_a_l, mu_d_l)
        mu_l_l, mu_b_l = a_l.equatorial_proper_motion2galactic_proper_motion()
        a_s = howfun.equatorial2galactic_coordinates(RA, Dec, mu_a_s, mu_d_s)
        mu_l_s, mu_b_s = a_s.equatorial_proper_motion2galactic_proper_motion()
        err_mu_l = max(abs(mu_l-mu_l_l), abs(mu_l-mu_l_s))
        err_mu_b = max(abs(mu_b-mu_b_l), abs(mu_b-mu_b_s))
        #print('mu_l = %f +- %f \n mu_b = %f +- %f' % (mu_l, err_mu_l, mu_b, err_mu_b))
        return mu_l, err_mu_l, mu_b, err_mu_b
    def load_table(s, tablename='Gaia_counterparts_for_10_PRE_bursters', tableformat='ascii'):
        tablename = s.path + '/' + tablename
        s.t = Table.read(tablename, format=tableformat)
        s.t['ra_x'] = s.t['ra_x']/15
    def write_out_new_table_summaring_l_b_mu_l_mu_d_and__v_t_pec(s, tableformat='ascii', D_step=0.1):
        s.load_table()
        ls = bs = mu_ls = err_mu_ls = mu_bs = err_mu_bs = np.array([])
        RAs, Decs, Dmins, Dmaxs = s.t['ra_x'], s.t['dec_x'], s.t['Dmin'], s.t['Dmax']
        mu_as, err_mu_as, mu_ds, err_mu_ds = s.t['pmra'], s.t['pmra_error'], s.t['pmdec'], s.t['pmdec_error']
        #Dmaxs, Dmins = s.t['Dmax'], s.t['Dmin']
        for i in range(len(s.t)):
            a = howfun.equatorial2galactic_coordinates(RAs[i], Decs[i])
            l, b = a.equatorial_position2galactic_position()
            ls, bs = np.append(ls, l), np.append(bs, b)
            mu_l, err_mu_l, mu_b, err_mu_b = s.calculate_proper_motion_in_Galactic_coordinates_with_uncertainties(\
                RAs[i], Decs[i], mu_as[i], mu_ds[i], err_mu_as[i], err_mu_ds[i])
            mu_ls, err_mu_ls, mu_bs, err_mu_bs = np.append(mu_ls, mu_l), np.append(err_mu_ls, err_mu_l), np.append(mu_bs, mu_b),\
                np.append(err_mu_bs, err_mu_b)
        #estimate v_t_pec and its l- and b-components 
        v_t_pec_mins = v_t_pec_maxs = v_b_pec_mins = v_b_pec_maxs = v_l_pec_mins = v_l_pec_maxs = np.array([])
        for i in range(len(s.t)):
            v_l_pec_min, v_l_pec_max = s.v_l_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(ls[i], bs[i],\
                mu_ls[i], err_mu_ls[i], Dmins[i], Dmaxs[i], D_step)   
            v_b_pec_min, v_b_pec_max = s.v_b_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(ls[i], bs[i],\
                mu_bs[i], err_mu_bs[i], Dmins[i], Dmaxs[i], D_step)   
            v_t_pec_min, v_t_pec_max = s.v_t_pec__and_uncertainty_for_the_first_and_fourth_Galactic_quadrants(ls[i], bs[i],\
                mu_ls[i], err_mu_ls[i], mu_bs[i], err_mu_bs[i], Dmins[i], Dmaxs[i], D_step)
            v_l_pec_mins = np.append(v_l_pec_mins, v_l_pec_min)
            v_l_pec_maxs = np.append(v_l_pec_maxs, v_l_pec_max)
            v_b_pec_mins = np.append(v_b_pec_mins, v_b_pec_min)
            v_b_pec_maxs = np.append(v_b_pec_maxs, v_b_pec_max)
            v_t_pec_mins = np.append(v_t_pec_mins, v_t_pec_min)
            v_t_pec_maxs = np.append(v_t_pec_maxs, v_t_pec_max)
        s.t3 = Table([s.t['srcname'], RAs, Decs, ls, bs, mu_ls, err_mu_ls, mu_bs, err_mu_bs,\
            v_l_pec_mins, v_l_pec_maxs, v_b_pec_mins, v_b_pec_maxs, v_t_pec_mins, v_t_pec_maxs], \
            names=['srcname', 'RA', 'Dec', 'l', 'b', 'mu_l', 'err_mu_l', 'mu_b', 'err_mu_b',\
            'v_l_pec_min', 'v_l_pec_max', 'v_b_pec_min', 'v_b_pec_max', 'v_t_pec_min', 'v_t_pec_max']) 
        output = s.path + '/space_velocities_of_PRE_bursters.' + tableformat
        s.t3.write(output, format=tableformat, overwrite=True)
    def write_out_table_summaring_l_b_mu_l_mu_d(s, tableformat='ascii', D_step=0.1):
        s.load_table()
        ls = bs = mu_ls = err_mu_ls = mu_bs = err_mu_bs = np.array([])
        RAs, Decs, Dmins, Dmaxs = s.t['ra_x'], s.t['dec_x'], s.t['Dmin'], s.t['Dmax']
        mu_as, err_mu_as, mu_ds, err_mu_ds = s.t['pmra'], s.t['pmra_error'], s.t['pmdec'], s.t['pmdec_error']
        #Dmaxs, Dmins = s.t['Dmax'], s.t['Dmin']
        for i in range(len(s.t)):
            a = howfun.equatorial2galactic_coordinates(RAs[i], Decs[i])
            l, b = a.equatorial_position2galactic_position()
            ls, bs = np.append(ls, l), np.append(bs, b)
            mu_l, err_mu_l, mu_b, err_mu_b = s.calculate_proper_motion_in_Galactic_coordinates_with_uncertainties(\
                RAs[i], Decs[i], mu_as[i], mu_ds[i], err_mu_as[i], err_mu_ds[i])
            mu_ls, err_mu_ls, mu_bs, err_mu_bs = np.append(mu_ls, mu_l), np.append(err_mu_ls, err_mu_l), np.append(mu_bs, mu_b),\
                np.append(err_mu_bs, err_mu_b)
        s.t1 = Table([s.t['srcname'], RAs, Decs, ls, bs, mu_ls, err_mu_ls, mu_bs, err_mu_bs], \
            names=['srcname', 'RA', 'Dec', 'l', 'b', 'mu_l', 'err_mu_l', 'mu_b', 'err_mu_b']) 
        output = s.path + '/transfer_to_Galactic_coordinate_proper_motions_and_positions.' + tableformat
        s.t1.write(output, format=tableformat, overwrite=True)
    def read__v_t_pec__and_estimate_weighted_average_and_weighted_std(s):
        tablename = s.path + '/space_velocities_of_PRE_bursters.ascii'
        t = Table.read(tablename, format='ascii')
        v_t_pec_mins, v_t_pec_maxs = t['v_t_pec_min'], t['v_t_pec_max']
        v_t_pecs = (v_t_pec_mins + v_t_pec_maxs)/2
        v_t_pec_errs = (v_t_pec_maxs - v_t_pec_mins)/2
        return howfun.weighted_avg_and_std(v_t_pecs, v_t_pec_errs), howfun.weightX(v_t_pecs, v_t_pec_errs)

class use_high_precision_Gaia_subset_to_estimate_parallax_zero_point:
    """
    This is an 'attempt' 
    to use Gaia DR2 sources with high-precision parallaxes to find sources with no detected proper motions (and parallaxes),
    then use these extragalactic-looking sources to constrain parallax zero point.
    It might not be of real use.
    """
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, srcname, radius=1, parallax_precision_threshold=1):
        """
        radius in arcmin, parallax_precision_threshold in mas;
        workflow is the one-stop function that chains all functions
        """
        s.srcname, s.radius, s.PPT = srcname, radius, parallax_precision_threshold
        s.cone_search_result = s.path + '/.' + s.srcname.replace(' ', '') + '_' + str(s.radius) + 'arcmin_radius_search_result'
    def workflow(s):
        import copy
        if not os.path.exists(s.cone_search_result):
            classA = plot_Gaia_sources_around_a_source(s.srcname, s.radius)
            RA_deg, Dec_deg = classA.Simbad_srcname_to_position(s.srcname)
            s.T = s.cone_search_Gaia_search_within_given_radius_in_deg(RA_deg, Dec_deg, s.radius)
            s.T.write(s.cone_search_result, format='ascii', overwrite=True)
        s.T = Table.read(s.cone_search_result, format='ascii')
        s.T0 = copy.deepcopy(s.T)
        classB = catalog_of_Gaia_counterparts_for_AGNs_to_zero_parallax_point(s.srcname)
        s.T = classB.delete_rows_where_the_required_parameter_is_absent(s.T, 'parallax')
        s.T1 = copy.deepcopy(s.T)
        s.T = s.filter_in_high_precision_subset(s.T, s.PPT)
        s.T2 = copy.deepcopy(s.T)
        s.T = s.filter_out_sources_with_detected_proper_motions(s.T, 1)
        s.T3 = copy.deepcopy(s.T)
        s.av_parallax, s.integral_error_parallax, s.std_parallax = classB.calculate_zero_parallax_and_its_sigma(s.T)
    def cone_search_Gaia_search_within_given_radius_in_deg(s, RA_deg, Dec_deg, radius):
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astroquery.gaia import Gaia
        coord = SkyCoord(ra=RA_deg, dec=Dec_deg, unit=(u.degree, u.degree), frame='icrs')
        radius_unit = u.Quantity(radius, u.arcmin)
        j = Gaia.cone_search_async(coord, radius_unit)
        return j.get_results()
    def filter_in_high_precision_subset(s, input_table, PPT):
        """
        PPT --> parallax_precision_threshold; PPT in mas
        """
        t = input_table
        index = t['parallax_error'] < PPT
        return t[index]
    def filter_out_sources_with_detected_proper_motions(s, input_table, threshold_SNR=1):
        """
        note that a parallax filter is also applied here
        """
        t = input_table
        index1 = abs(t['pmra'])/t['pmra_error'] < threshold_SNR
        index2 = abs(t['pmdec'])/t['pmdec_error'] < threshold_SNR
        index3 = t['parallax']/t['parallax_error'] < 3
        index = np.logical_and(np.array(index1), np.array(index2))
        index = np.logical_and(np.array(index), np.array(index3))
        return s.T[index]

class pulsars_based_GdotG_kD(object):
    def __init__(s):
        s.t = s.combine_pulsardats_to_table() 
        print s.t
        s.MCresults = s.foldername + '/GdotG_kD.monte.carlo'
        #s.MCresults_3paras = s.foldername + '/GdotG_kD_sp.monte.carlo'
        #s.MC_estimates_output=s.foldername+'/3para_monte_carlo.estimates.out'
        #[Xs, errXs, Ys, errYs, Zs, errZs] = s.prepare_pulsardats_for_abcfitting()
        #print s.prepare_pulsardats_for_abcfitting()
        fits = howfun.lsqfit()
        #print fits.abcfit3D(Xs, errXs, Ys, errYs, Zs, errZs)
        #print fits.abcfit3D1(Xs, errXs, Ys, errYs, Zs, errZs)
    def plot_pdfs_for_Gdot_kD_Sp(s, bins=100, HowManySigmaToDisplay=10):
        import matplotlib.pyplot as plt
        s.output_GdotG_kD_Sp_and_uncertainties()
        #plot pdfs now    
        for parameter in ['gdotg', 'kd', 'sp']:
            plt.rc('text', usetex=True)
            exec('value = s.value_%s' % parameter)
            exec('error = s.err_%s' % parameter)
            exec("minimum = %.16f" % (value-HowManySigmaToDisplay*error))
            exec("maximum = %.16f" % (value+HowManySigmaToDisplay*error))
            exec("parameters = s.%ss[s.%ss>minimum]" % (parameter,parameter))
            parameters = parameters[parameters<maximum]
            plt.hist(parameters, bins, density=True, facecolor='g', alpha=0.75)
            x = np.arange(minimum,maximum,error/4)
            exec('y = norm.pdf(x, %.16f, %.16f)' % (value, error))
            plt.plot(x,y)
            plt.ylabel('probability density')
            exec("plt.xlabel('%s')" % parameter)
            exec("plt.xlim(%.16f, %.16f)" % (minimum, maximum))
            exec("plt.savefig('%s/abcfit_%s_%sbins.eps')" % (s.foldername, parameter, bins))
            plt.clf()
        print "plots done."
    def covariance_2d_plots_with_chainconsumer(s, HowManySigma=9):
        from numpy.random import normal, multivariate_normal
        from chainconsumer import ChainConsumer
        from matplotlib import pyplot as plt
        plt.rc('text', usetex=True)
        #np.random.seed(2)
        #cov = 0.2 * normal(size=(3, 3)) + np.identity(3)
        truth = normal(size=3)
        #data = multivariate_normal(truth, np.dot(cov, cov.T), size=100000)
        s.output_GdotG_kD_Sp_and_uncertainties() 
        data = s.np2darray
        ## filter out the miss-fitted data falling into local max ############################
        for i in range(0,3):
            para_names = ['gdotg', 'kd', 'sp']
            exec("maximum = s.value_%s+HowManySigma*s.err_%s" % (para_names[i], para_names[i]))
            exec("minimum = s.value_%s-HowManySigma*s.err_%s" % (para_names[i], para_names[i]))
            index = data[:, i]<maximum
            data = data[index, :]
            index = data[:, i]>minimum
            data = data[index, :]
        c = ChainConsumer().add_chain(data, parameters=[r"$\dot G/G$", r"$\kappa_D$", r"$\alpha$"])
        fig = c.plotter.plot(truth=truth)
        fig.set_size_inches(3 + fig.get_size_inches())
        plt.savefig(s.foldername + '/abcfitting_covariance_2d_plots.eps')
        plt.clf()
        
    def output_GdotG_kD_Sp_and_uncertainties(s, inputfile='', MC_estimates_output=''):
        if inputfile == '':
            inputfile = s.MCresults_3paras
        else:
            inputfile = s.MCresults_3paras + '.new'
        if MC_estimates_output == '':
            MC_estimates_output = s.MC_estimates_output
        else:
            MC_estimates_output = s.MC_estimates_output + '.new'
        [Xs, errXs, Ys, errYs, Zs, errZs] = s.prepare_pulsardats_for_abcfitting()
        print s.prepare_pulsardats_for_abcfitting()
        fits = howfun.lsqfit()
        [s.gdotg0, s.kd0, s.sp0, junk1, junk2] = fits.abcfit3D(Xs, errXs, Ys, errYs, Zs, errZs)
        print fits.abcfit3D(Xs, errXs, Ys, errYs, Zs, errZs)
        if not os.path.exists(inputfile):
            print "no monte-carlo results; aborting"
            sys.exit()
        else:
            readfile = open(inputfile, 'r')
            s.np2darray = pickle.load(readfile)
            #print np2darray
            #print "\n\n\n"
            readfile.close()
        s.gdotgs = s.np2darray[:,0]
        s.kds = s.np2darray[:,1] 
        s.sps = s.np2darray[:,2]
        
        #MC_estimates_output = s.foldername + '/3para_monte_carlo.estimates.out'
        fileWrite = open(MC_estimates_output, 'w')
        fileWrite.write("marginalized estimates obtained with %d monte-carlo runs:\n" % len(s.np2darray))
        for i in [2,1]:
            [s.value_gdotg, s.err_gdotg] = howfun.sample2estimate(s.gdotgs, math.erf(i/2**0.5))
            #print s.err_gdotg*1e13
            #print "\n\n"
            [s.value_kd, s.err_kd] = howfun.sample2estimate(s.kds, math.erf(i/2**0.5))
            [s.value_sp, s.err_sp] = howfun.sample2estimate(s.sps, math.erf(i/2**0.5))
            print("confidencelevel = %f" % math.erf(i/2**0.5))
            fileWrite.write("confidencelevel = %f:\n" % math.erf(i/2**0.5))
            fileWrite.write("Gdot/G = %f + %f - %f (1e-13/yr)\n" % (s.gdotg0*1e13, (s.value_gdotg + s.err_gdotg - s.gdotg0)*1e13, (s.gdotg0 - s.value_gdotg + s.err_gdotg)*1e13))
            fileWrite.write("kappa_D = %f + %f -%f \n" % (s.kd0, s.value_kd+s.err_kd-s.kd0, s.kd0-s.value_kd+s.err_kd))
            fileWrite.write("s_p = %f + %f - %f\n" % (s.sp0, s.value_sp+s.err_sp-s.sp0, s.sp0-s.value_sp+s.err_sp))
        fileWrite.close()
        
    def monte_carlo_sampling_abc(s, monte_carlo_runs=100):
        [Xs, errXs, Ys, errYs, Zs, errZs] = s.prepare_pulsardats_for_abcfitting()
        fits = howfun.lsqfit()
        count = 0
        if os.path.exists(s.MCresults_3paras):
            readfile = open(s.MCresults_3paras, 'r')
            np2darray = pickle.load(readfile)
            readfile.close()

        while count < monte_carlo_runs:
            RXs = np.random.normal(Xs, errXs) #random Xs
            RYs = np.random.normal(Ys, errYs) #random Ys
            RZs = np.random.normal(Zs, errZs) #random Zs
            [GdotG, kD, Sp, junk1, junk2] = fits.abcfit3D(RXs, errXs, RYs, errYs, RZs, errZs)
            try:
                np2darray = np.row_stack((np2darray, [GdotG,kD,Sp]))
            except NameError:
                np2darray = np.array([GdotG,kD,Sp])
            count += 1
        writefile = open(s.MCresults_3paras, 'w')
        pickle.dump(np2darray, writefile)
        writefile.close()
        print len(np2darray)
        return np2darray
    def prepare_pulsardats_for_abcfitting(s): # Z=a*(c*X-1)+b*c**2*Y
        Pbs = s.Pbs/estimate_uncertainty.yr2d # in yr
        T = estimate_uncertainty.T/3600./24./estimate_uncertainty.yr2d # in yr
        Zs = -s.Pbdot_exs/2/Pbs # in 1/yr
        errZs = abs(s.err_Pbdot_exs/2/Pbs) #in 1/yr
        Xs = (1+0.5/(s.qs+1))*s.qs*s.m_cs 
        dXdMs = Xs/s.m_cs
        dXdQs = (1+0.5/(s.qs+1)**2)*s.m_cs
        errXs_sq = dXdMs**2*s.err_m_cs**2 + dXdQs**2*s.err_qs**2
        errXs = errXs_sq**0.5
        Ys = 2*math.pi**2*T/Pbs**2 # in 1/yr
        Ys *= s.qs**3/(s.qs+1)*s.m_cs**3 # 1/yr
        dYdMs = 3*Ys/s.m_cs
        dYdQs = Ys*(3./s.qs-1./(s.qs+1))
        errYs_sq = dYdMs**2*s.err_m_cs**2 + dYdQs**2*s.err_qs**2
        errYs = errYs_sq**0.5
        return Xs, errXs, Ys, errYs, Zs, errZs
    def plot_pdfs_for_Gdot_kD(s):
        s.output_GdotG_kD_and_uncertainties()
        #plot pdfs now    
        bins = 100
        for parameter in ['GdotG', 'kD']:
            plt.rc('text', usetex=True)
            exec("plt.hist(s.%ss, bins, density=True, facecolor='g', alpha=0.75)" % parameter)
            exec("minimum = min(s.%ss)" % parameter)
            exec("maximum = max(s.%ss)" % parameter)
            x = np.arange(minimum,maximum,(maximum-minimum)/1000)
            exec('y = norm.pdf(x, s.value_%s, s.err_%s)' % (parameter,parameter))
            plt.plot(x,y)
            plt.ylabel('probability density')
            exec("plt.xlabel('%s')" % parameter)
            #if parameter == 'GdotG':
            #    plt.xlim(-6e-12,3e-12)
            exec("plt.savefig('%s/%s_%sbins.eps' % (s.foldername, parameter, bins))")
            plt.clf()
        print "plots done." 
    def output_GdotG_kD_and_uncertainties(s):
        [Xs, errXs, Ys, errYs] = s.prepare_pulsardats_for_fitting()
        print s.prepare_pulsardats_for_fitting()
        fits = howfun.lsqfit()
        [s.GdotG0, s.kD0, junk1, junk2] = fits.linearfit2D(Xs, errXs, Ys, errYs)
        print fits.linearfit2D(Xs, errXs, Ys, errYs)
        if not os.path.exists(s.MCresults):
            print "no monte-carlo results; running monte-carlo..."
            s.monte_carlo_sampling(5000)
        readfile = open(s.MCresults, 'r')
        np2darray = pickle.load(readfile)
        readfile.close()
        s.GdotGs = np2darray[:,0]
        s.kDs = np2darray[:,1] 
        
        MC_estimates_output = s.foldername + '/monte_carlo.estimates.out'
        fileWrite = open(MC_estimates_output, 'w')
        fileWrite.write("marginalized estimates obtained with %d monte-carlo runs:\n" % len(np2darray))
        for i in [2,1]:
            #err_ranges = howfun.uncertainties_from_2Dsample(np2darray, i, 1)
            #[GdotG_lower, GdotG_upper] = err_ranges[:,0]
            #[kD_lower, kD_upper] = err_ranges[:,1]
            [GdotG1, errGdotG] = howfun.sample2estimate(s.GdotGs, math.erf(i/2**0.5))
            [kD1, errKD] = howfun.sample2estimate(s.kDs, math.erf(i/2**0.5))
            print("confidencelevel = %f" % math.erf(i/2**0.5))
            fileWrite.write("confidencelevel = %f:\n" % math.erf(i/2**0.5))
            fileWrite.write("Gdot/G = %f + %f - %f (1e-13/yr)\n" % (s.GdotG0*1e13, (GdotG1 + errGdotG - s.GdotG0)*1e13, (s.GdotG0 - GdotG1 + errGdotG)*1e13))
            fileWrite.write("kappa_D = %f + %f -%f \n" % (s.kD0, kD1+errKD-s.kD0, s.kD0-kD1+errKD))
        fileWrite.close()
        s.value_GdotG = GdotG1
        s.err_GdotG = errGdotG
        s.value_kD = kD1
        s.err_kD = errKD
        
    def monte_carlo_sampling(s, monte_carlo_runs=100):
        [Xs, errXs, Ys, errYs] = s.prepare_pulsardats_for_fitting()
        fits = howfun.lsqfit()
        count = 0
        if os.path.exists(s.MCresults):
            readfile = open(s.MCresults, 'r')
            np2darray = pickle.load(readfile)
            readfile.close()

        while count < monte_carlo_runs:
            RXs = np.random.normal(Xs, errXs) #random Xs
            RYs = np.random.normal(Ys, errYs) #random Ys
            [GdotG, kD, junk1, junk2] = fits.linearfit2D(RXs, errXs, RYs, errYs)
            try:
                np2darray = np.row_stack((np2darray, [GdotG,kD]))
            except NameError:
                np2darray = np.array([GdotG,kD])
            count += 1
        writefile = open(s.MCresults, 'w')
        pickle.dump(np2darray, writefile)
        writefile.close()
        print len(np2darray)
        return np2darray
    def prepare_pulsardats_for_fitting(s):
        Ys = -25*s.Pbs*s.Pbdot_exs/math.pi**2 #in d
        Ys *= 24*3600 #in s
        Ys /= estimate_uncertainty.T*s.qs**3/(s.qs+1)*s.m_cs**3 #in 1
        print Ys
        Cx = 50*s.Pbs**2/math.pi**2/estimate_uncertainty.T # in d**2/s
        Cx *= 24*3600 # in d
        Cx /= estimate_uncertainty.yr2d # in yr, so that Gdot/G is in 1/yr
        Xs = (s.qs+1)/s.qs**3/s.m_cs**3 - 0.05*(2*s.qs+3)/s.qs**2/s.m_cs**2
        Xs *= Cx
        print Xs
        
        dYdPs = Ys/s.Pbdot_exs #dY/dPbdot_ex
        dYdMs = -3*Ys/s.m_cs #dY/dm_c
        dYdQs = (1/(s.qs+1)-3/s.qs)*Ys
        errYs_sq = dYdPs**2*s.err_Pbdot_exs**2 + dYdMs**2*s.err_m_cs**2 + dYdQs**2*s.err_qs**2
        errYs = errYs_sq**0.5
        print errYs
        dXdMs = Cx
        dXdMs *= (s.qs+1)/s.qs**3/s.m_cs**3 - 0.05*(2*s.qs+3)/s.qs**2/s.m_cs**2
        dXdQs = Cx
        dXdQs *= (-2*s.qs-3)/s.qs**4/s.m_cs**3 + 0.1*(s.qs+3)/s.qs**3/s.m_cs**2
        errXs_sq = dXdQs**2*s.err_qs**2 + dXdMs**2*s.err_m_cs**2
        errXs = errXs_sq**0.5
        return Xs, errXs, Ys, errYs
    def combine_pulsardats_to_table(s):
        auxdir = os.environ['PSRVLBAUXDIR'] 
        s.foldername = auxdir + '/Gdot_kD' 
        compiled_table = s.foldername + '/psrtable.dat' 
        #if not os.path.exists(compiled_table):
        s.psrnames = np.array([])
        s.Pbs = np.array([])
        for estimate in ['Pbdot_ex', 'm_c', 'q']:
            exec('s.%ss = np.array([])' % estimate)
            exec('s.err_%ss = np.array([])' % estimate)
        
        pulsardats = glob.glob(r'%s/J*.dat' % s.foldername)
        pulsardats.sort()
        for pulsardat in pulsardats:
            [psrname, Pb, Pbdot_ex, err_Pbdot_ex, m_c, err_m_c, q, err_q] = s.read_pulsardat(pulsardat) 
            for estimate in ['psrname', 'Pb', 'Pbdot_ex', 'err_Pbdot_ex', 'm_c', 'err_m_c', 'q', 'err_q']:
                exec('s.%ss = np.append(s.%ss, %s)' % (estimate, estimate, estimate))
        t = Table([s.psrnames, s.Pbs, s.Pbdot_exs, s.err_Pbdot_exs, s.m_cs, s.err_m_cs, s.qs, s.err_qs], names=['psrname', 'Pb', 'Pbdot_ex', 'err_Pbdot_ex', 'm_c', 'err_m_c', 'q', 'err_q'])
        t.write(compiled_table, format='ascii', overwrite=True)
        #t1 = Table.read(compiled_table, format='ascii') 
        return t
    def read_pulsardat(s, pulsardat):
        for estimate in ['m_c', 'm_p', 'q']:
            exec("%s = 0" % estimate)
        lines = open(pulsardat).readlines()
        for line in lines:
            if 'psrname' in line:
                psrname = line.split('=')[-1].strip()
            if ('Pb' in line) and ('Pbdot' not in line):
                Pb = line.split('=')[-1].strip().split(' ')[0]
                Pb = float(Pb)
            for estimate in ['Pbdot_ex', 'm_c', 'm_p', 'q']:
                if estimate in line:
                    line1 = line.split('=')[-1].strip()
                    exec("%s = %s" % (estimate, line1.split('+-')[0].strip()))
                    exec("err_%s = %s" % (estimate, line1.split('+-')[-1].strip().split(' ')[0].strip()))
        if q == 0:
            q = m_p/m_c
            err_q_sq = 1/m_c**2*err_m_p**2 + m_p**2/m_c**4*err_m_c**2
            err_q = err_q_sq**0.5
        if m_c == 0:
            m_c = m_p/q
            err_m_c_sq = 1/q**2*err_m_p**2 + m_p**2/q**4*err_q**2
            err_m_c = err_m_c_sq**0.5
        return psrname, Pb, Pbdot_ex, err_Pbdot_ex, m_c, err_m_c, q, err_q

class pulsars_based_GdotG_kD_Sp_via_monte_carlo_from_base_observables(pulsars_based_GdotG_kD):
    def __init__(s):
        super(pulsars_based_GdotG_kD_Sp_via_monte_carlo_from_base_observables, s).__init__() #python2 way to use super
        s.MC_estimates_output_new = s.MC_estimates_output + '.new'
    def monte_carlo_base_observables(s):
        import copy
        t1 = copy.deepcopy(s.t)
        for parameter in ['Pbdot_ex', 'm_c', 'q']:
            err = 'err_' + parameter
            t1[parameter] = np.random.normal(s.t[parameter], s.t[err])
        return t1
    def get_X_Y_Z(s):
        t1 = s.monte_carlo_base_observables()
        for parameter in ['Pbdot_ex', 'm_c', 'q']:
            exec("%ss = t1['%s']" % (parameter, parameter))
        Pbs = s.Pbs/estimate_uncertainty.yr2d # in yr
        T = estimate_uncertainty.T/3600./24./estimate_uncertainty.yr2d # in yr
        Zs = -Pbdot_exs/2/Pbs # in 1/yr
        errZs = abs(s.err_Pbdot_exs/2/Pbs) #in 1/yr
        Xs = (1+0.5/(qs+1))*qs*m_cs 
        dXdMs = Xs/m_cs
        dXdQs = (1+0.5/(qs+1)**2)*m_cs
        errXs_sq = dXdMs**2*s.err_m_cs**2 + dXdQs**2*s.err_qs**2
        errXs = errXs_sq**0.5
        Ys = 2*math.pi**2*T/Pbs**2 # in 1/yr
        Ys *= qs**3/(qs+1)*m_cs**3 # 1/yr
        dYdMs = 3*Ys/m_cs
        dYdQs = Ys*(3./qs-1./(qs+1))
        errYs_sq = dYdMs**2*s.err_m_cs**2 + dYdQs**2*s.err_qs**2
        errYs = errYs_sq**0.5
        return Xs, errXs, Ys, errYs, Zs, errZs
    def monte_carlo_GdotG_kD_Sp(s, monte_carlo_runs=100):
        fits = howfun.lsqfit()
        s.MCresults_3paras_new = s.MCresults_3paras + '.new' 
        count = 0
        if os.path.exists(s.MCresults_3paras_new):
            readfile = open(s.MCresults_3paras_new, 'r')
            np2darray = pickle.load(readfile)
            readfile.close()

        while count < monte_carlo_runs:
            [Xs, errXs, Ys, errYs, Zs, errZs] = s.get_X_Y_Z()
            [GdotG, kD, Sp, junk1, junk2] = fits.abcfit3D(Xs, errXs, Ys, errYs, Zs, errZs)
            try:
                np2darray = np.row_stack((np2darray, [GdotG,kD,Sp]))
            except NameError:
                np2darray = np.array([GdotG,kD,Sp])
            count += 1
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round(count * 100 / monte_carlo_runs)) + " \r")
        writefile = open(s.MCresults_3paras_new, 'w')
        pickle.dump(np2darray, writefile)
        writefile.close()
        print len(np2darray)
        return np2darray

class simulate_fake_pulsars_for_abcfitting_check(pulsars_based_GdotG_kD):
    def __init__(s, number=100):
        super(simulate_fake_pulsars_for_abcfitting_check, s).__init__() #python2 way to use super
        s.number = number
    def monte_carlo_Pb_Mp_q(s, number=100):
        number = s.number
        m_ps = np.random.uniform(1.1, 2.1, number)
        qs = np.random.uniform(1, 11, number)
        Pbs = np.random.uniform(0.3, 70, number)
        return m_ps, qs, Pbs
    def calculate_Pbdot_ex_with_Pb_Mp_q(s, GdotG, kD, alpha):
        print GdotG, kD, alpha
        [m_ps, qs, Pbs1] = s.monte_carlo_Pb_Mp_q()
        Pbs = Pbs1/estimate_uncertainty.yr2d # in yr
        T = estimate_uncertainty.T/3600./24./estimate_uncertainty.yr2d # in yr
        Pbdot_dipoles = -4*math.pi**2*T*m_ps/Pbs/(qs+1)*kD*(alpha*m_ps)**2
        Pbdot_GdotGs = -2*GdotG*Pbs
        Pbdot_GdotGs *= 1 - alpha*m_ps*(1+0.5/(1+qs))
        Pbdot_exs = Pbdot_dipoles + Pbdot_GdotGs
        Pbdot_exs = np.random.normal(Pbdot_exs, abs(Pbdot_exs))
        return m_ps, qs, Pbs1, Pbdot_exs
    def write_out_fake_pulsar_inputfiles(s, GdotG=-73.25187e-12, kD=-9.042974297540898e-07, alpha=0.697557, number=100):
        number = s.number
        [m_ps, qs, Pbs, Pbdot_exs] = s.calculate_Pbdot_ex_with_Pb_Mp_q(GdotG, kD, alpha)
        err_Pbdot_exs = np.random.uniform(0.7, 5, number)
        err_Pbdot_exs *= abs(Pbdot_exs)
        err_qs = 0.1*qs
        err_m_ps = 0.1*m_ps
        for i in range(number):
            outfile = s.foldername + '/Jfake_pulsar' + str(i+1) + '.dat'
            print outfile
            fileWrite = open(outfile, 'w')
            fileWrite.write("psrname = fale_pulsar%s\n" % str(i+1))
            fileWrite.write("Pb = %.2f (d)\n" % Pbs[i])
            fileWrite.write("Pbdot_ex = %.60f +- %.60f\n" % (Pbdot_exs[i], err_Pbdot_exs[i]))
            fileWrite.write("q = %.2f +- %.2f\n" % (qs[i], err_qs[i]))
            fileWrite.write("m_p = %.2f +- %.2f\n" % (m_ps[i], err_m_ps[i]))
            fileWrite.close()
        print "All fake_pulsar*.dat generated."
    def estimate_simulation_uncertainties(s, GdotG=-73.25187e-12, kD=-9.042974297540898e-07, alpha=0.697557, runs=100):
        fits = howfun.lsqfit()
        count = 0
        gdotgs = np.array([])
        kds = np.array([])
        sps = np.array([])
        while count<runs:
            s.write_out_fake_pulsar_inputfiles(GdotG, kD, alpha)
            s.combine_pulsardats_to_table()
            [Xs, errXs, Ys, errYs, Zs, errZs] = s.prepare_pulsardats_for_abcfitting()
            [gdotg, kd, sp, funk1, funk2] = fits.abcfit3D(Xs, errXs, Ys, errYs, Zs, errZs)
            for parameter in ['gdotg', 'kd', 'sp']:
                exec('%ss = np.append(%ss, %s)' % (parameter, parameter, parameter))
            count += 1
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round(count * 100 / runs)) + " \r")
        
        print sps
        simulation_estimates_output = s.foldername + '/simulation.estimates.out'
        fileWrite = open(simulation_estimates_output, 'w')
        fileWrite.write("marginalized estimates obtained with %d simulation runs (where fake pulsars are used):\n" % runs)
        for i in [2,1]:
            for parameter in ['gdotg', 'kd', 'sp']:
                exec('[value_%s, err_%s] = howfun.sample2estimate(%ss, math.erf(i/2**0.5))' % (parameter, parameter, parameter))
            fileWrite.write("confidencelevel = %f:\n" % math.erf(i/2**0.5))
            fileWrite.write("Gdot/G = %.30f +- %.30f (1e-13/yr)\n" % (value_gdotg*1e13, err_gdotg*1e13))
            fileWrite.write("kappa_D = %.30f +- %.30f \n" % (value_kd, err_kd))
            fileWrite.write("s_p = %f +- %f\n" % (value_sp, err_sp))
        fileWrite.close()

    def bootstrapping_of_abcfitting_with_simulated_pulsars_plus_real_ones(s, bootstrapruns=100): #this code is unfinished!!!
   #this code is unfinished!!
        import random
        length = len(s.t)
        draws = int(length/2) #half of the sample volume makes the most variable combination
        gdotgs = np.array([])
        kds = np.array([])
        alphas = np.array([])
        count = 0
        while count<bootstrapruns:
            indexes = random.sample(range(length), draws)
            count += 1

class find_virtual_calibrator_position_with_colinear_calibrators:
    """
    Note
    ----
    1. In calibrator_search_mode, either srcname is in the form of 'HH:MM:SS.SSS,dd:mm:ss.ss'.
    2. When self.inverse_referencing==True, target and phscal2 are switched.
    
    Input parameters
    ----------------
    kwargs :
        1. inverse_referencing : bool (default : False)
            When inverse_referencing==True, target and phscal2 are switched.
    
    Return parameters
    -----------------
    return sep_V2T : float (in arcmin)
        Angular separation between the target and the virtual calibrator.
    inbeamsn_ratio : float
        The ratio by which the phase solutions are multiplied.
    """
    def __init__(s, targetname, phscal1name, phscal2name, **kwargs): ## phscal2 is the "inbeamcal"
        s.calibrator_search_mode = False
        if prepare_path_source(targetname) == False:
            s.calibrator_search_mode = True
        try:
            s.inverse_referencing = kwargs['inverse_referencing']
        except KeyError:
            s.inverse_referencing = False
        s.psrname = targetname
        s.phscal1name = phscal1name
        if s.inverse_referencing:
            s.targetname = phscal2name
            s.phscal2name = targetname
        else:
            s.targetname = targetname
            s.phscal2name = phscal2name
        #print s.targetname, s.phscal1name
        s.prepare_positions()
        s.quantify_the_plane1_paramters_defined_by_two_phscals_and_000point()
        s.get_the_plane2_perpendicular_to_plane1_and_pass_target_and_000point()
        s.solve_the_position_of_virtual_phscal()
        s.separation_between_virtual_phscal_and_sources()
    def prepare_positions(s):
        for source in ['target', 'phscal1', 'phscal2']:
            if not s.calibrator_search_mode:
                exec('[s.RA_%s, s.Dec_%s] = srcposition(s.psrname, s.%sname)' % (source, source, source))
            else:
                exec('s.RA_%s = s.%sname.split(",")[0].strip()' % (source, source))
                exec('s.Dec_%s = s.%sname.split(",")[1].strip()' % (source, source))
            exec('s.RA_deg_%s = 15*howfun.dms2deg(s.RA_%s)' % (source, source))
            exec('s.Dec_deg_%s = howfun.dms2deg(s.Dec_%s)' % (source, source))
            exec('s.RA_rad_%s = math.pi/180*s.RA_deg_%s' % (source, source))
            exec('s.Dec_rad_%s = math.pi/180*s.Dec_deg_%s' % (source, source))
    def quantify_the_plane1_paramters_defined_by_two_phscals_and_000point(s): #solve a system of two Ax+By=z eqs
        phis_phscal = np.array([s.RA_rad_phscal1, s.RA_rad_phscal2])
        thetas_phscal = math.pi/2 - np.array([s.Dec_rad_phscal1, s.Dec_rad_phscal2])
        Xs_phscal = np.sin(thetas_phscal) * np.cos(phis_phscal)
        Ys_phscal = np.sin(thetas_phscal) * np.sin(phis_phscal)
        Zs_phscal = np.cos(thetas_phscal)
        M1 = np.row_stack((Xs_phscal, Ys_phscal))
        M1 = np.mat(M1)
        M1 = M1.T
        M2 = np.mat(Zs_phscal)
        M2 = M2.T
        print M1
        print M2
        solution = M1.I * M2
        s.A = solution[0,0]
        s.B = solution[1,0] #the vector for plane1 is (A,B,1)
        print s.A, s.B
        return s.A, s.B
    def get_the_plane2_perpendicular_to_plane1_and_pass_target_and_000point(s):
        phi_target = s.RA_rad_target
        theta_target = math.pi/2 - s.Dec_rad_target
        X_target = math.sin(theta_target) * math.cos(phi_target)
        Y_target = math.sin(theta_target) * math.sin(phi_target)
        Z_target = math.cos(theta_target)
        s.Z_target = Z_target
        print X_target,Y_target,Z_target
        M1 = np.mat([[s.A,      s.B     ],
                     [X_target, Y_target]])
        M2 = np.mat([[-1      ],
                     [Z_target]])
        solution = M1.I * M2
        s.a = solution[0,0]
        s.b = solution[1,0]
        print s.a, s.b
        return s.a, s.b
    def solve_the_position_of_virtual_phscal(s):
        C1 = np.mat([[s.A, s.B],
                     [s.a, s.b]])
        c1 = np.linalg.det(C1)
        X = s.b - s.B
        Y = s.A - s.a
        Z = c1
        R = (X**2 + Y**2 + Z**2)**0.5
        X /= R
        Y /= R
        Z /= R
        if Z*s.Z_target<0: #this way might introduce bug in the vicinity of equator
            X *= -1
            Y *= -1
            Z *= -1
        theta_rad = math.acos(Z)
        if not Y/math.sin(theta_rad)<0:
            phi_rad = math.acos(X/math.sin(theta_rad))
        else:
            phi_rad = 2*math.pi - math.acos(X/math.sin(theta_rad))
   
        Dec_rad = math.pi/2 - theta_rad
        RA_h = phi_rad/math.pi*180/15
        Dec_deg = Dec_rad/math.pi*180
        print RA_h, Dec_deg
        s.RA_hms_V = howfun.deg2dms(RA_h)
        s.Dec_dms_V = howfun.deg2dms(Dec_deg)
        print s.RA_hms_V, s.Dec_dms_V
        return s.RA_hms_V, s.Dec_dms_V
    def separation_between_virtual_phscal_and_sources(s):
        ## P1 is the main phscal, P2 is the "inbeam"
        sep_V2T = howfun.separation(s.RA_target, s.Dec_target, s.RA_hms_V, s.Dec_dms_V)
        sep_V2P1 = howfun.separation(s.RA_phscal1, s.Dec_phscal1, s.RA_hms_V, s.Dec_dms_V) 
        sep_V2P2 = howfun.separation(s.RA_phscal2, s.Dec_phscal2, s.RA_hms_V, s.Dec_dms_V)
        sep_P1toP2 = howfun.separation(s.RA_phscal2, s.Dec_phscal2, s.RA_phscal1, s.Dec_phscal1)
        #print sep_V2T, sep_V2P1, sep_V2P2 #, sep_P1toP2, sep_P1toP2-sep_V2P1-sep_V2P2
        tolerance = 1e-5
        if abs(sep_V2P1 + sep_P1toP2 - sep_V2P2) < tolerance:
            inbeamsn_ratio = -sep_V2P1/sep_P1toP2
        elif abs(sep_V2P1 - sep_P1toP2 - sep_V2P2) < tolerance or abs(sep_P1toP2 - sep_V2P1 - sep_V2P2) < tolerance:
            inbeamsn_ratio = sep_V2P1/sep_P1toP2
        else:
            print("The calculation doesn't pass self-examination. Check the code!")
            sys.exit()
        print("inbeamsn_ratio = " + str(inbeamsn_ratio) + ", where the inbeamcal is phscal2.")
        print("The calibrator throw from the virtual calibrator is " + str(sep_V2T) + " arcmin.")
        return sep_V2T, inbeamsn_ratio
    
class solve_the_two_correction_factors_in_2D_interpolation:
    def __init__(s, target, src1, src2, src3):
        notice = "Note that the order of the three provided SRCs is important. That implies the first correction factor is calculated on the line of src1 and src2, and the second correction factor is calculated on the line of the target and src3"
        print notice
        for src in ['target', 'src1', 'src2', 'src3']:
            exec("s.%s = %s" % (src, src))
        s.prepare_positions_for_the_calculation()
        s.convert_to_cartesian_coordinate()
        s.calculate_the_cross_point_position()
        s.estimate_the_two_ratios()
    def prepare_positions_for_the_calculation(s):
        for src in ['target', 'src1', 'src2', 'src3']:
            exec("[s.RA_%s, s.Dec_%s] = srcposition(s.target, s.%s)" % (src, src, src))
            exec("phi_deg_%s = 15*howfun.dms2deg(s.RA_%s)" % (src, src))
            exec("theta_deg_%s = 90 - howfun.dms2deg(s.Dec_%s)" % (src, src))
            for parameter in ['phi', 'theta']:
                exec("s.%s_rad_%s = %s_deg_%s/180*math.pi" % (parameter, src, parameter, src))
        print s.phi_rad_src1, s.theta_rad_target
    def convert_to_cartesian_coordinate(s):
        for src in ['target', 'src1', 'src2', 'src3']:
            exec("s.x_%s = math.sin(s.theta_rad_%s) * math.cos(s.phi_rad_%s)" % (src, src, src))
            exec("s.y_%s = math.sin(s.theta_rad_%s) * math.sin(s.phi_rad_%s)" % (src, src, src))
            exec("s.z_%s = math.cos(s.theta_rad_%s)" % (src, src))
        print s.z_target, s.x_src2
    def calculate_the_cross_point_position(s):
        a11 = np.mat([[s.y_src1, s.z_src1],
                      [s.y_src2, s.z_src2]])
        a11 = np.linalg.det(a11)
        a12 = np.mat([[s.z_src1, s.x_src1],
                      [s.z_src2, s.x_src2]])
        a12 = np.linalg.det(a12)
        b1  = np.mat([[s.y_src1, s.x_src1],
                      [s.y_src2, s.x_src2]])
        b1  = np.linalg.det(b1)
        a21 = np.mat([[s.y_src3,   s.z_src3],
                      [s.y_target, s.z_target]])
        a21 = np.linalg.det(a21)
        a22 = np.mat([[s.z_src3,   s.x_src3],
                      [s.z_target, s.x_target]])
        a22 = np.linalg.det(a22)
        b2  = np.mat([[s.y_src3,   s.x_src3],
                      [s.y_target, s.x_target]])
        b2  = np.linalg.det(b2)

        A = np.mat([[a11, a12],
                    [a21, a22]])
        B = np.mat([[b1],
                    [b2]])
        solution = A.I * B
        XoverZ = solution[0,0]
        YoverZ = solution[1,0]
        Z = (XoverZ**2 + YoverZ**2 + 1)**(-0.5)
        if Z*s.z_target < 0: #this might fail in the vicinity of the equator
            Z *= -1
        X = Z * XoverZ
        Y = Z * YoverZ
        theta = math.acos(Z)
        phi = math.atan2(Y/math.sin(theta), X/math.sin(theta))
        if phi < 0:
            phi += 2*math.pi
        
        phi_deg = phi*180/math.pi
        theta_deg = theta*180/math.pi
        print phi_deg/15, 90-theta_deg
        s.RA = howfun.deg2dms(phi_deg/15)
        s.Dec = howfun.deg2dms(90-theta_deg)
        print s.RA, s.Dec
        return s.RA, s.Dec
    def estimate_the_two_ratios(s):
        convention = "The correction is done in two steps. In the first step, phase_cross = phase_src2+factor1*(phase_src1-phase_src2). In the second step, phase_ultimate = phase_src3+factor2*(phase_cross-phase_src3)."
        print convention
        sep1 = howfun.separation(s.RA, s.Dec, s.RA_src1, s.Dec_src1) 
        sep2 = howfun.separation(s.RA, s.Dec, s.RA_src2, s.Dec_src2)
        sep12 = howfun.separation(s.RA_src1, s.Dec_src1, s.RA_src2, s.Dec_src2)
        #print sep1, sep2, sep12
        if abs(sep1-sep2-sep12) < 1e-5:
            s.ratio1 = -sep2/sep12
        else:
            s.ratio1 = sep2/sep12
        
        sep3 = howfun.separation(s.RA, s.Dec, s.RA_src3, s.Dec_src3)
        sept = howfun.separation(s.RA, s.Dec, s.RA_target, s.Dec_target)
        sep3t = howfun.separation(s.RA_src3, s.Dec_src3, s.RA_target, s.Dec_target)
        #print sep3, sept, sep3t, sep3+sept-sep3t
        if abs(sept-sep3-sep3t) < 1e-5:
            s.ratio2 = -sep3t/sep3
        else:
            s.ratio2 = sep3t/sep3
        print s.ratio1, s.ratio2

class generatepmparin:
    """
    Functions:
    1) generate pmpar.in.preliminary or pmpar.in, using the 'write_out_preliminary/final_pmpar_in' function group
    2) while estimating systematics, you can choose dual-phscal mode
    3) then you can bootstrap with the pmparinfile
    4) in the end, you can make plots with the bootstrap instograms
    5) test potential outlying of a specific epoch, using the 'plot_to_justify_potential_outlier_using_bootstrapping_given_expno' function group
    6) make corner plots of the three or five astrometric parameters using the 'covariance_2d_plots_with_chainconsumer' function
    7) to make pmparin in the case of inverse referencing with respect to more than one IBCs, one needs to change 'primartytarget' (prIBC)
       through all IBCs; each time run and create a pmparin (either pmpar.in or pmpar.in.preliminary)

    Usage instructions:
    The functions are largely organized in the order of running: write_out_preliminary/final --> bootstrap_pmpar --> boostrapped_sample2estimates
    --> make plots

    Input:
    e.g. exceptions=['bh142','bh145a']

    Reference epoch:
    It will be automatically determined as the median of the epochs, and rounded to an integer.
    """
    def __init__(s, targetname, exceptions='', dualphscal=False, dualphscalratio=1, inverse_referencing=False):
        s.targetname = targetname
        s.exceptions = exceptions
        s.dualphscal = dualphscal
        s.dualphscalratio = dualphscalratio
        s.inverse_referencing = inverse_referencing
        #s.epoch = epoch
        [auxdir, s.configdir, s.targetdir, s.phscalname, s.prIBCname] = prepare_path_source(targetname, s.inverse_referencing)
        s.pmparesultsdir = s.targetdir + '/pmparesults/'
    def find_statsfiles(s, check_inbeam='', **kwargs):
        """
        Input parameters
        ----------------
        kwargs: key=value
            1) search_keyword : str
                A string used to look for statsfiles, e.g. 'preselfcal'
        """
        try:
            search_keyword = kwargs['search_keyword']
        except KeyError:
            search_keyword = ''
        usage = "When search_keyword != '', one needs to assign check_inbeam to search for statsfiles. In other words, in this case, check_inbeam does not necessarily refer to inbeamsrc."

        targetdir = s.targetdir
        if not os.path.exists(targetdir):
            print("%s doesn't exist; aborting\n" % targetdir)
            sys.exit()
        if search_keyword == '':
            if not s.inverse_referencing:
                if check_inbeam == '':
                    s.statsfiles=glob.glob(r'%s/*/*.gated.difmap.jmfit.stokesi.stats' % targetdir) #find statsfile for each epoch
                else:
                    s.statsfiles=glob.glob(r'%s/*/*_%s_preselfcal.difmap.jmfit*.stats' % (targetdir, check_inbeam)) #find statsfile for each epoch
            else:
                s.statsfiles = glob.glob(r'%s/*/*_%s_divided.difmap.jmfit.stokesi.stats' % (targetdir, s.prIBCname))
        elif search_keyword != '' and check_inbeam != '':
            s.statsfiles = glob.glob(r'%s/*/*%s*%s*.stats' % (targetdir, check_inbeam, search_keyword))
        else:
            print(usage)
            sys.exit()
        s.statsfiles.sort()
    def statsfile2expno_and_decyear(s, statsfile):
        b = plot_position_scatter(s.targetname, s.exceptions)
        [expno, decyear] = b.statsfile2expno_and_decyear(statsfile)
        return expno, decyear
    def statsfile2position_and_its_error(s, statsfile): 
        statsline = open(statsfile).readlines()[-7:]
        for line in statsline:
            if 'Actual RA' in line:
                RA = line.split('RA:')[-1].strip()
                print RA
            if 'Actual Dec' in line:
                Dec = line.split('Dec:')[-1].strip() 
                print Dec
            if 'RA error (hms)' in line:
                error0RA = float(line.split('(hms):')[-1].strip())*0.001 #unit: s
            if 'Dec error (mas)' in line:
                error0Dec = float(line.split('(mas):')[-1].strip())*0.001 #unit: "
        if s.exceptions != '':
            [expno, junk] = s.statsfile2expno_and_decyear(statsfile)
            if expno in s.exceptions:
                [RA, error0RA, Dec, error0Dec] = align_position_in_adhoc_experiment(RA, error0RA, Dec, error0Dec, s.targetname, s.exceptions, s.dualphscal, s.dualphscalratio)
        return RA, error0RA, Dec, error0Dec
    def statsfiles2median_decyear2epoch(s, statsfiles):
        from astropy.time import Time
        b = plot_position_scatter(s.targetname, s.exceptions)
        decyears = b.statsfiles2decyears(statsfiles).astype(float)
        median_decyear = howfun.sample2median(decyears)
        median = Time(median_decyear, format='decimalyear')
        return round(median.mjd)
    def write_out_pmpar_in(s, statsfiles, pulsitions):
        """
        Input parameters
        ----------------
        statsfiles : list of str
            A list of statsfile entailing the position information.
        pulsitions : str
            The output pmpar.in file contain pulsar positions (pulsitions).
        """
        if not os.path.exists(s.pmparesultsdir):
            os.system('mkdir %s' % s.pmparesultsdir)
        inverse_referenced_to = ''
        if s.inverse_referencing:
            inverse_referenced_to = '.to.' + s.prIBCname
        s.RAs = np.array([])
        s.Decs = np.array([])
        s.error0RAs = np.array([])
        s.error0Decs = np.array([])
        s.expnos = np.array([])
        s.decyears = np.array([])
        s.epoch = s.statsfiles2median_decyear2epoch(statsfiles)
        fileWrite = open(pulsitions, 'w')
        if s.inverse_referencing:
            fileWrite.write("### positions for %s inverse-referencing to %s\n" % (s.prIBCname, s.targetname))
        else:
            fileWrite.write("### positions for %s\n" % s.targetname)
        fileWrite.write("### pmpar format\n")
        fileWrite.write("#name = " + s.targetname + "\n")
        fileWrite.write("#ref = " + s.targetname + "\n")
        fileWrite.write("epoch = %f\n" % s.epoch)
        fileWrite.write("#pi = 0\n")
        fileWrite.write("# decimalyear RA +/- Dec +/-\n")
        for statsfile in statsfiles:
            [expno, decyear] = s.statsfile2expno_and_decyear(statsfile)
            s.expnos = np.append(s.expnos, expno)
            s.decyears = np.append(s.decyears, decyear)
            [RA, error0RA, Dec, error0Dec] = s.statsfile2position_and_its_error(statsfile)
            s.RAs = np.append(s.RAs, RA)
            s.Decs = np.append(s.Decs, Dec)
            s.error0RAs = np.append(s.error0RAs, error0RA)
            s.error0Decs = np.append(s.error0Decs, error0Dec)
            fileWrite.write("%s %s %.7f %s %.6f\n" % (decyear, RA, error0RA, Dec, error0Dec))
        fileWrite.close()
        s.nepoch = len(s.RAs)
        
    def write_out_preliminary_pmpar_in(s, check_inbeam=''):
        """
        for example, check_inbeam='IBC01433' means looking at bd1*IBC01433_preselfcal*stats and make the corresponding IBC01433.pmpar.in
        """
        if check_inbeam == '':
            pulsitions = s.pmparesultsdir + '/' + s.targetname + inverse_referenced_to + '.pmpar.in.preliminary'
        else:
            pulsitions = s.pmparesultsdir + '/' + check_inbeam + inverse_referenced_to + '.pmpar.in.preliminary'
        s.find_statsfiles(check_inbeam)
        s.write_out_pmpar_in(s.statsfiles, pulsitions)
    def write_out_pmpar_in_given_srcname_and_pmparinsuffix(s, srcname, search_keyword, pmparinsuffix=''):
        """
        Functions
        ---------
        Create a pmpar.in file by searching statsfiles containing srcname and search_keyword info.
        
        Input parameters
        ----------------
        srcname : str
            Source name info.
        search_keyword : str
            A keyword in statsfile used to pinpoint statsfiles.
        pmparinsuffix : str (default: '')
            Suffix for the output pmpar.in file; when it equals '', the search_keyword will be used as the suffix.
        """
        if pmparinsuffix == '':
            pmparinsuffix = search_keyword
        if (pmparinsuffix != '') and (not pmparinsuffix.startswith('.')):
            pmparinsuffix = '.' + pmparinsuffix    
        pulsitions = s.pmparesultsdir + '/' + srcname + '.pmpar.in' + pmparinsuffix
        s.find_statsfiles(srcname, search_keyword=search_keyword)
        s.write_out_pmpar_in(s.statsfiles, pulsitions)

    def write_out_pmparin_incl_sysErr(s, pmparesultsdir, targetname, pulsition_suffix, nepoch, epoch, decyears, expnos, RAs, error0RAs, Decs, error0Decs, dualphscal, paraA_rchsq, **kwargs):
        """
        Input parameters
        ----------------
        kwargs : 
            1. inverse_referencing : bool (default : False)
                This is not actually necessary since self.inverse_referencing can be called.
                But the introduction of this kwarg allow extracting this function from the class.
        """
        try:
            inverse_referencing = kwargs['inverse_referencing']
        except KeyError:
            inverse_referencing = False
        errorRAs   = np.array([])
        errorDecs  = np.array([])
        if pulsition_suffix != '':
            pulsition_suffix = '.' + pulsition_suffix
        inverse_referenced_to = ''
        if inverse_referencing:
            inverse_referenced_to = '.to.' + s.prIBCname
        pulsitions = pmparesultsdir + '/' + targetname + inverse_referenced_to + '.pmpar.in' + pulsition_suffix
        fileWrite = open(pulsitions, 'w')
        if s.inverse_referencing:
            fileWrite.write("### positions for %s inverse-referencing to %s\n" % (s.prIBCname, s.targetname))
        else:
            fileWrite.write("### positions for %s\n" % s.targetname)
        fileWrite.write("### pmpar format\n")
        fileWrite.write("#name = " + targetname + "\n")
        fileWrite.write("#ref = " + targetname + "\n")
        fileWrite.write("epoch = %f\n" % epoch)
        fileWrite.write("#pi = 0\n")
        fileWrite.write("# decimalyear RA +/- Dec +/-\n")
        for i in range(nepoch):
            sysError = expno_sysErr(expnos[i], dualphscal, paraA_rchsq, inverse_referencing)
            [sysErrRA_mas, sysErrDec_mas, sysErrRA_ms] = sysError.sysErr()
            sysErrRA_s = sysErrRA_ms/1000
            sysErrDec_as = sysErrDec_mas/1000
            errorRA = (error0RAs[i]**2 + sysErrRA_s**2)**0.5 # in s
            errorRAs = np.append(errorRAs, errorRA)
            errorDec = (error0Decs[i]**2 + sysErrDec_as**2)**0.5 # in "
            errorDecs = np.append(errorDecs, errorDec)
            fileWrite.write("%s %s %.7f %s %.6f\n" % (decyears[i], RAs[i], errorRA, Decs[i], errorDec))
        fileWrite.close()
        return pulsitions, errorRAs, errorDecs
    def write_out_pmparin_incl_sysErr_two_paraA_rchsq(s, pmparesultsdir, targetname, exceptions, pulsition_suffix, nepoch, epoch, decyears, expnos, RAs, error0RAs, Decs, error0Decs, dualphscal, paraA1_rchsq, paraA2_rchsq):
        errorRAs   = np.array([])
        errorDecs  = np.array([])
        if pulsition_suffix != '':
            pulsition_suffix = '.' + pulsition_suffix
        pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in' + pulsition_suffix
        fileWrite = open(pulsitions, 'w')
        fileWrite.write("### pmpar format\n")
        fileWrite.write("#name = " + targetname + "\n")
        fileWrite.write("#ref = " + targetname + "\n")
        fileWrite.write("epoch = %f\n" % epoch)
        fileWrite.write("#pi = 0\n")
        fileWrite.write("# decimalyear RA +/- Dec +/-\n")
        for i in range(nepoch):
            if expnos[i] in exceptions:
                paraA_rchsq = paraA1_rchsq
            else:
                paraA_rchsq = paraA2_rchsq
            sysError = expno_sysErr(expnos[i], dualphscal, paraA_rchsq)
            [sysErrRA_mas, sysErrDec_mas, sysErrRA_ms] = sysError.sysErr()
            sysErrRA_s = sysErrRA_ms/1000
            sysErrDec_as = sysErrDec_mas/1000
            errorRA = (error0RAs[i]**2 + sysErrRA_s**2)**0.5 # in s
            errorRAs = np.append(errorRAs, errorRA)
            errorDec = (error0Decs[i]**2 + sysErrDec_as**2)**0.5 # in "
            errorDecs = np.append(errorDecs, errorDec)
            fileWrite.write("%s %s %.7f %s %.6f\n" % (decyears[i], RAs[i], errorRA, Decs[i], errorDec))
        fileWrite.close()
        return pulsitions, errorRAs, errorDecs
    ## get parallax, mu_a and mu_d
    def pulsitions2paras(s, pulsitions, pmparesultsdir, targetname):
        pmparout = pmparesultsdir + '/' + targetname + '.pmpar.out'
        os.system("pmpar %s > %s" % (pulsitions, pmparout))
        [RA0, Dec0, junk1, PI0, mu_a0, mu_d0, junk2, junk3, junk4, junk5, junk6, rchsq] = readpmparout(pmparout)
        D0 = 1/PI0
        return D0, PI0, mu_a0, mu_d0, RA0, Dec0, rchsq
    def write_out_final_pmpar_in(s, paraA_rchsq=1e-3, paraA_rchsq_step=1e-3, paraA1_rchsq=3.102e-4):
        """
        Note
        ----
        1. inverse_referencing is directly feeded into self.write_out_pmparin_incl_sysErr() as a kwarg.
        """
        s.write_out_preliminary_pmpar_in()
        if not s.dualphscal:
            [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, '', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, 1e-3, inverse_referencing=s.inverse_referencing)
            [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
        if s.dualphscal:
            [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, 'use.A1', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA_rchsq, inverse_referencing=s.inverse_referencing)
            #[pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr_two_paraA_rchsq(s.pmparesultsdir, s.targetname, s.exceptions, 'unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA1_rchsq, paraA_rchsq)
            [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
            while False: #iteratively getting paraA_rchsq is turned off
                if rchsq < 1:
                    print("Reduced chi-square already less than unity without systematics; aborting")
                    sys.exit()
                while rchsq > 1:
                    paraA_rchsq += paraA_rchsq_step
                    [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, '.unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA_rchsq, inverse_referencing=s.inverse_referencing)
                    #[pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr_two_paraA_rchsq(s.pmparesultsdir, s.targetname, s.exceptions, 'unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA1_rchsq, paraA_rchsq)
                    [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
        return D0, PI0, mu_a0, mu_d0, RA0, Dec0, rchsq, paraA_rchsq
    def parse_prior_for_bootstrap_pmpar(s, prior, bootstrapruns):
        variable, value = prior.split('=')
        if '+-' in value:
            value, error = value.split('+-')
        else:
            error = 0
        values = np.random.normal(float(value), float(error), bootstrapruns)
        str_values = map(str, values)
        prior_randoms = [variable + '=' + str_value for str_value in str_values]
        return prior_randoms
    def parse_priors_for_bootstrap_pmpar(s, priors, bootstrapruns):
        """
        example of priors: 'mu_a=3.1+-0.2 mu_d=-1.2+-0.1' or 'mu_a=3.1 mu_d=-1.2'
        the output is a 2D string array
        """
        prior_strings = []
        for prior in priors.split(' '):
            prior_randoms = s.parse_prior_for_bootstrap_pmpar(prior, bootstrapruns)
            prior_strings.append(prior_randoms)
        return prior_strings
        
    def bootstrap_pmpar(s, pmparinfile, bootstrapruns, priors='', overwrite_table=False):
        """
        priors can be 'mu_a=3.1 mu_d=-2.3' or 'mu_a=3.1+-0.2 mu_d=-2.3+-0.3' or 'mu_a=3.1 mu_d=-2.3+-0.3' or just 'mu_a=2.3+-0.1'
        """
        from astropy.table import vstack
        bootstrapped_relative_positions_table = s.pmparesultsdir + '/.' + s.targetname + '_relative_positions.dat'
        if os.path.exists(bootstrapped_relative_positions_table):
            os.remove(bootstrapped_relative_positions_table)
        saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_five_histograms_plot_parameters.pickle' 
        if os.path.exists(saved_plot_parameters):
            os.remove(saved_plot_parameters)
        pulsitions = s.pmparesultsdir + '/.' + s.targetname + '.pmpar.in.bootstrap'
        pmparinfile = s.pmparesultsdir + '/' + pmparinfile
        if not os.path.exists(pmparinfile):
            print("%s does not exists; aborting\n" % pmparinfile)
            sys.exit()
        #pmparout_bootstrap = pmparesultsdir + '/.' + targetname + '.pmpar.out.bootstrap'
        positions = []
        if priors != '':
            prior_strings = s.parse_priors_for_bootstrap_pmpar(priors, bootstrapruns)
        lines = open(pmparinfile).readlines()
        for line in lines:
            if 'epoch' in line:
                epochline = line
            if line.count(':') == 4 and (not line.strip().startswith('#')):
                positions.append(line)
        nepoch = len(positions)
        PIs = np.array([])
        mu_as = np.array([])
        mu_ds = np.array([])
        RAs = np.array([])
        Decs = np.array([])
        count = 0
        while count < bootstrapruns:
            random_indices = np.random.randint(0, nepoch, nepoch)
            if len(np.unique(random_indices)) < 3: #use at least 3 different positions for astrometric fit
                continue
            fileWrite = open(pulsitions, 'w')
            fileWrite.write(epochline)
            if priors != '':
                for prior_string in prior_strings:
                    fileWrite.write(prior_string[count] + '\n')
            #fileWrite.write("mu_a = -3.82\nmu_d = -16.09\n")
            for i in random_indices:
                fileWrite.write(positions[i])
            fileWrite.close()
            [D, PI, mu_a, mu_d, RA, Dec, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname) 
            PIs = np.append(PIs, PI)
            mu_as = np.append(mu_as, mu_a)
            mu_ds = np.append(mu_ds, mu_d)
            RAs = np.append(RAs, RA)
            Decs = np.append(Decs, Dec)
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 1000 / bootstrapruns)/10) + " \r")
            count += 1
        t = Table([PIs, mu_as, mu_ds, RAs, Decs], names=['PI', 'mu_a', 'mu_d', 'RA', 'Dec']) 
        print t
        bootstrapped_five_parameters_table = s.pmparesultsdir + '/.' + s.targetname + '_five_parameters.dat'
        if not overwrite_table:
            if os.path.exists(bootstrapped_five_parameters_table):
                t0 = Table.read(bootstrapped_five_parameters_table, format='ascii')
                t = vstack([t0, t])
        t.write(bootstrapped_five_parameters_table, format='ascii', overwrite=True)
    
    def calculate_v_t(s, PI, mu_a, mu_d):
        a = estimate_uncertainty(s.targetname)
        v_t = (mu_a**2+mu_d**2)**0.5 / PI * a.A
        return v_t
    def bootstrapped_sample2measured_value(s):
        bootstrapped_five_parameters_table = s.pmparesultsdir + '/.' + s.targetname + '_five_parameters.dat'
        if not os.path.exists(bootstrapped_five_parameters_table):
            print("%s (bootstrap results) does not exist; aborting\n" % bootstrapped_five_parameters_table)
            sys.exit()
        t = Table.read(bootstrapped_five_parameters_table, format='ascii')
        for estimate in ['PI', 'mu_a', 'mu_d', 'RA', 'Dec']:
            exec("%ss = np.array(t['%s'])" % (estimate, estimate))
            #exec("%ss.sort()" % estimate)
        ## relative positions #############
        bootstrapped_relative_positions_table = s.pmparesultsdir + '/.' + s.targetname + '_relative_positions.dat'
        if os.path.exists(bootstrapped_relative_positions_table):
            t_RP = Table.read(bootstrapped_relative_positions_table, format='ascii')
            for estimate in ['rRA', 'rDec']:
                exec("%ss = np.array(t_RP['%s'])" % (estimate, estimate))
        ## plot parameters ################
        saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_five_histograms_plot_parameters.pickle' 
        if os.path.exists(saved_plot_parameters):
            readfile = open(saved_plot_parameters, 'r')
            plot_parameter_dictionary = pickle.load(readfile)
            readfile.close()
            for key in plot_parameter_dictionary:
                exec("%s = plot_parameter_dictionary[key]" % key)
            for estimate in ['PI', 'mu_a', 'mu_d', 'RA', 'Dec']:
                exec("s.most_probable_%s = howfun.sample2most_probable_value(%ss, binno_%s)" % (estimate, estimate, estimate))
            most_probable_v_t = s.calculate_v_t(s.most_probable_PI, s.most_probable_mu_a, s.most_probable_mu_d)
            #[rRAs_MP, rDecs_MP] = diffposition(howfun.deg2dms(RAs), howfun.deg2dms(s.most_probable_RA), howfun.deg2dms(Decs), howfun.deg2dms(s.most_probable_Dec))
        ## write out estimates #########################################################
        bootstrap_estimates_output = s.pmparesultsdir + '/' + s.targetname + '.bootstrap.estimates.out'
        fileWrite = open(bootstrap_estimates_output, 'w')
        fileWrite.write("estimates obtained with %d bootstrap runs:\n" % len(PIs))
        for how_many_sigma in [3, 2, 1]: #only saves s.value_* and so on for the 1sigma
            CL = math.erf(how_many_sigma/2**0.5) 
            for estimate in ['PI', 'mu_a', 'mu_d', 'RA', 'Dec']:
                exec("[s.value_%s, s.error_%s] = howfun.sample2estimate(%ss, CL)" % (estimate, 
                estimate, estimate))
                exec("s.median_%s = howfun.sample2median(%ss)" % (estimate, estimate))
                exec("s.error_%s_symm = howfun.sample2uncertainty(%ss, s.median_%s, CL)" % (
                estimate, estimate, estimate))
                exec("[s.median_low_end_%s, s.median_high_end_%s] = howfun.sample2median_range(%ss, CL)" % (estimate, estimate, estimate)) #low and high end of the central 68% (or else) of the sample
            if os.path.exists(bootstrapped_relative_positions_table):
                for estimate in ['rRA', 'rDec']:
                    exec("s.median_low_end_%s, s.median_high_end_%s = howfun.sample2median_range(%ss, CL)" % (estimate, estimate, estimate))
            ## error in mas relative to most+probable_RA/Dec, get the symmetric error form
            if os.path.exists(saved_plot_parameters):
                [error_RA_mas1, error_Dec_mas1] = diffposition(howfun.deg2dms(s.value_RA+s.error_RA), howfun.deg2dms(s.most_probable_RA),
                                                               howfun.deg2dms(s.value_Dec+s.error_Dec), howfun.deg2dms(s.most_probable_Dec))
                [error_RA_mas2, error_Dec_mas2] = diffposition(howfun.deg2dms(s.value_RA-s.error_RA), howfun.deg2dms(s.most_probable_RA),
                                                               howfun.deg2dms(s.value_Dec-s.error_Dec), howfun.deg2dms(s.most_probable_Dec))
                error_RA_mas = max(abs(error_RA_mas1), abs(error_RA_mas2))
                error_Dec_mas = max(abs(error_Dec_mas1), abs(error_Dec_mas2))
            ## transverse velocity ########################################################
            min_v_t = s.calculate_v_t(s.value_PI+s.error_PI, min(abs(s.value_mu_a-s.error_mu_a), abs(s.value_mu_a+s.error_mu_a)), 
                                                             min(abs(s.value_mu_d-s.error_mu_d), abs(s.value_mu_d+s.error_mu_d)))
            max_v_t = s.calculate_v_t(s.value_PI-s.error_PI, max(abs(s.value_mu_a-s.error_mu_a), abs(s.value_mu_a+s.error_mu_a)), 
                                                             max(abs(s.value_mu_d-s.error_mu_d), abs(s.value_mu_d+s.error_mu_d)))
            median_min_v_t = s.calculate_v_t(s.median_high_end_PI,  min(abs(s.median_low_end_mu_a), abs(s.median_high_end_mu_a)),
                                                                   min(abs(s.median_low_end_mu_d), abs(s.median_high_end_mu_d))) 
            median_max_v_t = s.calculate_v_t(s.median_low_end_PI, max(abs(s.median_low_end_mu_a), abs(s.median_high_end_mu_a)),
                                                                   max(abs(s.median_low_end_mu_d), abs(s.median_high_end_mu_d))) 
            median_v_t = s.calculate_v_t(s.median_PI, s.median_mu_a, s.median_mu_d) #unnecessarily run three times, unsatisfied with this
            fileWrite.write(70*"=" + "\n")
            fileWrite.write("confidencelevel = %f:\n" % CL)
            try:
                fileWrite.write("epoch = %f\n" % s.epoch)
            except AttributeError:
                pass
            fileWrite.write("pi = %f +- %f (mas)\n" % (s.value_PI, s.error_PI))
            fileWrite.write("mu_a = %f +- %f (mas/yr) #mu_a=mu_ra*cos(dec)\n" % (s.value_mu_a, s.error_mu_a))
            fileWrite.write("mu_d = %f +- %f (mas/yr)\n" % (s.value_mu_d, s.error_mu_d))
            fileWrite.write("PI_symm = %f +- %f (mas) # symmetric uncertainty interval around median\n" % (s.median_PI, s.error_PI_symm))
            fileWrite.write("mu_a_symm = %f +- %f (mas/yr)\n" % (s.median_mu_a, s.error_mu_a_symm))
            fileWrite.write("mu_d_symm = %f +- %f (mas/yr)\n" % (s.median_mu_d, s.error_mu_d_symm))
            s.error_Dec_symm *= 3600*1000 #in mas
            s.error_RA_symm *= 3600*1000 * 15 * np.cos(s.median_Dec/180*np.pi) #in mas
            fileWrite.write("RA_symm = %s +- %f (mas)\n" % (howfun.deg2dms(s.median_RA), s.error_RA_symm))
            fileWrite.write("Dec_symm = %s +- %f (mas)\n" % (howfun.deg2dms(s.median_Dec), s.error_Dec_symm))
            if os.path.exists(bootstrapped_relative_positions_table):
                fileWrite.write("median_RA = %s + %f - %f (mas) # central 68 percent around median\n" % (howfun.deg2dms(s.median_RA), s.median_high_end_rRA, abs(s.median_low_end_rRA)))
                fileWrite.write("median_Dec = %s + %f - %f (mas)\n" % (howfun.deg2dms(s.median_Dec), s.median_high_end_rDec, abs(s.median_low_end_rDec)))
            fileWrite.write("median_PI = %f + %f - %f (mas)\n" % (s.median_PI, s.median_high_end_PI-s.median_PI, s.median_PI-s.median_low_end_PI))
            fileWrite.write("median_mu_a = %f + %f - %f (mas)\n" % (s.median_mu_a, s.median_high_end_mu_a-s.median_mu_a, s.median_mu_a-s.median_low_end_mu_a))
            fileWrite.write("median_mu_d = %f + %f - %f (mas)\n" % (s.median_mu_d, s.median_high_end_mu_d-s.median_mu_d, s.median_mu_d-s.median_low_end_mu_d))
            fileWrite.write("median_D = %f + %f - %f (kpc) \n" % (1/s.median_PI, 1/s.median_low_end_PI-1/s.median_PI, 1/s.median_PI-1/s.median_high_end_PI))
            fileWrite.write("median_v_t = %f + %f - %f (km/s)\n" % (median_v_t, median_max_v_t-median_v_t, median_v_t-median_min_v_t))
            if os.path.exists(saved_plot_parameters):
                fileWrite.write("most_probable_PI = %f + %f - %f (mas)\n" % (s.most_probable_PI, s.value_PI+s.error_PI-s.most_probable_PI, s.most_probable_PI-s.value_PI+s.error_PI))
                fileWrite.write("most_probable_mu_a = %f + %f - %f (mas/yr)\n" % (s.most_probable_mu_a, s.value_mu_a+s.error_mu_a-s.most_probable_mu_a, s.most_probable_mu_a-s.value_mu_a+s.error_mu_a))
                fileWrite.write("most_probable_mu_d = %f + %f - %f (mas/yr)\n" % (s.most_probable_mu_d, s.value_mu_d+s.error_mu_d-s.most_probable_mu_d, s.most_probable_mu_d-s.value_mu_d+s.error_mu_d))
                #[value_RA_MP, error_RA_MP] = howfun.sample2estimate(rRAs_MP, CL)
                #[value_Dec_MP, error_Dec_MP] = howfun.sample2estimate(rDecs_MP, CL)
                fileWrite.write("most_probable_RA = %s +- %f (mas)\n" % (howfun.deg2dms(s.most_probable_RA), error_RA_mas))
                fileWrite.write("most_probable_Dec = %s +- %f (mas)\n" % (howfun.deg2dms(s.most_probable_Dec), error_Dec_mas))
                fileWrite.write("most_probable_D = %f + %f - %f (kpc)\n" % (1/s.most_probable_PI, 1/(s.value_PI-s.error_PI)-1/s.most_probable_PI,
                                                                            1/s.most_probable_PI-1/(s.value_PI+s.error_PI)))
                fileWrite.write("most_probable_v_t = %f + %f - %f (km/s)\n" % (most_probable_v_t, max_v_t-most_probable_v_t, most_probable_v_t-min_v_t))
        fileWrite.close()
        os.system("cat %s" % bootstrap_estimates_output)
        # recover the original data, since they are apparently sorted when running the functions 
        #for estimate in ['PI', 'mu_a', 'mu_d', 'RA', 'Dec']:
        #    exec("%ss = np.array(t['%s'])" % (estimate, estimate))
        return PIs, mu_as, mu_ds, RAs, Decs
    
    def plot_three_histograms_with_bootstrap_results(s, use_saved_plot_parameters=False, binno_PI=500, binno_mu_a=300, binno_mu_d=150, xlim_PI=[], xlim_mu_a=[], xlim_mu_d=[], save_plot_parameters=False):
        from matplotlib import pyplot as plt
        import matplotlib.gridspec as gridspec
        plt.rc('text', usetex=True)
        [PIs, mu_as, mu_ds, RAs, Decs] = s.bootstrapped_sample2measured_value()
        #use saved plot parameters (or not)
        saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_three_histograms_plot_parameters.pickle' 
        if use_saved_plot_parameters and os.path.exists(saved_plot_parameters):
            readfile = open(saved_plot_parameters, 'r')
            plot_parameter_dictionary = pickle.load(readfile)
            readfile.close()
            for key in plot_parameter_dictionary:
                exec("%s = plot_parameter_dictionary[key]" % key)
        most_probable_PI = howfun.sample2most_probable_value(PIs, binno_PI)
        most_probable_mu_a = howfun.sample2most_probable_value(mu_as, binno_mu_a)
        most_probable_mu_d = howfun.sample2most_probable_value(mu_ds, binno_mu_d)
        ## start the plot
        fig = plt.figure(figsize=[8,3])
        gs = gridspec.GridSpec(2, 6)
        #subplot1
        ax1 = fig.add_subplot(gs[:2, :2])
        ax1.hist(PIs, binno_PI, density=True, facecolor='g', alpha=0.75)
        ax1.set_xlabel('parallax (mas)')
        if xlim_PI != []:
            ax1.set_xlim(xlim_PI[0], xlim_PI[1])
        ax1.set_ylabel('probability density')
        ax1.axvline(x=most_probable_PI, c='black', linestyle='--', linewidth=0.5)
        ax1.axvline(x=s.value_PI+s.error_PI, c='black', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.value_PI-s.error_PI, c='black', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.median_PI, c='blue', linestyle='--', linewidth=0.5)
        ax1.axvline(x=s.median_low_end_PI, c='blue', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.median_high_end_PI, c='blue', linestyle='-.', linewidth=0.5)
        #subplot2
        ax2 = fig.add_subplot(gs[:2, 2:4])
        ax2.hist(mu_as, binno_mu_a, density=True, facecolor='g', alpha=0.75)
        if xlim_mu_a != []:
            ax2.set_xlim(xlim_mu_a[0], xlim_mu_a[1])
        ax2.set_xlabel(r'$\rm \mu_{\alpha}~(mas~yr^{-1})$')
        ax2.axvline(x=most_probable_mu_a, c='black', linestyle='dashed', linewidth=0.5)
        ax2.axvline(x=s.value_mu_a+s.error_mu_a, c='black', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.value_mu_a-s.error_mu_a, c='black', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.median_mu_a, c='blue', linestyle='dashed', linewidth=0.5)
        ax2.axvline(x=s.median_low_end_mu_a, c='blue', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.median_high_end_mu_a, c='blue', linestyle='-.', linewidth=0.5)
        #subplot3
        ax3 = fig.add_subplot(gs[:2, 4:6])
        ax3.hist(mu_ds, binno_mu_d, density=True, facecolor='g', alpha=0.75)
        if xlim_mu_d != []:
            ax3.set_xlim(xlim_mu_d[0], xlim_mu_d[1])
        ax3.set_xlabel(r'$\rm \mu_{\delta}~(mas~yr^{-1})$')
        ax3.axvline(x=most_probable_mu_d, c='black', linestyle='dashed', linewidth=0.5)
        ax3.axvline(x=s.value_mu_d+s.error_mu_d, c='black', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.value_mu_d-s.error_mu_d, c='black', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.median_mu_d, c='blue', linestyle='dashed', linewidth=0.5)
        ax3.axvline(x=s.median_low_end_mu_d, c='blue', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.median_high_end_mu_d, c='blue', linestyle='-.', linewidth=0.5)
        #whole setup
        gs.tight_layout(fig) #rect=[0, 0.1, 1, 1])
        plt.savefig('%s/three_histograms.pdf' % s.pmparesultsdir)
        print("\nplots have been made.\n")
        #save plot parameters if required
        if save_plot_parameters:
            plot_parameter_dictionary = {}
            for newkey in ['binno_PI', 'binno_mu_a', 'binno_mu_d', 'xlim_PI', 'xlim_mu_a', 'xlim_mu_d']:
                exec("plot_parameter_dictionary[newkey] = %s" % newkey)
            writefile = open(saved_plot_parameters, 'w')
            pickle.dump(plot_parameter_dictionary, writefile)
            writefile.close()

    def plot_five_histograms_with_bootstrap_results(s, use_saved_plot_parameters=False, binno_PI=500, binno_mu_a=300, binno_mu_d=150, binno_RA=300, binno_Dec=300, xlim_PI=[], xlim_mu_a=[], xlim_mu_d=[], xlim_rRA=[], xlim_rDec=[], save_plot_parameters=False):
        from matplotlib import pyplot as plt
        import matplotlib.gridspec as gridspec
        plt.rc('text', usetex=True)
        [PIs, mu_as, mu_ds, RAs, Decs] = s.bootstrapped_sample2measured_value()
        for parameter in ['PI', 'mu_a', 'mu_d']:
            exec("index = %ss > s.value_%s-9*s.error_%s" % (parameter, parameter, parameter))
            exec("%s1s = %ss[index]" % (parameter, parameter))
            exec("index = %s1s < s.value_%s+9*s.error_%s" % (parameter, parameter, parameter))
            exec("%ss = %s1s[index]" % (parameter, parameter))
        saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_five_histograms_plot_parameters.pickle' 
        if use_saved_plot_parameters and os.path.exists(saved_plot_parameters):
            readfile = open(saved_plot_parameters, 'r')
            plot_parameter_dictionary = pickle.load(readfile)
            readfile.close()
            for key in plot_parameter_dictionary:
                exec("%s = plot_parameter_dictionary[key]" % key)
        most_probable_PI = howfun.sample2most_probable_value(PIs, binno_PI)
        most_probable_mu_a = howfun.sample2most_probable_value(mu_as, binno_mu_a)
        most_probable_mu_d = howfun.sample2most_probable_value(mu_ds, binno_mu_d)
        most_probable_RA = howfun.sample2most_probable_value(RAs, binno_RA)
        most_probable_Dec = howfun.sample2most_probable_value(Decs, binno_Dec)
        [most_probable_rRA, most_probable_rDec] = diffposition(howfun.deg2dms(most_probable_RA), howfun.deg2dms(s.median_RA), 
                                                              howfun.deg2dms(most_probable_Dec), howfun.deg2dms(s.median_Dec))
        ##start plotting
        fig = plt.figure()
        gs = gridspec.GridSpec(4, 6)
        #plot PIs
        ax1 = fig.add_subplot(gs[:2, :2])
        ax1.hist(PIs, binno_PI, density=True, facecolor='g', alpha=0.75)
        ax1.set_xlabel('parallax (mas)')
        if xlim_PI != []:
            ax1.set_xlim(xlim_PI[0], xlim_PI[1])
        ax1.set_ylabel('probability density')
        #ax1.axvline(x=most_probable_PI, c='black', linestyle='--', linewidth=0.5)
        ax1.axvline(x=s.value_PI+s.error_PI, c='black', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.value_PI-s.error_PI, c='black', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.median_PI, c='blue', linestyle='--', linewidth=0.5)
        ax1.axvline(x=s.median_low_end_PI, c='blue', linestyle='-.', linewidth=0.5)
        ax1.axvline(x=s.median_high_end_PI, c='blue', linestyle='-.', linewidth=0.5)
        #plot mu_as
        ax2 = fig.add_subplot(gs[:2, 2:4])
        ax2.hist(mu_as, binno_mu_a, density=True, facecolor='g', alpha=0.75)
        if xlim_mu_a != []:
            ax2.set_xlim(xlim_mu_a[0], xlim_mu_a[1])
        ax2.set_xlabel(r'$\rm \mu_{\alpha}~(mas~yr^{-1})$')
        #ax2.axvline(x=most_probable_mu_a, c='black', linestyle='dashed', linewidth=0.5)
        ax2.axvline(x=s.value_mu_a+s.error_mu_a, c='black', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.value_mu_a-s.error_mu_a, c='black', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.median_mu_a, c='blue', linestyle='dashed', linewidth=0.5)
        ax2.axvline(x=s.median_low_end_mu_a, c='blue', linestyle='-.', linewidth=0.5)
        ax2.axvline(x=s.median_high_end_mu_a, c='blue', linestyle='-.', linewidth=0.5)
        #plot mu_ds
        ax3 = fig.add_subplot(gs[:2, 4:6])
        ax3.hist(mu_ds, binno_mu_d, density=True, facecolor='g', alpha=0.75)
        if xlim_mu_d != []:
            ax3.set_xlim(xlim_mu_d[0], xlim_mu_d[1])
        ax3.set_xlabel(r'$\rm \mu_{\delta}~(mas~yr^{-1})$')
        #ax3.axvline(x=most_probable_mu_d, c='black', linestyle='dashed', linewidth=0.5)
        ax3.axvline(x=s.value_mu_d+s.error_mu_d, c='black', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.value_mu_d-s.error_mu_d, c='black', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.median_mu_d, c='blue', linestyle='dashed', linewidth=0.5)
        ax3.axvline(x=s.median_low_end_mu_d, c='blue', linestyle='-.', linewidth=0.5)
        ax3.axvline(x=s.median_high_end_mu_d, c='blue', linestyle='-.', linewidth=0.5) 
        #prepare relative positions for plotting 
        bootstrapped_relative_positions_table = s.pmparesultsdir + '/.' + s.targetname + '_relative_positions.dat'
        if os.path.exists(bootstrapped_relative_positions_table):
            t_RP = Table.read(bootstrapped_relative_positions_table, format='ascii')
            for estimate in ['rRA', 'rDec']:
                exec("%ss = np.array(t_RP['%s'])" % (estimate, estimate))
                exec("%ss.sort()" % estimate)
        else:
            RAs_dms = howfun.deg2dms(RAs)
            Decs_dms = howfun.deg2dms(Decs)
            median_RA_dms = howfun.deg2dms(s.median_RA)
            median_Dec_dms = howfun.deg2dms(s.median_Dec)
            [rRAs, rDecs] = diffposition(RAs_dms, median_RA_dms, Decs_dms, median_Dec_dms)
            t_RP = Table([rRAs, rDecs], names=['rRA', 'rDec'])
            t_RP.write(bootstrapped_relative_positions_table, format='ascii', overwrite=True)
        [value_rRA, error_rRA] = howfun.sample2estimate(rRAs, 1)
        [value_rDec, error_rDec] = howfun.sample2estimate(rDecs, 1)
        #trim rRAs and rDecs to within +-9 sigma to ease plotting
        index = rRAs > value_rRA - 9*error_rRA
        rRA1s = rRAs[index]
        index = rRA1s < value_rRA + 9*error_rRA
        rRAs = rRA1s[index]
        index = rDecs > value_rDec - 9*error_rDec
        rDec1s = rDecs[index]
        index = rDec1s < value_rDec + 9*error_rDec
        rDecs = rDec1s[index]
        #plot rRAs
        ax4 = fig.add_subplot(gs[2:4, 1:3])
        ax4.hist(rRAs, binno_RA, density=True, facecolor='g', alpha=0.75)
        if xlim_rRA != []:
            ax4.set_xlim(xlim_rRA[0], xlim_rRA[1])
        ax4.set_xlabel(r'relative RA.\,(mas)')
        #ax4.axvline(x=most_probable_rRA, c='black', linestyle='dashed', linewidth=0.5)
        ax4.axvline(x=value_rRA+error_rRA, c='black', linestyle='-.', linewidth=0.5)
        ax4.axvline(x=value_rRA-error_rRA, c='black', linestyle='-.', linewidth=0.5)
        ax4.axvline(x=0, c='blue', linestyle='dashed', linewidth=0.5)
        #plot rDecs
        ax5 = fig.add_subplot(gs[2:4, 3:5])
        ax5.hist(rDecs, binno_Dec, density=True, facecolor='g', alpha=0.75)
        if xlim_rDec != []:
            ax5.set_xlim(xlim_rDec[0], xlim_rDec[1])
        ax5.set_xlabel(r'relative Decl.\,(mas)')
        #ax5.axvline(x=most_probable_rDec, c='black', linestyle='dashed', linewidth=0.5)
        ax5.axvline(x=value_rDec+error_rDec, c='black', linestyle='-.', linewidth=0.5)
        ax5.axvline(x=value_rDec-error_rDec, c='black', linestyle='-.', linewidth=0.5)
        ax5.axvline(x=0, c='blue', linestyle='dashed', linewidth=0.5)
        #finalize plotting
        gs.tight_layout(fig) #rect=[0, 0.1, 1, 1])
        plt.savefig('%s/five_histograms.eps' % s.pmparesultsdir)
        print("\nplots have been made.\n")
        #save plot parameters if required
        if save_plot_parameters:
            plot_parameter_dictionary = {}
            for newkey in ['binno_PI', 'binno_mu_a', 'binno_mu_d', 'binno_RA', 'binno_Dec', 'xlim_PI', 'xlim_mu_a', 'xlim_mu_d', 'xlim_rRA', 'xlim_rDec']:
                exec("plot_parameter_dictionary[newkey] = %s" % newkey)
            writefile = open(saved_plot_parameters, 'w')
            pickle.dump(plot_parameter_dictionary, writefile)
            writefile.close()
    
    def covariance_2d_plots_with_chainconsumer(s, plot_bins=130, HowManyParameters=3, HowManySigma=11, plot_extents=[(),(),(),(),()], mark_median_instead_of_most_probable_value=False):
        """
        corner plot for pi/mu_a/mu_d or pi/mu_a/mu_d/RA/Dec, indicating covariance between the parameters.
        due to different binning scheme (chainconsumer use one parameter to change binning, while my previous code use a separate bin_no for each), need to separately determine the binno in other plot functions, in order to align the truth value to the peak of the histograms.
        """
        from numpy.random import normal, multivariate_normal
        from chainconsumer import ChainConsumer
        from matplotlib import pyplot as plt
        plt.rc('text', usetex=True)
        [PIs, mu_as, mu_ds, RAs, Decs] = s.bootstrapped_sample2measured_value()
        bootstrapped_five_parameters_table = s.pmparesultsdir + '/.' + s.targetname + '_five_parameters.dat'
        #print(type(HowManyParameters))
        data = np.array([PIs, mu_as, mu_ds])
        if HowManyParameters > 3:
            bootstrapped_relative_positions_table = s.pmparesultsdir + '/.' + s.targetname + '_relative_positions.dat'
            if not os.path.exists(bootstrapped_relative_positions_table):
                print("%s does not exist; aborting" % bootstrapped_relative_positions_table)
                sys.exit()
            t_RP = Table.read(bootstrapped_relative_positions_table, format='ascii')
            data = np.array([PIs, mu_as, mu_ds, t_RP['rRA'], t_RP['rDec']])
            [s.value_rRA, s.error_rRA] = howfun.sample2estimate(np.array(t_RP['rRA']), 1)
            [s.value_rDec, s.error_rDec] = howfun.sample2estimate(np.array(t_RP['rDec']), 1)
        para_names = ['PI', 'mu_a', 'mu_d', 'rRA', 'rDec']
        labels = [r"$\rm parallax\,(mas)$", r"$\rm \mu_{\alpha}~(mas~yr^{-1})$", r"$\rm \mu_{\delta}~(mas~yr^{-1})$", r'relative RA.\,(mas)', r'relative Decl.\,(mas)']
        #data = data[:HowManyParameters, :]
        #para_names = para_names[:HowManyParameters]
        #labels = labels[:HowManyParameters]
        truths = [s.median_PI, s.median_mu_a, s.median_mu_d, 0, 0]
        if not mark_median_instead_of_most_probable_value:
            saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_five_histograms_plot_parameters.pickle' #the following passage to re-calculate most_probable_values can be omitted because they have been calculated in s.bootstrapped_sample2measured_value
            if not os.path.exists(saved_plot_parameters):
                print('%s does not exist; aborting' % saved_plot_parameters)
                sys.exit()
            readfile = open(saved_plot_parameters, 'r')
            plot_parameter_dictionary = pickle.load(readfile)
            readfile.close()
            for key in plot_parameter_dictionary:
                exec("%s = plot_parameter_dictionary[key]" % key)
            most_probable_PI = howfun.sample2most_probable_value(PIs, binno_PI)
            most_probable_mu_a = howfun.sample2most_probable_value(mu_as, binno_mu_a)
            most_probable_mu_d = howfun.sample2most_probable_value(mu_ds, binno_mu_d)
            most_probable_rRA = howfun.sample2most_probable_value(RAs, binno_RA) - s.median_RA
            most_probable_rDec = howfun.sample2most_probable_value(Decs, binno_Dec) - s.median_Dec
            truths = [most_probable_PI, most_probable_mu_a, most_probable_mu_d, most_probable_rRA, most_probable_rDec]
        #truths = []
        ## filter out the miss-fitted data falling into local max ############################
        for i in range(HowManyParameters):
            exec("maximum = s.value_%s+HowManySigma*s.error_%s" % (para_names[i], para_names[i]))
            exec("minimum = s.value_%s-HowManySigma*s.error_%s" % (para_names[i], para_names[i]))
            index = data[i, :]<maximum
            data = data[:, index]
            index = data[i, :]>minimum
            data = data[:, index]
        #for i in range(HowManyParameters):
        #    binno = math.floor(plot_bins*0.1*(np.size(data,1))**0.5)
        #    truths.append(howfun.sample2most_probable_value(data[i,:], binno))
        ## make corner plot now ##############################################################
        c = ChainConsumer().add_chain(np.transpose(data), parameters=labels[:HowManyParameters])
        #fig = c.plotter.plot(truth=truth)
        #outputfigure = s.pmparesultsdir + '/covariance_2d_plots_' + str(HowManyParameters) + '_parameters.pdf'
        if not mark_median_instead_of_most_probable_value:
            plot_extents_from_saved = [tuple(xlim_PI), tuple(xlim_mu_a), tuple(xlim_mu_d), tuple(xlim_rRA), tuple(xlim_rDec)]
            c.configure(summary=False, colors="#388E3C", kde=False, bins=plot_bins, sigma2d=False)
            fig = c.plotter.plot(parameters=labels[:HowManyParameters], extents = plot_extents_from_saved[:HowManyParameters], figsize='page', truth=truths[:HowManyParameters])
        if plot_extents != [(),(),(),(),()]:
            fig = c.plotter.plot(parameters=labels[:HowManyParameters], extents=plot_extents[:HowManyParameters], figsize='page', truth=truths[:HowManyParameters])
        #fig.set_size_inches(3 + fig.get_size_inches())
        plt.savefig('%s/covariance_2d_plots_%d_parameters.pdf' % (s.pmparesultsdir, HowManyParameters))
        plt.clf()
    
    def bootstrap_and_save_predicted_positions_at_specific_expno(s, pmparinfile, expno, bootstrapruns=10000):
        pulsitions = s.pmparesultsdir + '/.' + s.targetname + '.pmpar.in.bootstrap'
        pmparinfile = s.pmparesultsdir + '/' + pmparinfile
        b = plot_position_scatter(s.targetname, s.exceptions)
        decyear = b.expno2decyear(expno)
        if not os.path.exists(pmparinfile):
            print("%s does not exists; aborting\n" % pmparinfile)
            sys.exit()
        bootstrapped_predicted_position_file = s.pmparesultsdir + '/.' + s.targetname + '.pmpar.predicted.position.bootstrap.at.' + expno
        positions = []
        lines = open(pmparinfile).readlines()
        for line in lines:
            if 'epoch' in line:
                epochline = line
            if line.count(':') == 4 and (not line.strip().startswith('#')):
                positions.append(line)
                if decyear in line:
                    [junk, RA, errorRA, Dec, errorDec] = line.strip().split(' ')
                    errorRA = float(errorRA)*1000 #s to ms
                    errorRA *= 15*math.cos(howfun.dms2deg(Dec)/180*math.pi) #ms to mas
                    errorDec = float(errorDec)*1000 #s to mas
        nepoch = len(positions)
        count = 0
        while count < bootstrapruns:
            random_indices = np.random.randint(0, nepoch, nepoch)
            if len(np.unique(random_indices)) < 3: #use at least 3 different positions for astrometric fit
                continue
            fileWrite = open(pulsitions, 'w')
            fileWrite.write(epochline)
            #fileWrite.write(priors)
            for i in random_indices:
                fileWrite.write(positions[i])
            fileWrite.close()
            os.system('pmpar %s -z %s >> %s' % (pulsitions, decyear, bootstrapped_predicted_position_file)) 
            print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 1000 / bootstrapruns)/10) + " \r")
            count += 1
        return RA, errorRA, Dec, errorDec
    def prepare_for_plot_to_justify_potential_outlier_using_bootstrapping_given_expno(s, pmparinfile, expno, bootstrapruns=10000, overwrite=False):
        #pmparinfile = s.pmparesultsdir + '/' + pmparinfile
        bootstrapped_predicted_position_file = s.pmparesultsdir + '/.' + s.targetname + '.pmpar.predicted.position.bootstrap.at.' + expno
        if bootstrapruns == 0:
            if not os.path.exists(bootstrapped_predicted_position_file):
                print("bootstrapruns==0 and %s does not exist; aborting" % bootstrapped_predicted_position_file)
        else:
            if overwrite:
                os.system("rm %s" % bootstrapped_predicted_position_file)
        [refRA, errorRA, refDec, errorDec] = s.bootstrap_and_save_predicted_positions_at_specific_expno(pmparinfile, expno, bootstrapruns)
        [RAs, Decs] = s.read_bootstrapped_predicted_position_file(bootstrapped_predicted_position_file)
        [diffRAs, diffDecs] = diffposition(RAs, refRA, Decs, refDec)
        return diffRAs, diffDecs, errorRA, errorDec
    def plot_to_justify_potential_outlier_using_bootstrapping_given_expno(s, pmparinfile, expno, bootstrapruns=10000, overwrite=False, binno_diffRA=500, binno_diffDec=500, xlim_diffRA=[], xlim_diffDec=[]):
        from matplotlib import pyplot as plt
        import matplotlib.gridspec as gridspec
        from scipy.stats import norm
        plt.rc('text', usetex=True)
        [diffRAs, diffDecs, errorRA, errorDec] = s.prepare_for_plot_to_justify_potential_outlier_using_bootstrapping_given_expno(pmparinfile, expno, bootstrapruns, overwrite)
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 4)
        #first subplot
        ax1 = fig.add_subplot(gs[:2, :2])
        ax1.hist(diffRAs, binno_diffRA, density=True, facecolor='g', alpha=0.75)
        ax1.set_ylabel('probability density')
        ax1.set_xlabel(r'relative RA.\,(mas)')
        maximum = max(diffRAs)
        minimum = min(diffRAs)
        if xlim_diffRA != []:
            [minimum, maximum] = xlim_diffRA
            ax1.set_xlim(minimum, maximum)
        x1 = np.arange(minimum, maximum, (maximum-minimum)/1000)
        y1 = norm.pdf(x1, 0, errorRA)
        ax1.plot(x1,y1)
        #second subplot
        ax2 = fig.add_subplot(gs[:2, 2:4])
        ax2.hist(diffDecs, binno_diffDec, density=True, facecolor='g', alpha=0.75)
        ax2.set_xlabel(r'relative Decl.\,(mas)')
        maximum = max(diffDecs)
        minimum = min(diffDecs)
        if xlim_diffDec != []:
            [minimum, maximum] = xlim_diffDec
            ax2.set_xlim(minimum, maximum)
        x2 = np.arange(minimum, maximum, (maximum-minimum)/1000)
        y2 = norm.pdf(x2, 0, errorDec)
        ax2.plot(x2,y2)
        gs.tight_layout(fig)
        plt.savefig('%s/test_outlying_at_%s.eps' % (s.pmparesultsdir, expno))
    def read_bootstrapped_predicted_position_file(s, bootstrapped_predicted_position_file):
        RAs = np.array([])
        Decs = np.array([])
        lines = open(bootstrapped_predicted_position_file).readlines()
        for line in lines:
            RA = line.strip().split('  ')[1]
            Dec = line.strip().split('  ')[2]
            RAs = np.append(RAs, RA)
            Decs = np.append(Decs, Dec)
        return RAs, Decs
    def estimate_outlying_level_given_expno(s, pmparinfile, expno, bootstrapruns=10000, overwrite=False):
        import scipy.special as sp       
        [diffRAs, diffDecs, errorRA, errorDec] = s.prepare_for_plot_to_justify_potential_outlier_using_bootstrapping_given_expno(pmparinfile, expno, bootstrapruns, overwrite)
        for i in ['RA', 'Dec']:
            exec('[value_diff%s, error_diff%s] = howfun.sample2estimate(diff%ss, 1)' % (i, i, i))
            exec('error%s_total = (error%s**2 + error_diff%s**2)**0.5' % (i, i, i))
            exec('SNR_diff%s = abs(value_diff%s)/error%s_total' % (i, i, i))
            exec('possibility_of_consistency_%s = 1 - math.erf(SNR_diff%s/2**0.5)' % (i, i))
        possibility_of_consistency = possibility_of_consistency_RA * possibility_of_consistency_Dec
        print(SNR_diffRA, SNR_diffDec)
        print("The %s position is at %f (or %f sigma) confidence outlying." % (expno, 1-possibility_of_consistency, 2**0.5*sp.erfinv(1-possibility_of_consistency)))
        return 1-possibility_of_consistency
    def read_most_probable_psr_position_obtained_with_bootstrap(s, targetname, HowManySigma=1):
        """
        Function:
        Read target position and its uncertainty from bootstrap.estimates.out file.

        Notice for use:
        a) only works for HowManySigma=1, 2, or 3
        """
        if int(HowManySigma) not in [1,2,3]:
            print('the second parameter should be either 1, 2 or 3 (sigma); aborting')
            sys.exit()
        [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
        BSresults = targetdir + '/pmparesults/' + targetname + '.bootstrap.estimates.out'
        lines = open(BSresults).readlines()
        for estimate in ['RA', 'Dec']:
            exec('%ss = np.array([])' % estimate)
            exec('err%ss = np.array([])' % estimate)
            for line in lines:
                if 'most_probable_' in line and estimate in line:
                    phrase = line.split('=')[-1].split('(mas)')[0].strip()
                    exec("%s = phrase.split('+-')[0].strip()" % estimate)
                    exec("err%s = phrase.split('+-')[-1].strip()" % estimate) 
                    exec("%ss = np.append(%ss, %s)" % (estimate, estimate, estimate))
                    exec("err%ss = np.append(err%ss, float(err%s))" % (estimate, estimate, estimate))
            exec("%s = %ss[3-int(HowManySigma)]" % (estimate, estimate))
            exec("err%s = err%ss[3-int(HowManySigma)]" % (estimate, estimate))
        position = np.array([RA, Dec])
        err = np.array([errRA, errDec])
        return position, err
    def abspsrposition_enhanced(s, targetname, dualphscal=False, dualphscalratio=1, 
            prIBC_has_defined_position_in_catalog=False, exceptive_epochs='', HowManySigma=1):
        """
        Notice for use
        --------------
        a) The function is designed to estimate the absolute position for the target anchored to a primary in-beam calibrator, 
        and indirectly to a main phase calibrator. This works for all PSRPI and MSRPI targets.
        b) It can also work for dualphscal setup.
        c) It would also work for the simplest phscal-target setup.

        Measured value
        --------------
        This function adopts the most probable RA/Dec that is kept by the boostrap.estimates.out file.
        Then it will be shifted when the target is tied to the main phscal, and sebsequently aligned to the catalog phscal position.
        
        Uncertainty
        -----------
        a) For phscal-prIBC-target case: This function reports two sets of results (scatter of positions for prIBC or phscal). 
        The uncertainties include the bootstrap uncertainty, prIBC position uncertainty and catalog uncertainty for the main phscal.
        b) For dualphscal setup: Assuming both calibrators have well defined positions, then no posterior position as well as its
        scatter is needed. So the uncertainties come from the uncertainties of absolute positions of each phscal, added in quadrature 
        by the bootstrap uncertainties.
        """
        dualphscalratio = float(dualphscalratio)
        epoch = readpulsition(targetname)[2]
        [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
        [phscalstats, prIBCstats] = target2positionscatter(targetname, exceptive_epochs)
        phscalRADEC = [phscalstats[0][0], phscalstats[1][0]]
        phscalRADECstd = np.array([phscalstats[0][1], phscalstats[1][1]])
        [psrRADEC0, err_psr0] = s.read_most_probable_psr_position_obtained_with_bootstrap(targetname, HowManySigma)
        rfcinfos_phscal = rfcposition(phscalname)
        phscal_absRADEC = [rfcinfos_phscal[0], rfcinfos_phscal[2]]
        errRAphscal = float(rfcinfos_phscal[1])
        errDECphscal = float(rfcinfos_phscal[3])
        err_abs_phscal = np.array([errRAphscal, errDECphscal])
        phscal_refRADEC = srcposition(targetname, phscalname)
        if phscalname != prIBCname:
            prIBC_RADEC = [prIBCstats[0][0], prIBCstats[1][0]]
            prIBC_RADECstd = np.array([prIBCstats[0][1], prIBCstats[1][1]])
            prIBC_refRADEC = srcposition(targetname, prIBCname)
            if prIBC_has_defined_position_in_catalog:
                rfcinfos_prIBC = rfcposition(prIBCname)
                prIBC_absRADEC = [rfcinfos_prIBC[0], rfcinfos_prIBC[2]]
                errRAprIBC = float(rfcinfos_prIBC[1])
                errDECprIBC = float(rfcinfos_prIBC[3])
                err_abs_prIBC = np.array([errRAprIBC, errDECprIBC])
        if phscalname != prIBCname and not dualphscal and not prIBC_has_defined_position_in_catalog:
            ## 1) using prIBC shift ####################
            psrRADEC1 = howfun.dms2deg(psrRADEC0) + howfun.dms2deg(prIBC_RADEC) - howfun.dms2deg(prIBC_refRADEC)
            ## 2) using phscal shift ###################
            psrRADEC2 = howfun.dms2deg(psrRADEC0) + howfun.dms2deg(phscal_refRADEC) - howfun.dms2deg(phscalRADEC)
            ## 3) update absolute position for phscal #####
            psrRADEC2 = psrRADEC2 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
            psrRADEC2 = howfun.deg2dms(psrRADEC2)
            psrRADEC1 = psrRADEC1 + howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC)
            psrRADEC1 = howfun.deg2dms(psrRADEC1)
            ## error estimation #####################################################
            err1 = (err_psr0**2 + err_abs_phscal**2 + prIBC_RADECstd**2)**0.5
            err2 = (err_psr0**2 + err_abs_phscal**2 + phscalRADECstd**2)**0.5
            return psrRADEC1, err1, psrRADEC2, err2, epoch
        elif phscalname != prIBCname and dualphscal and prIBC_has_defined_position_in_catalog:
            ## update absolute position for phscal and prIBC ########################
            psrRADEC = howfun.dms2deg(psrRADEC0)
            psrRADEC += (1-dualphscalratio) * (howfun.dms2deg(phscal_absRADEC) - howfun.dms2deg(phscal_refRADEC))
            psrRADEC += dualphscalratio * (howfun.dms2deg(prIBC_absRADEC) - howfun.dms2deg(prIBC_refRADEC))
            psrRADEC = howfun.deg2dms(psrRADEC)
            ## error estimation #####################################################
            err = (err_psr0**2 + dualphscalratio**2*err_abs_prIBC**2 + (1-dualphscalratio)**2*err_abs_phscal**2)**0.5
            return psrRADEC, err, epoch
        else:
            pass
    def calculate_mu_b_and_mu_l(s):
        """
        estimating the uncertainties is not easy as mu_a and mu_d are strongly associated
        """
        pmparout = s.pmparesultsdir + s.targetname + '.pmpar.out'
        RA, Dec, epoch, pi, mu_a, mu_d, error_pi, error_mu_a, error_mu_d, l, b, rchsq = readpmparout(pmparout)
        a = howfun.equatorial2galactic_coordinates(RA, Dec, mu_a, mu_d, error_mu_a, error_mu_d)
        return a.equatorial_proper_motion2galactic_proper_motion()
    def calculate_local_transverse_velocity_in_the_Galaxy(s, d, l, b, mu_l, mu_b, v_t): #in kpc, deg, deg, mas/yr, mas/yr and km/s
        """
        It is originally designed for J1810 falling into the first quadrant of the Galactic coordinates.
        It might need extra inspection for the validity in other qudrants.
        This function will be improved later.
        Non-circular motion of the Sun is neglected.
        """
        R0 = 8.15 #(+-0.15kpc) Reid et al. 2019
        V_Sun2GC = 247 #(+-4km/s)
        V0 = 236 #(+-7km/s)
        l *= math.pi/180
        b *= math.pi/180
        v_b_obs = v_t * mu_b / (mu_b**2 + mu_l**2)**0.5
        v_l_obs = v_t * mu_l / (mu_b**2 + mu_l**2)**0.5
        v_b_gsr = v_b_obs - V_Sun2GC*math.sin(l)*math.sin(b)
        v_l_gsr = v_l_obs + V_Sun2GC*math.cos(l)
        alpha = math.atan((R0*math.cos(l)-d)/(R0*math.sin(l)))
        v_l_loc = v_l_gsr - V0*math.sin(alpha)
        v_b_loc = v_b_gsr + V0*math.cos(alpha)*math.sin(b)
        return v_b_loc, v_l_loc         
    def calculate_apparent_proper_motion_in_LSR(s, RA, Dec, d, v_l_loc=0, v_b_loc=0): #in h, deg
        """
        calculate the apparent proper motion of LSR (v_b_loc=v_l_loc=0) observed from the solar system

        Input parameters
        ----------------
        RA : str
            in HH:MM:SS.SSS.
        Dec : str
            in dd:mm:ss.ssss.
        d : float
            distance in kpc.

        Return parameters
        -----------------
        mu_a : float
            in mas/yr.
        mu_d : float
            in mas/yr.
        """
        R0 = 8.15 #(+-0.15kpc) Reid et al. 2019
        V_Sun2GC = 247 #(+-4km/s)
        V0 = 236 #(+-7km/s)
        aa = howfun.equatorial2galactic_coordinates(RA, Dec)
        l, b = aa.equatorial_position2galactic_position() #deg, deg
        l *= math.pi/180
        b *= math.pi/180
        alpha = math.atan((R0*math.cos(l)-d)/(R0*math.sin(l)))
        v_l_gsr = V0*math.sin(alpha) + v_l_loc
        v_b_gsr = -V0*math.cos(alpha)*math.sin(b) + v_b_loc
        v_b_obs = v_b_gsr + V_Sun2GC*math.sin(l)*math.sin(b) #in km/s
        v_l_obs = v_l_gsr - V_Sun2GC*math.cos(l) #in km/s
        ab = estimate_uncertainty('J1012+5307') # targetname irrelevant
        mu_l, mu_b = v_l_obs/d/ab.A, v_b_obs/d/ab.A #in mas/yr
        mu_a, mu_d = aa.galactic_proper_motion2equatorial_proper_motion(mu_l, mu_b)
        return mu_a, mu_d

class look_for_indirect_SNR_associations:
    """
    Sometimes a pulsar is not directly associated with an SNR, but its divorced companion is. That we term indirect SNR association.
    This class deals with and studies the SNRs from the Green SNR catalog and pulsars recorded in PSRCAT.
    We are looking for pulsars not directly associated with their immediately nearby SNRs (beta=sep/R_SNR>3). 
    We focus on pulsars with precise proper motion measurements and search within their characteristic ages for SNRs on its reverse trajectory.
    Due to the apparent proper motion of the progenitors of SNRs (~2'/100kyr), we start from a pulsar sample with age_i<100kyr (age_i is defined in PSRCAT).
    """
    SNR_dir = '/fred/oz002/hding/SNR/'
    def __init__(s, age_threshold=1e6, filter_maximum_distance=10, years_forward=False, SNR_catalog='SNR_catalog.txt'): #age in yr, distance in deg
        """
        years_forward is for chance-alignment testing, see the docstring for the 'search...' funtion
        """
        s.SNR_infile, s.age_threshold = s.SNR_dir + SNR_catalog, age_threshold
        s.parse_SNR_catalog_to_table(s.SNR_infile)
        s.gather_pulsar_side_infos_from_PSRCAT(s.age_threshold)
        s.search_for_closest_SNR_for_a_pulsar_trajectory(filter_maximum_distance, years_forward)
        s.recognize_potential_directly_and_indirectly_associated_SNR_pulsar_pairs()
    def parse_SNR_catalog_to_table(s, SNR_infile):
        SNR_table = s.SNR_dir + 'SNR_table'
        if os.path.exists(SNR_table):
            s.t_SNR = Table.read(SNR_table, format='ascii')
        else:
            SNRnames = RAs = Decs = sizes = shapes = np.array([])
            lines = open(SNR_infile).readlines()
            for line in lines:
                line = line.strip()
                if not line[0].isdigit():
                    continue
                infos = line.split('  ')
                SNRname = 'G' + infos[0].strip() + infos[1].strip()
                RA = infos[2].strip().replace(' ',':')
                Dec = infos[3].strip().replace(' ',':')
                SNRnames, RAs, Decs = np.append(SNRnames, SNRname), np.append(RAs, RA), np.append(Decs, Dec)
                size1 = infos[4].replace('?','').strip().split('x')
                size = np.array(size1).astype(np.float)
                size = max(size)
                shape = infos[5].strip()
                sizes, shapes = np.append(sizes, size), np.append(shapes, shape)
            RAs_h, Decs_deg = howfun.dms2deg(RAs), howfun.dms2deg(Decs)
            s.t_SNR = Table([SNRnames, RAs, Decs, RAs_h, Decs_deg, sizes, shapes], \
                names=['SNRname', 'RA', 'Dec', 'RA_h', 'Dec_deg', 'size', 'shape']) 
            s.t_SNR.write(SNR_table, format='ascii', overwrite=True)
    def gather_pulsar_side_infos_from_PSRCAT(s, age_threshold):
        """
        run psrcat command -> write to a file -> parse the file -> produce a table
        """
        psrcat_command = ("psrcat -l 'age_i<%f && age_i>0 pmra!=0 && pmdec!=0' -s age_i -c\
            'jname age_i pmra pmdec raj decj posepoch' -o short_csv" % age_threshold)
        psrcat_output = s.SNR_dir + 'psrcat_output'
        os.system(psrcat_command + '>' + psrcat_output) #inquire and output to csv format
        PSRJnames = age_Is = mu_as = mu_ds = RAs = Decs = refepochs = np.array([])
        lines = open(psrcat_output).readlines()
        for line in lines:
            line = line.strip()
            if not line[0].isdigit():
                continue
            infos = line.split(';')
            junk1, PSRJname, age_I, mu_a, mu_d, RA, Dec, refepoch, junk2 = infos
            PSRJnames, age_Is, mu_as, mu_ds, RAs, Decs, refepochs = np.append(PSRJnames, PSRJname), np.append(age_Is, age_I),\
                np.append(mu_as, mu_a), np.append(mu_ds, mu_d), np.append(RAs, RA), np.append(Decs, Dec), np.append(refepochs, refepoch)
        RAs_h, Decs_deg = howfun.dms2deg(RAs), howfun.dms2deg(Decs)
        age_Is, mu_as, mu_ds, refepochs  = np.array(age_Is).astype(np.float), np.array(mu_as).astype(np.float),\
            np.array(mu_ds).astype(np.float), np.array(refepochs).astype(np.float)
        s.t_psr = Table([PSRJnames, age_Is, mu_as, mu_ds, RAs, Decs, RAs_h, Decs_deg, refepochs],\
            names=['PSRJname', 'age_I', 'mu_a', 'mu_d', 'RA', 'Dec', 'RA_h', 'Dec_deg', 'refepoch'])
        psr_table = s.SNR_dir + 'psr_table'
        s.t_psr.write(psr_table, format='ascii', overwrite=True)
    def search_for_closest_SNR_for_a_pulsar_trajectory1(s, filter_maximum_distance): #in deg
        """
        This function is an Euclidean approximation, and has been deprecated!!!
        """
        PSRJnames = age_Is = age_Ks = SNRs = minDs = SNRsizes = nowDs = np.array([])
        X1s_SNR, Y1s_SNR = s.t_SNR['RA_deg'], s.t_SNR['Dec_deg']
        for i in range(len(s.t_psr)):
            ## narrow SNR with a radius
            x2, y2, mu_a, mu_d, age_I = s.t_psr[i]['RA_deg'], s.t_psr[i]['Dec_deg'], s.t_psr[i]['mu_a'],\
                s.t_psr[i]['mu_d'], s.t_psr[i]['age_I']
            D1s = howfun.distance_from_an_array_of_positions_to_one_specific_position(x2,y2, X1s_SNR, Y1s_SNR)
            t_SNR = s.t_SNR[D1s<filter_maximum_distance]
            if len(t_SNR) == 0:
                continue
            ## solve the closest distance from an SNR center to the line segment defined by two endpoints of the puslar
            Xs_SNR, Ys_SNR = t_SNR['RA_deg'], t_SNR['Dec_deg']
            x1 = x2 - mu_a/1000./3600.*age_I
            y1 = y2 - mu_d/1000./3600.*age_I
            Ds, lamdas = howfun.solve_the_distance_from_a_point_to_a_line_segment(Xs_SNR,Ys_SNR,x1,y1,x2,y2,True)
            if min(Ds) == float('inf'):
                continue
            minD = min(Ds)
            t_SNR1 = t_SNR[Ds==minD]
            lamda = lamdas[Ds==minD][0]
            print(min(Ds)*60, lamda)
            print(t_SNR1)
            x_SNR, y_SNR = t_SNR1['RA_deg'][0], t_SNR1['Dec_deg'][0]
            print(x_SNR, y_SNR, x2, y2)
            nowD = 60*((x_SNR-x2)**2+(y_SNR-y2)**2)**0.5 #in arcmin
            ## print a table summarizing the results
            PSRJnames = np.append(PSRJnames, s.t_psr[i]['PSRJname'])
            age_Is    = np.append(age_Is, s.t_psr[i]['age_I'])
            age_Ks    = np.append(age_Ks, age_I*lamda)
            SNRs      = np.append(SNRs, t_SNR1['SNRname'][0])
            minDs     = np.append(minDs, 60*minD) #in arcmin
            SNRsizes  = np.append(SNRsizes, t_SNR1['size'][0]) #in arcmin
            nowDs     = np.append(nowDs, nowD)
        s.t_paired = Table([PSRJnames, age_Is, age_Ks, SNRs, minDs, SNRsizes, nowDs],\
            names=['PSRJnames', 'age_I', 'age_K', 'SNR', 'minD', 'SNRsize', 'nowD'])
        indirectly_paired_table = s.SNR_dir + 'indirectly_paired_table'
        s.t_paired.write(indirectly_paired_table, format='ascii', overwrite=True)
    def search_for_closest_SNR_for_a_pulsar_trajectory(s, filter_maximum_distance, years_forward=False): #in deg
        """
        compatible with spherical geomemtry, yet might run into some bugs near Decl=+/-90deg
        years_forward=True is for chance alignment testing because indirect assoc made with years_foward=True must be fake.
        """
        years_forward = -2 * int(years_forward) + 1 #False -> 1, True -> -1
        PSRJnames = age_Is = age_Ks = SNRs = minDs = SNRsizes = nowDs = np.array([])
        RAs_SNR, Decs_SNR = s.t_SNR['RA'], s.t_SNR['Dec']
        print(RAs_SNR, Decs_SNR)
        for i in range(len(s.t_psr)):
            ## narrow SNR with a radius
            RA2, Dec2, mu_a, mu_d, age_I = s.t_psr[i]['RA'], s.t_psr[i]['Dec'], s.t_psr[i]['mu_a'],\
                s.t_psr[i]['mu_d'], s.t_psr[i]['age_I']
            Ds = 1./60. * howfun.separations(str(RA2), str(Dec2), RAs_SNR, Decs_SNR) #in deg
            t_SNR = s.t_SNR[Ds<filter_maximum_distance]
            if len(t_SNR) == 0:
                continue
            ## solve the closest distance from an SNR center to the line segment defined by two endpoints of the puslar
            RA1s_SNR, Dec1s_SNR = t_SNR['RA'], t_SNR['Dec']
            years_back = age_I * 70./31 #the ratio of interest particularly for J1810-197
            aa = howfun.spherical_astrometry()
            RA1, Dec1 = aa.calculate_positions_at_another_time_with_intial_position_and_proper_motion(RA2, Dec2, mu_a, mu_d, -years_back*years_forward)
            D12 = howfun.separation(RA1, Dec1, RA2, Dec2) #in arcmin
            position2 = RA2 + ',' + Dec2
            position1 = RA1 + ',' + Dec1
            D1s = lamdas = np.array([])
            for j in range(len(t_SNR)):
                position_SNR = RA1s_SNR[j] + ',' + Dec1s_SNR[j]
                ab = find_virtual_calibrator_position_with_colinear_calibrators(position_SNR, position1, position2)
                RA_v, Dec_v = ab.solve_the_position_of_virtual_phscal()
                D_v2 = howfun.separation(RA_v, Dec_v, RA2, Dec2) #in arcmin
                if howfun.separation(RA_v, Dec_v, RA1, Dec1) > D12:
                    #RA_v = Dec_v = 'inf'
                    D1s = np.append(D1s, float('inf'))
                    lamdas = np.append(lamdas, float('inf'))
                    continue
                elif D_v2 > D12:
                    RA_v, Dec_v = RA1, Dec1
                D_VS = howfun.separation(RA_v, Dec_v, RA1s_SNR[j], Dec1s_SNR[j])
                D1s = np.append(D1s, D_VS) #in arcmin
                D1_v2 = howfun.separation(RA_v, Dec_v, RA2, Dec2) #in arcmin
                lamda = D1_v2/D12
                lamdas = np.append(lamdas, lamda)
            minD = min(D1s)
            if minD == float('inf'):
                continue
            t_SNR1 = t_SNR[D1s==minD]
            lamda = lamdas[D1s==minD]
            print(minD, lamda)
            print(t_SNR1)
            RA_SNR, Dec_SNR = t_SNR1['RA'][0], t_SNR1['Dec'][0]
            nowD = howfun.separation(RA_SNR, Dec_SNR, RA2, Dec2)
            ## print a table summarizing the results
            PSRJnames = np.append(PSRJnames, s.t_psr[i]['PSRJname'])
            age_Is    = np.append(age_Is, s.t_psr[i]['age_I'])
            age_Ks    = np.append(age_Ks, years_back*lamda)
            SNRs      = np.append(SNRs, t_SNR1['SNRname'][0])
            minDs     = np.append(minDs, minD) #in arcmin
            SNRsizes  = np.append(SNRsizes, t_SNR1['size'][0]) #in arcmin
            nowDs     = np.append(nowDs, nowD)
        s.t_paired = Table([PSRJnames, age_Is, age_Ks, SNRs, minDs, SNRsizes, nowDs],\
            names=['PSRJnames', 'age_I', 'age_K', 'SNR', 'minD', 'SNRsize', 'nowD'])
        indirectly_paired_table = s.SNR_dir + 'indirectly_paired_table'
        s.t_paired.write(indirectly_paired_table, format='ascii', overwrite=True)
    def recognize_potential_directly_and_indirectly_associated_SNR_pulsar_pairs(s):
        s.t_paired['direct_assoc'] = s.t_paired['nowD'] < s.t_paired['SNRsize'] #SNRsize is 2 times the radius for shell SNRs
        s.t_paired['indirect_assoc'] = (s.t_paired['nowD']>s.t_paired['SNRsize']) &\
            (s.t_paired['minD']<1./2.*s.t_paired['SNRsize']) #note that we compare minD to half of SNRsize (roughly the radius), stricter than direct assoc.
        s.t_indirect_assoc = s.t_paired[s.t_paired['indirect_assoc']]

class measure_the_angular_broadened_size_of_the_target:
    def __init__(s, targetname, exception_epochs=['']):
        """
        e.g. exception_epochs = ['bd179h0']
        you need to run 1) generate_statsfiles*jmfitfromfile(), then run
                        2) compile_the_deconvolved_*statistics()
        uncomment the '#' in __init__() will do the trick.
        """
        s.targetname = targetname
        [auxdir, configdir, s.targetdir, phscalname, prIBCname] = prepare_path_source(s.targetname)
        s.expnos = s.targetname2expnos(targetname, exception_epochs)
        print(s.expnos)
        #s.generate_statsfiles_including_deconvolved_size_info_with_jmfitfromfile()
        #s.compile_the_deconvolved_size_of_the_target_and_calculate_statistics()
    def targetname2expnos(s, targetname, exception_epochs): 
        [auxdir, configdir, targetdir, phscalname, prIBCname] = prepare_path_source(targetname)
        expnos = []
        vexfiles = glob.glob(r'%s/*/*.vex' % targetdir)
        vexfiles.sort()
        for vexfile in vexfiles:
            expno = vexfile.split('/')[-2].strip()
            if expno not in exception_epochs:
                expnos.append(expno)
        return expnos
    def generate_statsfiles_including_deconvolved_size_info_with_jmfitfromfile(s):
        for expno in s.expnos:
            expdir = s.targetdir + '/' + expno
            os.system('jmfitfromfile.py %s/*difmap.gated.fits.ii.a %s/junk > %s/%s_jmfit_%s.stats' % (expdir, expdir, expdir, expno, s.targetname))
            print('jmfitfromfile.py %s/*difmap.gated.fits.ii.a %s/junk > %s/%s_jmfit_%s.stats' % (expdir, expdir, expdir, expno, s.targetname))
    def compile_the_deconvolved_size_of_the_target_and_calculate_statistics(s):
        """
        return av_major_ax, std_major_ax, av_minor_ax, std_minor_ax, av_size, std_size (in arcsecond)
        """
        s.statsfiles=glob.glob(r'%s/*/*_jmfit_%s.stats' % (s.targetdir, s.targetname)) #find statsfile for each epoch
        # here, the statsfiles were generated separately with "jmfitfromfile.py *ii.a junk > *.stats" 
        s.statsfiles.sort()
        print(s.statsfiles)
        major_axs = np.array([])
        minor_axs = np.array([])
        for statsfile in s.statsfiles:
            lines = open(statsfile).readlines()[-10:]
            for line in lines:
                if 'Major ax' in line:
                    line = line.replace('      ','')
                    major_ax = line.split('    ')[2].strip()
                    major_axs = np.append(major_axs, major_ax)
                if 'Minor ax' in line:
                    line = line.replace('      ','')
                    minor_ax = line.split('    ')[2].strip()
                    minor_axs = np.append(minor_axs, minor_ax)
        major_axs = major_axs.astype(np.float)
        minor_axs = minor_axs.astype(np.float)
        sizes = (major_axs + minor_axs) / 2.
        av_major_ax = np.average(major_axs)
        std_major_ax = np.std(major_axs)
        av_minor_ax = np.average(minor_axs)
        std_minor_ax = np.std(minor_axs)
        av_size = np.average(sizes)
        std_size = np.std(sizes)
        return av_major_ax, std_major_ax, av_minor_ax, std_minor_ax, av_size, std_size
      



        
