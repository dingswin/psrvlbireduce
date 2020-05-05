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
def expno2sources(expno):
    auxdir    = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    expconfigfile = configdir + '/' + expno + '.yaml'
    expconfig = yaml.load(open(expconfigfile))
    targetname = expconfig['targets'][0]
    configfile = configdir + '/' + targetname + '.yaml'
    config = yaml.load(open(configfile))
    prIBCname = config['primaryinbeam']
    try:
        othertargetname = config['othertargetname']
    except KeyError:
        othertargetname = ''
    if othertargetname != '':
        foldername = othertargetname
        phscalname = target2cals(foldername, expno)[0]
    else:
        phscalname = target2cals(targetname, expno)[0]
        foldername = targetname
    return targetname, foldername, phscalname, prIBCname

def prepare_path_source(targetname):
    auxdir    = os.environ['PSRVLBAUXDIR']
    configdir = auxdir + '/configs/'
    expconfigfile = configdir + '/' + targetname + '.yaml'
    targetdir = auxdir + '/processing/' + targetname
    if not os.path.exists(targetdir):
        print("%s doesn't exist; aborting\n" % targetdir)
        return False
        sys.exit()
    expconfig = yaml.load(open(expconfigfile))
    prIBCname = expconfig['primaryinbeam']
    phscalname = target2cals(targetname)[0]
    print phscalname, prIBCname
    return auxdir, configdir, targetdir, phscalname, prIBCname

def target2cals(targetname, expno=''): #get phscal and bandpass cal given targetname
    auxdir = os.environ['PSRVLBAUXDIR']
    targetdir = auxdir + '/processing/' + targetname
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
    auxdir = os.environ['PSRVLBAUXDIR']
    targetdir = auxdir + '/processing/' + targetname
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
    
    if expno != '':
        for sourcefile1 in sourcefiles:
            if expno in sourcefile1:
                sourcefile = sourcefile1        
    return cals

def print_out_calibrator_plan(targetname, expno=''):
    import matplotlib.pyplot as plt
    cals = target2phscals_and_inbeamcals(targetname, expno)
    RAs = np.array([])
    Decs = np.array([])
    for cal in cals:
        [RA, Dec] = srcposition(targetname, cal)
        RAs = np.append(RAs, RA)
        Decs = np.append(Decs, Dec)
    RAs_h = howfun.dms2deg(RAs)
    Decs_deg = howfun.dms2deg(Decs)
    [RA_t, Dec_t] = srcposition(targetname, targetname)
    RA_t_h = howfun.dms2deg(RA_t)
    Dec_t_deg = howfun.dms2deg(Dec_t)

    plt.scatter(RAs_h, Decs_deg)
    plt.plot(RA_t_h, Dec_t_deg, 'rs')
    plt.xlabel('Right Ascension (h)')
    plt.gca().invert_xaxis()
    plt.ylabel('Declination (deg)')
    for i, cal in enumerate(cals):
        plt.annotate(cal, (RAs_h[i], Decs_deg[i]))
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
    Function:
    Consume calibrator positions (relative to another calibrator) and calculate the average position and the position scatter.
    Outlier exclusion:
    3-sigma theshold is used to exclude outliers in an iterative way.
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
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.difmap.jmfit.stats' % targetdir)
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
    print phscalRAs, phscalDecs
    prIBCstatsfiles = targetdir2prIBCstatsfiles(targetdir, exceptions)
    [prIBC_RAs, prIBC_Decs] = nonpulsar_statsfiles2positions(prIBCstatsfiles)
    print prIBC_RAs, prIBC_Decs
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
        expnos = statsfiles2expnos(s.prIBCstatsfiles)
        [RAs, Decs] = s.get_prIBC_preselfcal_positions()
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
                       std_RA_mas*2, std_Dec_mas*2, facecolor="0.91", alpha=0.1, zorder=0)
        ax1.add_patch(rect)
        
        ## rfc2018 position of J1819-2036 for comparison
        [RA_rfc2018, Dec_rfc2018] = ['18:19:36.895534', '-20:36:31.57089']
        [diffRA_rfc2018, diffDec_rfc2018] = diffposition(RA_rfc2018, refRA_prIBC, Dec_rfc2018, refDec_prIBC)
        plt.plot(diffRA_rfc2018, diffDec_rfc2018, 'r*', zorder=1)
        plt.savefig("%s/pmparesults/prIBC_preselfcal_scatter.eps" % s.targetdir)
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
        startime_search_key = 'source=' + targetname            
        ## read obs_time from vex file
        vexlines = open(vexfile).readlines()
        for line in vexlines:
            if 'MJD' in line:
                MJD = int(line.split(':')[-1])
            if startime_search_key in line:
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
    prIBCstatsfiles = glob.glob(r'%s/*/*_preselfcal.difmap.jmfit.stats' % targetdir)
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
    [psrRA, psrDec, epoch, junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8] = readpulsition(targetname)
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
    sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
    sourcefiles.sort()
    sourcefile = sourcefiles[0]
    expno = sourcefile.split('/')[-1].split('.')[0].strip().lower()
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
    def __init__(s, expno, dualphscal='', paraA_rchsq=1e-4):
        [s.targetname, s.foldername, s.phsrefname, s.prIBCname] = expno2sources(expno)
        s.expno = expno
        auxdir    = os.environ['PSRVLBAUXDIR']
        targetdir = auxdir + '/processing/' + s.foldername
        s.expdir = targetdir + '/' + expno
        s.paraA_rchsq = paraA_rchsq
        s.dualphscal = dualphscal
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
    def targetbeam(s): #also get SNprIBCs, using prIBC statsfile,
        prIBCstatsfiles = glob.glob(r'%s/*%s.difmap.jmfit.stokesi.stats' % (s.expdir, s.prIBCname))
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
        lines = open(sumfile).readlines()
        for line in lines:
            if s.targetname in line:
                #if ':' in line.split(s.targetname)[0] and howfun.no_alphabet(line.split(s.targetname)[-1]):
                if ':' in line.split(s.targetname)[0]:
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
    RA2 = howfun.dms2deg(RA) + diffRA_ms/1000/3600
    RA2 = howfun.deg2dms(RA2)
    Dec2 = howfun.dms2deg(Dec) + diffDec_mas/1000/3600
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
    RA2 = howfun.dms2deg(RA) + diffRA_ms/1000/3600
    RA2 = howfun.deg2dms(RA2)
    Dec2 = howfun.dms2deg(Dec) + diffDec_mas/1000/3600
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

    def __init__(s, targetname, timingfityear=''): #s -> self
        print readpulsition(targetname)
        s.targetname = targetname
        [s.RA, s.Dec, s.epoch, s.pi, s.mu_a, s.mu_d, s.error0_pi, s.error0_mu_a, 
        s.error0_mu_d, s.l, s.b, funk1] = readpulsition(targetname)
        if timingfityear != '':
            if type(timingfityear) != str:
                timingfityear = str(timingfityear)
            s.timingfityear = timingfityear
            print '\n' + timingfityear
            [s.epochTm, s.DM, s.Pb, s.estimatesTm, s.errorsTm] = s.readtimingfit(timingfityear)
            [s.RATm, s.DecTm, s.piTm, s.mu_aTm, s.mu_dTm] = s.estimatesTm
            [s.error_RATm, s.error_DecTm, s.error_piTm, s.error_mu_aTm, s.error_mu_dTm] = s.errorsTm
            print s.readtimingfit(timingfityear)
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
            muTm = muTm/1000/3600/180*math.pi #mas/yr -> rad/yr
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
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, catalog='catalog1', tableformat='ascii'):
        s.catalog = s.path + catalog
        s.srcnames = s.read_srcname_from_catalog1()
        s.convert_simbad_coordinates_to_topcat_friendly_file(tableformat='ascii')
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
    def read_srcname_from_catalog1(s):
        srcnames = np.array([])
        lines = open(s.catalog).readlines()
        for line in lines:
            srcname = line.strip().split('  ')[-1]
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

class Simbad_source_to_Gaia_count_in_1deg_radius_to_AGNs_crossmatch_radius:
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, srcname, AGN_number, count_radius=1, howmanysigma=3):
        s.srcname = srcname
        s.AGN_number = AGN_number
        howmanysigma1 = 1
        CL1 = 100*math.erf(howmanysigma1/2**0.5)
        [RA_deg, Dec_deg] = s.Simbad_srcname_to_position(s.srcname)
        print RA_deg, Dec_deg
        [Gaia_count, CpSqAs] = s.position_to_Gaia_count_in_1deg_radius(RA_deg, Dec_deg, count_radius)
        print "%d Gaia sources counted within %f deg radius around %s, equivalent to %f count/arcsecond^2" % (Gaia_count, count_radius, srcname, CpSqAs)
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
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astroquery.gaia import Gaia
        coord = SkyCoord(ra=RA_deg, dec=Dec_deg, unit=(u.degree, u.degree), frame='icrs')
        radius_unit = u.Quantity(radius, u.deg)
        j = Gaia.cone_search_async(coord, radius_unit)
        r = j.get_results()
        Gaia_count = len(r)
        count_per_sq_as = Gaia_count/(radius**2)/(3600**2)
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

class catalog_of_Gaia_counterparts_for_AGNs_to_zero_parallax_point(object):
    path = "/fred/oz002/hding/AQLX-1/PREBursters_catalog/"
    def __init__(s, srcname, magnitude_match_offset=0.02):
        import copy
        s.srcname = srcname
        s.MMO = magnitude_match_offset
        s.T0 = s.read_catalog_of_Gaia_counterparts_for_AGNs(srcname)
        print s.T0
        s.T2 = s.delete_rows_where_the_required_parameter_is_absent(s.T0, 'parallax')
        s.T1 = copy.deepcopy(s.T2)
        
        s.mag_g = s.read_phot_g_mean_mag_of_target(s.srcname)
        s.T3 = s.get_like_magnitude_AGNs(s.T2, s.mag_g, s.MMO)
        print s.T3
        [zero_parallax_point, err_ZPP, std_parallax] = s.calculate_zero_parallax_and_its_sigma(s.T3)
        print zero_parallax_point, err_ZPP, std_parallax
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
    def get_like_magnitude_AGNs(s, catalogtable, mag_g, magnitude_match_offset):
        MMO = magnitude_match_offset
        i = catalogtable['phot_g_mean_mag']<mag_g*(1+MMO)
        T1 = catalogtable[i]
        i = T1['phot_g_mean_mag']>mag_g*(1-MMO)
        T2 = T1[i]
        return T2
    def calculate_zero_parallax_and_its_sigma(s, catalogtable):
        parallaxes = catalogtable['parallax']
        errs_parallax = catalogtable['parallax_error']
        [av_parallax, std_parallax] = howfun.weighted_avg_and_std(parallaxes, errs_parallax)
        [average_parallax, integral_err] = howfun.weightX(parallaxes, errs_parallax)
        return average_parallax, integral_err, std_parallax

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
        T = estimate_uncertainty.T/3600/24/estimate_uncertainty.yr2d # in yr
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
        T = estimate_uncertainty.T/3600/24/estimate_uncertainty.yr2d # in yr
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
        T = estimate_uncertainty.T/3600/24/estimate_uncertainty.yr2d # in yr
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
    in calibrator_search_mode, either srcname is in the form of 'HH:MM:SS.SSS,dd:mm:ss.ss'
    """
    def __init__(s, targetname, phscal1name, phscal2name): ## phscal2 is the "inbeamcal"
        s.targetname = targetname
        s.phscal1name = phscal1name
        s.phscal2name = phscal2name
        print s.targetname, s.phscal1name
        s.calibrator_search_mode = False
        if prepare_path_source(s.targetname) == False:
            s.calibrator_search_mode = True
        s.prepare_positions()
        s.quantify_the_plane1_paramters_defined_by_two_phscals_and_000point()
        s.get_the_plane2_perpendicular_to_plane1_and_pass_target_and_000point()
        s.solve_the_position_of_virtual_phscal()
        s.separation_between_virtual_phscal_and_sources()
    def prepare_positions(s):
        for source in ['target', 'phscal1', 'phscal2']:
            if not s.calibrator_search_mode:
                exec('[s.RA_%s, s.Dec_%s] = srcposition(s.targetname, s.%sname)' % (source, source, source))
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

    Usage instructions:
    The functions are largely organized in the order of running: write_out_preliminary/final --> bootstrap_pmpar --> boostrapped_sample2estimates
    --> make plots

    Input:
    e.g. exceptions=['bh142','bh145a']
    """
    def __init__(s, targetname, exceptions='', dualphscal=False, dualphscalratio=1, epoch=57700):
        s.targetname = targetname
        s.exceptions = exceptions
        s.dualphscal = dualphscal
        s.dualphscalratio = dualphscalratio
        s.epoch = epoch
        [auxdir, s.configdir, s.targetdir, s.phscalname, s.prIBCname] = prepare_path_source(targetname)
        s.pmparesultsdir = s.targetdir + '/pmparesults/'
    def find_statsfiles(s):
        targetdir = s.targetdir
        if not os.path.exists(targetdir):
            print("%s doesn't exist; aborting\n" % targetdir)
            sys.exit()
        s.statsfiles=glob.glob(r'%s/*/*.gated.difmap.jmfit.stokesi.stats' % targetdir) #find statsfile for each epoch
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
    def write_out_preliminary_pmpar_in(s):
        s.find_statsfiles()
        if not os.path.exists(s.pmparesultsdir):
            os.system('mkdir %s' % s.pmparesultsdir)
        pulsitions = s.pmparesultsdir + '/' + s.targetname + '.pmpar.in.preliminary'
        s.RAs = np.array([])
        s.Decs = np.array([])
        s.error0RAs = np.array([])
        s.error0Decs = np.array([])
        s.expnos = np.array([])
        s.decyears = np.array([])
        fileWrite = open(pulsitions, 'w')
        fileWrite.write("### pmpar format\n")
        fileWrite.write("#name = " + s.targetname + "\n")
        fileWrite.write("#ref = " + s.targetname + "\n")
        fileWrite.write("epoch = %f\n" % s.epoch)
        fileWrite.write("#pi = 0\n")
        fileWrite.write("# decimalyear RA +/- Dec +/-\n")
        for statsfile in s.statsfiles:
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

    def write_out_pmparin_incl_sysErr(s, pmparesultsdir, targetname, pulsition_suffix, nepoch, epoch, decyears, expnos, RAs, error0RAs, Decs, error0Decs, dualphscal, paraA_rchsq):
        errorRAs   = np.array([])
        errorDecs  = np.array([])
        pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in.' + pulsition_suffix
        fileWrite = open(pulsitions, 'w')
        fileWrite.write("### pmpar format\n")
        fileWrite.write("#name = " + targetname + "\n")
        fileWrite.write("#ref = " + targetname + "\n")
        fileWrite.write("epoch = %f\n" % epoch)
        fileWrite.write("#pi = 0\n")
        fileWrite.write("# decimalyear RA +/- Dec +/-\n")
        for i in range(nepoch):
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
    def write_out_pmparin_incl_sysErr_two_paraA_rchsq(s, pmparesultsdir, targetname, exceptions, pulsition_suffix, nepoch, epoch, decyears, expnos, RAs, error0RAs, Decs, error0Decs, dualphscal, paraA1_rchsq, paraA2_rchsq):
        errorRAs   = np.array([])
        errorDecs  = np.array([])
        pulsitions = pmparesultsdir + '/' + targetname + '.pmpar.in.' + pulsition_suffix
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
    def write_out_final_pmpar_in(s, paraA_rchsq=3.102e-4, paraA_rchsq_step=1e-3, paraA1_rchsq=3.102e-4):
        s.write_out_preliminary_pmpar_in()
        if not s.dualphscal:
            [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, '', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, 1e-3)
            [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
        if s.dualphscal:
            [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, 'use.A1', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA_rchsq)
            #[pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr_two_paraA_rchsq(s.pmparesultsdir, s.targetname, s.exceptions, 'unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA1_rchsq, paraA_rchsq)
            [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
            while False: #iteratively getting paraA_rchsq is turned off
                if rchsq < 1:
                    print("Reduced chi-square already less than unity without systematics; aborting")
                    sys.exit()
                while rchsq > 1:
                    paraA_rchsq += paraA_rchsq_step
                    [pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr(s.pmparesultsdir, s.targetname, '.unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA_rchsq)
                    #[pulsitions, errorRAs, errorDecs] = s.write_out_pmparin_incl_sysErr_two_paraA_rchsq(s.pmparesultsdir, s.targetname, s.exceptions, 'unity.rchsq', s.nepoch, s.epoch, s.decyears, s.expnos, s.RAs, s.error0RAs, s.Decs, s.error0Decs, s.dualphscal, paraA1_rchsq, paraA_rchsq)
                    [D0, PI0,mu_a0,mu_d0,RA0,Dec0, rchsq] = s.pulsitions2paras(pulsitions, s.pmparesultsdir, s.targetname)
        return D0, PI0, mu_a0, mu_d0, RA0, Dec0, rchsq, paraA_rchsq

    def bootstrap_pmpar(s, pmparinfile, bootstrapruns, priors='', overwrite_table=False):
        from astropy.table import vstack
        pulsitions = s.pmparesultsdir + '/.' + s.targetname + '.pmpar.in.bootstrap'
        pmparinfile = s.pmparesultsdir + '/' + pmparinfile
        if not os.path.exists(pmparinfile):
            print("%s does not exists; aborting\n" % pmparinfile)
            sys.exit()
        #pmparout_bootstrap = pmparesultsdir + '/.' + targetname + '.pmpar.out.bootstrap'
        positions = []
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
            fileWrite.write(priors)
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
        saved_plot_parameters = s.pmparesultsdir + '/.' + s.targetname + '_five_histograms_plot_parameters.pickle' 
        if os.path.exists(saved_plot_parameters):
            readfile = open(saved_plot_parameters, 'r')
            plot_parameter_dictionary = pickle.load(readfile)
            readfile.close()
            for key in plot_parameter_dictionary:
                exec("%s = plot_parameter_dictionary[key]" % key)
            for estimate in ['PI', 'mu_a', 'mu_d']:
                exec("s.most_probable_%s = howfun.sample2most_probable_value(%ss, binno_%s)" % (estimate, estimate, estimate))
            most_probable_v_t = s.calculate_v_t(s.most_probable_PI, s.most_probable_mu_a, s.most_probable_mu_d)
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
                ## transverse velocity ########################################################
            min_v_t = s.calculate_v_t(s.value_PI+s.error_PI, min(abs(s.value_mu_a-s.error_mu_a), abs(s.value_mu_a+s.error_mu_a)), 
                                                             min(abs(s.value_mu_d-s.error_mu_d), abs(s.value_mu_d+s.error_mu_d)))
            max_v_t = s.calculate_v_t(s.value_PI-s.error_PI, max(abs(s.value_mu_a-s.error_mu_a), abs(s.value_mu_a+s.error_mu_a)), 
                                                             max(abs(s.value_mu_d-s.error_mu_d), abs(s.value_mu_d+s.error_mu_d)))
            fileWrite.write(70*"=" + "\n")
            fileWrite.write("confidencelevel = %f:\n" % CL)
            fileWrite.write("epoch = %f\n" % s.epoch)
            fileWrite.write("pi = %f +- %f (mas)\n" % (s.value_PI, s.error_PI))
            fileWrite.write("mu_a = %f +- %f (mas/yr) #mu_a=mu_ra*cos(dec)\n" % (s.value_mu_a, s.error_mu_a))
            fileWrite.write("mu_d = %f +- %f (mas/yr)\n" % (s.value_mu_d, s.error_mu_d))
            fileWrite.write("PI_symm = %f +- %f (mas) # symmetric uncertainty interval around median\n" % (s.median_PI, s.error_PI_symm))
            fileWrite.write("mu_a_symm = %f +- %f (mas/yr)\n" % (s.median_mu_a, s.error_mu_a_symm))
            fileWrite.write("mu_d_symm = %f +- %f (mas/yr)\n" % (s.median_mu_d, s.error_mu_d_symm))
            if os.path.exists(saved_plot_parameters):
                fileWrite.write("most_probable_PI = %f + %f - %f (mas)\n" % (s.most_probable_PI, s.value_PI+s.error_PI-s.most_probable_PI, s.most_probable_PI-s.value_PI+s.error_PI))
                fileWrite.write("most_probable_mu_a = %f + %f - %f (mas/yr)\n" % (s.most_probable_mu_a, s.value_mu_a+s.error_mu_a-s.most_probable_mu_a, s.most_probable_mu_a-s.value_mu_a+s.error_mu_a))
                fileWrite.write("most_probable_mu_d = %f + %f - %f (mas/yr)\n" % (s.most_probable_mu_d, s.value_mu_d+s.error_mu_d-s.most_probable_mu_d, s.most_probable_mu_d-s.value_mu_d+s.error_mu_d))
                fileWrite.write("most_probable_D = %f + %f - %f (kpc)\n" % (1/s.most_probable_PI, 1/(s.value_PI-s.error_PI)-1/s.most_probable_PI,
                                                                            1/s.most_probable_PI-1/(s.value_PI+s.error_PI)))
                fileWrite.write("most_probable_v_t = %f + %f - %f (km/s)\n" % (most_probable_v_t, max_v_t - most_probable_v_t, most_probable_v_t-min_v_t))
        fileWrite.close()
        os.system("cat %s" % bootstrap_estimates_output)
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
        plt.savefig('%s/three_histograms.eps' % s.pmparesultsdir)
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
        most_probable_rRA = howfun.sample2most_probable_value(RAs, binno_RA) - s.median_RA
        most_probable_rDec = howfun.sample2most_probable_value(Decs, binno_Dec) - s.median_Dec
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
        ax1.axvline(x=most_probable_PI, c='black', linestyle='--', linewidth=0.5)
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
        ax2.axvline(x=most_probable_mu_a, c='black', linestyle='dashed', linewidth=0.5)
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
        ax3.axvline(x=most_probable_mu_d, c='black', linestyle='dashed', linewidth=0.5)
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
        #plot rRAs
        ax4 = fig.add_subplot(gs[2:4, 1:3])
        ax4.hist(rRAs, binno_RA, density=True, facecolor='g', alpha=0.75)
        if xlim_rRA != []:
            ax4.set_xlim(xlim_rRA[0], xlim_rRA[1])
        ax4.set_xlabel(r'relative RA.\,(mas)')
        ax4.axvline(x=most_probable_rRA, c='black', linestyle='dashed', linewidth=0.5)
        ax4.axvline(x=value_rRA+error_rRA, c='black', linestyle='-.', linewidth=0.5)
        ax4.axvline(x=value_rRA-error_rRA, c='black', linestyle='-.', linewidth=0.5)
        ax4.axvline(x=0, c='blue', linestyle='dashed', linewidth=0.5)
        #plot rDecs
        ax5 = fig.add_subplot(gs[2:4, 3:5])
        ax5.hist(rDecs, binno_Dec, density=True, facecolor='g', alpha=0.75)
        if xlim_rDec != []:
            ax5.set_xlim(xlim_rDec[0], xlim_rDec[1])
        ax5.set_xlabel(r'relative Decl.\,(mas)')
        ax5.axvline(x=most_probable_rDec, c='black', linestyle='dashed', linewidth=0.5)
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
    
    def covariance_2d_plots_with_chainconsumer(s, HowManyParameters=3, HowManySigma=11, plot_extents=[(),(),(),(),()], plot_bins=1):
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
        c.configure(summary=False, colors="#388E3C", kde=False, bins=plot_bins, sigma2d=False)
        fig = c.plotter.plot(parameters=labels[:HowManyParameters], figsize='page', truth=truths[:HowManyParameters])
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
