####################################################################
## Codes written for general use
## Hao Ding
####################################################################

import math,os,sys,random
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def table_str(string_in_table):
    a = str(string_in_table)
    b = a.split('\n')[-1].strip()
    return b

def dms_str2deg(string): #-- convert dd:mm:ss.ssss to dd.dddd
    a = string.split(':')
    b = []
    for item in a:
        b.append(float(item))
    if b[0] < 0:
        sign = -1
    else:
        sign = 1
    b[0] = abs(b[0])
    i=len(b)
    degree = b[i-1]
    while i != 1:
        degree /= 60
        i -= 1
        degree += b[i-1]
    degree *= sign
    #degree = ((b[2]/60+b[1])/60+b[0])*sign
    return degree
def dms2deg(array):
    if type(array)==str:
        degrees = dms_str2deg(array)
    else:
        degrees = np.array([])
        for string in array:
            degree = dms_str2deg(string)
            degrees = np.append(degrees, degree)
    return degrees

def separation(RA1,Dec1,RA2,Dec2): #-- calculate angular separation (in min) given RAs/Decs in dd:mm:ss.ssss format
    RA0 = np.array([dms2deg(str(RA1)), dms2deg(str(RA2))])
    Dec = np.array([dms2deg(str(Dec1)), dms2deg(str(Dec2))])
    Dec_rad = Dec*math.pi/180
    #RA = RA0*np.cos(Dec_rad)*15 # hour to deg
    RA = RA0*15 #hour to deg
    RA_rad = RA*math.pi/180
    cos_sep = math.cos(Dec_rad[0])*math.cos(Dec_rad[1])*math.cos(RA_rad[0]-RA_rad[1]) + math.sin(Dec_rad[0])*math.sin(Dec_rad[1])
    #diffRA = abs(RA[0]-RA[1])*60 #deg to min
    #diffDec = abs(Dec[0]-Dec[1])*60 
    sep = math.acos(cos_sep)
    sep = sep*180/math.pi*60 # rad to arcmin
    #sep = math.sqrt(diffRA**2+diffDec**2)
    return sep
def separations(RA, Dec, RAs, Decs):
    """
    like the 'separation' function, but working for array as well.
    """
    if type(RAs) == str:
        sep = separation(str(RA), str(Dec), RAs, Decs)
        return sep
    else:
        seps = np.array([])
        for i in range(len(RAs)):
            sep = separation(str(RA), str(Dec), str(RAs[i]), str(Decs[i]))
            seps = np.append(seps, sep)
        return seps

def distance_from_an_array_of_positions_to_one_specific_position(X,Y,xs,ys):
    if type(xs)==int or type(xs)==float:
        D = ((X-xs)**2+(Y-ys)**2)**0.5
        return D
    else:
        Ds = np.array([])
        for i in range(len(xs)):
            D = ((X-xs[i])**2+(Y-ys[i])**2)**0.5
            Ds = np.append(Ds, D)
        return Ds
def solve_the_distance_from_a_point_to_a_line_segment(xs,ys,x1,y1,x2,y2, SNR_assoc_mode=False):
    """
    work for an array of xs and ys;
    see solve_the_distance_from_a_point_to_a_line_segment1 for more explanation
    Note that in SNR_assoc_mode, x2, y2 is the current position of pulsar
    """
    if type(xs)==int or type(xs)==float:
        D = solve_the_distance_from_a_point_to_a_line_segment1(xs,ys,x1,y1,x2,y2, SNR_assoc_mode)
        return D, lamda
    else:
        Ds = lamdas = np.array([])
        for i in range(len(xs)):
            D, lamda = solve_the_distance_from_a_point_to_a_line_segment1(xs[i],ys[i],x1,y1,x2,y2, SNR_assoc_mode)
            Ds, lamdas = np.append(Ds,D), np.append(lamdas, lamda)
        return Ds, lamdas
def solve_the_distance_from_a_point_to_a_line_segment1(x,y,x1,y1,x2,y2, SNR_assoc_mode):
    """
    (x1,y1) and (x2,y2) are the endpoints of a continuous line segment; 
    (x2,y2) for current position, (x1,y1) stands for the earlier end;
    the function is made to solve the distance from (x,y) to the line segment;
    output: distance D (from X,Y to x,y), lamda - tracing the proportional position of X,Y
    """
    if (x2-x1)*(x-x1)+(y2-y1)*(y-y1) < 0: 
        X, Y = x1, y1
    elif (x1-x2)*(x-x2)+(y1-y2)*(y-y2) < 0:
        X, Y = x2, y2
        if SNR_assoc_mode:
            X = Y = float('inf')
    else:
        A = np.mat([[y2-y1, x1-x2],
                    [x2-x1, y2-y1]])
        B = np.mat([[x1*y2-y1*x2],
                    [x*(x2-x1)+y*(y2-y1)]])
        solution = A.I * B
        X, Y = solution[0,0], solution[1,0]
    lamda = (X-x2)/(x1-x2)
    D = ((X-x)**2+(Y-y)**2)**0.5
    return D, lamda

class spherical_astrometry:
    """
    This class does astrometric calculations on a spherical 2-D surface. 
    Some functions that belong here but have been developed won't be included.
    """
    def __init__(s):
        pass
    def calculate_positions_at_another_time_with_intial_position_and_proper_motion(s, RA, Dec, mu_a, mu_d, T):
        """
        RA in 'HH:MM:SS.SS', Dec in 'dd:mm:ss.s', mu_a and mu_d in mas/yr, t in yr and in the positive direction.
        Caveat: the bevahior around the two poles (Dec=+/-90deg) is not well modeled.
        """
        a0, d0 = dms2deg(str(RA)), dms2deg(str(Dec)) #a0 in h, d0 in deg,
        mu_a, mu_d = mu_a/1000./3600., mu_d/1000./3600. #both in deg/yr
        d = mu_d*T + d0
        d_rad, d0_rad = d*math.pi/180, d0*math.pi/180
        a = np.log(abs(1./np.cos(d_rad) + np.tan(d_rad)))
        a -= np.log(abs(1./np.cos(d0_rad) + np.tan(d0_rad)))
        a *= mu_a/mu_d/15.*180/math.pi
        a += a0
        return deg2dms(a), deg2dms(d) #a in h-str, d in deg-str; (a,d) is the position extrapolated at time t (t can be negative, which is the past)
    def plot_the_trajectory_of_a_pulsar_with_a_constant_proper_motion_on_an_Euclidean_canvas(s, RA, Dec, mu_a, mu_d, T):
        """
        RA in 'HH:MM:SS.SS', Dec in 'dd:mm:ss.s', mu_a and mu_d in mas/yr, t in yr and in the positive direction.
        the plot is made by sampling time from 0 to t (t can be negative, which means the past)
        """
        if T < 0:
            ts = np.arange(T, 0, -T/1000)
        else:
            ts = np.arange(0, T, T/1000)
        RA1s = Dec1s = np.array([])
        for t in ts:
            RA1, Dec1 = s.calculate_positions_at_another_time_with_intial_position_and_proper_motion(RA, Dec, mu_a, mu_d, t)
            RA1, Dec1 = dms2deg(RA1), dms2deg(Dec1)
            RA1s, Dec1s = np.append(RA1s, RA1), np.append(Dec1s, Dec1)
        plt.plot(RA1s, Dec1s)
        plt.xlabel('RA. (h)')
        plt.ylabel('Decl. (deg)')
        plt.gca().invert_xaxis()
        plt.show()

def colonizedms(string):
    import re
    string = string.strip()
    sign = ''
    if string.startswith('-') or string.startswith('+'):
        sign = string[0]
        string = string[1:]
    #    sign = -1
    #else:
    #    sign = 1
    #string = filter(str.isdigit,string)
    newstr = re.sub('\s+', ':', string.strip())
    #newstr = string[:2] + ':' + string[2:4] + ':' + string[4:6] + '.' + string[6:]
    #if sign == -1:
    #    newstr = "-" + newstr
    newstr = sign + newstr
    return newstr

def deg2dms(array):
    if type(array)==float or type(array)==np.float64:
        dmss = deg_float2dms(array)
    else:
        dmss = np.array([])
        for number in array:
            dms = deg_float2dms(number)
            dmss = np.append(dmss, dms)
    return dmss
def deg_float2dms(degree): #-- degree to dd:mm:ss.sssssss
    sign = np.sign(degree)
    degree = float(abs(degree)) 
    d = math.floor(degree)
    m = math.floor((degree-d)*60)
    s = ((degree-d)*60-m)*60
    dms="%02d:%02d:%010.7f" % (d,m,s)
    if sign == -1:
        dms = '-' + dms
    return dms 

## get extension in x and y direction with analytic solutions
def ellipse2xy(LA,SA,PA): #-- LA ->full long axis; SA ->full short axis; PA ->position angle
    b = 0.5*LA
    a = 0.5*SA
    d = PA*math.pi/180
    # for y_max:
    lamda = (a**2-b**2)*math.sin(d)*math.cos(d)/((b*math.cos(d))**2+(a*math.sin(d))**2) #lamda=x/y
    y_max = ((lamda*math.cos(d)+math.sin(d))/a)**2 + ((math.cos(d)-lamda*math.sin(d))/b)**2
    y_max = y_max**(-0.5)
    # for x_max:
    mu = (a**2-b**2)*math.sin(d)*math.cos(d)/((a*math.cos(d))**2+(b*math.sin(d))**2)
    x_max = ((math.cos(d)+mu*math.sin(d))/a)**2 + ((mu*math.cos(d)-math.sin(d))/b)**2
    x_max = x_max**(-0.5)
    return x_max,y_max
def deprojectbeam2xy(LA,SA,PA):
    [x_max,y_max] = ellipse2xy(LA,SA,PA)
    xyratio = x_max/y_max
    x = (math.pi*xyratio*LA/2*SA/2)**0.5
    y = (math.pi*LA/2*SA/2/xyratio)**0.5
    return x,y

def mas2ms(ErrArray_mas,Dec):
    d = dms2deg(Dec)*math.pi/180
    ErrArray_ms = ErrArray_mas/15/math.cos(d)
    return ErrArray_ms
def ms2mas(ErrArray_ms, Dec):
    d = dms2deg(Dec)*math.pi/180
    ErrArray_mas = ErrArray_ms*15*math.cos(d)
    return ErrArray_mas

def sample2estimate(array1,confidencelevel):
    array = sorted(array1)
    CL = confidencelevel
    if CL<1:
        SV = int(CL*len(array)) #SV -> significant volume
    elif CL>=1 and CL<10:
        CL = math.erf(CL/2**0.5)
        SV = int(CL*len(array))
    elif CL>=10:
        SV = CL    
    delta = float('inf')
    for i in range(len(array)-SV-1):
        diff = array[i+SV] - array[i]
        if diff < delta:
            j=i
            delta = diff
    confidence_min = array[j]
    confidence_max = array[j+SV]
    value = 0.5*(confidence_min+confidence_max)
    error = 0.5*(confidence_max-confidence_min)
    return value,error

def sample2most_probable_value(array1, bins=1000):
    array = sorted(array1)
    bins = int(bins)
    [counts, values] = np.histogram(array, bins)
    index_max_count = np.argmax(counts)
    most_probable_value = 0.5*(values[index_max_count] + values[index_max_count+1])
    return most_probable_value

def sample2median(array1):
    array = sorted(array1)
    length = len(array)
    if length % 2 == 0:
        median = 0.5*(array[length/2-1] + array[length/2])
    else:
        median = array[(length-1)/2]
    return median
def sample2median_range(array1, confidencelevel):
    array = sorted(array1)
    CL = confidencelevel
    if CL<1:
        SV = int(CL*len(array)) #SV -> significant volume
    elif CL>=1 and CL<10:
        CL = math.erf(CL/2**0.5)
        SV = int(CL*len(array))
    elif CL>=10:
        SV = CL
    index_start = (len(array)-SV)/2-1
    index_end = index_start + SV
    return array[index_start], array[index_end]

def uncertainties_from_2Dsample(table, howmanysigma=1, step=1): #ascii-format table where only FLOAT allowed!!
# this only gives the most compact pi, or t[:,0], not globally; if we want to have global minimum, we have to look at MCMC
    CL = math.erf(howmanysigma/2**0.5) 
    if type(table) == Table:
        t = astropytable2nparray(table)
    else: #if not, it has to be an numpy array, or would go wrong
        t = table
    t = t[np.argsort(t[:,0])] 
    L = len(t)
    n_para = np.size(t,1)
    SV = int(CL*L)
    err_volume = float('inf')
    for i in range(0,L-SV-1,step):
        err = t[i+SV,0] - t[i,0] #err for 0th column
        for j in range(1,np.size(t,1)): # non-0th column
            maxj = max(t[i:i+SV+1,j])
            minj = min(t[i:i+SV+1,j])
            errj = maxj - minj
            err *= errj # err volume
        if err < err_volume:
            err_volume = err
            k = i
        print("\x1B[1A\x1B[Kprogress:{0}%".format(round((i + 1) * 100 / (L-SV-1))) + " \r")
    starts = np.array([t[k,0]])
    ends = np.array([t[k+SV,0]])
    for j in range(1,np.size(t,1)):
        starts = np.append(starts, min(t[k:k+SV+1,j]))
        ends = np.append(ends, max(t[k:k+SV+1,j]))
    err_ranges = np.row_stack((starts, ends))
    return err_ranges
    
def astropytable2nparray(table): #only for 2d table, not workable for 1d data!
    t1 = Table.read(table, format='ascii')
    for item in t1.colnames:
        t1_row = np.array(t1[item])
        try:
            np2darray = np.row_stack((np2darray, t1_row))
        except NameError:
            np2darray = t1_row
    t3 = np2darray.transpose()
    t = t3[np.argsort(t3[:,0])] 
    #a=a[np.argsort(a[:,0])]
    return t

def nD_sample2estimates(howmanysigma=1, *arrays):
    #instructions below
    CL = math.erf(howmanysigma/2**0.5)
    N = len(arrays)
    if N == 0:
        print "There is no sample; aborting..."
        sys.exit()
    if N == 1:
        return sample2estimate(arrays[0], CL)
        sys.exit()

    margin = 0.01
    start_ratios = np.array([])
    end_ratios = np.array([])
    l = []
    for array1 in arrays:
        array = sorted(array1)
        start_ratios = np.append(start_ratios, CL)
        end_ratios = np.append(end_ratios, 1-margin)
    #start_ratios = np.array([0.83,0.90,0.82])
    #end_ratios = np.array([0.86,0.93,0.85])

    #use langsam mode (400,40,3), then use prior start/end_ratios and switch to fast search (5,40,30)
    iterations = 30
    j = 0
    while j < iterations: 
        steps = 40
        count = 0
        [junk1, ratios_stack, junk2, junk3] = nD_ratios(arrays, start_ratios, end_ratios, CL)
        while count < steps-1:
            [junk1, ratios, junk2, junk3] = nD_ratios(arrays, start_ratios, end_ratios, CL)
            ratios_stack = np.row_stack((ratios_stack, ratios))
            count += 1
        start_ratios = np.array([])
        end_ratios = np.array([])
        for i in range(N):
            start_ratios = np.append(start_ratios, min(ratios_stack[:,i]))
            end_ratios = np.append(end_ratios, max(ratios_stack[:,i]))
        print start_ratios, end_ratios
        print start_ratios, end_ratios
        j += 1
    return nD_ratios(arrays, start_ratios, end_ratios, CL, 400)

def nD_ratios(arrays, start_ratios, end_ratios, CL, trials=5):
    count = 0
    N = len(arrays)
    V_compare = float('inf')
    while count < trials:
        Ve = 1
        values = np.array([])
        errors = np.array([])
        ratios = np.array([])
        p = CL
        for i in range(N):
            ratio = random.uniform(start_ratios[i], end_ratios[i])
            if i != N-1:
                [value, error] = sample2estimate(arrays[i], ratio)
                ratios = np.append(ratios, ratio)
                p /= ratio
                start_ratios[i+1] = max(p, start_ratios[i+1])
            else:
                if p > 0.99:
                    p = 0.99
                [value, error] = sample2estimate(arrays[i], p)
                ratios = np.append(ratios, p)
            Ve *= error
            values = np.append(values, value)
            errors = np.append(errors, error)
        if Ve < V_compare:
            V_compare = Ve
            ratios_optimal = ratios
            values_optimal = values
            errors_optimal = errors
        print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 100 / trials)) + " \r")
        count += 1
    return V_compare, ratios_optimal, values_optimal, errors_optimal

def sample2uncertainty(array1,estimate,confidencelevel): #offered with an estimate and present symmetric format
    array = sorted(array1)
    CL = confidencelevel
    SV = int(CL*len(array))
    delta = float('inf')
    for i in range(len(array)-SV-1):
        diff = array[i+SV] + array[i] - 2*estimate
        diff = abs(diff)
        if diff < delta:
            j=i
            delta = diff
    uncertainty = estimate - array[j]
    return uncertainty

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass
            
def pmparout2paras(infile):                        
    for parameter in ['mu_a', 'mu_d', 'PI', 'RA', 'Dec', 'rchsq']:
        exec('%ss = np.array([])' % parameter)
    lines = open(infile).readlines()
    for line in lines:
        if '+-' in line:
            line = line.split('+-')[0].strip()
            line1 = line.split('=')[-1].strip()
            if 'RA' in line:
                RAs = np.append(RAs, dms2deg(line1))
            if 'Dec' in line:
                Decs = np.append(Decs, dms2deg(line1))
            if 'mu_a' in line:
                mu_as = np.append(mu_as, float(line1))
            if 'mu_d' in line:
                mu_ds = np.append(mu_ds, float(line1))
            if 'pi' in line:
                PIs = np.append(PIs, float(line1))
        if 'Reduced' in line:
            rchsq = line.split('=')[-1].strip()
            rchsqs = np.append(rchsqs, float(rchsq))
    for parameter in ['mu_a', 'mu_d', 'PI', 'RA', 'Dec', 'rchsq']:
        exec('length = len(%ss)' % parameter)
        if length == 1:
            exec('%ss = %ss[0]' % (parameter, parameter))
    return PIs,mu_as,mu_ds,RAs,Decs, rchsqs

## written by Adam: run pmpar and read pmpar_e/t
def pmpar(pmparexec, filename, nopmsubtract=False):
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

def weightX(Xs,errs_X): #both numpy array
    errs_X = 1./errs_X**2
    sum_err = sum(errs_X)
    Xbar = sum(errs_X/sum_err*Xs)
    err = 1./sum_err**0.5
    return Xbar, err

def weighted_avg_and_std(Xs, errs_X):
    errs_X = 1./errs_X**2
    sum_err = sum(errs_X)
    w = errs_X/sum_err 
    Xbar = np.average(Xs, weights=w)
    variance = np.average((Xs-Xbar)**2, weights=w)
    return (Xbar, math.sqrt(variance))

def is_pure_number_or_space(str):
    is_pure_number_or_space = True
    for letter in str:
        if not letter.isdigit() and letter != ' ':
            is_pure_number_or_space = False
            break
    return is_pure_number_or_space

def no_alphabet(str):
    no_alphabet = True
    for letter in str:
        if letter.isalpha():
            no_alphabet = False
            break
    return no_alphabet

class lsqfit: #y=ax+b
    def __init__(s):
        #Xs = np.array([1,2,4,7,10])
        #errXs = np.array([0.1,0.3,0.3,0.5,0.2])
        #Ys = np.array([1,3,6,2,4])
        #errYs = errXs
        #print s.linearfit2D(Xs, errXs, Ys, errYs)
        #[a, b, r_chi_sq, run] = s.linearfit2D(Xs, errXs, Ys, errYs)
        #s.plot_linearfit2D(a, b, Xs, Ys, errXs, errYs)
        pass
    def abcfit3D(s, Xs, errXs, Ys, errYs, Zs, errZs): #Z=a*(c*X-1)+b*c**2*Y
        DoF = len(Xs) - 3
        X1s = Xs
        Y1s = Ys
        chi_sq_old = float('inf') 
        chi_sq = 1e10
        run = 0
        while abs(chi_sq/chi_sq_old)<1-1e-5:
            run += 1
            chi_sq_old = chi_sq
            # calculate a,b,c based on X1s, Y1s and Zs (and errX/Y/Zs)
            A1 = np.mat([[sum(X1s/errZs**2),     sum(-1./errZs**2),  sum(Y1s/errZs**2)], 
                         [sum(X1s**2/errZs**2),  sum(-X1s/errZs**2), sum(Y1s*X1s/errZs**2)], 
                         [sum(X1s*Y1s/errZs**2), sum(-Y1s/errZs**2), sum(Y1s**2/errZs**2)]])
            B1 = np.mat([[sum(Zs/errZs**2)], [sum(Zs*X1s/errZs**2)], [sum(Zs*Y1s/errZs**2)]])
            abc_solutions = A1.I*B1
            ac = abc_solutions[0,0]
            a = abc_solutions[1,0]
            bc2 = abc_solutions[2,0]
            c = ac/a
            b = bc2/c**2
            # calculate chi_sq
            chi_sq_i = (X1s-Xs)**2/errXs**2 + (Y1s-Ys)**2/errYs**2
            chi_sq_i += (Zs-a*(c*X1s-1)-b*c**2*Y1s)**2/errZs**2
            chi_sq = sum(chi_sq_i)
            # calculate X1s and Y1s, based on a,b,c,Xs,Ys
            A2 = np.array([[1./errXs**2 + a**2*c**2/errZs**2, a*b*c**3/errZs**2],
                         [a*b*c**3/errZs**2, b**2*c**4/errZs**2 + 1./errYs**2]])
            B2 = np.array([[Xs/errXs**2 + a*c*(Zs+a)/errZs**2], [Ys/errYs**2 + b*c**2*(Zs+a)/errZs**2]])
            X1sY1s_solutions = np.zeros((2,1,len(Xs)))
            for i in range(len(Xs)):
                X1sY1s_solutions[:,:,i] = np.linalg.solve(A2[:,:,i], B2[:,:,i])
            X1s = X1sY1s_solutions[0,0,:]
            Y1s = X1sY1s_solutions[1,0,:]
        r_chi_sq = chi_sq/DoF
        return a, b, c, r_chi_sq, run
    def abcfit3D1(s, Xs, errXs, Ys, errYs, Zs, errZs): #Z=a*(c*X-1)+b*c**2*Y; this function is replaced by abcfit3D
        DoF = len(Xs) - 3
        X1s = Xs
        Y1s = Ys
        chi_sq_old = float('inf') 
        chi_sq = 1e10
        run = 0
        while abs(chi_sq/chi_sq_old)<1-1e-5:
            run += 1
            chi_sq_old = chi_sq
            # calculate a,b,c based on X1s, Y1s and Zs (and errX/Y/Zs)
            C1 = sum(X1s**2/errZs**2)*sum(Y1s/errZs**2) - sum(X1s/errZs**2)*sum(X1s*Y1s/errZs**2)
            C2 = sum(X1s/errZs**2)*sum(Y1s/errZs**2) - sum(1./errZs**2)*sum(X1s*Y1s/errZs**2)
            C3 = sum(Zs*X1s/errZs**2)*sum(Y1s/errZs**2) - sum(Zs/errZs**2)*sum(X1s*Y1s/errZs**2) 
            C4 = sum(X1s*Y1s/errZs**2)*sum(Y1s/errZs**2) - sum(X1s/errZs**2)*sum(Y1s**2/errZs**2)
            C5 = sum(Y1s/errZs**2)*sum(Y1s/errZs**2) - sum(1./errZs**2)*sum(Y1s**2/errZs**2)
            C6 = sum(Zs*Y1s/errZs**2)*sum(Y1s/errZs**2) - sum(Zs/errZs**2)*sum(Y1s**2/errZs**2)
            C7 = C5*C1 - C2*C4
            C8 = C6*C1 - C3*C4
            a = -C8/C7
            c = (a*C2+C3)/a/C1
            b = sum(Zs/errZs**2) - a*c*sum(X1s/errZs**2) + a*sum(1./errZs**2)
            b /= c**2*sum(Y1s/errZs**2)
            # calculate chi_sq
            chi_sq_i = (X1s-Xs)**2/errXs**2 + (Y1s-Ys)**2/errYs**2
            chi_sq_i += (Zs-a*(c*X1s-1)-b*c**2*Y1s)**2/errZs**2
            chi_sq = sum(chi_sq_i)
            # calculate X1s and Y1s, based on a,b,c,Xs,Ys
            D1s = 1./errXs**2 + a**2*c**2/errZs**2
            D2s = a*b*c**3/errZs**2
            D3s = Xs/errXs**2 + a*c*(Zs+a)/errZs**2
            D4s = D2s
            D5s = b**2*c**4/errZs**2 + 1./errYs**2
            D6s = Ys/errYs**2 + b*c**2*(Zs+a)/errZs**2
            Deltas = D1s*D5s - D2s*D4s
            X1s = D3s*D5s - D2s*D6s
            X1s /= Deltas
            Y1s = D1s*D6s - D3s*D4s
            Y1s /= Deltas
        r_chi_sq = chi_sq/DoF
        return a, b, c, r_chi_sq, run

    def linearfit2D(s, Xs, errXs, Ys, errYs):
        DoF = len(Xs) - 2 #degree of freedom
        X1s = Xs
        chi_sq_old = float('inf') 
        chi_sq = 1e10
        run = 0
        while abs(chi_sq/chi_sq_old)<1-1e-5:
            run += 1
            chi_sq_old = chi_sq
            Delta = sum(1./errYs**2)*sum(X1s**2/errYs**2) 
            Delta -= (sum(X1s/errYs**2))**2
            b = sum(X1s**2/errYs**2)*sum(Ys/errYs**2) - sum(X1s/errYs**2)*sum(X1s*Ys/errYs**2)
            b /= Delta
            a = sum(1/errYs**2)*sum(X1s*Ys/errYs**2) - sum(X1s/errYs**2)*sum(Ys/errYs**2)
            a /= Delta
            chi_sq = sum((X1s-Xs)**2/errXs**2 + (Ys-a*X1s-b)**2/errYs**2)
            X1s = Xs/errXs**2 + a*(Ys-b)/errYs**2
            X1s /= 1/errXs**2 + a**2/errYs**2
        r_chi_sq = chi_sq/DoF
        return a, b, r_chi_sq, run
    def plot_linearfit2D(s, a, b, Xs, Ys, errXs, errYs):
        x = np.arange(min(Xs), max(Xs), float(max(Xs)-min(Xs))/100)
        plt.plot(x, a*x+b)
        plt.errorbar(Xs, Ys, xerr=errXs, yerr=errYs, fmt='.', markersize=3)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('linear fitting')
        #ax = plt.gca()
        #ax.invert_xaxis()
        plt.show()

class equatorial2galactic_coordinates:
    """
    transform position and proper motion in equatorial to galactic coordinates,
    uncertainties propagation is included.

    input: RA in h, Dec in deg, mu_a and mu_d in mas/yr
    """
    a0 = 192.85948 #in deg; (a0,d0) is the equatorial coordinates of the North Galactic Pole
    d0 = 27.12825 #in deg;
    l_NCP = 122.93192 #in deg; the Galactic longitude of the North Celestial Pole
    def __init__(s, RA, Dec, mu_a=0, mu_d=0, err_mu_a=1, err_mu_d=1):
        s.RA, s.Dec, s.mu_a, s.mu_d, s.err_mu_a, s.err_mu_d = RA, Dec, mu_a, mu_d, err_mu_a, err_mu_d
        if type(s.RA) != float:
            s.RA = dms2deg(s.RA)
        if type(s.Dec) != float:
            s.Dec = dms2deg(s.Dec)
    def equatorial_position2galactic_position(s):
        """
        output: l and b, both in deg
        """
        a = 15*s.RA - s.a0 #in deg
        a *= math.pi/180 #in rad
        d, d0 = math.pi/180*s.Dec, math.pi/180*s.d0 #in rad
        ## calculate b
        b = math.sin(d)*math.sin(d0) + math.cos(d)*math.cos(d0)*math.cos(a)
        b = math.asin(b) #in rad
        ## calculate l now
        cosl1 = math.sin(d)*math.cos(d0) - math.cos(d)*math.sin(d0)*math.cos(a)
        cosl1 /= math.cos(b)
        sinl1 = math.cos(d)*math.sin(a)/math.cos(b)
        l1 = math.atan2(sinl1, cosl1) #in rad
        l = s.l_NCP - 180/math.pi*l1 #in deg
        return l, b*180/math.pi
    def equatorial_proper_motion2galactic_proper_motion(s):
        l, b = s.equatorial_position2galactic_position()
        a = 15*s.RA - s.a0 #in deg
        a, d, d0, l_N, l, b = math.pi/180*np.array([a, s.Dec, s.d0, s.l_NCP, l, b]) #in rad
        A_G = np.mat([[0,             np.cos(b)              ], #galactic matrix
                      [np.cos(l_N-l), np.sin(b)*np.sin(l_N-l)]]) 
        A_E = np.mat([[-np.cos(d0)*np.sin(a), np.cos(d)*np.sin(d0)-np.sin(d)*np.cos(d0)*np.cos(a)],
                      [-np.cos(a),            np.sin(d)*np.sin(a)                                ]])
        A = A_G.I * A_E
        B = np.mat([[s.mu_a],
                    [s.mu_d]])
        solution = A * B
        mu_l, mu_b = solution[0,0], solution[1,0]
        return mu_l, mu_b
    def galactic_proper_motion2equatorial_proper_motion(s, mu_l, mu_b):
        l, b = s.equatorial_position2galactic_position()
        a = 15*s.RA - s.a0 #in deg
        a, d, d0, l_N, l, b = math.pi/180*np.array([a, s.Dec, s.d0, s.l_NCP, l, b]) #in rad
        A_G = np.mat([[0,             np.cos(b)              ], #galactic matrix
                      [np.cos(l_N-l), np.sin(b)*np.sin(l_N-l)]]) 
        A_E = np.mat([[-np.cos(d0)*np.sin(a), np.cos(d)*np.sin(d0)-np.sin(d)*np.cos(d0)*np.cos(a)],
                      [-np.cos(a),            np.sin(d)*np.sin(a)                                ]])
        A = A_G.I * A_E
        B = np.mat([[mu_l],
                    [mu_b]])
        solution = A.I * B
        mu_a, mu_d = solution[0,0], solution[1,0]
        return mu_a, mu_d
