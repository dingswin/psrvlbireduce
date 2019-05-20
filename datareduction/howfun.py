import math,os,sys
import numpy as np

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
    degree = ((b[2]/60+b[1])/60+b[0])*sign
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
    RA0 = np.array([dms2deg(RA1), dms2deg(RA2)])
    Dec = np.array([dms2deg(Dec1), dms2deg(Dec2)])
    Dec_rad = Dec*math.pi/180
    RA = RA0*np.cos(Dec_rad)*15 # hour to deg
    diffRA = abs(RA[0]-RA[1])*60 #deg to min
    diffDec = abs(Dec[0]-Dec[1])*60 
    sep = math.sqrt(diffRA**2+diffDec**2)
    return sep

def colonizedms(string):
    if string.startswith('-'):
        sign = -1
    else:
        sign = 1
    string = filter(str.isdigit,string)
    newstr = string[:2] + ':' + string[2:4] + ':' + string[4:6] + '.' + string[6:]
    if sign == -1:
        newstr = "-" + newstr
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
    dms="%02d:%02d:%09.7f" % (d,m,s)
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

def sample2estimate(array,confidencelevel):
    array.sort()
    CL = confidencelevel
    SV = int(CL*len(array)) #SV -> significant volume
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

def sample2uncertainty(array,estimate,confidencelevel): #offered with an estimate and present symmetric format
    array.sort()
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
    mu_as = np.array([])
    mu_ds = np.array([])
    PIs   = np.array([])
    RAs   = np.array([])
    Decs  = np.array([])
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
    if len(mu_as) == 1:
        mu_as = mu_as[0]
    if len(mu_ds) == 1:
        mu_ds = mu_ds[0]
    if len(PIs) == 1:
        PIs = PIs[0]
    if len(RAs) == 1:
        RAs = RAs[0]
    if len(Decs) == 1:
        Decs = Decs[0]
    return PIs,mu_as,mu_ds,RAs,Decs

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
