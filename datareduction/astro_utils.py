#!/usr/bin/env python

import math, os, sys

def getParallaxExtrema(ra, dec):
    os.system("rm -f pxextremajunk*.txt")
    output = open("pxextremajunk1.txt", "w")
    output.write("pi = 1\n")
    output.write("mu_a = 0\n")
    output.write("mu_d = 0\n")
    output.write("55000 %s 0.01 %s 0.01\n" % (ra, dec))
    output.write("55180 %s 0.01 %s 0.01\n" % (ra, dec))
    output.close()
    os.system("pmpar pxextremajunk1.txt -m 55555 > pxextremajunk2.txt")
    input = open("pxextremajunk2.txt")
    lines = input.readlines()
    input.close()
    ra1splitline = lines[0].split()
    ra2splitline = lines[1].split()
    totalsplitline = lines[4].split()
    radate1 = ra1splitline[11] + " " + ra1splitline[12]
    radate2 = ra2splitline[11] + " " + ra2splitline[12]
    totaldate1 = totalsplitline[12] + " " + totalsplitline[13][:-1]
    totaldate2 = totalsplitline[15] + " " + totalsplitline[16]
    ratio = float(ra1splitline[5])/float(totalsplitline[4])
    return [radate1 + ", " + radate2, totaldate1 + ", " + totaldate2, ratio]

def stringToRad(posstr, was_hms):
    posstr = posstr.strip()
    if was_hms and "h" in posstr:
        posstr = posstr.replace("h", ":")
        posstr = posstr.replace("m", ":")
        posstr = posstr.replace("s", "")
    if not was_hms and "d" in posstr:
        posstr = posstr.replace("d", ":")
        posstr = posstr.replace("'", ":")
        posstr = posstr.replace("\"", "")
    splits = posstr.split(':')
    try:
        d = int(splits[0])
        m = 0
        s = 0
        negmult = 1.0
        if len(splits) > 1 and len(splits[1].strip()) > 0:
            m = int(splits[1])
        if len(splits) > 2 and len(splits[2].strip()) > 0:
            s = float(splits[2])
        if posstr[0] == "-":
            d = -d
            negmult = -1.0
        radval = d + m/60.0 + s/3600.0
        radval *= math.pi/180.0
        radval *= negmult
        if was_hms:
            radval *= 180.0/12.0
    except ValueError:
        print "Bad position string", posstr
        radval = -999
    return radval

def getPosError(rastr, raerrstr, decstr, decerrstr):
    splitra = rastr.split(':')
    raerr = 1.0
    if len(splitra) == 1:
        raerr = 3600.0
    elif len(splitra) == 2:
        raerr = 60.0
    elif len(splitra[2]) > 2:
        for i in range(len(splitra[2])-3):
            raerr /= 10.0
    raerr *= 15*math.cos(stringToRad(decstr, False))*float(raerrstr)
    splitdec = decstr.split(':')
    decerr = 1.0
    if len(splitdec) == 1:
        decerr = 3600.0
    elif len(splitdec) == 2:
        decerr = 60.0
    elif len(splitdec[2]) > 2:
        for i in range(len(splitdec[2])-3):
            decerr /= 10.0
    decerr *= float(decerrstr)
    return math.sqrt(raerr*raerr + decerr*decerr)

def posdiff(targetrarad, targetdecrad, calrarad, caldecrad):
    sinsqdecdiff = math.sin((targetdecrad-caldecrad)/2.0)
    sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
    sinsqradiff  = math.sin((targetrarad-calrarad)/2.0)
    sinsqradiff  = sinsqradiff*sinsqradiff

    return 2*math.asin(math.sqrt(sinsqdecdiff +
                       math.cos(targetdecrad)*math.cos(caldecrad)*sinsqradiff))

def posradians2string(rarad, decrad):
    rah = rarad * 12 / math.pi
    rhh = int(rah)
    rmm = int(60*(rah - rhh))
    rss = 3600*rah - (3600*rhh + 60*rmm)
    decd = decrad * 180 / math.pi
    decformat = "+%02d:%02d:%010.7f"
    if decd < 0:
        decd = -decd
        decformat = '-' + decformat[1:]
    ddd = int(decd)
    dmm = int(60*(decd - ddd))
    dss = 3600*decd - (3600*ddd + 60*dmm)
    rastring  = "%02d:%02d:%011.8f" % (rhh,rmm,rss)
    decstring = decformat % (ddd, dmm, dss)
    return rastring, decstring

def name2approxradians(name):
    strippedname = name
    while len(strippedname) > 0 and strippedname[0].isalpha():
        strippedname = strippedname[1:]
    if len(strippedname) == 0:
        print "Unparseable name ", name
    decmultiplier = 1.0
    if '+' in strippedname:
        splitname = strippedname.split('+')
    elif '-' in strippedname:
        splitname = strippedname.split('-')
        decmultiplier = -1.0
    else:
        print "Unparseable name ", name
        sys.exit()
    
    hh = float(splitname[0][0:2])
    if len(splitname[0]) > 2:
        mm = float(splitname[0][2:4])
        if len(splitname[0]) > 4:
            ss = float(splitname[0][4:])
            if len(splitname[0][4:]) <= 2:
                ss += 0.5
        else:
            mm += 0.5
            ss = 0.0
    else:
        hh += 0.5
        mm = 0.0
        ss = 0.0
    raradapprox = (hh/12.0 + mm/(12*60.0) + ss/(12*60*60.0))*math.pi

    dd = float(splitname[1][0:2])
    if len(splitname[1]) > 2:
        mm = float(splitname[1][2:4])
        if len(splitname[1]) > 4:
            ss = float(splitname[1][4:])
            if len(splitname[1][4:]) <= 2:
                ss += 0.5
        else:
            mm += 0.5
            ss = 0.0
    else:
        dd += 0.5
        mm = 0.0
        ss = 0.0
    decradapprox = decmultiplier*(dd/180.0 + mm/(180*60.0) + ss/(180*60*60.0))*math.pi

    return raradapprox, decradapprox

def posradians2name(rarad, decrad, radecimals=-1, decdecimals=-1):
    absdecrad = math.fabs(decrad)
    rah = rarad * 12 / math.pi
    rhh = int(rah)
    rmm = int(60*(rah - rhh))
    rss = 3600*rah - (3600*rhh + 60*rmm)
    decd = absdecrad * 180 / math.pi
    ddd = int(decd)
    dmm = int(60*(decd - ddd))
    dss = 3600*decd - (3600*ddd + 60*dmm)

    if radecimals < 0:
       rastring = "J%02d%02" % (rhh, rmm)
    elif radecimals == 0:
       rastring = "J%02d%02d%02d" % (rhh, rmm, int(rss))
    else:
       raformat = "J%%02d%%02d%%02f.%d" % radecimals
       rastring = raformat % (rhh, rmm, rss)
    if decrad < 0:
       sign = '-'
    else:
       sign = "+"
    if decdecimals < 0:
       decstring = "%02d%02d" % (ddd, dmm)
    elif decdecimals == 0:
       decstring = "%02d%02d%02d" % (ddd, dmm, int(dss))
    else:
       decformat = "%%02d%%02d%%02f.%d" % decdecimals
       decstring = decformat % (rhh, rmm, rss)

    return rastring + sign + decstring

def ymd2mjd(year, month, day):
    return year*367 - int(7*(year + int((month + 9)/12))/4) + \
           int(275*month/9) + day - 678987

def mjd2ymdhms(mjd):
    imjd = int(mjd)
    fmjd = mjd - imjd

    j = imjd + 32044 + 2400001
    g = j / 146097
    dg = j % 146097
    c = ((dg/36524 + 1)*3)/4
    dc = dg - c*36524
    b = dc / 1461
    db = dc % 1461
    a = ((db/365 + 1)*3)/4
    da = db - a*365
    y = g*400 + c*100 + b*4 + a
    m = (da*5 + 308)/153 - 2
    d = da - ((m + 4)*153)/5 + 122

    year = y - 4800 + (m + 2)/12;
    month = (m + 2)%12 + 1;
    day = d + 1;

    hour = int(24*fmjd)
    minute = int(1440*fmjd - 60*hour)
    second = 86400*fmjd - (3600*hour + 60*minute)

    return year, month, day, hour, minute, second 
