#!/usr/bin/env python
import os, sys
from astropy.time import Time

usage = "prepare_astrometric_epoch.py <.vex file>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

vexin = open(sys.argv[1])
vexlines = vexin.readlines()
vexin.close()

startfound = False
stopfound = False
for line in vexlines:
    if 'MJD' in line:
        MJD = int(line.split(':')[-1])
    if "exper_nominal_start" in line:
        splitline = line.split('=')
        syear     = int(splitline[1][:4])
        sdoy      = int(splitline[1][5:8])
        syy       = syear - 100*(syear/100)
        startfound = True
    elif "exper_nominal_stop" in line:
        splitline = line.split('=')
        eyear     = int(splitline[1][:4])
        edoy      = int(splitline[1][5:8])
        eyy       = eyear - 100*(eyear/100)
        stopfound = True
    if startfound and stopfound:
        break

if not (startfound and stopfound):
    print "Couldn't find start and/or stop date! Aborting."
    sys.exit()

os.system("mkdir logs")
os.system("mkdir tables")
os.system("mkdir images")
os.chdir("logs")
os.system("rm *")
os.system("wget -T 12 https://vlbi.gsfc.nasa.gov/apriori/usno_finals.erp")
# in case the erp file service was down...
if not os.path.exists("usno_finals.erp"):
    print "\ndownload usno_finals.erp from another route\n"
    today=Time.now()
    if today.mjd-MJD<30:
        print "\nchoose the newer erp file\n"
        os.system("wget ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/usno500_finals.erp")
        os.rename("usno500_finals.erp","usno_finals.erp")
    else:
        os.system("wget ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/usno_finals.erp")

if edoy != sdoy:
    sdoy  = edoy
    syear = eyear
os.system("wget ftp://cddis.gsfc.nasa.gov/gps/products/ionex/%04d/%03d/*.Z" % (syear, sdoy))
if os.path.exists("igsg%03d0.%02di.Z" % (sdoy,syy)):
    print ("\ngunzip igsg%03d0.%02di.Z...\n" % (sdoy,syy))
    os.system("gunzip igsg%03d0.%02di.Z" % (sdoy,syy))
else:
    print ("\ngunzip jgsg%03d0.%02di.Z...\n" % (sdoy,syy))
    os.system("gunzip jplg%03d0.%02di.Z" % (sdoy, syy))
