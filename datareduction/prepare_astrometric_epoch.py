#!/usr/bin/python
import os, sys

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
os.system("wget https://vlbi.gsfc.nasa.gov/apriori/usno_finals.erp")
os.system("wget ftp://cddis.gsfc.nasa.gov/gps/products/ionex/%04d/%03d/*.Z" % (syear, sdoy))
os.system("gunzip jplg%03d0.%02di.Z" % (sdoy, syy))
if edoy != sdoy:
    os.system("wget ftp://cddis.gsfc.nasa.gov/gps/products/ionex/%04d/%03d/*.Z" % (eyear, edoy))
    os.system("gunzip jplg%03d0.%02di.Z" % (edoy, eyy))
