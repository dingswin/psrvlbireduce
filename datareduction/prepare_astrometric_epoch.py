#!/usr/bin/python2
import os, sys, ftplib

def ftpget(url, directory, filename):
    """Return contents of a file on an ftp-ssl site"""
    contents = []
    ftps = ftplib.FTP_TLS(url)
    # login and encrypt connection
    ftps.login()
    ftps.prot_p()
    ftps.cwd(directory)
    ftps.retrlines("RETR {:s}".format(filename), contents.append)

    return contents

usage = "prepare_astrometric_epoch.py <.vex file>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

gsfc_url = "gdc.cddis.eosdis.nasa.gov"
eop_dir = "vlbi/gsfc/ancillary/solve_apriori/"
eop_filename = "usno_finals.erp"

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
#os.system("wget https://vlbi.gsfc.nasa.gov/apriori/usno_finals.erp")
#os.system("wget -4 ftp://cddis.gsfc.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp")
eop_page = ftpget(gsfc_url, eop_dir, eop_filename)
#os.system('wget --auth-no-challenge "https://cddis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp"')
#os.system("wget -4 ftp://cddis.gsfc.nasa.gov/gps/products/ionex/%04d/%03d/*.Z" % (syear, sdoy))
os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/jplg%03d0.%02di.Z > jplg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/igsg%03d0.%02di.Z > igsg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/esag%03d0.%02di.Z > esag%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/codg%03d0.%02di.Z > codg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
os.system("gunzip jplg%03d0.%02di.Z" % (sdoy, syy))
if edoy != sdoy:
    os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/jplg%03d0.%02di.Z > jplg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
    os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/igsg%03d0.%02di.Z > igsg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
    os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/esag%03d0.%02di.Z > esag%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
    os.system("curl -u anonymous:adeller@astro.swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/codg%03d0.%02di.Z > codg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
    os.system("gunzip jplg%03d0.%02di.Z" % (edoy, eyy))
