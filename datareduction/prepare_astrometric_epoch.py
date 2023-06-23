#!/usr/bin/env python
import os, sys, ftplib, glob
from astropy.time import Time
from datetime import datetime

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

def find_obs_date_from_idifits(idifitsfile):
    header = open(idifitsfile, 'rb').readline()[:1000]
    msgs = header.decode().split('   ')
    for msg in msgs:
        if 'DATE-OBS' in msg:
            date_obs = msg.strip().split('=')[-1].strip()
            #[year, month, day] = date_obs.split("'")[1].split('-')
    return date_obs.split("'")[1].split('-')
def idifits2obs_month(idifitsfile):
    [year, month, day] = find_obs_date_from_idifits(idifitsfile)
    obs_date = datetime(int(year.strip()), int(month.strip()), int(day.strip()))
    output_str = obs_date.strftime("%b%y").lower() #e.g. Mar16
    return output_str
def find_idifits_file__and__get_obsmonth():
    idifitsfiles = glob.glob(r'*.idifits')
    print(idifitsfiles)
    obs_month = idifits2obs_month(idifitsfiles[0])
    return obs_month

def main():
    """
    main program
    """
    usage = "prepare_astrometric_epoch.py <.vex file>"
    if len(sys.argv) != 2:
        print(usage)
        sys.exit()

    gsfc_url = "gdc.cddis.eosdis.nasa.gov"
    eop_dir = "vlbi/gsfc/ancillary/solve_apriori/"
    eop_filename = "usno_finals.erp"

    vexfile = sys.argv[1]
    experiment = vexfile.split('.')[0].strip()
    VLBAkeyfileftproot = 'www.vlba.nrao.edu/astro/VOBS/astronomy/'
    sumfile = experiment + '.sum'

    if not os.path.exists(vexfile):
        print(("\n%s doesn't exists. trying to download one...\n" % vexfile))
        obsmonth = find_idifits_file__and__get_obsmonth()
        #obsmonth = raw_input("What's the observation month (in mmmyy format, e.g. mar19)\n")
        VLBAkeyfileftpdir = VLBAkeyfileftproot + obsmonth + '/' + experiment + '/'
        os.system('wget -T 5 %s%s' % (VLBAkeyfileftpdir, vexfile))
    if not os.path.exists(vexfile):
        print(("%s not found on ftp server; aborting\n" % vexfile))
        sys.exit()

    if os.path.exists(vexfile):
        print("Vex file found. Moving on.")
    
    vexin = open(sys.argv[1])
    vexlines = vexin.readlines()
    vexin.close()

    startfound = False
    stopfound = False
    for line in vexlines:
        if 'date' in line:
            obsdate  = line.split(':')[-1].split(' ')
            obsmonth = obsdate[-2].strip().lower()
            obsmonth = obsmonth + obsdate[-1].strip()[2:4]
            print(obsmonth)
        if 'MJD' in line:
            MJD = int(line.split(':')[-1])
        # setting start date parameters
        if "exper_nominal_start" in line:
            splitline = line.split('=')
            syear     = int(splitline[1][:4])
            sdoy      = int(splitline[1][5:8])
            #syy       = syear - 100*(syear/100)
            syy       = syear % 100 
            startfound = True
        # setting stop date parameters
        elif "exper_nominal_stop" in line:
            splitline = line.split('=')
            eyear     = int(splitline[1][:4])
            edoy      = int(splitline[1][5:8])
            #eyy       = eyear - 100*(eyear/100)
            eyy       = eyear % 100
            stopfound = True
        if startfound and stopfound:
            break

    if not (startfound and stopfound):
        print("Couldn't find start and/or stop date! Aborting.")
        sys.exit()

    if not os.path.exists(sumfile):
        print(("try to download %s.sum...\n" % experiment))
        VLBAkeyfileftpdir = VLBAkeyfileftproot + obsmonth + '/' + experiment + '/'
        os.system('wget -T 5 %s%s.sum' % (VLBAkeyfileftpdir, experiment))

    os.system("mkdir logs")
    os.system("mkdir tables")
    os.system("mkdir images")
    os.chdir("logs")
    print("\ndeleting logfiles...\n")
    try:
        os.system("rm *")
    except OSError:
        pass
    """DOWNLOADING ERP FILES"""
    #eop_page = ftpget(gsfc_url, eop_dir, eop_filename)
    os.system("curl -u anonymous:daip@nrao.edu --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp > usno_finals.erp")
    if not os.path.exists("usno_finals.erp"):
        print("\ndownload usno_finals.erp from another route\n")
        today=Time.now()
        if today.mjd-MJD<30:
            print("\nchoose the newer erp file\n")
            os.system("wget ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/usno500_finals.erp")
            os.rename("usno500_finals.erp","usno_finals.erp")
        else:
            os.system("wget ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/usno_finals.erp")
    #os.system('wget --auth-no-challenge "https://cddis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp"')
    #os.system("wget -4 ftp://cddis.gsfc.nasa.gov/gps/products/ionex/%04d/%03d/*.Z" % (syear, sdoy))
    """DOWNLOADING IONEX FILES"""
    os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/jplg%03d0.%02di.Z\
        > jplg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
    os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/igsg%03d0.%02di.Z\
        > igsg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
    os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/esag%03d0.%02di.Z\
        > esag%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
    os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/codg%03d0.%02di.Z\
        > codg%03d0.%02di.Z" % (syear, sdoy, sdoy, syy, sdoy, syy))
    os.system("gunzip igsg%03d0.%02di.Z" % (sdoy, syy))
    # if end date is different than start date, need to download two sets of files
    if edoy != sdoy:
        os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/jplg%03d0.%02di.Z\
            > jplg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
        os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/igsg%03d0.%02di.Z\
            > igsg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
        os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/esag%03d0.%02di.Z\
            > esag%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
        os.system("curl -u anonymous:haoding@swin.edu.au -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/gps/products/ionex/%04d/%03d/codg%03d0.%02di.Z\
            > codg%03d0.%02di.Z" % (eyear, edoy, edoy, eyy, edoy, eyy))
        os.system("gunzip igsg%03d0.%02di.Z" % (edoy, eyy))


if __name__ == "__main__":
    main()
