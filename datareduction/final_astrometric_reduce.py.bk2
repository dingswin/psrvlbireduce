#!/usr/bin/env ParselTongue
################################################################################
## final_astrometric_reduce.py: A ParselTongue script for phase referencing VLBA 
## data and getting positions of pulsars and inbeam calibrators
## Adam Deller, 08 Jan 2013
################################################################################

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

################################################################################
# General python imports
################################################################################
import sys, os, string, math, warnings, subprocess, yaml, glob
import interaction, vlbatasks
from time import gmtime, strftime
from optparse import OptionParser
warnings.defaultaction = "always"


################################################################################
# Little logger class to put print statements to the log file
################################################################################
class Logger(object):
    def __init__(self, loghandle):
        self.terminal = sys.stdout
        self.log = loghandle

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def isatty(self):
        return False

################################################################################
# .source file parser
################################################################################
def parsesourcefile(sourcefile):
    sourcein = open(sourcefile)
    sourcelines = sourcein.readlines()
    sourcein.close()

    if not "GATED" in sourcelines[0]:
        print "Error parsing source file, no GATED on first line"
        sys.exit()
    gateduvfile = sourcelines[0].split(':')[1].strip()
    if not "UNGATED" in sourcelines[1]:
        print "Error parsing source file, no UNGATED on second line"
        sys.exit()
    ungateduvfile = sourcelines[1].split(':')[1].strip()
    numinbeams = 0
    inbeamfiles = []
    targetnames = []
    inbeamnames = []
    phscalnames = []
    inbeamuvdatas = []
    atline = 2
    while "INBEAM" in sourcelines[atline]:
        inbeamfiles.append(sourcelines[atline].split(':')[1].strip())
        inbeamuvdatas.append(AIPSUVData(experiment.upper() + "_I" + str(numinbeams+1), klass, 1, uvsequence))
        numinbeams += 1
        atline += 1
    if not "BANDPASS" in sourcelines[atline]:
        print "Error parsing source file, no BANDPASS on line %d" % atline
        sys.exit()
    ampcalsrc = sourcelines[atline].split(':')[-1].strip()
    atline += 1
    numtargets = int(sourcelines[atline].split(':')[1])
    atline += 1
    for i in range(numtargets):
        targetnames.append(sourcelines[atline].split(':')[1].strip())
        atline += 1
        phscalnames.append(sourcelines[atline].split(':')[1].strip())
        atline += 1
        #print "PHSREF: " + phscalnames[-1]
        inbeamnames.append([])
        while atline < len(sourcelines) and "INBEAM" in sourcelines[atline]:
            inbeamnames[-1].append(sourcelines[atline].split(':')[1].strip())
            #print sourcelines[atline].split(':')[1].strip()
            atline += 1
    return gateduvfile, ungateduvfile, numinbeams, inbeamfiles, inbeamuvdatas, \
           targetnames, inbeamnames, phscalnames, ampcalsrc

################################################################################
# Inbeam selfcal (phase only or amplitude and phase)
################################################################################
def inbeamselfcal(doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata,
                  expconfig, targetconfigs, modeldir, modeltype, targetonly,
                  calonly, beginif, endif, doampcal, dosecondary, sumifs,
                  clversion, targetnames, leakagedopol = 0):
    normed = []
    tocaluvdata = []
    tocalimagedata = []
    tocalnames = []
    tocalindices = []
    tocalconfigs = []
    for i in range(numtargets):
        normed.append([])
        ######### Need to figure out a better way of inserting the "not simple" information...
        #tocaluvdata.append(None()
        #tocalimagedata.append(None)
        #tocalnames.append(targetnames[i])
        #tocalindices.append(i)
        #tocalconfigs.append(targetconfigs[i])

    if not targetonly:
        simplecal = True
        print doneinbeams
        for (inbeamsrc,targetfilenum) in zip(doneinbeams,inbeamfilenums):
            parenttarget = -1
            for i in range(numtargets):
                 if parenttarget >= 0: break
                 for s in inbeamnames[i]:
                     if inbeamsrc == s:
                         parenttarget = i
                         break
                 if inbeamsrc == targetnames[i]:
                     parenttarget = i
                     break
            if parenttarget < 0:
                print "Couldn't find a parent target for " + inbeamsrc + \
                      " - this should never happen, aborting!"
                sys.exit()
            config = targetconfigs[parenttarget]
            print "Parenttarget is", parenttarget
            print "Sumifs is", sumifs
            print "Doampcal is", doampcal
            print "Desecondary is", dosecondary
            solmins = -1
            if doampcal:
                if not sumifs:
                    try:
                        solmins = config['inbeamcalibapnmins']
                    except KeyError:
                        print "No key for inbeamcalibapnmins, will skip inbeamapn calib for this source"
                else:
                    try:
                        solmins = config['inbeamcalibap1mins']
                    except KeyError: pass
            else:
                if dosecondary:
                    if not sumifs:
                        print "Can't do separate IFs secondary!"
                        sys.exit()
                    try:
                        solmins = config['inbeamcalibsp1mins']
                    except KeyError: pass
                else:
                    if sumifs:
                        solmins = config['inbeamcalibp1mins']
                    else:
                        try:
                            solmins = config['inbeamcalibpnmins']
                        except KeyError: pass
            if solmins < 0:
                print "Skipping " + inbeamsrc
                if sumifs:
                    if doampcal:
                        sntablepath = tabledir + inbeamsrc + '.icalib.ap1.sn'
                        pspath  = tabledir + inbeamsrc + '.icalib.ap1.ps'
                    elif dosecondary:
                        sntablepath = tabledir + inbeamsrc + '.icalib.sp1.sn'
                        pspath  = tabledir + inbeamsrc + '.icalib.sp1.ps'
                    else:
                        sntablepath = tabledir + inbeamsrc + '.icalib.p1.sn'
                        pspath      = tabledir + inbeamsrc + '.icalib.p1.ps'
                else:
                    if dosecondary:
                        print "Can't do separate IFs secondary!"
                        sys.exit()
                    if doampcal:
                        sntablepath = tabledir + inbeamsrc + '.icalib.apn.sn'
                        pspath  = tabledir + inbeamsrc + '.icalib.apn.ps'
                    else:
                        sntablepath = tabledir + inbeamsrc + '.icalib.pn.sn'
                        pspath  = tabledir + inbeamsrc + '.icalib.pn.ps'
                print "Removing", sntablepath, " and ", pspath
                os.system("rm -f " + sntablepath)
                os.system("rm -f " + pspath)
                continue
            if len(config['primaryinbeam'].split(','))  == 1 or config['primaryinbeam'].split(',')[0] == inbeamsrc:  
                tocalindices.append(parenttarget)
            subtractif = 0
            if config['skiplastif']:
                subtractif = 1
            shortname = inbeamsrc
            if len(inbeamsrc) > 12:
                shortname = inbeamsrc[:12]
            inbeam_image_data = None
            for j in range(10,0,-1):
                inbeam_uv_data = AIPSUVData(shortname, 'CALIB', 1, j)
                if inbeam_uv_data.exists():
                    inbeam_uv_data.zap()
            if targetfilenum >= 0:
                uvdata = inbeamuvdatas[targetfilenum]
            else:
                uvdata = gateduvdata
            domulti = False
            combineifs = False
            haveampcal = False
            if expconfig['ampcalscan'] > 0:
                haveampcal = True
            # split out inbeam_uv_data=inbeamsrc.CALIB.1 from inbeamuvdata
            vlbatasks.splittoseq(uvdata, clversion, 'CALIB', inbeamsrc,
                                 1, domulti, haveampcal, beginif, 
                                 endif-subtractif, combineifs, leakagedopol)
            inbeam_uv_data.table('NX', 1).zap()
            if targetfilenum >= 0:
                inbeam_image_file = modeldir + inbeamsrc + ".clean.fits"
                rawuvoutputfile = directory + inbeamsrc + ".formodeling.uv.fits"
                if not os.path.exists(inbeam_image_file):
                    print "Can't find " + modeltype + " inbeam model  " + inbeam_image_file
                    if modeltype == "preliminary":
                        print "I will write out a data file for this inbeam to " + rawuvoutputfile
                        print "Please image it with your favourite tool"
                        print "(clean only if using difmap, no modelfitting)"
                        print "When complete, copy the image fits file to " + inbeam_image_file
                        vlbatasks.writedata(inbeam_uv_data, rawuvoutputfile, True)
                    else:
                        print "Aborting!!"
                    sys.exit(1)
                inbeam_image_data = AIPSImage(shortname, "CLEAN", 1, 1)
                if inbeam_image_data.exists():
                    inbeam_image_data.zap()
                vlbatasks.fitld_image(inbeam_image_file, inbeam_image_data)
                if not dosecondary and len(config['primaryinbeam'].split(',')) > 1:
                    simplecal = False
                    normdata = AIPSUVData(shortname, 'NORMUV', 1, 1)
                    if normdata.exists():
                        normdata.zap()
                    vlbatasks.normaliseUVData(inbeam_uv_data, inbeam_image_data, normdata)
                    normed[parenttarget].append(normdata)
                    if config['primaryinbeam'].split(',')[0] == inbeamsrc:
                        tocaluvdata.append(AIPSUVData('CONCAT' + str(parenttarget), 'UVSRT', 1, 1))
                        tocalimagedata.append(None)
                        tocalnames.append('CONCAT' + str(parenttarget))
                        tocalconfigs.append(config)
                elif config['separateifmodel']:
                    simplecal = False
                    for i in range(beginif,endif+1-subtractif):
                        normdata = AIPSUVData(shortname, 'NORMUV', 1, i)
                        if normdata.exists():
                            normdata.zap()
                        inbeam_image_file = "%s%s.IF%d.clean.fits" % (modeldir, inbeamsrc, i)
                        if not os.path.exists(inbeam_image_file):
                            print "Can't find inbeam model file (multi-IF %d) %s" % (i, inbeam_image_file)
                            sys.exit()
                        inbeam_image_data = AIPSImage(shortname, "CLEAN", 1, i)
                        if inbeam_image_data.exists():
                            inbeam_image_data.zap()
                        vlbatasks.fitld_image(inbeam_image_file, inbeam_image_data)
                        vlbatasks.normaliseUVData(inbeam_uv_data, inbeam_image_data, normdata, 
                                                  i, i)
                        normed[parenttarget].append(normdata)
                        normdata.table('NX', 1).zap()
                    tocaluvdata.append(AIPSUVData('CONCAT' + str(parenttarget), 'UVSRT', 1, 1))
                    tocalimagedata.append(None)
                    tocalnames.append('CONCAT' + str(parenttarget))
                    tocalconfigs.append(config)
                else:
                    tocaluvdata.append(inbeam_uv_data)
                    tocalimagedata.append(inbeam_image_data)
                    tocalnames.append(inbeamsrc)
                    tocalconfigs.append(config)
            else:
                tocaluvdata.append(inbeam_uv_data)
                tocalimagedata.append(None)
                tocalnames.append(inbeamsrc)
                tocalconfigs.append(config)
        if not simplecal:
            for i in range(numtargets):
                solmins = -1
                if doampcal:
                    if sumifs:
                        try:
                            solmins = targetconfigs[i]['inbeamcalibap1mins']
                        except KeyError: pass
                    else:
                        try:
                            solmins = targetconfigs[i]['inbeamcalibapnmins']
                        except KeyError: pass
                else:
                    if dosecondary:
                        try:
                            solmins = targetconfigs[i]['inbeamcalibsp1mins']
                        except KeyError: pass
                    else:
                        if sumifs:
                            solmins = config['inbeamcalibp1mins']
                        else:
                            try:
                                solmins = config['inbeamcalibpnmins']
                            except KeyError: pass
                if solmins < 0:
                    continue
                if len(normed[i]) == 1:
                    #tocaluvdata.append(normed[i][0])
                    normed[i][0].rename('CONCAT' + str(i), 'UVSRT', 1)
                elif len(normed[i]) == 0:
                    print "No need to do anything for " + targetnames[i]
                    continue
                else:
                    for j in range(10,0,-1):
                        concatuvdata = AIPSUVData('CONCAT' + str(i), 'DBCON', 1, j)
                        if concatuvdata.exists():
                            concatuvdata.zap()
                        calibuvdata = AIPSUVData('CONCAT' + str(i), 'CALIB', 1, j)
                        if calibuvdata.exists():
                            calibuvdata.zap()
                        uvsrtuvdata = AIPSUVData('CONCAT' + str(i), 'UVSRT', 1, j)
                        if uvsrtuvdata.exists():
                            uvsrtuvdata.zap()
                    for n in normed[i][1:]:
                        vlbatasks.match_headersource(normed[i][0], n)
                    vlbatasks.dbcon(normed[i], uvsrtuvdata)
                    #tocaluvdata.append(uvsrtuvdata)
                    #vlbatasks.uvsrt(concatuvdata, uvsrtuvdata)
                    #tocalnames.append('CONCAT' + str(i))
                    for n in normed[i]:
                        n.zap()
        for (inbeamsrc,config,inbeam_uv_data,inbeam_image_data) in \
             zip(tocalnames,tocalconfigs,tocaluvdata,tocalimagedata):
            dostokesi   = True
            try:
                inbeamuvrange = config['inbeamuvrange']
            except KeyError:
                inbeamuvrange = [0, 0]
            if sumifs:
                dostokesi = True
                if doampcal:
                    sntablepath = tabledir + inbeamsrc + '.icalib.ap1.sn'
                    pspath  = tabledir + inbeamsrc + '.icalib.ap1.ps'
                    solmins = config['inbeamcalibap1mins']
                    solsnr  = config['inbeamcalibap1snr']
                    soltype = config['inbeamcalibap1type']
                    try:
                        weightit = config['inbeamcalibap1weightit']
                    except KeyError:
                        weightit = 0
                elif dosecondary:
                    sntablepath = tabledir + inbeamsrc + '.icalib.sp1.sn'
                    pspath  = tabledir + inbeamsrc + '.icalib.sp1.ps'
                    solmins = config['inbeamcalibsp1mins']
                    solsnr  = config['inbeamcalibsp1snr']
                    soltype = config['inbeamcalibsp1type']
                    try:
                        weightit = config['inbeamcalibsp1weightit']
                    except KeyError:
                        weightit = 0
                    try:
                        inbeamuvrange = config['secondaryinbeamuvrange']
                    except KeyError:
                        inbeamuvrange = [0, 0]
                else:
                    sntablepath = tabledir + inbeamsrc + '.icalib.p1.sn'
                    pspath      = tabledir + inbeamsrc + '.icalib.p1.ps'
                    solmins     = config['inbeamcalibp1mins']
                    solsnr      = config['inbeamcalibp1snr']
                    soltype     = config['inbeamcalibp1type']
                    try:
                        weightit = config['inbeamcalibp1weightit']
                    except KeyError:
                        weightit = 0
            else:
                if dosecondary:
                    print "Can't do separate IFs secondary!"
                    sys.exit()
                if doampcal:
                    sntablepath = tabledir + inbeamsrc + '.icalib.apn.sn'
                    pspath  = tabledir + inbeamsrc + '.icalib.apn.ps'
                    solmins = config['inbeamcalibapnmins']
                    solsnr  = config['inbeamcalibapnsnr']
                    soltype = config['inbeamcalibapntype']
                    try:
                        weightit = config['inbeamcalibapnweightit']
                    except KeyError:
                        weightit = 0
                    try:
                        dostokesi = config['inbeamcalibapnstokesi']
                    except KeyError:
                        dostokesi = False
                else:
                    sntablepath = tabledir + inbeamsrc + '.icalib.pn.sn'
                    pspath      = tabledir + inbeamsrc + '.icalib.pn.ps'
                    solmins     = config['inbeamcalibpnmins']
                    solsnr      = config['inbeamcalibpnsnr']
                    soltype     = config['inbeamcalibpntype']
                    try:
                        weightit = config['inbeamcalibpnweightit']
                    except KeyError:
                        weightit = 0
                    try:
                        dostokesi = config['inbeamcalibpnstokesi']
                    except KeyError:
                        dostokesi = False
            print "Using solution interval " + str(solmins) + ', requiring S/N ' + str(solsnr)
            try:
                flagwheremodelbelow = expconfig['inbeamminmodelflux'] #Jy
            except KeyError:
                flagwheremodelbelow = -1
            if not alwayssaved or not os.path.exists(sntablepath):
                vlbatasks.singlesource_calib(inbeam_uv_data, inbeam_image_data,
                                             1, expconfig['refant'], doampcal, 
                                             solmins, dostokesi, soltype, 
                                             solsnr, sumifs, inbeamuvrange, weightit,
                                             flagwheremodelbelow)
                snoutver = 1
                if not expconfig['skipsnedt']:
                    vlbatasks.snedt(inbeam_uv_data, 1)
                    snoutver = 2
                if os.path.exists(sntablepath):
                    os.remove(sntablepath)
                vlbatasks.writetable(inbeam_uv_data, 'SN', snoutver, sntablepath)
                #if not expconfig['skipsnedt']:
                #     if os.path.exists(sntablepath + '.unedited'):
                #         os.remove(sntablepath + '.unedited')
                #     vlbatasks.writetable(inbeam_uv_data, 'SN', 1, sntablepath + '.unedited')
                type = 'PHAS'
                if doampcal:
                    type = 'AMP'
                if sumifs:
                    vlbatasks.plottops(inbeam_uv_data, 'SN', snoutver, type, 1, 1, 1,
                                       pspath)
                else:
                    vlbatasks.plottops(inbeam_uv_data, 'SN', snoutver, type, 0, 2, 4,
                                       pspath)
            inbeam_uv_data.zap()
    return tocalnames, tocalindices

def target2cals(targetname): #get phscal and bandpass cal given targetname
    auxdir = os.environ['PSRVLBAUXDIR']
    targetdir = auxdir + '/processing/' + targetname
    sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
    if sourcefiles == []:
        targetname = raw_input("What's the name of the root directory for this target?\n")
        targetdir = auxdir + '/processing/' + targetname
        sourcefiles = glob.glob(r'%s/*/*.source' % targetdir)
        if sourcefiles == []:
            print("source files not found; abort")
            sys.exit()
        
    sourcefiles.sort()
    sourcefile = sourcefiles[0] # this might cause problem in adhoc observations
    lines = open(sourcefile).readlines()
    for line in lines:
        if 'BANDPASS' in line:
            bpcal = line.split(':')[-1].strip()
        if 'PHSREF' in line:
            phscal = line.split(':')[-1].strip()
            cals = [phscal, bpcal]
    return cals

#def find_alternative_targetname(targetname):
#    expconfigfile = configdir + experiment + '.yaml'
#    if not os.path.exists(expconfigfile):
#        parser.error("Experiment config file %s does not exist!" % expconfigfile)
#    expconfig     = yaml.load(open(expconfigfile))
#    rootdir       = expconfig['rootdir']

def applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, 
                     expconfig, targetconfigs, targetonly, calonly, doampcal,
                     dosecondary, sumifs, clversion, snversion, inbeamnames, 
                     targetnames):
    phsrefnames = []
    for targetname in targetnames:
        [phscalname, junk] = target2cals(targetname)
        phsrefnames.append(phscalname)
    sncount = 0
    numinbeams = len(inbeamuvdatas)
    #for i in range(numtargets):
    for i in range(20):
        if not targetonly:
            for j in range(numinbeams):
                vlbatasks.deletetable(inbeamuvdatas[j], 'SN', snversion+sncount)
        if not calonly:
            vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
            if haveungated:
                vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
        sncount += 1
    sncount = 0
    for inbeamsrc in tocalnames:
        if sumifs:
            if doampcal:
                calibstring = 'ap1'
            elif dosecondary:
                calibstring = 'sp1'
            else:
                calibstring = 'p1'
                if int(dualphscal_setup[0])>0:
                    calibstring += '.dualphscal'
        else:
            if dosecondary:
                print "Can't do separate IFs secondary cal!"
                sys.exit()
            if doampcal:
                calibstring = 'apn'
            else:
                calibstring = 'pn'
                if int(dualphscal_setup[0])>0:
                    calibstring += '.dualphscal'
        print inbeamsrc
        if "CONCAT" in inbeamsrc:
            calibtablepath = "%sCONCAT%d.icalib.%s.sn" % (tabledir, sncount, calibstring)
        else:
            calibtablepath = "%s%s.icalib.%s.sn" % (tabledir, inbeamsrc, calibstring)
        ## for the first imbeamselfcalp1 or inbeamselfcalpn apply
        if not os.path.exists(calibtablepath):
            calibtablepath = calibtablepath.replace('.dualphscal', '')

        print "Applying table " + calibtablepath
        if targetonly and (not os.path.exists(calibtablepath)):
            print "For target-only, the SN file must exist already - aborting!"
            print calibtablepath
            sys.exit(1)
        if not targetonly:
            for i in range(numinbeams):
                vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath,
                                    snversion+sncount)
        if not calonly:
            vlbatasks.loadtable(gateduvdata, calibtablepath,
                                snversion+sncount)
            if haveungated:
                vlbatasks.loadtable(ungateduvdata, calibtablepath,
                                    snversion+sncount)
        sncount += 1
    if not targetonly:
        for i in range(numinbeams):
            #vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+sncount)
            vlbatasks.mergesntables(inbeamuvdatas[i], snversion, sncount, 
                                    expconfig['refant'])
    if not calonly:
        #vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
        #print "Deleted SN table version " + str(snversion+sncount)
        print "Merging SN tables between " + str(snversion) + " and " + str(snversion + sncount -1)
        vlbatasks.mergesntables(gateduvdata, snversion, sncount, expconfig['refant'])
        if haveungated:
            #vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
            vlbatasks.mergesntables(ungateduvdata, snversion, sncount, expconfig['refant'])
    if not targetonly:
        for i in range(numinbeams):
            for j in range(10):
                vlbatasks.deletetable(inbeamuvdatas[i], 'CL', clversion+j+1)
            sourcelist = []
            for j in tocalindices:
                if i < len(inbeamnames[j]):
                    sourcelist.append(inbeamnames[j][i]) 
            print "Applying inbeamsn for ", sourcelist, " to file ", i
            vlbatasks.applysntable(inbeamuvdatas[i], snversion+sncount, '2PT', 
                                   clversion, expconfig['refant'], sourcelist, 'CALP') #does not necessarily go through the sourcelist
        sourcelist = []
        for i in tocalindices: 
            sourcelist.append(phsrefnames[i]) 
        vlbatasks.applysntable(inbeamuvdatas[0], snversion+sncount, '2PT', 
                               clversion, expconfig['refant'], sourcelist, 'CALP')
    if not calonly:
        sourcelist = []
        for i in tocalindices:
            sourcelist.append(targetnames[i])
        for j in range(10):
            vlbatasks.deletetable(gateduvdata, 'CL', clversion+j+1)
        print sourcelist
        vlbatasks.applysntable(gateduvdata, snversion+sncount, '2PT', clversion, 
                               expconfig['refant'], sourcelist, 'CALP')
        if haveungated:
            for j in range(10):
                vlbatasks.deletetable(ungateduvdata, 'CL', clversion+j+1)
            vlbatasks.applysntable(ungateduvdata, snversion+sncount, '2PT', 
                                   clversion, expconfig['refant'], sourcelist, 'CALP')
    return sncount + 1

################################################################################
# Print a summary of CL versions, SN versions, and runlevel
################################################################################
def printTableAndRunlevel(runlevel, snversion, clversion, uvdata):
    maxclversion = -1
    maxsnversion = -1
    if uvdata.exists():
        tables = uvdata.tables
        for table in tables:
            if table[1][-2:] == "CL" and table[0] >= maxclversion:
                maxclversion = table[0]
            if table[1][-2:] == "SN" and table[0] >= maxsnversion:
                maxsnversion = table[0]
    prtstr = "Runlevel is %d, clversion is %d (from uvdata, maxclversion is %d)" % \
             (runlevel, clversion, maxclversion)
    prtstr2 = ", snversion is %d (from uvdata, maxsnversion is %d)" % \
              (snversion, maxsnversion)
    print prtstr + prtstr2

################################################################################
# Option parsing and defaulted global variables
################################################################################
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    try:
        aipsver = os.environ['AIPS_VERSION'].split('/')[-1]
    except KeyError:
        aipsver = '31DEC18'
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-e", "--experiment", dest="experiment", default="",
                  help="Experiment name")
parser.add_option("--targetonly", dest="targetonly", default=False,
                  action="store_true", help="Reduce target only")
parser.add_option("--calonly", dest="calonly", default=False,
                  action="store_true", help="Reduce calibrators only")
parser.add_option("-r", "--runlevel", dest="runlevel", default="1",
                  help="Runlevel to start[,stop] at")
parser.add_option("--alwayssaved", dest="alwayssaved", default=False,
                  action="store_true", help="Always used saved cal tables")
parser.add_option("--clearcatalog", dest="clearcatalog", default=False,
                  action="store_true", help="Zap existing catalog data")
parser.add_option("--zapallcaltables", dest="zapallcaltables", default=False,
                  action="store_true", help="Zap all SN and CL tables")
parser.add_option("--noimageoutofbeam", dest="noimageoutofbeam", default=False,
                  action="store_true", help="Don't split the amp cal source")
parser.add_option("--startlocaltv", dest="startlocaltv", default=False,
                  action="store_true", help="Start a local AIPS TV " + \
                                            "(for X11, VNC etc")
(options, junk) = parser.parse_args()
auxdir          = ""
rootdir         = ""
targetonly      = options.targetonly
calonly         = options.calonly
experiment      = options.experiment
zapallcaltables = options.zapallcaltables
imageoutofbeam  = not options.noimageoutofbeam
alwayssaved     = options.alwayssaved
startlocaltv    = options.startlocaltv
runsplit        = options.runlevel.split(',')
runfromlevel    = int(runsplit[0])
runtolevel      = 999
if len(runsplit) > 1:
    runtolevel = int(runsplit[1])


################################################################################
# Check validity of inputs, load up the config files
################################################################################
if experiment == "":
    parser.error("You must supply an experiment name")
if runfromlevel > runtolevel:
    parser.error("Runtolevel (" + str(runtolevel) + ") must be greater " \
                 "than or equal to runfromlevel (" + str(runfromlevel) + ")")
try:
    auxdir = os.environ['PSRVLBAUXDIR']
except KeyError:
    print "PSRVLBAUXDIR is not defined - aborting!"
    sys.exit(1)
try:
    codedir = os.environ['PSRVLBICODEDIR']
except KeyError:
    print "PSRVLBICODEDIR is not defined - aborting!"
    sys.exit(1)
configdir     = auxdir + '/configs/'
finalmodeldir = auxdir + '/sourcemodels/final/'
prelimmodeldir= auxdir + '/sourcemodels/preliminary/'
expconfigfile = configdir + experiment + '.yaml'
if not os.path.exists(expconfigfile):
    parser.error("Experiment config file %s does not exist!" % expconfigfile)
expconfig     = yaml.load(open(expconfigfile))
rootdir       = expconfig['rootdir']
try:
    if expconfig['uselocalmodels']:
        finalmodeldir = rootdir + '/models/final/'
        prelimmodeldir = rootdir + '/models/preliminary/'
except KeyError:
    pass
directory     = rootdir + '/' + experiment.lower() + '/'
tabledir      = directory + "/tables/"
logdir        = directory + "/logs/"
clversion     = 1
snversion     = 1
klass         = 'UVDATA'
uvsequence    = 1
logfile       = directory + "/" + experiment.lower() + ".datacheck.log"
logmode       = "a"
AIPS.userno   = expconfig['userno']
numtargets    = len(expconfig["targets"])
modeldir      = finalmodeldir
modeltype     = "final"
if options.clearcatalog:
    expconfig['clearcatalog'] = True
if expconfig['useprelimmodels']:
    modeldir  = prelimmodeldir
    modeltype = "preliminary"
targetconfigs = []
for i in range(numtargets):
    targetconfigfile = configdir + expconfig["targets"][i] + '.yaml'
    if not os.path.exists(targetconfigfile):
        parser.error("Target config file %s does not exist!" % targetconfigfile)
    targetconfigs.append(yaml.load(open(targetconfigfile)))
skiplastif = targetconfigs[0]['skiplastif']
if numtargets > 1:
    for config in targetconfigs[1:]:
        if skiplastif != targetconfigs[i]['skiplastif']:
            print "All targets must have the same value for skiplastif! Aborting"
            sys.exit()
if os.path.exists(logfile) and runfromlevel == 1:
    os.system("rm -f " + logfile)
    logmode = "w"
if startlocaltv:
    tv = AIPSTV("local")
    tv.start()
if alwayssaved:
    interaction.setalwaysuse(True)
try:
    dualphscal_setup = targetconfigs[0]['dualphscal'].split(',')
except KeyError:
    dualphscal_setup = ['-1','0']

try:
    triphscal_setup = targetconfigs[0]['triphscal'].split(';')
    dotriphscal = True
except KeyError:
    dotriphscal = False

if float(dualphscal_setup[0]) > 0 and dotriphscal:
    print("Can't do both dualphscal and triphscal; aborting...")
    sys.exit()


################################################################################
# Parse the source file and set some more variables
################################################################################
sourcefile = experiment.lower() + ".source"
if sourcefile[0] != '/':
    sourcefile = directory + sourcefile
gateduvdata    = AIPSUVData(experiment.upper() + "_G", klass, 1, uvsequence)
ungateduvdata  = AIPSUVData(experiment.upper() + "_U", klass, 1, uvsequence)
alluvdatas     = [gateduvdata, ungateduvdata]
gateduvfile, ungateduvfile, numinbeams, inbeamfiles, inbeamuvdatas, \
targetnames, inbeamnames, phscalnames, ampcalsrc = parsesourcefile(sourcefile)
haveungated = True
if ungateduvfile == "":
    haveungated = False
    alluvdatas = [gateduvdata]
if calonly:
    alluvdatas = []
for uvdata in inbeamuvdatas:
    alluvdatas.append(uvdata)

################################################################################
# Save a record of the time and the arguments of this run, and start the log
################################################################################
AIPS.log = open(logfile, logmode)
sys.stdout = Logger(AIPS.log)
logout = open(directory + 'runlog.txt', 'a')
logout.write(strftime("%a, %d %b %Y %H:%M:%S +0000\n", gmtime()))
for a in sys.argv:
    logout.write(a + ' ')
logout.write('\n\n')
logout.close()

################################################################################
# Zap the existing cal tables if requested
################################################################################
if zapallcaltables and runfromlevel > 1:
    print "Zapping all SN, BP, and CL tables (except CL 1)!"
    for uvdata in alluvdatas:
        tables = uvdata.tables
        for table in tables:
            if table[1][-2:] == "SN":
                print "Zapping SN " + str(table[0]) + " from " + uvdata.name
                uvdata.table('SN', table[0]).zap()
            elif table[1][-2:] == "BP":
                print "Zapping BP " + str(table[0]) + " from " + uvdata.name
                uvdata.table('BP', table[0]).zap()
            elif table[1][-2:] == "CL" and table[0] > 1:
                print "Zapping CL " + str(table[0]) + " from " + uvdata.name
                uvdata.table('CL', table[0]).zap()

################################################################################
# Begin the actual data processing
################################################################################
runlevel = 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load the uv data ############################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Loading UVDATA from disk"
    vexfile  = directory + '/' + experiment.lower() + '.vex'
    refdate = vlbatasks.getrefdate(vexfile)
    try:
        cltablemins = expconfig['cltablemins']
    except KeyError:
        cltablemins = 0.166667
    if not targetonly:
        for i in range(numinbeams):
            os.system("rm -f %s/templink.fits" % directory)
            os.system("ln -s %s %s/templink.fits" % \
                      (inbeamfiles[i], directory))
            if expconfig['clearcatalog'] and inbeamuvdatas[i].exists():
                inbeamuvdatas[i].zap()
            vlbatasks.fitld_vlba("%s/templink.fits" % directory, inbeamuvdatas[i], [], 0.05, refdate, cltablemins)
    if not calonly:
        if (expconfig['clearcatalog'] or targetonly) and gateduvdata.exists():
            gateduvdata.zap(True)
        if (expconfig['clearcatalog'] or targetonly) and ungateduvdata.exists():
            ungateduvdata.zap(True)
        os.system("rm -f %s/templink.fits" % directory)
        os.system("ln -s %s %s/templink.fits" % (gateduvfile, directory))
        vlbatasks.fitld_vlba("%s/templink.fits" % directory, gateduvdata, targetnames, 0.05, refdate, cltablemins)
        if haveungated:
            os.system("rm -f %s/templink.fits" % directory)
            os.system("ln -s %s %s/templink.fits" % (ungateduvfile, directory))
            vlbatasks.fitld_vlba("%s/templink.fits" % directory, ungateduvdata, targetnames, 0.05, refdate, cltablemins)
else:
    print "Skipping UVDATA load"

runlevel = runlevel + 1
try:
    moriffactor = expconfig['moriffactor']
except KeyError:
    moriffactor = 0
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Increase the number of IFs if requested #####################################
if runfromlevel <= runlevel and runtolevel >= runlevel and moriffactor > 1:
    print "Runlevel " + str(runlevel) + ": Splitting each IF " + str(moriffactor) + " ways"
    if not targetonly:
        for i in range(numinbeams):
            newinbeamdata = AIPSUVData(experiment.upper() + "_I" + str(i+1), "MORIF", 1, uvsequence)
            if newinbeamdata.exists():
                newinbeamdata.zap()
            vlbatasks.morif(inbeamuvdatas[i], newinbeamdata, moriffactor)
            inbeamuvdatas[i].zap()
            newinbeamdata.rename(experiment.upper() + "_I" + str(i+1), klass, uvsequence)
    if not calonly:
        newgateddata = AIPSUVData(experiment.upper() + "_G", "MORIF", 1, uvsequence)
        if newgateddata.exists():
            newgateddata.zap()
        vlbatasks.morif(gateduvdata, newgateddata, moriffactor)
        gateduvdata.zap()
        newgateddata.rename(experiment.upper() + "_G", klass, uvsequence)
        if haveungated:
            newungateddata = AIPSUVData(experiment.upper() + "_U", "MORIF", 1, uvsequence)
            if newungateddata.exists():
                newungateddata.zap()
            vlbatasks.morif(ungateduvdata, newungateddata, moriffactor)
            ungateduvdata.zap()
            newungateddata.rename(experiment.upper() + "_U", klass, uvsequence)
else:
    print "Not running MORIF"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Unflag Pie Town if requested ################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and expconfig['unflagpietown']:
    print "Runlevel " + str(runlevel) + ": Unflagging Pie Town"
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.unflagantenna(inbeamuvdatas[i], 1, 1, 9)
    if not calonly:
        if haveungated:
            vlbatasks.unflagantenna(ungateduvdata, 1, 1, 9)
        vlbatasks.unflagantenna(gateduvdata, 1, 1, 9)
else:
    print "Not unflagging Pie Town"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load the user flags (if any) ################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Loading user flags"
    userflagfiles = glob.glob(tabledir + '/*.flag')
    if not targetonly:
        userflagfile = tabledir + "additionaledit.flag"
        if os.path.exists(userflagfile):
            for i in range(numinbeams):
                vlbatasks.userflag(inbeamuvdatas[i], 1, userflagfile)
        for i in range(numtargets):
            for j in range(len(inbeamnames[i])):
                #for inbeamname in inbeamnames[i]:
                extraflagfile = tabledir + "additionaledit." + inbeamnames[i][j] + ".flag"
                if os.path.exists(extraflagfile):
                    vlbatasks.userflag(inbeamuvdatas[j], 1, extraflagfile)
    if not calonly:
        userflagfile = tabledir + "additionaledit.flag"
        if os.path.exists(userflagfile):
            vlbatasks.userflag(gateduvdata, 1, userflagfile)
            if haveungated:
                vlbatasks.userflag(ungateduvdata, 1, userflagfile)
        for targetname in targetnames:
            extraflagfile = tabledir + "additionaledit." + targetname + ".flag"
            if os.path.exists(extraflagfile):
                vlbatasks.userflag(gateduvdata, 1, extraflagfile)
                if haveungated:
                    vlbatasks.userflag(ungateduvdata, 1, extraflagfile)
    if len(userflagfiles) == 0:
        print "No user flag file - skipping"
    for userflagfile in userflagfiles:
        if "target" in userflagfile:
            print "Warning: old-style additionaledit.target.flag file no longer supported!"
else:
    print "Skipping flagging from user flag file"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run the flagging for zero-fringe rate effects ###############################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Running UVFLG for zero fringe rates"
    try:
        suppressionfactor = targetconfigs[0]['fringerateflagsuppressionfactor']
    except KeyError:
        suppressionfactor = 7
    if not targetonly:
        for inbeamuvdata in inbeamuvdatas:
            vlbatasks.fringerateflag(inbeamuvdata, 1, suppressionfactor)
    if not calonly:
        if haveungated:
            vlbatasks.fringerateflag(ungateduvdata, 1, suppressionfactor)
        vlbatasks.fringerateflag(gateduvdata, 1, suppressionfactor)
else:
    print "Skipping flagging times of zero fringe rate"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run the autoflagger (if requested) ##########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and expconfig['autoflag']:
    print "Runlevel " + str(runlevel) + ": Running autoflagger"
    autoflagfile = tabledir + "auto.rawuvdata.flag"
    if targetonly and (not os.path.exists(autoflagfile)):
        print "For target-only, autoflag file %s must already exist!" % autoflagfile
        print "Aborting!"
        sys.exit(1)
    if not targetonly:
        if not os.path.exists(autoflagfile) or \
           not interaction.yesno("Do you wish to used saved autoflag table?"):
            runline =  "./autoflaggerplotter.py "
            runline += "--i=%d.%s_I1.UVDATA.1.1 " % (AIPS.userno, experiment.upper())
            runline += "--makeplot "
            runline += "--usepickle=%s/%s.aipsdata.pickle " % (tabledir, experiment.lower())
            runline += "--noint --flagfile=%s" % (autoflagfile)
            os.system(runline)
        for i in range(numinbeams):
            vlbatasks.userflag(inbeamuvdatas[i], 1, autoflagfile)
    if not calonly:
        vlbatasks.userflag(gateduvdata, 1, autoflagfile)
        if haveungated:
            vlbatasks.userflag(ungateduvdata, 1, autoflagfile)
else:
    print "Skipping autoflagging"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Duplicate the initial CL table ##############################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Duplicating CL table 1"
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.tacop(inbeamuvdatas[i], 'CL', 1, inbeamuvdatas[i], 2)
    if not calonly:
        vlbatasks.tacop(gateduvdata, 'CL', 1, gateduvdata, 2)
        if haveungated:
            vlbatasks.tacop(ungateduvdata, 'CL', 1, ungateduvdata, 2)
else:
    print "Skipping duplication of initial CL table"

runlevel = runlevel + 1
clversion = clversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Correct positions ###########################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Running CLCOR to shift sources"
    uvdatas = alluvdatas
    if targetonly:
        uvdatas = [gateduvdata]
        if haveungated:
            uvdatas = [gateduvdata, ungateduvdata]
    elif calonly:
        uvdatas = inbeamuvdatas
    try:
        shifts = expconfig["shifts"]
    except KeyError:
        shifts = []
        print "No sources to shift"
    for s in shifts:
        if s == "": continue
        srcname  = s.split(',')[0]
        rashift  = float(s.split(',')[1])
        decshift = float(s.split(',')[2])
        print "Shifting " + srcname + " by " + str(rashift) + "," + str(decshift)
        foundone = False
        for dataset in uvdatas:
            present = False
            for row in dataset.table('SU', 1):
                print row.source.strip()
                if row.source.strip() == srcname:
                    present = True
                    foundone = True
            if present:
                print "About to operate on ", dataset, srcname, rashift, decshift, clversion
                vlbatasks.shift_source(dataset, srcname, rashift, decshift,
                                       clversion)
        if not foundone:
            print "Didn't find source " + srcname + " in any datasets!"
            sys.exit()
else:
    print "Skipping source shifting"

runlevel = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run TECOR ###################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not expconfig['skiptecor']:
    print "Runlevel " + str(runlevel) + ": Running TECOR to correct ionosphere"
    try:
        follow = expconfig['tecorfollow']
    except KeyError:
        print "Follow not specified in expconfig file, defaulting to 0.2"
        follow = 0.2
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.correct_iono(inbeamuvdatas[i], logdir, clversion, follow)
    if not calonly:
        vlbatasks.correct_iono(gateduvdata, logdir, clversion, follow)
        if haveungated:
            vlbatasks.correct_iono(ungateduvdata, logdir, clversion, follow)
else:
    print "Skipping ionospheric corrections"

runlevel  = runlevel + 1
if not expconfig['skiptecor']:
    clversion = clversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run CLCOR to correct EOPs ###################################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Running CLCOR to correct for EOPs"
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.correct_eops(inbeamuvdatas[i], logdir, clversion)
    if not calonly:
        vlbatasks.correct_eops(gateduvdata, logdir, clversion)
        if haveungated:
            vlbatasks.correct_eops(ungateduvdata, logdir, clversion)
else:
    print "Skipping EOP corrections"

runlevel  = runlevel + 1
clversion = clversion + 1
leakagedopol = 0
try:
    xpolscan = expconfig['xpolscan']
except KeyError:
    xpolscan = -1
try:
    xpolsource = expconfig['xpolsource']
except KeyError:
    xpolsource = ampcalsrc
try:
    leakagescan = expconfig['leakagescan']
    leakagedopol = 2
except KeyError:
    leakagescan = -1
try:
    leakagesource = expconfig['leakagesource']
except KeyError:
    leakagesource = ampcalsrc
try:
    leakageuvrange = expconfig['leakageuvrange']
except KeyError:
    leakageuvrange = [0,0]
try:
    leakageweightit = expconfig['leakageweightit']
except KeyError:
    leakageweightit = 0
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run CLCOR to correct PANG ###################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and leakagescan > 0:
    print "Runlevel " + str(runlevel) + ": Running CLCOR to correct PANG"
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.clcor_pang(inbeamuvdatas[i], clversion)
    if not calonly:
        vlbatasks.clcor_pang(gateduvdata, clversion)
        if haveungated:
            vlbatasks.clcor_pang(ungateduvdata, clversion)
else:
    print "Skipping PANG corrections"

runlevel  = runlevel + 1
if leakagescan > 0:
    clversion = clversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
adelayfile = tabledir + 'aprioridelays.txt'
gatedadelayfile = tabledir + 'aprioridelays.gated.txt'
## Correct for a priori delays if available using CLCOR ########################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    if os.path.exists(adelayfile):
        print "Runlevel " + str(runlevel) + ": Using CLCOR to correct measured delays"
    else:
        print "Runlevel " + str(runlevel) + ": Not correcting a priori delays, just duplicating CL table"
    uvdatas = alluvdatas
    if targetonly:
        uvdatas = [gateduvdata]
        if haveungated:
            uvdatas = [gateduvdata, ungateduvdata]
    elif calonly:
        uvdatas = inbeamuvdatas
    for uvdata in uvdatas:
        vlbatasks.tacop(uvdata, 'CL', clversion, uvdata, clversion+1)
        thisfile = adelayfile
        if uvdata == gateduvdata and os.path.exists(gatedadelayfile):
            useroverride = interaction.yesno("There is an override file for a priori delays for the gated dataset! Do you wish to use this file? ONLY IF THERE WAS A MISTAKE WITH THE CLOCKS AT CORRELATION TIME!!! : ")
            if useroverride:
                thisfile = gatedadelayfile
        if os.path.exists(thisfile):
            vlbatasks.clcordelaysfromfile(uvdata, thisfile, clversion+1)
else:
    print "Skipping correction of a priori clock delays"

runlevel  = runlevel + 1
clversion = clversion + 1
try:
    skipapcalaccor = expconfig['skipapcalaccor']
except KeyError:
    skipapcalaccor = False
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Do the amplitude calibration (either load existing table or do and then edit) and inspect
if runfromlevel <= runlevel and runtolevel >= runlevel and not skipapcalaccor:
    print "Runlevel " + str(runlevel) + ": Doing amplitude calibration"
    if targetonly and (not os.path.exists(tabledir + 'accor.sn') or \
                       not os.path.exists(tabledir + 'apcal.sn')):
        print "For target-only, the SN files must exist already - aborting!"
        sys.exit(1)
    uvdatas = alluvdatas
    if targetonly:
        uvdatas = [gateduvdata]
        if haveungated:
            uvdatas = [gateduvdata, ungateduvdata]
    elif calonly:
        uvdatas = inbeamuvdatas
    try:
        accorinterpol = expconfig['accorinterpol']
    except KeyError:
        accorinterpol = '2PT'
    if not targetonly and ((not os.path.exists(tabledir + 'accor.sn') and not\
        os.path.exists(tabledir + 'apcal.sn')) or not \
       interaction.yesno("Do you wish to used saved SN table of amplitude calibration?")):
        for i in range(7):
            vlbatasks.deletetable(inbeamuvdatas[0], 'SN', i+snversion)
        dotsysfix = False
        try:
            if expconfig['allowtsysfix'] and 'HN' in inbeamuvdatas[0].antennas:
                dotsysfix = True
        except KeyError:
            pass
        try:
            accorampsmo = expconfig['accorampsmo']
        except KeyError:
            accorampsmo = 0.03
        try:
            apcalampsmo = expconfig['apcalampsmo']
        except KeyError:
            apcalampsmo = 1.5
        try:
            accorinterval = expconfig['accorinterval']
        except KeyError:
            accorinterval = -2
        try:
            copytsys = expconfig['copytsys']
        except KeyError:
            copytsys = ""
        if copytsys != "" and dotsysfix:
            print "Can't do tsys fix and copytsys!"
            sys.exit()
        vlbatasks.accor(inbeamuvdatas[0], accorinterval)
        if accorampsmo > 0:
            vlbatasks.snsmo(inbeamuvdatas[0], 'BOTH', 20, accorampsmo, 0.0, 0.0, snversion, expconfig['refant'], True)
        else:
            vlbatasks.tacop(inbeamuvdatas[0], 'SN', snversion, inbeamuvdatas[0], snversion+1)
        ### load .antab file if available
        antabfile = logdir + '/' + experiment + '.antab'
        if os.path.exists(antabfile):
            for uvdata in alluvdatas:
                vlbatasks.antab(uvdata, antabfile, 1, 1)
        vlbatasks.apcal(inbeamuvdatas[0], snversion+2)
        if apcalampsmo > 0:
            vlbatasks.snsmo(inbeamuvdatas[0], 'BOTH', 20, apcalampsmo, 0.0, 0.0, snversion+2, expconfig['refant'], True)
        else:
            vlbatasks.tacop(inbeamuvdatas[0], 'SN', snversion+2, inbeamuvdatas[0], snversion+3)
        if dotsysfix:
            print "Doing tsys fix"
            vlbatasks.fixtsys(inbeamuvdatas[0], snversion+3)
        elif copytsys != "":
            groups = copytsys.split(':')
            for g in groups:
                copytsyssplit = g.split(',')
                if not len(copytsyssplit) == 3:
                    print "Copytsys was " + copytsys + ", must be --copytsys=AN1,AN2,scalefactor"
                    sys.exit()
                if copytsyssplit[0] == copytsyssplit[1]:
                    ifsplit = copytsyssplit[2].split('@')
                    if len(ifsplit) != 2 or not ifsplit[0][:2] == "if" or not ifsplit[1][:2] == "if":
                        print "If copying between IFs of the same antenna, must use"
                        print "--copytsys=AN,AN,ifX@ifY where X and Y are input and output IFs respectively (-1 for all)"
                    vlbatasks.copytsys(inbeamuvdatas[0], snversion+3, copytsyssplit[0],
                                       copytsyssplit[1], 1.0, int(ifsplit[0][2:]), int(ifsplit[1][2:]))
                else:
                    vlbatasks.copytsys(inbeamuvdatas[0], snversion+3, copytsyssplit[0],
                                       copytsyssplit[1], float(copytsyssplit[2]))
                if not g==groups[-1]:
                    vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+3)
                    vlbatasks.tacop(inbeamuvdatas[0], 'SN', snversion+4, inbeamuvdatas[0], snversion+3)
                    vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+4)
        else:
            vlbatasks.tacop(inbeamuvdatas[0], 'SN', snversion+3, inbeamuvdatas[0], snversion+4)
        snoutver1 = snversion + 1
        snoutver2 = snversion + 4
        if not expconfig['skipsnedt']:
            vlbatasks.snedt(inbeamuvdatas[0], snversion+1)
            vlbatasks.snedt(inbeamuvdatas[0], snversion+4)
            snoutver1 = snversion+5
            snoutver2 = snversion+6
        if os.path.exists(tabledir + 'accor.sn'):
            os.remove(tabledir + 'accor.sn')
        if os.path.exists(tabledir + 'apcal.sn'):
            os.remove(tabledir + 'apcal.sn')
        vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver1, tabledir + 'accor.sn')
        vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver2, tabledir + 'apcal.sn')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver1, 'AMP', 0, 2, 4, tabledir + 'accor.ps')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver2, 'AMP', 0, 2, 4, tabledir + 'apcal.ps')
        for i in range(7):
            vlbatasks.deletetable(inbeamuvdatas[0], 'SN', i+snversion)
    for uvdata in uvdatas:
        vlbatasks.loadtable(uvdata, tabledir + 'accor.sn', snversion)
        vlbatasks.loadtable(uvdata, tabledir + 'apcal.sn', snversion+1)
        vlbatasks.applysntable(uvdata, snversion, accorinterpol, clversion, expconfig['refant'])
        vlbatasks.applysntable(uvdata, snversion+1, '2PT', clversion+1, expconfig['refant'])
else:
    print "Skipping amplitude calibration"

runlevel  = runlevel + 1
if not skipapcalaccor:
    clversion = clversion + 2
    snversion = snversion + 2
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Correct for primary beam attenuation ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not expconfig['skippbcor']:
    print "Runlevel " + str(runlevel) + ": Doing primary beam correction"
    vexfile  = directory + '/' + experiment.lower() + '.vex'
    scanlist = vlbatasks.getvexscaninfo(vexfile)
    fieldsourcenames = {}
    fieldsourcenames[ampcalsrc] = ampcalsrc
    if not targetonly:
        for i in range(numinbeams):
            for j in range(numtargets):
                if i==0:
                    fieldsourcenames[phscalnames[j]] = phscalnames[j]
                if len(inbeamnames[j]) > i:
                    if expconfig['dodefaultnames']:
                        fieldsourcenames["TARGETPT"] = inbeamnames[j][i]
                    elif "174" in experiment:
                        fieldsourcenames[targetnames[j][:5] + "PT"] = inbeamnames[j][i]
                    else:
                        fieldsourcenames[targetnames[j] + "PT"] = inbeamnames[j][i]
            pbsntable = tabledir + 'pbcor.cal' + str(i) + '.sn'
            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion)
            vlbatasks.correct_primarybeam(inbeamuvdatas[i], snversion-1, i, scanlist, fieldsourcenames, False, False)
            if os.path.exists(pbsntable):
                os.remove(pbsntable)
            vlbatasks.writetable(inbeamuvdatas[i], 'SN', snversion, pbsntable)
            vlbatasks.plottops(inbeamuvdatas[i], 'SN', snversion, 'AMP', 0, 2, 4, tabledir + 'pbcor.cal' + str(i) + '.ps')
            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion)
            vlbatasks.loadtable(inbeamuvdatas[i], pbsntable, snversion)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion, '2PT', clversion, expconfig['refant'])
    if not calonly:
        for i in range(numtargets):
            if expconfig['dodefaultnames']:
                fieldsourcenames["TARGETPT"] = targetnames[i]
            elif "174" in experiment:
                fieldsourcenames[targetnames[j][:5] + "PT"] = targetnames[i]
            else:
                fieldsourcenames[targetnames[i] + "PT"] = targetnames[i]
        uvdatas = [gateduvdata]
        if haveungated:
            uvdatas = [gateduvdata, ungateduvdata]
        for uvdata in uvdatas:
            pbsntable = tabledir + 'pbcor.target.sn'
            vlbatasks.deletetable(uvdata, 'SN', snversion)
            vlbatasks.correct_primarybeam(uvdata, snversion-1, 0, scanlist, fieldsourcenames, False, False)
            if os.path.exists(pbsntable):
                os.remove(pbsntable)
            vlbatasks.writetable(uvdata, 'SN', snversion, pbsntable)
            vlbatasks.deletetable(uvdata, 'SN', snversion)
            vlbatasks.loadtable(uvdata, pbsntable, snversion)
            vlbatasks.applysntable(uvdata, snversion, 'SELN', clversion, expconfig['refant'])
        vlbatasks.plottops(gateduvdata, 'SN', snversion, 'AMP', 0, 2, 4, tabledir + 'pbcor.target.ps')
else:
    print "Skipping primary beam correction"

runlevel  = runlevel + 1
if not expconfig['skippbcor']:
    clversion = clversion + 1
    snversion = snversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Do PCAL correction and inspect ##############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    not expconfig['skippcal'] and expconfig['ampcalscan'] > 0:
    print "Runlevel " + str(runlevel) + ": Doing pulse cal calibration"
    if targetonly and (not os.path.exists(tabledir + 'pccor.sn')):
        print "For target-only, the PC file must exist already - aborting!"
        sys.exit(1)
    if  not (targetonly or (os.path.exists(tabledir + 'pccor.sn') and \
                            interaction.yesno("Do you wish to used saved SN table " + \
                                  "for pulse cal?"))):
        try:
            inbeamuvdatas[0].table('SN', snversion).zap()
        except IOError:
            print "No need to delete old SN table"
        vlbatasks.pccor(inbeamuvdatas[0], ampcalsrc, snversion, expconfig['ampcalscan'], expconfig['refant'])
        if os.path.exists(tabledir + 'pccor.sn'):
            os.remove(tabledir + 'pccor.sn')
        vlbatasks.writetable(inbeamuvdatas[0], 'SN', snversion, tabledir + 'pccor.sn')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snversion, 'PHAS', 0, 2, 4, tabledir + 'pccor.ps')
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion)
    if not calonly:
        vlbatasks.loadtable(gateduvdata, tabledir + 'pccor.sn', snversion)
        vlbatasks.applysntable(gateduvdata, snversion, '2PT', clversion, expconfig['refant'])
        if haveungated:
            vlbatasks.loadtable(ungateduvdata, tabledir + 'pccor.sn', snversion)
            vlbatasks.applysntable(ungateduvdata, snversion, '2PT', clversion, expconfig['refant'])
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'pccor.sn', snversion)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion, '2PT', clversion, expconfig['refant'])
else:
    print "Skipping pulse cal corrections"

runlevel  = runlevel + 1
if not expconfig['skippcal'] and expconfig['ampcalscan'] > 0:
    snversion = snversion + 1
    clversion = clversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run FRING (bandpass calibrator) #############################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    expconfig['ampcalscan'] > 0:
    if not os.path.exists(tabledir + 'ampcalfring.sn') or \
           not interaction.yesno("Do you wish to used saved SN table for ampcal FRING?"):
        print "Runlevel " + str(runlevel) + ": FRING'ing amplitude calibrator"
        sumifs = False
        sumrrll = False
        zerorates = True
        try:
            bandpassuvrange = expconfig['bandpassuvrange']
        except KeyError:
            bandpassuvrange = [0, 0]
        try:
            inttimesecs = expconfig['inttimesecs']
        except KeyError:
            inttimesecs = 2
        ampcalmodeldata = None
        ampcalmodelfile = modeldir + ampcalsrc + '.clean.fits'
        if os.path.exists(ampcalmodelfile):
            aipscalname = ampcalsrc
            if len(ampcalsrc) > 12:
                aipscalname = ampcalsrc[:12]
            print "Using " + modeltype + " model for " + ampcalsrc
            ampcalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
            if ampcalmodeldata.exists():
                ampcalmodeldata.zap()
            vlbatasks.fitld_image(ampcalmodelfile, ampcalmodeldata)
        try:
            delaywin = expconfig['delaywindowns']
        except KeyError:
            delaywin = 400
        try:
            ratewin = expconfig['ratewindowmhz']
        except KeyError:
            ratewin = 30
        try:
            halfbandwidth = expconfig['ampcalfringhalfbandwidth']
        except KeyError:
            halfbandwidth = False
        try:
            dispersivefit = expconfig['ampcalfringdispersivefit']
        except KeyError:
            dispersivefit = False
        try:
            doexhaustive = expconfig['exhaustivefring']
        except KeyError:
            doexhaustive = False
        vlbatasks.fring(inbeamuvdatas[0], snversion, clversion, 
                        expconfig['bandpassfringmins'], inttimesecs, ampcalsrc, expconfig['refant'],
                        False, expconfig['bandpassfringsnr'], sumifs, ampcalmodeldata,
                        sumrrll, bandpassuvrange, zerorates, delaywin, ratewin,
                        doexhaustive, halfbandwidth, dispersivefit)
        print "Still at runlevel " + str(runlevel) + ": editing FRING for ampcal calibrator"
        vlbatasks.snsmo(inbeamuvdatas[0], 'DELA', 20, 0.0, 0.0, 10.0, snversion,
                        expconfig['refant'])
        snoutver = snversion+1
        if not expconfig['skipsnedt']:
            vlbatasks.snedt(inbeamuvdatas[0], snversion+1)
            snoutver = snversion+2
        if os.path.exists(tabledir + 'ampcalfring.sn'):
            os.remove(tabledir + 'ampcalfring.sn')
        vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver, tabledir + \
                             'ampcalfring.sn')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver, 'DELA', 0, 2, 4,
                           tabledir + 'ampcalfring.ps')
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion)
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+1)
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+2)
else:
    print "Skipping FRING and SNSMO/SNEDT of ampcal FRING results"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Copy the FRING SN table around and apply it #################################
if runfromlevel <= runlevel and runtolevel >= runlevel and expconfig['ampcalscan'] > 0:
    print "Runlevel " + str(runlevel) + ": Loading ampcal FRING SN table & calibrating"
    if targetonly and (not os.path.exists(tabledir + 'ampcalfring.sn')):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'ampcalfring.sn', snversion)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion, '2PT', clversion,
                                   expconfig['refant'], [], 'CALP')
    if not calonly:
        vlbatasks.loadtable(gateduvdata, tabledir + 'ampcalfring.sn', snversion)
        vlbatasks.applysntable(gateduvdata, snversion, '2PT', clversion, 
                               expconfig['refant'], [], 'CALP')
        if haveungated:
            vlbatasks.loadtable(ungateduvdata, tabledir + 'ampcalfring.sn', snversion)
            vlbatasks.applysntable(ungateduvdata, snversion, '2PT', clversion, 
                                   expconfig['refant'], [], 'CALP')
else:
    print "Skipping calibration of ampcal FRING results"

runlevel  = runlevel + 1
if expconfig['ampcalscan'] > 0:
    snversion = snversion + 1
    clversion = clversion + 1
xpolsnfilename = tabledir + '/xpolfring.sn'
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run xpoldelaycal ############################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    xpolscan > 0:
    if not os.path.exists(xpolsnfilename) or \
       not interaction.yesno("Do you wish to used saved SN table for Xpol cal?"):
        print "Runlevel " + str(runlevel) + ": Calculating XPOL delays"
        try:
            inttimesecs = expconfig['inttimesecs']
        except KeyError:
            inttimesecs = 2
        try:
            xpolsolintmins = expconfig['xpolsolintmins']
        except KeyError:
            xpolsolintmins = 3
        try:
            delaywin = expconfig['delaywindowns']
        except KeyError:
            delaywin = 400
        try:
            ratewin = expconfig['ratewindowmhz']
        except KeyError:
            ratewin = 30
        if os.path.exists(xpolsnfilename):
            os.remove(xpolsnfilename)
        xpolmodel = AIPSImage("XPOLSRC","CLNMOD",1,1)
        if xpolmodel.exists():
            xpolmodel.zap()
        xpolmodelfile = modeldir + '/' + xpolsource + ".clean.fits"
        if not os.path.exists(xpolmodelfile):
            print xpolmodelfile, "doesn't exit: must have model for xpol calibration"
            sys.exit()
        vlbatasks.fitld_image(xpolmodelfile, xpolmodel)
        print "int time in seconds is ", inttimesecs
        vlbatasks.xpoldelaycal(inbeamuvdatas[0], clversion, expconfig['refant'], 
                               xpolsource, xpolscan, xpolmodel, xpolsolintmins, 
                               inttimesecs, xpolsnfilename, delaywin, ratewin)
else:
    print "Skipping determining xpoldelays"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load and apply xpoldelaycal #################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and xpolscan > 0:
    print "Runlevel " + str(runlevel) + ": Loading XPOL CAL FRING SN table & calibrating"
    if targetonly and (not os.path.exists(xpolsnfilename)):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.loadtable(inbeamuvdatas[i], xpolsnfilename, snversion)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion, '2PT', clversion,
                                   expconfig['refant'])
    if not calonly:
        vlbatasks.loadtable(gateduvdata, xpolsnfilename, snversion)
        vlbatasks.applysntable(gateduvdata, snversion, '2PT', clversion,
                               expconfig['refant'])
        if haveungated:
            vlbatasks.loadtable(ungateduvdata, xpolsnfilename, snversion)
            vlbatasks.applysntable(ungateduvdata, snversion, '2PT', clversion,
                                   expconfig['refant'])
else:
    print "Skipping calibration of xpolcal cross-pol FRING results"

## Run leakagecal on the leakage calibrator and save AN file ###################
runlevel  = runlevel + 1
if xpolscan > 0:
    snversion += 1 
    clversion += 1
leakagefilename = tabledir + "/leakage.an"
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    leakagescan > 0:
    if not os.path.exists(leakagefilename) or \
           not interaction.yesno("Do you wish to used saved AN table for leakage cal?"):
        print "Runlevel " + str(runlevel) + ": Running leakage cal"
        if os.path.exists(leakagefilename):
            os.remove(leakagefilename)
        leakagemodel = AIPSImage("LEAKSRC","CLNMOD",1,1)
        if leakagemodel.exists():
            leakagemodel.zap()
        leakagemodelfile = modeldir + '/' + leakagesource + ".clean.fits"
        leakageoutputfile = directory + '/' + experiment + "_" + leakagesource + "_leakagecal_uv.fits"
        if not os.path.exists(leakagemodelfile):
            print leakagemodelfile, "doesn't exit: must have model for leakage calibration"
            sys.exit()
        vlbatasks.fitld_image(leakagemodelfile, leakagemodel)
        try:
            leakageacalmins = expconfig['leakageampcalmins']
        except KeyError:
            leakageacalmins = 0.3333333
        try:
            leakagepcalmins = expconfig['leakagephscalmins']
        except KeyError:
            leakagepcalmins = 0.02
        hasbptable = False
        vlbatasks.leakagecalc(inbeamuvdatas[0], leakagesource, leakagemodel, leakagefilename, 
                    expconfig['refant'], leakageacalmins, leakagepcalmins, leakagescan, clversion,
                    hasbptable, leakageoutputfile, leakageuvrange, leakageweightit)
else:
    print "Skipping calibration of leakage"

## Load up the leakage-calibrated AN table #####################################
runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
if runfromlevel <= runlevel and runtolevel >= runlevel and leakagescan > 0:
    print "Runlevel " + str(runlevel) + ": Loading leakage-calculated AN table"
    if targetonly and (not os.path.exists(xpolsnfilename)):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.deletetable(inbeamuvdatas[i], "AN", 1)
            vlbatasks.loadtable(inbeamuvdatas[i], leakagefilename, 1)
    if not calonly:
        vlbatasks.deletetable(gateduvdata, "AN", 1)
        vlbatasks.loadtable(gateduvdata, leakagefilename, 1)
        if haveungated:
            vlbatasks.deletetable(ungateduvdata, "AN", 1)
            vlbatasks.loadtable(ungateduvdata, leakagefilename, 1)
else:
    print "Skipping loading the leakage calibrated AN table"

## Remember to set dopol=2 on everything hereon!!!!

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run BPASS ###################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
    expconfig['ampcalscan'] > 0:
    print "Runlevel " + str(runlevel) + ": Generating bandpass corrections"
    if not os.path.exists(tabledir + 'bpass.bp') or not \
       interaction.yesno("Do you wish to used saved BP table for bandpass?"):
        ampcalmodeldata = None
        ampcalmodelfile = modeldir + ampcalsrc + '.clean.fits'
        if os.path.exists(ampcalmodelfile):
            aipscalname = ampcalsrc
            if len(ampcalsrc) > 12:
                aipscalname = ampcalsrc[:12]
            print "Using " + modeltype + " model for " + ampcalsrc
            ampcalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
            if ampcalmodeldata.exists():
                ampcalmodeldata.zap()
            vlbatasks.fitld_image(ampcalmodelfile, ampcalmodeldata)
        vlbatasks.bpass(inbeamuvdatas[0], ampcalsrc, clversion, 
                        expconfig['ampcalscan'], ampcalmodeldata, leakagedopol)
        if os.path.exists(tabledir + 'bpass.bp'):
            os.remove(tabledir + 'bpass.bp')
        vlbatasks.writetable(inbeamuvdatas[0], 'BP', 1, tabledir + 'bpass.bp')
        vlbatasks.deletetable(inbeamuvdatas[0], 'BP', 1)
else:
    print "Skipping bandpass corrections"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load BPASS ##################################################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    expconfig['ampcalscan'] > 0:
    print "Runlevel " + str(runlevel) + ": Loading bandpass corrections"
    if not os.path.exists(tabledir + 'bpass.bp'):
        print "Error - bandpass table " + tabledir + 'bpass.bp does not exist'
        sys.exit(1)
    if not calonly:
        vlbatasks.loadtable(gateduvdata, tabledir + 'bpass.bp', 1)
        if haveungated:
            vlbatasks.loadtable(ungateduvdata, tabledir + 'bpass.bp', 1)
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'bpass.bp', 1)
else:
    print "Skipping loading of bandpass corrections"

bandpassclversion = clversion

runlevel  = runlevel + 1
maxphsreffringmins = -1
maxphsrefcalibapnmins = -1
maxphsrefcalibpnmins = -1
for config in targetconfigs:
    if config['phsreffringmins'] > maxphsreffringmins:
        maxphsreffringmins = config['phsreffringmins']
    if config['phsrefcalibapnmins'] > maxphsrefcalibapnmins:
        maxphsrefcalibapnmins = config['phsrefcalibapnmins']
    try:
        if config['phsrefcalibpnmins'] > maxphsrefcalibpnmins:
            maxphsrefcalibpnmins = config['phsrefcalibpnmins']
    except KeyError:
        pass
donephscalnames = []
doneconfigs = []
for i in range(len(phscalnames)):
    phscal = phscalnames[i]
    if phscal in donephscalnames:
        continue
    donephscalnames.append(phscal)
    doneconfigs.append(targetconfigs[i])
dophscaldump = False
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run FRING (phase reference calibrators) #####################################
if runfromlevel <= runlevel and runtolevel >= runlevel and not targetonly and \
   maxphsreffringmins > 0:
    if not os.path.exists(tabledir + 'phsreffring.sn') or \
           not interaction.yesno("Do you wish to used saved SN table for phsref FRING?"):
        print "Runlevel " + str(runlevel) + ": FRING'ing phase calibrator"
        for phscal, config in zip(donephscalnames, doneconfigs):
            phscalmodeldata = None
            phscalmodelfile = modeldir + phscal + '.clean.fits'
            if os.path.exists(phscalmodelfile):
                aipscalname = phscal
                if len(phscal) > 12:
                    aipscalname = phscal[:12]
                print "Using " + modeltype + " model for " + phscal
                phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                if phscalmodeldata.exists():
                    phscalmodeldata.zap()
                vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
            else:
                print "Currently no " + modeltype + " model for " + phscal
                if modeltype == "preliminary":
                    try:
                        allowphscalpointmodel = expconfig['allowphscalpointmodel']
                        if not allowphscalpointmodel:
                            print "Using a point source, will dump output after FRING"
                            dophscaldump = True
                    except KeyError:
                        print "Using a point source, will dump output after FRING"
                        dophscaldump = True
                else:
                    print "Aborting!"
                    sys.exit()
            doband = True
            if expconfig['ampcalscan'] <= 0:
                doband = False
            try:
                phsrefuvrange = config['phsrefuvrange']
            except KeyError:
                phsrefuvrange = [0, 0]
            try:
                delaywin = expconfig['delaywindowns']
            except KeyError:
                delaywin = 400
            try:
                ratewin = expconfig['ratewindowmhz']
            except KeyError:
                ratewin = 30
            try:
                sumifs = config['phsreffringsumifs']
                if expconfig['ampcalscan'] <= 0 and "bd152" in experiment:
                    print "Overriding sumifs because no amp cal fring was done"
                    sumifs = False
            except KeyError:
                sumifs = False
            try:
                ratemwf = config['phsreffringratemwfhz']
            except KeyError:
                ratemwf = 10
            try:
                delaymwf = config['phsreffringdelaymwfns']
            except KeyError:
                delaymwf = 10
            try:
                sumrrll = config['phsreffringsumpols']
            except KeyError:
                sumrrll = False
            try:
                halfbandwidth = config['phsreffringhalfbandwidth']
            except KeyError:
                halfbandwidth = False
            try:
                dispersivefit = config['phsreffringdispersivefit']
            except KeyError:
                dispersivefit = False
            try:
                inttimesecs = expconfig['inttimesecs']
            except KeyError:
                inttimesecs = 2
            try:
                doexhaustive = expconfig['exhaustivefring']
            except KeyError:
                doexhaustive = False
            print "Delay and rate windows: ", delaywin, ratewin
            print "UV range: ", phsrefuvrange
            vlbatasks.fring(inbeamuvdatas[0], snversion, clversion,  
                            config['phsreffringmins'], inttimesecs, phscal, 
                            expconfig['refant'], doband, 
                            config['phsreffringsnr'], sumifs, phscalmodeldata,
                            sumrrll, phsrefuvrange,False,delaywin,ratewin,
                            doexhaustive, halfbandwidth, dispersivefit, leakagedopol)
        try:
            dosnsmo = config['dosnsmo']
        except KeyError:
            dosnsmo = True
        try:
            smoothedrefant = expconfig['smoothedrefant']
        except KeyError:
            smoothedrefant = expconfig['refant']
        try:
            mwfminutes = config['mwfminutes']
        except:
            mwfminutes = 20
        snoutver = snversion
        if dosnsmo:
            #vlbatasks.snsmo(inbeamuvdatas[0], 'DELA', 20, 0.0, 0.0, 10.0, snversion,
            #                smoothedrefant, True)
            if sumifs:
                doratesmoothing = False
                ratesmoothingmins = 0.0
            else:
                doratesmoothing = True
                ratesmoothingmins = 3.0
            vlbatasks.snmwfclip(inbeamuvdatas[0], mwfminutes, 0.0, 0.0, delaymwf, ratemwf, snversion,
                                smoothedrefant, doratesmoothing, ratesmoothingmins)
            snoutver = snversion+1
        if not expconfig['skipsnedt']:
            vlbatasks.snedt(inbeamuvdatas[0], snoutver)
            snoutver += 1
        if os.path.exists(tabledir + 'phsreffring.sn'):
            os.remove(tabledir + 'phsreffring.sn')
        vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver, tabledir + \
                             'phsreffring.sn')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver, 'DELA', 0, 2, 4,
                           tabledir + 'phsreffring.delay.ps')
        vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver, 'RATE', 0, 2, 4,
                           tabledir + 'phsreffring.rate.ps')
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion)
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+1)
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+2)
else:
    print "Skipping FRING and SNSMO/SNEDT of FRING results"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Copy the phsref FRING SN table around and apply it ##########################
if runfromlevel <= runlevel and runtolevel >= runlevel and maxphsreffringmins > 0:
    print "Runlevel " + str(runlevel) + ": Loading phsref FRING SN table & calibrating"
    if targetonly and (not os.path.exists(tabledir + 'phsreffring.sn')):
        print "For target-only, the SN file must exist already - aborting!"
        sys.exit(1)
    try:
        interptype = expconfig['phsreffringinterp']
    except KeyError:
        interptype = 'AMBG'
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'phsreffring.sn', snversion)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion, interptype, clversion, 
                                   expconfig['refant'])
    if not calonly:
        vlbatasks.loadtable(gateduvdata, tabledir + 'phsreffring.sn', snversion)
        vlbatasks.applysntable(gateduvdata, snversion, interptype, clversion, 
                               expconfig['refant'])
        if haveungated:
            vlbatasks.loadtable(ungateduvdata, tabledir + 'phsreffring.sn', snversion)
            vlbatasks.applysntable(ungateduvdata, snversion, interptype, clversion, 
                                   expconfig['refant'])
else:
    print "Skipping calibration of FRING results"

runlevel  = runlevel + 1
if maxphsreffringmins > 0:
    snversion = snversion + 1
    clversion = clversion + 1
print donephscalnames
#donephscalnames = []
#doneconfigs = []
#for i in range(len(phscalnames)):
#    phscal = phscalnames[i]
#    try:
#        if phscal in donephscalnames or (targetconfigs[i]['phsrefcalibpnmins'] < 0 and not "IBC" in phscal):
#            continue
#        donephscalnames.append(phscal)
#        doneconfigs.append(targetconfigs[i])
#    except KeyError:
#        continue
if dophscaldump: # Need to dump out the phs cal sources so we can make models of them
    for phscal in donephscalnames:
        for i in range(20): #Clear any old CALIB split catalog entries
            phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, i)
            if phscal_uv_data.exists():
                phscal_uv_data.zap()
        phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, 1)
        rawuvoutputfile = directory + experiment.upper() + '_' + \
                                   phscal + ".formodeling.uv.fits"
        doband = False
        domulti = False
        if expconfig['ampcalscan'] > 0:
            doband = True
        combineifs = False
        beginif = -1
        endif = -1
        vlbatasks.splittoseq(inbeamuvdatas[0], clversion, 'CALIB', phscal, 1, domulti,
                             doband, beginif, endif, combineifs, leakagedopol)
        vlbatasks.writedata(phscal_uv_data, rawuvoutputfile, True)
    print "UV datasets of the phase reference sources have been written out to model"
    sys.exit()

printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run phase CALIB on the phase reference sources ##############################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxphsrefcalibpnmins > 0:
    haveall = True
    for phscal in donephscalnames:
        calibtablepath = tabledir + phscal + '.calibpn.sn'
        if not os.path.exists(calibtablepath):
            haveall = False
    if not targetonly and (not haveall or \
           not interaction.yesno("Do you wish to use saved SN table for phsref phase CALIB?")):
        print "Runlevel " + str(runlevel) + ": Running CALIB on phsref sources"
        for phscal, config in zip(donephscalnames, doneconfigs):
            for i in range(20): #Clear any old CALIB split catalog entries
                phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, i)
                if phscal_uv_data.exists():
                    phscal_uv_data.zap()
            phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, 1)
            doband = False
            domulti = False
            if expconfig['ampcalscan'] > 0:
                doband = True
            combineifs = False
            beginif = -1
            endif = -1
            vlbatasks.splittoseq(inbeamuvdatas[0], clversion, 'CALIB', phscal, 1, domulti,
                                 doband, beginif, endif, combineifs, leakagedopol)
            phscal_uv_data.table('NX', 1).zap()
            phscalmodeldata = None
            phscalmodelfile = modeldir + phscal + '.clean.fits'
            try:
                phsrefuvrange = config['phsrefuvrange']
            except KeyError:
                phsrefuvrange = [0, 0]
            try:
                phsrefweightit = expconfig['phsrefcalibpnweightit']
            except KeyError:
                phsrefweightit = 0
            if os.path.exists(phscalmodelfile):
                print "Using " + modeltype + " model for " + phscal
                aipscalname = phscal
                if len(phscal) > 12:
                    aipscalname = phscal[:12]
                phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                if phscalmodeldata.exists():
                    phscalmodeldata.zap()
                vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
            else:
                print "Currently no " + modeltype + " model for " + phscal + '; aborting!'
                sys.exit()
            calibsoltype = config['phsrefcalibpntype']
            calibtablepath = tabledir + phscal + '.calibpn.sn'
            try:
                flagwheremodelbelow = expconfig['phsrefminmodelflux'] #Jy
            except KeyError:
                flagwheremodelbelow = -1
            vlbatasks.singlesource_calib(phscal_uv_data, phscalmodeldata,
                                         1, expconfig['refant'], False, 
                                         config['phsrefcalibpnmins'],
                                         False, calibsoltype, 
                                         config['phsrefcalibpnsnr'], False,
                                         phsrefuvrange, phsrefweightit,
                                         flagwheremodelbelow)
            snoutver = 1
            if not expconfig['skipsnedt']:
                snoutver = 2
                vlbatasks.snedt(phscal_uv_data, 1)
            else:
                snoutver = 2
                vlbatasks.snsmo(phscal_uv_data, 'AMP', 300, 0.20, 0.0, 0.0, 1, expconfig['refant'])
            if os.path.exists(calibtablepath):
                os.remove(calibtablepath)
            vlbatasks.writetable(phscal_uv_data, 'SN', snoutver, calibtablepath)
            vlbatasks.plottops(phscal_uv_data, 'SN', snoutver, 'PHAS', 0, 2, 4,
                               tabledir + phscal + '.calibpn.phs.ps')
            phscal_uv_data.zap()
else:
    print "Skipping phase cal CALIB on the phase reference sources"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load all the phase CALIB solutions #########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxphsrefcalibpnmins > 0:
    print "Runlevel " + str(runlevel) + ": Loading phs ref PN CALIB " + \
          "solutions and applying"
    sncount = 0
    for phscal in phscalnames:
        if not targetonly:
            for i in range(numinbeams):
                vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+sncount)
        if not calonly:
            vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
            if haveungated:
                vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
        sncount += 1
    sncount = 0
    for phscal in donephscalnames:
        calibtablepath = tabledir + phscal + '.calibpn.sn'
        if targetonly and (not os.path.exists(calibtablepath)):
            print "For target-only, the SN file must exist already - aborting!"
            sys.exit(1)
        if not targetonly:
            for i in range(numinbeams):
                vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath,
                                    snversion+sncount)
        if not calonly:
            vlbatasks.loadtable(gateduvdata, calibtablepath,
                                snversion+sncount)
            if haveungated:
                vlbatasks.loadtable(ungateduvdata, calibtablepath,
                                    snversion+sncount)
        sncount += 1
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+sncount)
            vlbatasks.mergesntables(inbeamuvdatas[i], snversion, sncount, 
                                    expconfig['refant'])
    if not calonly:
        vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
        vlbatasks.mergesntables(gateduvdata, snversion, sncount, 
                                expconfig['refant'])
        if haveungated:
            vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
            vlbatasks.mergesntables(ungateduvdata, snversion, sncount, 
                                    expconfig['refant'])
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.deletetable(inbeamuvdatas[i], 'CL', clversion+1)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion+sncount, 'SELF', 
                                   clversion, expconfig['refant'])
    if not calonly:
        vlbatasks.deletetable(gateduvdata, 'CL', clversion+1)
        vlbatasks.applysntable(gateduvdata, snversion+sncount, 'SELF', 
                               clversion, expconfig['refant'])
        if haveungated:
            vlbatasks.deletetable(ungateduvdata, 'CL', clversion+1)
            vlbatasks.applysntable(ungateduvdata, snversion+sncount, 'SELF', 
                                   clversion, expconfig['refant'])
else:
    print "Skipping loading/application of phase reference source phase CALIB solutions"
    sncount = 0
    if maxphsrefcalibpnmins > 0:
        sncount += len(donephscalnames)

runlevel += 1
if maxphsrefcalibpnmins > 0:
    snversion = snversion + sncount + 1
    clversion = clversion + 1
donephscalnames = []
doneconfigs = []
for i in range(len(phscalnames)):
    phscal = phscalnames[i]
    try:
        if phscal in donephscalnames or targetconfigs[i]['phsrefcalibapnmins'] < 0:
            continue
        donephscalnames.append(phscal)
        doneconfigs.append(targetconfigs[i])
    except KeyError:
        continue
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Run amp CALIB on the phase reference sources ################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxphsrefcalibapnmins > 0:
    haveall = True
    for phscal in donephscalnames:
        calibtablepath = tabledir + phscal + '.calibapn.sn'
        if not os.path.exists(calibtablepath):
            haveall = False
    if not targetonly and (not haveall or \
           not interaction.yesno("Do you wish to used saved SN table for phsref amp CALIB?")):
        print "Runlevel " + str(runlevel) + ": Running CALIB on phsref sources"
        for phscal, config in zip(donephscalnames, doneconfigs):
            for i in range(20): #Clear any old CALIB split catalog entries
                phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, i)
                if phscal_uv_data.exists():
                    phscal_uv_data.zap()
            phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, 1)
            doband = False
            domulti = False
            if expconfig['ampcalscan'] > 0:
                doband = True
            combineifs = False
            beginif = -1
            endif = -1
            vlbatasks.splittoseq(inbeamuvdatas[0], clversion, 'CALIB', phscal, 1, domulti,
                                 doband, beginif, endif, combineifs, leakagedopol)
            phscal_uv_data.table('NX', 1).zap()
            phscalmodeldata = None
            phscalmodelfile = modeldir + phscal + '.clean.fits'
            try:
                phsrefuvrange = config['phsrefuvrange']
            except KeyError:
                phsrefuvrange = [0, 0]
            try:
                phsrefweightit = expconfig['phsrefcalibapnweightit']
            except KeyError:
                phsrefweightit = 0
            try:
                flagwheremodelbelow = expconfig['phsrefminmodelflux'] #Jy
            except KeyError:
                flagwheremodelbelow = -1
            if os.path.exists(phscalmodelfile):
                print "Using " + modeltype + " model for " + phscal
                aipscalname = phscal
                if len(phscal) > 12:
                    aipscalname = phscal[:12]
                phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                if phscalmodeldata.exists():
                    phscalmodeldata.zap()
                vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
            else:
                print "Currently no " + modeltype + " model for " + phscal + '; aborting!'
                sys.exit()
            calibsoltype = config['phsrefcalibapntype']
            calibtablepath = tabledir + phscal + '.calibapn.sn'
            try:
                normalise = config['phsrefcalibapnnorm']
            except KeyError:
                normalise = False
            phsrefcalibapnclip = -1
            try:
                phsrefcalibapnclip = config['phsrefcalibapnclip']
            except KeyError:
                phsrefcalibapnclip = 0.2
            try: 
                tmpphsrefcalibapnclip = expconfig['phsrefcalibapnclip']
                if phsrefcalibapnclip > 0:
                   print "Overriding default value of phsrefcalibapnclip (", phsrefcalibapnclip, ") for this source with", tmpphsrefcalibapnclip
                phsrefcalibapnclip = tmpphsrefcalibapnclip
            except KeyError:
                pass
            vlbatasks.singlesource_calib(phscal_uv_data, phscalmodeldata,
                                         1, expconfig['refant'], True,
                                         config['phsrefcalibapnmins'],
                                         False, calibsoltype,
                                         config['phsrefcalibapnsnr'], False,
                                         phsrefuvrange, phsrefweightit,
                                         flagwheremodelbelow, normalise)
            snoutver = 1
            if not expconfig['skipsnedt']:
                snoutver = 2
                vlbatasks.snedt(phscal_uv_data, 1)
            else:
                snoutver = 2
                vlbatasks.snsmo(phscal_uv_data, 'AMP', 300, phsrefcalibapnclip, 0.0, 0.0, 1, expconfig['refant'])
            if os.path.exists(calibtablepath):
                os.remove(calibtablepath)
            vlbatasks.writetable(phscal_uv_data, 'SN', snoutver, calibtablepath)
            vlbatasks.plottops(phscal_uv_data, 'SN', snoutver, 'AMP', 0, 2, 4,
                               tabledir + phscal + '.calibapn.amp.ps')
            vlbatasks.plottops(phscal_uv_data, 'SN', snoutver, 'PHAS', 0, 2, 4,
                               tabledir + phscal + '.calibapn.phs.ps')
            phscal_uv_data.zap()
else:
    print "Skipping amp cal CALIB on the phase reference sources"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load all the amp CALIB solutions ###########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxphsrefcalibapnmins > 0:
    print "Runlevel " + str(runlevel) + ": Loading phs ref APN CALIB " + \
          "solutions and applying"
    sncount = 0
    for phscal in phscalnames:
        if not targetonly:
            for i in range(numinbeams):
                vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+sncount)
        if not calonly:
            vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
            if haveungated:
                vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
        sncount += 1
    sncount = 0
    for phscal in donephscalnames:
        calibtablepath = tabledir + phscal + '.calibapn.sn'
        if targetonly and (not os.path.exists(calibtablepath)):
            print "For target-only, the SN file must exist already - aborting!"
            sys.exit(1)
        if not targetonly:
            for i in range(numinbeams):
                vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath,
                                    snversion+sncount)
        if not calonly:
            vlbatasks.loadtable(gateduvdata, calibtablepath,
                                snversion+sncount)
            if haveungated:
                vlbatasks.loadtable(ungateduvdata, calibtablepath,
                                    snversion+sncount)
        sncount += 1
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+sncount)
            vlbatasks.mergesntables(inbeamuvdatas[i], snversion, sncount,
                                    expconfig['refant'])
    if not calonly:
        vlbatasks.deletetable(gateduvdata, 'SN', snversion+sncount)
        vlbatasks.mergesntables(gateduvdata, snversion, sncount,
                                expconfig['refant'])
        if haveungated:
            vlbatasks.deletetable(ungateduvdata, 'SN', snversion+sncount)
            vlbatasks.mergesntables(ungateduvdata, snversion, sncount,
                                    expconfig['refant'])
    if not targetonly:
        for i in range(numinbeams):
            vlbatasks.deletetable(inbeamuvdatas[i], 'CL', clversion+1)
            vlbatasks.applysntable(inbeamuvdatas[i], snversion+sncount, 'SELF',
                                   clversion, expconfig['refant'])
    if not calonly:
        vlbatasks.deletetable(gateduvdata, 'CL', clversion+1)
        vlbatasks.applysntable(gateduvdata, snversion+sncount, 'SELF',
                               clversion, expconfig['refant'])
        if haveungated:
            vlbatasks.deletetable(ungateduvdata, 'CL', clversion+1)
            vlbatasks.applysntable(ungateduvdata, snversion+sncount, 'SELF',
                                   clversion, expconfig['refant'])
else:
    print "Skipping loading/application of phase reference source amp CALIB solutions"
    sncount = 0
    if maxphsrefcalibapnmins > 0:
        sncount += len(donephscalnames)
        
doneinbeams = []
inbeamfilenums = []
secondaryinbeams = []
secondaryfilenums = []
numifs = 4
maxinbeamcalibp1mins = -1 # Phase-only, primary inbeam, summed IFs
maxinbeamcalibpnmins = -1 # Phase-only, primary inbeam, separate IFs
maxinbeamcalibsp1mins = -1 # Phase-only, secondary inbeam (rarely used)
maxinbeamcalibap1mins = -1 # Amplitude and phase, primary inbeam, combined IFs
maxinbeamcalibapnmins = -1 # Amplitude and phase, primary inbeam, separate IFs
if inbeamuvdatas[0].exists():
    numifs = vlbatasks.getNumIFs(inbeamuvdatas[0])
beginif = 1
endif = numifs
for config in targetconfigs:
    if config['inbeamcalibp1mins'] > maxinbeamcalibp1mins:
        maxinbeamcalibp1mins = config['inbeamcalibp1mins']
    if config['inbeamcalibap1mins'] > maxinbeamcalibap1mins:
        maxinbeamcalibap1mins = config['inbeamcalibap1mins']
    try:
        if config['inbeamcalibsp1mins'] > maxinbeamcalibsp1mins:
            maxinbeamcalibsp1mins = config['inbeamcalibsp1mins']
    except KeyError:
        print "No secondary inbeam calibrator..."
    try:
        if config['inbeamcalibpnmins'] > maxinbeamcalibpnmins:
            maxinbeamcalibpnmins = config['inbeamcalibpnmins']
    except KeyError:
        print "No separate IF inbeam selfcal..."
    try:
        if config['inbeamcalibapnmins'] > maxinbeamcalibapnmins:
            maxinbeamcalibapnmins = config['inbeamcalibapnmins']
    except KeyError:
        print "No separate IF inbeam amp selfcal..."
print "maxinbeamcalibp1mins", maxinbeamcalibp1mins
print "maxinbeamcalibpnmins", maxinbeamcalibpnmins
print "maxinbeamcalibap1mins", maxinbeamcalibap1mins
print "maxinbeamcalibapnmins", maxinbeamcalibapnmins
print "maxinbeamcalibsp1mins", maxinbeamcalibsp1mins
for i in range(numtargets):
    primaryinbeams = targetconfigs[i]['primaryinbeam'].split(',')
    try:
        secondaryinbeam = targetconfigs[i]["secondaryinbeam"]
        secondaryinbeam = secondaryinbeam.split(',')
    except KeyError:
        secondaryinbeam = "XXXXX"
    for j in range(len(inbeamnames[i])):
        for inbeamsrc in secondaryinbeam:
            if not secondaryinbeam == None and (inbeamsrc.strip() == inbeamnames[i][j].strip()):
                secondaryinbeams.append(inbeamsrc.strip())
                secondaryfilenums.append(j)
        for primaryinbeam in primaryinbeams:
            if primaryinbeam.strip() == inbeamnames[i][j].strip():
                if not primaryinbeam.strip() in doneinbeams:
                    doneinbeams.append(primaryinbeam.strip())
                    inbeamfilenums.append(j)
            if primaryinbeam.strip() == targetnames[i]:
                if not targetnames[i] in doneinbeams:
                    doneinbeams.append(targetnames[i])
                    inbeamfilenums.append(-1)
for primaryinbeam in primaryinbeams:
    if not primaryinbeam.strip() in doneinbeams:
        print "Didn't find primary inbeam " + primaryinbeam.strip() + " amongst the data! Check the inbeam name in your config file!"
        print doneinbeams
        sys.exit()
 
runlevel  = runlevel + 1
if maxphsrefcalibapnmins > 0:
    snversion = snversion + sncount + 1
    clversion = clversion + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Generate a raw (no phase selfcal) inbeam dataset if requested ##############
if runfromlevel <= runlevel and runtolevel >= runlevel and \
   expconfig['writerawinbeam']:
    print "Runlevel " + str(runlevel) + ": Writing raw inbeam outputs"
    if not targetonly: # Make a "raw" inbeam fileset if required
        for i in range(numtargets):
            count = 0
            subtractif = 0
            if targetconfigs[i]['skiplastif']:
                subtractif = 1
            for inbeamsrc in inbeamnames[i]:
                aipssrcname = inbeamsrc
                if len(inbeamsrc) > 12:
                    aipssrcname = inbeamsrc[:12]
                rawuvoutputfile =  directory + experiment.upper() + '_' + \
                                   inbeamsrc + ".formodeling.uv.fits"
                splitdata = AIPSUVData(aipssrcname, 'NOIB', 1, 1)
                if splitdata.exists():
                    splitdata.zap()
                uvdata = inbeamuvdatas[count]
                doband = False
                domulti = False
                if expconfig['ampcalscan'] > 0:
                    doband = True
                combineifs = False
                vlbatasks.splittoseq(inbeamuvdatas[0], clversion, 'CALIB', phscal, 1,
                                     domulti, doband, beginif, endif-subtractif, 
                                     combineifs, leakagedopol)
                vlbatasks.writedata(splitdata, rawuvoutputfile, True)
                count += 1
else:
    print "Skipping dump of raw inbeam data (pre-selfcal)"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Do a combined IF (and pol) phase selfcal on the inbeams if requested #################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibp1mins > 0:
    print "Runlevel " + str(runlevel) + ": Doing phase-only inbeam selfcal (combined IFs)"
    tocalnames, tocalindices = inbeamselfcal(doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata, 
                               expconfig, targetconfigs, modeldir, modeltype, targetonly, 
                               calonly, beginif, endif, False, False, True, clversion, targetnames,
                               leakagedopol)
else:
    print "Skipping inbeam phase-only selfcal (combined IFs)"
    tocalnames = []
    tocalindices = []
    for i in range(numtargets):
        config = targetconfigs[i]
        if config['separateifmodel'] or \
           len(config['primaryinbeam'].split(',')) > 1:
            tocalnames.append('CONCAT' + str(i))
        else:
            tocalnames.append(targetconfigs[i]['primaryinbeam'])
        tocalindices.append(i)

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion, inbeamuvdatas[0])
## Load all the inbeam CALIB solutions ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibp1mins > 0:
    print "Runlevel " + str(runlevel) + ": Applying inbeam CALIB p1 sols"
    sncount = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, 
                               targetconfigs, targetonly, calonly, False, False, True,
                               clversion, snversion, inbeamnames, targetnames)
    #sncount is used to point at SN table in post-phscal stage
else:
    print "Skipping application of inbeam phase-only selfcal (combined IFs)"
    if maxinbeamcalibp1mins > 0:
        sncount = len(tocalnames) + 1
    else:
        sncount = 0

runlevel = runlevel + 1
targetcl = 0 #used for pointing at new CL table in post-phscal stage.
if maxinbeamcalibp1mins > 0:
    snversion = snversion + sncount
    targetcl += 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Do a separate IF phase selfcal on the inbeams if requested #################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibpnmins > 0:
    print "Runlevel " + str(runlevel) + ": Doing phase-only inbeam selfcal (separate IFs)"
    tocalnames, tocalindices = inbeamselfcal(doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata,
                               expconfig, targetconfigs, modeldir, modeltype, targetonly,
                               calonly, beginif, endif, False, False, False, clversion+targetcl, 
                               targetnames, leakagedopol)
else:
    print "Skipping inbeam phase-only selfcal (separate IFs)"
    tocalnames = []
    tocalindices = []
    for i in range(numtargets):
        config = targetconfigs[i]
        if config['separateifmodel'] or \
           len(config['primaryinbeam'].split(',')) > 1:
            tocalnames.append('CONCAT' + str(i))
        else:
            tocalnames.append(targetconfigs[i]['primaryinbeam'])
        tocalindices.append(i)

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Do dual-phscal calibration if requested: 1). correct INBEAM.icalib.p1.sn ###################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
   int(dualphscal_setup[0].strip()) > 0:
    print("Adopt dual-phscal mode now...")
    inbeamselfcal_phase_time_folder = directory+'/inbeamselfcal_phase_time_evolution'
    if not os.path.exists(inbeamselfcal_phase_time_folder):
        os.system('mkdir %s' % inbeamselfcal_phase_time_folder)
    
    inbeamselfcalp1sntable = tabledir + config['primaryinbeam'].split(',')[0].strip() + '.icalib.p1.sn'
    for i in range(20):
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+i)
        ## reverse the last applyinbeamcalib only on target data
        if not calonly:
            vlbatasks.deletetable(gateduvdata, 'CL', clversion+1+i) 
            if haveungated:
                vlbatasks.deletetable(ungateduvdata, 'CL', clversion+1+i)
    vlbatasks.loadtable(inbeamuvdatas[0], inbeamselfcalp1sntable, snversion)
    dualphscalp1 = vlbatasks.calibrate_target_phase_with_two_colinear_phscals(inbeamuvdatas[0])
    dualphscalp1.compile_into_table()
    #originalinbeamselfcalp1sntable = dualphscal.copy_inbeamselfcal_sntable(inbeamselfcalp1sntable) 
    
    final_inbeamselfcal_phase_edit = inbeamselfcal_phase_time_folder + '/.corrected_phases_inbeam_selfcal.final'
    if not os.path.exists(final_inbeamselfcal_phase_edit):
        print("the final saved_inbeamselfcal_phase_edit not found, now heading to interactive phase correction. When you finalize the edit, make a copy of the output file, rename it to .corrected_phases_inbeam_selfcal.final and rerun the pipeline.")
        dualphscalp1.interactively_solve_phase_ambiguity(inbeamselfcal_phase_time_folder)
        sys.exit()
    else:
        phase_correction_factor = float(dualphscal_setup[1].strip())
        dualphscalp1.load_final_inbeamselfcal_phase_edit_and_prepare_for_edit_in_AIPS(final_inbeamselfcal_phase_edit,
                                                                                    phase_correction_factor)
        dualphscaloutputsn = inbeamselfcalp1sntable.replace('.sn', '.dualphscal.sn')
        dualphscalp1.edit_AIPS_sntable_and_write_out(snversion, dualphscaloutputsn)
        #apply the corrected p1.sn only to the target
        junk = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, 
                               targetconfigs, True, calonly, False, False, True,
                               clversion, snversion, inbeamnames, targetnames)
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion)

runlevel += 1

## Do dual-phscal calibration if requested: 2). correct INBEAM.icalib.pn.sn ###################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
   int(dualphscal_setup[0].strip()) > 0 and maxinbeamcalibpnmins > 0:
    inbeamselfcalpnsntable = tabledir + config['primaryinbeam'].split(',')[0].strip() + '.icalib.pn.sn'
    for i in range(20):
        vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion+i)
    vlbatasks.loadtable(inbeamuvdatas[0], inbeamselfcalpnsntable, snversion)
    dualphscalpn = vlbatasks.calibrate_target_phase_with_two_colinear_phscals(inbeamuvdatas[0])
    dualphscalpn.read_inbeamselfcalpn_solutions()
    dualphscaloutputsn = inbeamselfcalpnsntable.replace('.sn', '.dualphscal.sn')
    dualphscalpn.edit_inbeamselfcalpn_in_AIPS_and_write_out(phase_correction_factor, snversion, dualphscaloutputsn)
    vlbatasks.deletetable(inbeamuvdatas[0], 'SN', snversion)

runlevel += 1
## Load all the inbeam CALIB solutions ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibpnmins > 0:
    print "Runlevel " + str(runlevel) + ": Applying inbeam CALIB pn sols"
    sncount = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                               targetconfigs, targetonly, calonly, False, False, False,
                               clversion+targetcl, snversion, inbeamnames, targetnames)
else:
    print "Skipping application of inbeam phase-only selfcal (separate IFs)"
    if maxinbeamcalibpnmins > 0:
        sncount = len(tocalnames) + 1
    else:
        sncount = 0

runlevel = runlevel + 1
if maxinbeamcalibpnmins > 0:
    snversion = snversion + sncount
    targetcl += 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Do a combined IF amp + phase selfcal on the inbeams if requested ###########
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibap1mins > 0:
    print "Runlevel " + str(runlevel) + ": Doing amp+phase inbeam selfcal"
    tocalnames, tocalindices = inbeamselfcal(doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata,
                               expconfig, targetconfigs, modeldir, modeltype, targetonly,
                               calonly, beginif, endif, True, False, True, clversion+targetcl, 
                               targetnames, leakagedopol)
else:
    print "Skipping inbeam amp+phase selfcal (combined IFs)"
    tocalnames = []
    tocalindices = []
    for i in range(numtargets):
        config = targetconfigs[i]
        if config['inbeamcalibap1mins'] > 0:
            if config['separateifmodel'] or \
                len(config['primaryinbeam'].split(',')) > 1:
                tocalnames.append('CONCAT' + str(i))
            else:
                tocalnames.append(targetconfigs[i]['primaryinbeam'])
            tocalindices.append(i)

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Load all the inbeam CALIB solutions ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibap1mins > 0:
    print "Runlevel " + str(runlevel) + ": Applying inbeam CALIB ap1 sols"
    sncount = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                               targetconfigs, targetonly, calonly, True, False, True,
                               clversion+targetcl, snversion, inbeamnames, targetnames)
else:
    print "Skipping application of inbeam amp+phase selfcal (combined IFs)"
    if maxinbeamcalibap1mins > 0:
        sncount = len(tocalnames) + 1
    else:
        sncount = 0

runlevel = runlevel + 1
if maxinbeamcalibap1mins > 0:
    snversion = snversion + sncount
    targetcl += 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Do a separate IF amp + phase selfcal on the inbeams if requested ###########
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibapnmins > 0:
    print "Runlevel " + str(runlevel) + ": Doing amp+phase inbeam selfcal"
    tocalnames, tocalindices = inbeamselfcal(doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata,
                               expconfig, targetconfigs, modeldir, modeltype, targetonly,
                               calonly, beginif, endif, True, False, False, clversion+targetcl,
                               targetnames, leakagedopol)
else:
    print "Skipping inbeam amp+phase selfcal (separate IFs)"
    tocalnames = []
    tocalindices = []
    for i in range(numtargets):
        config = targetconfigs[i]
        try:
            if config['inbeamcalibapnmins'] > 0:
                if config['separateifmodel'] or \
                    len(config['primaryinbeam'].split(',')) > 1:
                    tocalnames.append('CONCAT' + str(i))
                else:
                    tocalnames.append(targetconfigs[i]['primaryinbeam'])
                tocalindices.append(i)
        except KeyError:
            pass

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Load all the inbeam CALIB solutions ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibapnmins > 0:
    print "Runlevel " + str(runlevel) + ": Applying inbeam CALIB apn sols"
    sncount = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                               targetconfigs, targetonly, calonly, True, False, False,
                               clversion+targetcl, snversion, inbeamnames, targetnames)
else:
    print "Skipping application of inbeam amp+phase selfcal (separate IFs)"
    if maxinbeamcalibapnmins > 0:
        sncount = len(tocalnames) + 1
    else:
        sncount = 0

runlevel = runlevel + 1
if maxinbeamcalibapnmins > 0:
    snversion = snversion + sncount
    targetcl += 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Do a secondary phase selfcal on the inbeam(s) if requested #################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibsp1mins > 0:
    print "Runlevel " + str(runlevel) + ": Doing phase-only inbeam selfcal on secondary inbeam"
    tocalnames, tocalindices = inbeamselfcal(secondaryinbeams, secondaryfilenums, 
                               inbeamuvdatas, gateduvdata, expconfig, targetconfigs, 
                               modeldir, modeltype, targetonly,
                               calonly, beginif, endif, False, True, True, 
                               clversion+targetcl, targetnames, leakagedopol)
else:
    print "Skipping secondary inbeam phase-only selfcal (combined IFs)"
    tocalnames = []
    tocalindices = []
    for i in range(numtargets):
        config = targetconfigs[i]
        try:
            if config['inbeamcalibsp1mins'] > 0:
                tocalnames.append(targetconfigs[i]['secondaryinbeam'])
                tocalindices.append(i)
        except KeyError:
            continue

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Load all the inbeam CALIB solutions ########################################
if runfromlevel <= runlevel and runtolevel >= runlevel and \
    maxinbeamcalibsp1mins > 0:
    print "Runlevel " + str(runlevel) + ": Applying secondary inbeam CALIB sp1 sols"
    sncount = applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                               targetconfigs, targetonly, calonly, False, True, True,
                               clversion+targetcl, snversion, inbeamnames, targetnames)
else:
    print "Skipping application of secondary inbeam phase-only selfcal (combined IFs)"
    if maxinbeamcalibsp1mins > 0:
        sncount = len(tocalnames) + 1
    else:
        sncount = 0

runlevel = runlevel + 1
if maxinbeamcalibsp1mins > 0:
    snversion = snversion + sncount
    targetcl += 1
scinttablepaths = []
maxscintcorrmins = -1
for i in range(numtargets):
    scinttablepaths.append(tabledir + targetnames[i] + '.scintcorrect.sn')
    if targetconfigs[i]['scintcorrmins'] > maxscintcorrmins:
        maxscintcorrmins = targetconfigs[i]['scintcorrmins']
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Calculate the scintillation correction #####################################
if runfromlevel <= runlevel and runtolevel >= runlevel and maxscintcorrmins > 0:
    print "Doing scintillation correction"
    for i in range(numtargets):
        pspath = tabledir + targetnames[i]+ '.scintcorr.ps'
        try:
            noscintcorr = expconfig['noscintcorr']
        except KeyError:
            noscintcorr = False
        if targetconfigs[i]['scintcorrmins'] > 0 and not noscintcorr:
            splituvdata = AIPSUVData(targetnames[i], "SCINT", 1, 1)
            if splituvdata.exists():
                splituvdata.zap()
            doband = False
            domulti = False
            if expconfig['ampcalscan'] > 0:
                doband = True
            combineifs = False
            beginif = -1
            endif = -1
            vlbatasks.splittoseq(gateduvdata, clversion+targetcl, 'SCINT', 
                                 targetnames[i], 1, domulti, doband, 
                                 beginif, endif, combineifs, leakagedopol)
            vlbatasks.wizCorrectScint(gateduvdata, 1, snversion, splituvdata,
                                      targetconfigs[i]['scintcorrmins'], 
                                      scinttablepaths[i])
            vlbatasks.plottops(gateduvdata, 'SN', snversion, 'AMP', numifs, 1, 
                               numifs, pspath)
            vlbatasks.deletetable(gateduvdata, 'SN', snversion)
        else:
            os.system("rm -f " + scinttablepaths[i])
            os.system("rm -f " + pspath)

else:
    print "Skipping scintillation correction calculation"

## Scintillation correction is applied later at the stage of the final split

runlevel += 1
for i in range(20): #Clear any old FINAL split catalog entries
    for j in range(numtargets):
        for inbeamsrc in inbeamnames[j]:
            aipssrcname = inbeamsrc
            if len(inbeamsrc) > 12:
                aipssrcname = inbeamsrc[:12]
            splitdata = AIPSUVData(aipssrcname, 'FINAL', 1, 1)
            if splitdata.exists():
                splitdata.zap()
        splitdata = AIPSUVData(phscalnames[j], 'FINAL', 1, 1)
        if splitdata.exists():
            splitdata.zap()
        splitdata = AIPSUVData(targetnames[j], 'UFINL', 1, 1)
        if splitdata.exists():
            splitdata.zap()
        splitdata = AIPSUVData(targetnames[j], 'GFINL', 1, 1)
        if splitdata.exists():
            splitdata.zap()
phscaluvfiles = []
inbeamuvfiles = []
gateduvfiles = []
dividedinbeamuvfiles = []
ungateduvfiles = []
ungatedpresent = []
dividesourcelist = []
ibshiftdivphscaluvfiles = []
divinbeampreselfcaluvfiles = []
try:
    if not config['dividesources'] == None:
        dividesourcelist = [x.strip(' ') for x in config['dividesources'].split(',')]
except KeyError:
    pass
for phscalsrc in phscalnames:
    if len(phscalsrc) > 12:
        phscalsrc = phscalsrc[:12]
    phscaluvfiles.append(directory + '/' + experiment + "_" + phscalsrc + \
                         "_pipeline_uv.fits")
    ibshiftdivphscaluvfiles.append(directory + '/' + experiment + '_' + phscalsrc + '_ibshiftdiv_uv.fits')
    inbeamuvfiles.append([])
    divinbeampreselfcaluvfiles.append([])
    dividedinbeamuvfiles.append([])
for targetoutname in targetnames:
    gfile = directory + '/' + experiment + "_" + targetoutname + \
            "_pipeline_uv.gated.fits"
    ufile = directory + '/' + experiment + "_" + targetoutname + \
            "_pipeline_uv.ungated.fits"
    if expconfig['dodefaultnames']:
        gfile = directory + '/' + experiment + "_pulsar_pipeline_uv.gated.fits"
        ufile = directory + '/' + experiment + "_pulsar_pipeline_uv.ungated.fits"
    gateduvfiles.append(gfile)
    ungateduvfiles.append(ufile)
for i in range(numtargets):
    icount = 0
    for inbeamoutname in inbeamnames[i]:
        ifile = directory + '/' + experiment + "_" + \
                inbeamoutname + "_pipeline_uv.fits"
        dfile = directory + '/' + experiment + "_" + \
                inbeamoutname + "_pipeline_divided_uv.fits"
        ifile_preselfcal = directory + '/' + experiment + "_preselfcal_" + \
                inbeamoutname + "_pipeline_divided_uv.fits"
        if expconfig['dodefaultnames']:
            ifile = directory + '/' + experiment + "_inbeam-" + str(i) + \
                    "_" + str(icount) + "_pipeline_uv.fits"
            ifile_preselfcal = directory + '/' + experiment + "_preselfcal_inbeam-" + str(i) + \
                    "_" + str(icount) + "_pipeline_uv.fits"
            dfile = directory + '/' + experiment + "_inbeam-" + str(i) + \
                    "_" + str(icount) + "_pipeline_divided_uv.fits"
        inbeamuvfiles[i].append(ifile)
        divinbeampreselfcaluvfiles[i].append(ifile_preselfcal)
        dividedinbeamuvfiles[i].append(dfile)
        icount += 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Split, image and write all three ############################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Runlevel " + str(runlevel) + ": Splitting and writing final images"
    for i in range(numtargets):
        config = targetconfigs[i]
        try:
            combineifs = config['combinefinalifs']
        except KeyError:
            combineifs = False
        splitmulti = False
        splitseqno = 1
        splitbeginif = -1
        splitendif = -1
        if expconfig['ampcalscan'] > 0:
            splitband = True
        else:
            splitband = False
        if not targetonly:
            count = 0
            if imageoutofbeam:
                # First the ampcal source, only do it once
                if i==0:
                    splitdata = AIPSUVData(ampcalsrc, 'FINAL', 1, 1)
                    if splitdata.exists():
                        splitdata.zap()
                    if expconfig['ampcalscan'] > 0:
                        vlbatasks.splittoseq(inbeamuvdatas[0], bandpassclversion, 'FINAL', 
                                             ampcalsrc, splitseqno, splitmulti, splitband, splitbeginif, 
                                             splitendif, combineifs, leakagedopol)
                        ampcaluvfile = directory + '/' + experiment + "_" + ampcalsrc + "_pipeline_uv.fits"
                        vlbatasks.writedata(splitdata, ampcaluvfile, True)
                    #else:
                    #    vlbatasks.splittoseqnoband(inbeamuvdatas[0], bandpassclversion, 'FINAL', 
                    #                               ampcalsrc, 1)
                    #ampcaluvfile = directory + '/' + experiment + "_" + ampcalsrc + "_pipeline_uv.fits"
                    #vlbatasks.writedata(splitdata, ampcaluvfile, True)
                # Then the phase reference source
                aipssrcname = phscalnames[i]
                if len(aipssrcname) > 12:
                    aipssrcname = aipssrcname[:12]
                for split_phscal_option in ['FINAL','IBS']: #ibshift
                    splitdata = AIPSUVData(aipssrcname, split_phscal_option, 1, 1)
                    if splitdata.exists():
                        splitdata.zap()
                    if split_phscal_option == 'FINAL':
                        vlbatasks.splittoseq(inbeamuvdatas[0], clversion, split_phscal_option, aipssrcname, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, leakagedopol)
                        vlbatasks.writedata(splitdata, phscaluvfiles[i], True)
                    else:
                        vlbatasks.splittoseq(inbeamuvdatas[0], clversion+targetcl, split_phscal_option, aipssrcname, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, leakagedopol)
                        phscal_image_file = modeldir + aipssrcname + ".clean.fits"
                        if not os.path.exists(phscal_image_file):
                            print "Need a model for " + aipssrcname + " since --divideinbeammodel=True"
                            print "But " + phscal_image_file + " was not found."
                            sys.exit()
                        modeldata = AIPSImage(aipssrcname, "CLEAN", 1, 1)
                        if modeldata.exists():
                            modeldata.zap()
                        vlbatasks.fitld_image(phscal_image_file, modeldata)
                        divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                        if divideddata.exists():
                            divideddata.zap()
                        vlbatasks.normaliseUVData(splitdata, modeldata,  divideddata)
                        vlbatasks.writedata(divideddata, ibshiftdivphscaluvfiles[i], True)
                    #plotfile = directory + '/' + experiment + '_' + aipssrcname + '.clean.ps'
                    #if not skipplots:
                    #    vlbatasks.image(splitdata, 0.5, 512, 75, 0.5, phscalnames[i], plotfile, False,
                    #                    fullauto, stokesi)
                
            # Then the inbeams
            for inbeamsrc in inbeamnames[i]:
                aipssrcname = inbeamsrc
                if len(inbeamsrc) > 12:
                    aipssrcname = inbeamsrc[:12]
                splitdata = AIPSUVData(aipssrcname, 'FINAL', 1, 1)
                if splitdata.exists():
                    splitdata.zap()
                vlbatasks.splittoseq(inbeamuvdatas[count], clversion+targetcl, 'FINAL', inbeamsrc, 
                                     splitseqno, splitmulti, splitband, splitbeginif, splitendif, 
                                     combineifs, leakagedopol)
                #plotfile = directory + '/' + experiment + '_' + aipssrcname + '.clean.ps'
                #if not skipplots:
                #    vlbatasks.image(splitdata, 0.5, 512, 75, 0.5, inbeamsrc, plotfile, False,
                #                    fullauto, stokesi)
                vlbatasks.writedata(splitdata, inbeamuvfiles[i][count], True)
                #split primary in-beam calibrator referenced to phscal (pre-inbeamselfcal)
                if inbeamsrc in config['primaryinbeam']:
                    splitdata_PS = AIPSUVData(aipssrcname, 'PRESEL', 1, 1) #pre-selfcalibration
                    if splitdata_PS.exists():
                        splitdata_PS.zap()
                    vlbatasks.splittoseq(inbeamuvdatas[count], clversion, 'PRESEL', inbeamsrc, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, leakagedopol)
                #divide inbeamuvfiles by model
                tempdivfile = directory + '/temp.fits'
                if inbeamsrc in dividesourcelist:
                    if config['separateifmodel'] and inbeamsrc in config['primaryinbeam']:
                        ifdivideddatas = []
                        for j in range(beginif,endif+1):
                            inbeam_image_file = modeldir + inbeamsrc + ".IF" + str(j) + ".clean.fits"
                            if not os.path.exists(inbeam_image_file):
                                print "Need a frequency-dependent model for " + inbeamsrc + " since --divideinbeammodel=True and separateifmodel=True"
                                print "But " + inbeam_image_file + " was not found."
                                sys.exit()
                            shortname = inbeamsrc
                            if len(inbeamsrc) > 12:
                                shortname = inbeamsrc[:12]
                            modeldata = AIPSImage(shortname, "CLEAN", 1, 1)
                            if modeldata.exists():
                                modeldata.zap()
                            vlbatasks.fitld_image(inbeam_image_file, modeldata)
                            ifdivideddata = AIPSUVData(aipssrcname, 'IFDIV', 1, j)
                            if ifdivideddata.exists():
                                ifdivideddata.zap()
                            vlbatasks.normaliseUVData(splitdata, modeldata, ifdivideddata, j, j)
                            #uvflg = AIPSTask('uvflg', version = aipsver)
                            #uvflg.indata = ifdivideddata
                            #for f in range(beginif,endif+1):
                            #    if not f == j:
                            #        uvflg.bif = f
                            #        uvflg.eif = f
                            #        uvflg()
                            #flaggeddivideddata = AIPSUVData(aipssrcname, 'FGDIV', 1, j)
                            #if flaggeddivideddata.exists():
                            #    flaggeddivideddata.zap()
                            #splat = AIPSTask('splat', version = aipsver)
                            #splat.indata = ifdivideddata
                            #splat.flagver = 1
                            #splat.outdata = flaggeddivideddata
                            #splat()
                            #ifdivideddatas.append(flaggeddivideddata)
                            ifdivideddatas.append(ifdivideddata)
                        divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                        if divideddata.exists():
                            divideddata.zap()
                        vlbatasks.dbcon(ifdivideddatas, divideddata)
                    else:
                        inbeam_image_file = modeldir + inbeamsrc + ".clean.fits"
                        if not os.path.exists(inbeam_image_file):
                            print "Need a model for " + inbeamsrc + " since --divideinbeammodel=True"
                            print "But " + inbeam_image_file + " was not found."
                            sys.exit()
                        shortname = inbeamsrc
                        if len(inbeamsrc) > 12:
                            shortname = inbeamsrc[:12]
                        modeldata = AIPSImage(shortname, "CLEAN", 1, 1)
                        if modeldata.exists():
                            modeldata.zap()
                        vlbatasks.fitld_image(inbeam_image_file, modeldata)
                        divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                        if divideddata.exists():
                            divideddata.zap()
                        vlbatasks.normaliseUVData(splitdata, modeldata, divideddata)
                    os.system("rm -f " + tempdivfile)
                    vlbatasks.writedata(divideddata, tempdivfile, True)
                    os.system("mv -f " + tempdivfile + " " + dividedinbeamuvfiles[i][count])
                #divide pre-selfcal primary inbeam uvfile by model
                if inbeamsrc in config['primaryinbeam']:
                    divideddata_PS = AIPSUVData(aipssrcname, 'PSDIV', 1, 1) #pre-selfcal div
                    if divideddata_PS.exists():
                        divideddata_PS.zap()
                    inbeam_image_file = modeldir + inbeamsrc + ".clean.fits"
                    if not os.path.exists(inbeam_image_file):
                        print "Need a model for " + inbeamsrc + " since --divideinbeammodel=True"
                        print "But " + inbeam_image_file + " was not found."
                        sys.exit()
                    shortname = inbeamsrc
                    if len(inbeamsrc) > 12:
                        shortname = inbeamsrc[:12]
                    modeldata = AIPSImage(shortname, "CLEAN", 1, 1)
                    if modeldata.exists():
                        modeldata.zap()
                    vlbatasks.fitld_image(inbeam_image_file, modeldata)
                    vlbatasks.normaliseUVData(splitdata_PS, modeldata,  divideddata_PS)
                    #os.system("rm -f " + tempdivfile)
                    vlbatasks.writedata(divideddata_PS, tempdivfile, True)
                    os.system("mv -f " + tempdivfile + " " + divinbeampreselfcaluvfiles[i][count])
                count += 1
        if not calonly:
            splitdata1 = AIPSUVData(targetnames[i], 'UFINL', 1, 1)
            if splitdata1.exists():
                splitdata1.zap()
            if haveungated:
                try:
                    vlbatasks.splittoseq(ungateduvdata, clversion+targetcl, 'UFINL', targetnames[i], 
                                         splitseqno, splitmulti, splitband, splitbeginif, splitendif, 
                                         combineifs, leakagedopol)
                    if os.path.exists(scinttablepaths[i]) and config['scintcorrmins'] > 0:
                        vlbatasks.loadtable(splitdata1, scinttablepaths[i], 1)
                        splitdataS = AIPSUVData(targetnames[i], 'UFINS', 1, 1)
                        if splitdataS.exists():
                            splitdataS.zap()
                        vlbatasks.splat(splitdata1, 1, [0,0,0,0,0,0,0,0], splitdataS, 1)
                        splitdata1.zap()
                        splitdata1 = splitdataS
                    ungatedpresent.append(True)
                except RuntimeError:
                    print "Guess there was no ungated for " + targetnames[i]
                    ungatedpresent.append(False)
            splitdata2 = AIPSUVData(targetnames[i], 'GFINL', 1, 1)
            if splitdata2.exists():
                splitdata2.zap()
            vlbatasks.splittoseq(gateduvdata, clversion+targetcl, 'GFINL', targetnames[i], splitseqno,
                                 splitmulti, splitband, splitbeginif, splitendif, combineifs, leakagedopol)

            if os.path.exists(scinttablepaths[i]) and config['scintcorrmins'] > 0:
                vlbatasks.loadtable(splitdata2, scinttablepaths[i], 1)
                splitdataS = AIPSUVData(targetnames[i], 'GFINS', 1, 1)
                if splitdataS.exists():
                    splitdataS.zap()
                vlbatasks.splat(splitdata2, 1, [0,0,0,0,0,0,0,0], splitdataS, 1)
                splitdata2.zap()
                splitdata2 = splitdataS
            #plotfile1 = directory + '/' + experiment + '_' + targetnames[i] + '.ungated.clean.ps'
            #plotfile2 = directory + '/' + experiment + '_' + targetnames[i] + '.gated.clean.ps'
            #if fullauto or skipplots:
            #    print "Skipping the tedious imaging step..."
            #else:
            #    vlbatasks.image(splitdata1, 0.5, 512, 1, 0.00002, targetnames[i], plotfile1, True,
            #                    fullauto, stokesi)
            #    vlbatasks.image(splitdata2, 0.5, 512, 1, 0.00002, targetnames[i], plotfile2, True,
            #                    fullauto, stokesi)
            if haveungated and ungatedpresent[i]:
                vlbatasks.writedata(splitdata1, ungateduvfiles[i], True)
            vlbatasks.writedata(splitdata2, gateduvfiles[i], True)
else:
    print "Skipping splitting/writing of final images"
    for i in range(numtargets):
        if os.path.exists(ungateduvfiles[i]):
            ungatedpresent.append(True)
        else:
            ungatedpresent.append(False)

runlevel  = runlevel + 1
fullauto  = True
stokesi   = True
gaussiantarget = False
gaussianinbeam = True
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
## Image targets using Difmap and fit for position #############################
if runfromlevel <= runlevel and runtolevel >= runlevel and not calonly:
    print "Runlevel " + str(runlevel) + ": Fitting target positions using difmap"
    for i in range(numtargets):
        config = targetconfigs[i]
        icount = 0
        subtractif = 0
        try:
            fullauto = not expconfig['manualdifmapediting']
        except KeyError:
            fullauto = True
        if targetconfigs[i]['skiplastif']:
            subtractif = 1
        targetimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                          '_difmap.gated.fits'
        jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                    '.gated.difmap.jmfit'
        if expconfig['dodefaultnames']:
            targetimagefile = directory + '/' + experiment + '_pulsar' + \
                              '_difmap.gated.fits'
            jmfitfile = directory + '/' + experiment + '_pulsar' + \
                            '.gated.difmap.jmfit'
        vlbatasks.difmap_maptarget(gateduvfiles[i], targetimagefile, fullauto, stokesi,
                                   config['difmappixelmas'], config['difmapnpixels'],
                                   config['difmapweightstring'], expconfig['difmaptargetuvaverstring'], config['usegaussiantarget'],
                                   beginif, endif-subtractif)
        vlbatasks.jmfit(targetimagefile, jmfitfile, targetnames[i], stokesi, endif-subtractif)
        targetimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                          '_difmap.ungated.fits'
        jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                    '.ungated.difmap.jmfit'
        if expconfig['dodefaultnames']:
            targetimagefile = directory + '/' + experiment + '_pulsar' + \
                              '_difmap.ungated.fits'
            jmfitfile = directory + '/' + experiment + '_pulsar' + \
                            '.ungated.difmap.jmfit'
        if haveungated and ungatedpresent[i]:
            vlbatasks.difmap_maptarget(ungateduvfiles[i], targetimagefile, fullauto, 
                                       stokesi, config['difmappixelmas'], config['difmapnpixels'], 
                                       config['difmapweightstring'], expconfig['difmaptargetuvaverstring'], config['usegaussiantarget'], 
                                       beginif, endif-subtractif)
            vlbatasks.jmfit(targetimagefile, jmfitfile, targetnames[i], stokesi, endif-subtractif)
        
        #difmap-image and jmfit ibshiftdivphscal
        #targetimagefile = directory + '/' + experiment + '_' + phscalnames[i] + '_ibshiftdiv_difmap.fits'
        aipssrcname = phscalnames[i]
        if len(aipssrcname) > 12:
            aipssrcname = aipssrcname[:12]
        ibshiftedimage = AIPSImage(aipssrcname, "ICL001", 1, 1)
        if ibshiftedimage.exists():
            ibshiftedimage.zap()
        divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
        jmfitfile = directory + '/' + experiment + '_' + aipssrcname + '_ibshiftdiv_difmap.jmfit'
        #vlbatasks.difmap_maptarget(ibshiftdivphscaluvfiles[i], targetimagefile, fullauto, stokesi, config['difmappixelmas'], config['difmapnpixels'], config['difmapweightstring'], config['usegaussiantarget'], beginif, endif-subtractif)
        #vlbatasks.jmfit(targetimagefile, jmfitfile, phscalnames[i], stokesi, endif-subtractif)
        vlbatasks.widefieldimage(divideddata, aipssrcname, 256, 0.75, True, 0.050,
                   0, 0, 0, 100, 20)
        vlbatasks.nonpulsarjmfit("", jmfitfile, aipssrcname, -1, -1, True,
                             False,ibshiftedimage,48)

        #then in-beam cals    
        for j in range(len(inbeamnames[i])):
            targetimagefile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                              '_difmap.fits'
            jmfitfile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                         '.difmap.jmfit'
            if expconfig['dodefaultnames']:
                targetimagefile = directory + '/' + experiment + '_inbeam' + str(i) + \
                                  "_" + str(icount) + '_difmap.fits'
                jmfitfile = directory + '/' + experiment + '_inbeam' + str(i) + \
                            "_" + str(icount) + '.difmap.jmfit'
            vlbatasks.difmap_maptarget(inbeamuvfiles[i][j], targetimagefile, 
                                       fullauto, stokesi, config['difmappixelmas'],
                                       config['difmapnpixels'], config['difmapweightstring'], '20, False',
                                       config['usegaussianinbeam'], beginif,
                                       endif-subtractif)
            vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                            stokesi, endif-subtractif)
            if inbeamnames[i][j] in dividesourcelist:
                targetimagefile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                                  '_divided_difmap.fits'
                jmfitfile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                            '_divided.difmap.jmfit'
                if expconfig['dodefaultnames']:
                    targetimagefile = directory + '/' + experiment + '_inbeam' + str(i) + \
                                      "_" + str(icount) + '_pipeline_divided_uv.fits'
                    jmfitfile = directory + '/' + experiment + '_inbeam' + str(i) + \
                                "_" + str(icount) + '_divided.difmap.jmfit'
                vlbatasks.difmap_maptarget(dividedinbeamuvfiles[i][j], targetimagefile, 
                                           fullauto, stokesi, config['difmappixelmas'], 
                                           config['difmapnpixels'], config['difmapweightstring'], '20, False', 
                                           config['usegaussianinbeam'], beginif, 
                                           endif-subtractif)
                vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                                stokesi, endif-subtractif)
            if inbeamnames[i][j] in config['primaryinbeam']:
                #targetimagefile = directory + '/' + experiment + '_' + inbeamnames[i][j] + '_preselfcal.divided_uv.fits'
                #vlbatasks.difmap_maptarget(divinbeampreselfcaluvfiles[i][j], targetimagefile, 
                #                           fullauto, stokesi, config['difmappixelmas'],
                #                           config['difmapnpixels'], config['difmapweightstring'],
                #                           config['usegaussianinbeam'], beginif,
                #                           endif-subtractif)
                #vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                #                stokesi, endif-subtractif)
                aipssrcname = inbeamnames[i][j]
                if len(aipssrcname) > 12:
                    aipssrcname = aipssrcname[:12]
                jmfitfile = directory + '/' + experiment + '_' + aipssrcname + '_preselfcal.divided.difmap.jmfit' 
                ibshiftedimage = AIPSImage(aipssrcname, "ICL001", 1, 1)
                if ibshiftedimage.exists():
                    ibshiftedimage.zap()
                #divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                vlbatasks.widefieldimage(divideddata_PS, aipssrcname, 256, 0.75, True, 0.050,
                           0, 0, 0, 100, 20)
                vlbatasks.nonpulsarjmfit("", jmfitfile, aipssrcname, -1, -1, True,
                                   False,ibshiftedimage,48)

            icount += 1
else:
    print "Skipping imaging/position fitting of target"

runlevel  = runlevel + 1
printTableAndRunlevel(runlevel, snversion, clversion+targetcl, inbeamuvdatas[0])
"""
## Pass inbeamsn back to phscal ################################################
expconfig     = yaml.load(open(expconfigfile))
expdir        = expconfig['rootdir'] + '/' + experiment + '/'
numtargets    = len(expconfig['targets'])
targetconfigs = []
for i in range(numtargets):
    targetconfigfile = configdir + expconfig["targets"][i] + '.yaml'
    if not os.path.exists(targetconfigfile):
        parser.error("Target config file %s does not exist!" % targetconfigfile)
    targetconfigs.append(yaml.load(open(targetconfigfile)))

for i in range(numtargets):
    targetname = expconfig["targets"][i]
    targetconfig = targetconfigs[i]
    sourcefile = expdir + experiment + '.source'
    sourcefilelines = open(sourcefile).readlines()
    phsrefname = ""
    modeltypestr = 'final'
    for j in range(len(sourcefilelines)):
        if targetname in sourcefilelines[j]:
            phsrefname = sourcefilelines[j+1].split(':')[-1].strip()
    if phsrefname == "":
        print "Couldn't find the phase reference source for " + targetname
        sys.exit()
    if expconfig['useprelimmodels']:
        modeltypestr = 'preliminary'
    aips8phsrefname = phsrefname
    aips12phsrefname = phsrefname
    if len(phsrefname) > 8:
        aips8phsrefname = phsrefname[:8]
    if len(phsrefname) > 12:
        aips12phsrefname = phsrefname[:12]
    phsrefdatafile = expdir + experiment + '_' + phsrefname + '_pipeline_uv.fits'
    phsrefimagefile = auxdir + 'sourcemodels/' + modeltypestr + '/' + phsrefname + '.clean.fits'
    shiftedfile = expdir + experiment + '_' + phsrefname + '_ibshiftdiv_uv.fits'
    jmfitfile = expdir + experiment + '_' + phsrefname + '_ibshiftdiv.jmfit'
    phsrefdata = AIPSUVData(aips8phsrefname, "TEMP", 1, 1)
    phsrefimage = AIPSImage(aips12phsrefname, "TMPCLN", 1, 1)
    splatdata = AIPSUVData(aips8phsrefname, "SPLAT", 1, 1)
    normdata = AIPSUVData(aips8phsrefname, "NORM", 1, 1)
    shiftedimage = AIPSImage(aips12phsrefname, "ICL001", 1, 1)
    shiftedbeam = AIPSImage(aips12phsrefname, "IBM001", 1, 1)
    if phsrefdata.exists():
        phsrefdata.zap()
    if splatdata.exists():
        splatdata.zap()
    if phsrefimage.exists():
        phsrefimage.zap()
    if normdata.exists():
        normdata.zap()
    if shiftedimage.exists():
        shiftedimage.zap()
    if shiftedbeam.exists():
        shiftedbeam.zap()
    vlbatasks.fitld_uvfits(phsrefdatafile, phsrefdata, "")
    vlbatasks.fitld_image(phsrefimagefile, phsrefimage)
    if "," in targetconfig['primaryinbeam']:
        caltable = expdir + 'tables/CONCAT0.icalib.p1.sn'
    else:
        caltable = expdir + 'tables/' + targetconfig['primaryinbeam'] + '.icalib.p1.sn'
    if not os.path.exists(caltable):
        print caltable + ' does not exist'
        sys.exit()
    vlbatasks.loadtable(phsrefdata, caltable, 1)
    vlbatasks.splat(phsrefdata, 1, [0], splatdata, 1, 0)
    vlbatasks.normaliseUVData(splatdata, phsrefimage, normdata)
    vlbatasks.writedata(normdata, shiftedfile, True)
    vlbatasks.widefieldimage(normdata, aips12phsrefname, 256, 0.75, True, 0.050,
                   0, 0, 0, 100, 20)
    vlbatasks.nonpulsarjmfit("", jmfitfile, aips12phsrefname, -1, -1, True,
                             False,shiftedimage,48)
"""
## Make some nice diagnostic plots #############################################
if runfromlevel <= runlevel and runtolevel >= runlevel:
    print "Making final diagnostic plots"
    #os.chdir(directory)
    #os.system("%s/make_final_diagnostic.py" % codedir)
else:
    print "Skipping making of diagnostic plots"
