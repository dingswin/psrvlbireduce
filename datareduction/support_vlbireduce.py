#!/usr/bin/env ParselTongue
"""
see the docstring in vlbi_astrometry.py
"""
################################################################################
## AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
AIPSTask.isbatch = 0
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

################################################################################
## General python imports
################################################################################
import sys, os, string, math, warnings, subprocess, yaml, glob
import interaction, vlbatasks
from time import gmtime, strftime
from optparse import OptionParser
warnings.defaultaction = "always"

################################################################################
## vlbireduce class
################################################################################
class support_vlbireduce(object):
    """
    Purpose
    -------
    This is a module serving as the core of the psrvlbireduce project. It is based upon Adam's well 
    tested and versatile code. The introduction of this module would make main code short and flexible. 
    The vision is to not only focus on VLBI data reduction for PSRPI/MSPSRPI projects, but also for 
    other compact radio sources observed in various observing setups.
    
    Structure
    ---------
    The module is functionally split into two classes. The 'vlbireduce' class is the one to be 
    called, which includes all the calibration steps. The 'support_vlbireduce' class provides the data-preparation 
    and supporting (repetively used) functions.
    
    Note
    ----
    Don't use the same function name in the two classes! The class shall only be called once!
    """
    ## initiate the variables that iterates during the run #################
    runlevel = 1
    clversion = 1
    snversion = 1
    targetcl  = 0
    def __init__(self, runfromlevel, runtolevel): 
        self.runfromlevel = runfromlevel
        self.runtolevel = runtolevel
    def printTableAndRunlevel(self, runlevel, snversion, clversion, uvdata):
        """
        Print a summary of CL versions, SN versions, and runlevel
        """
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

    def parsesourcefile(self, sourcefile, experiment, klass, uvsequence):
        """
        source file parser
        Note that the source file feature will be replaced by yaml configuration from 2021.
        """
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

    def inbeamselfcal(self, doneinbeams, inbeamfilenums, inbeamuvdatas, gateduvdata,
                      expconfig, targetconfigs, modeldir, modeltype, targetonly,
                      calonly, beginif, endif, doampcal, dosecondary, sumifs,
                      clversion, targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved, 
                      leakagedopol=0):
        """
        Functionality
        -------------
        1. Run CALIB (self-calibration) on in-beam calibrators.
        2. Output sntables and delete others generated in the process of this function. 
        
        Note
        ----
        1. Need to be very cautious when removing the input variables of this function 
            because replace them with self.XX won't work out. As an example, self.doneinbeams
            and self.secondaryinbeams share the same input (doneinbeams).
        2. This function does not change snversion or clversion.
        
        Input parameters
        ----------------
        doneinbeams : list of str
            Includes inbeam calibrator (can be the pulsar in the case of inverse referencing)
            for each target to do self-calibration.
        inbeamfilenums : list of int
            used to control which uvdata to do self-calibration, essential for inverse referencing.
        clversion : int
            clversion of uvdata (see local parameters) to run CALIB. Since uvdata can be for inbeam or
            for target (when applying inverse referencing), clversion need to be aligned to avoid
            complaints.
        
        Local parameters
        ----------------
        parenttarget : int
            the index for the target.
        uvdata : str
            physical uvdata to run self-calibration.
        simplecal : bool 
            If True, then either config['separateifmodel']==True, or 
            ((not dosecondary) and len(config['primaryinbeam'].split(',')) > 1)==True.
            Plainly put, if True, either a separate-IF model is used, or the concatenation of
            multiple divided inbeamuvdatas is used to do self-calibration.

        Return parameters
        -----------------
        tocalnames : list of str
        tocalindices : list of int
            indice of targets to do inbeam self-calibration (incl. inverse referencing).
        """
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
            print(doneinbeams)
            ## >>> for inverse referencing, inbeamsrc is the pulsar!
            for (inbeamsrc, targetfilenum) in zip(doneinbeams, inbeamfilenums): 
            ## <<<
                ## >>> get parenttarget, the index for the target
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
                    print("Couldn't find a parent target for " + inbeamsrc + \
                          " - this should never happen, aborting!")
                    sys.exit()
                ## <<<
                config = targetconfigs[parenttarget]
                print("Parenttarget is", parenttarget)
                print("Sumifs is", sumifs)
                print("Doampcal is", doampcal)
                print("Desecondary is", dosecondary)
                solmins = -1
                if doampcal:
                    if not sumifs:
                        try:
                            solmins = config['inbeamcalibapnmins']
                        except KeyError:
                            print("No key for inbeamcalibapnmins, will skip inbeamapn calib for this source")
                    else:
                        try:
                            solmins = config['inbeamcalibap1mins']
                        except KeyError: pass
                else:
                    if dosecondary:
                        if not sumifs:
                            print("Can't do separate IFs secondary!")
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
                    print("Skipping " + inbeamsrc)
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
                            print("Can't do separate IFs secondary!")
                            sys.exit()
                        if doampcal:
                            sntablepath = tabledir + inbeamsrc + '.icalib.apn.sn'
                            pspath  = tabledir + inbeamsrc + '.icalib.apn.ps'
                        else:
                            sntablepath = tabledir + inbeamsrc + '.icalib.pn.sn'
                            pspath  = tabledir + inbeamsrc + '.icalib.pn.ps'
                    print("Removing", sntablepath, " and ", pspath)
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
                ## >>> here is the second ingredient for inverse referencing to happen
                if targetfilenum >= 0: ## see prepare_variables_for_inbeamselfcal()
                    uvdata = inbeamuvdatas[targetfilenum]
                else:
                    uvdata = gateduvdata
                ## <<<
                domulti = False
                combineifs = False
                haveampcal = False
                if expconfig['ampcalscan'] > 0:
                    haveampcal = True
                ## >>> split out inbeam_uv_data=inbeamsrc.CALIB.1 from uvdata
                vlbatasks.splittoseq(uvdata, clversion, 'CALIB', inbeamsrc, ## split to inbeam_uv_data
                                     1, domulti, haveampcal, beginif, 
                                     endif-subtractif, combineifs, leakagedopol)
                ## <<<
                inbeam_uv_data.table('NX', 1).zap()
                if targetfilenum >= 0:
                    inbeam_image_file = modeldir + inbeamsrc + self.cmband + ".clean.fits"
                    rawuvoutputfile = directory + inbeamsrc + self.cmband + ".formodeling.uv.fits"
                    if not os.path.exists(inbeam_image_file):
                        print("Can't find " + modeltype + " inbeam model  " + inbeam_image_file)
                        if modeltype == "preliminary":
                            print("I will write out a data file for this inbeam to " + rawuvoutputfile)
                            print("Please image it with your favourite tool")
                            print("(clean only if using difmap, no modelfitting)")
                            print("When complete, copy the image fits file to " + inbeam_image_file)
                            vlbatasks.writedata(inbeam_uv_data, rawuvoutputfile, True)
                        else:
                            print("Aborting!!")
                        sys.exit(1)
                    inbeam_image_data = AIPSImage(shortname, "CLEAN", 1, 1)
                    if inbeam_image_data.exists():
                        inbeam_image_data.zap()
                    vlbatasks.fitld_image(inbeam_image_file, inbeam_image_data)
                    ## >>> use the concatenation of more than one divided inbeamuvdatas for calibration
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
                    ## <<<
                    elif config['separateifmodel']:
                        simplecal = False
                        for i in range(beginif,endif+1-subtractif):
                            normdata = AIPSUVData(shortname, 'NORMUV', 1, i)
                            if normdata.exists():
                                normdata.zap()
                            inbeam_image_file = "%s%s.IF%d%s.clean.fits" % (modeldir, inbeamsrc, i, self.cmband)
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
                ## >>> targetfilenum<0  --> inverse referencing
                else:
                    tocaluvdata.append(inbeam_uv_data)
                    tocalimagedata.append(None) ## no image data used for pulsar
                    tocalnames.append(inbeamsrc)
                    tocalconfigs.append(config)
                ## <<<
            ## >>> use the concatenation of more than one divided inbeamuvdatas for calibration
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
            ## <<< use the concatenation of more than one divided inbeamuvdatas for calibration
            
            ######################################################
            ## the actual self-calibration starts here
            ######################################################
            for (inbeamsrc,config,inbeam_uv_data,inbeam_image_data) in \
                 zip(tocalnames,tocalconfigs,tocaluvdata,tocalimagedata):
                ## >>> some of the local parameters have been defined previously;
                ## >>> this redundancy increase robustness against future extensions.
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
                print("Using solution interval " + str(solmins) + ', requiring S/N ' + str(solsnr))
                try:
                    flagwheremodelbelow = expconfig['inbeamminmodelflux'] ## Jy
                except KeyError:
                    flagwheremodelbelow = -1
                ## <<< 
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
    
    def normalise_UVData_with_separate_IF_model_and_concatenate(self, srcname, config, expconfig, uvdata, modeldir, clversion):
        """
        used when config['phscalseparateifmodel']==True (for phscal) or config['separateifmodel']==True (for inbeams), where config=targetconfig[i]
        """
        print('Normalise uvdata of %s against its separate_IF_models, then concatenate back together.' % srcname)
        normed_IF_datas = []
        subtractif = 0
        if config['skiplastif']:
            subtractif = 1
        shortname = srcname
        if len(srcname) > 12:
            shortname = srcname[:12]
        beginif = 1
        endif = vlbatasks.getNumIFs(uvdata)
        for j in range(10,0,-1):
            splitted_uv_data = AIPSUVData(shortname, 'CALIB', 1, j)
            if splitted_uv_data.exists():
                splitted_uv_data.zap()
        # split out inbeam_uv_data=inbeamsrc.CALIB.1 from inbeamuvdata
        domulti = False
        combineifs = False
        haveampcal = False
        if expconfig['ampcalscan'] > 0:
            haveampcal = True
        vlbatasks.splittoseq(uvdata, clversion, 'CALIB', shortname,
                             1, domulti, haveampcal, beginif, 
                             endif-subtractif, combineifs, self.leakagedopol)
        splitted_uv_data.table('NX', 1).zap()

        for i in range(beginif, endif+1-subtractif):
            normdata = AIPSUVData(shortname, 'NORMUV', 1, i)
            if normdata.exists():
                normdata.zap()
            IF_image_file = "%s%s.IF%d%s.clean.fits" % (modeldir, srcname, i, self.cmband)
            if not os.path.exists(IF_image_file):
                print "Can't find inbeam model file (multi-IF %d) %s" % (i, IF_image_file)
                sys.exit()
            IF_image_data = AIPSImage(shortname, "CLEAN", 1, i)
            if IF_image_data.exists():
                IF_image_data.zap()
            vlbatasks.fitld_image(IF_image_file, IF_image_data)
            vlbatasks.normaliseUVData(splitted_uv_data, IF_image_data, normdata, i, i)
            normed_IF_datas.append(normdata)
            normdata.table('NX', 1).zap()
        for j in range(10,0,-1):
            concatuvdata = AIPSUVData('CONCAT' , 'UVSRT', 1, j)
            if concatuvdata.exists():
                concatuvdata.zap()
        for normed_IF_data in normed_IF_datas:
            vlbatasks.match_headersource(normed_IF_datas[0], normed_IF_data)
        vlbatasks.dbcon(normed_IF_datas, concatuvdata)
        for normed_IF_data in normed_IF_datas:
            normed_IF_data.zap()
        return concatuvdata
                

    def target2cals(self, targetname, expno=''): #get phscal and bandpass cal given targetname
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
        return expconfig2cals(expconfig, expno)
    
    def expconfig2cals(self, expconfig, expno=''):
        targetdir = expconfig['rootdir']
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
    
    def applyinbeamcalib(self, tocalnames, tocalindices, inbeamuvdatas, gateduvdata, 
                         expconfig, targetconfigs, targetonly, calonly, doampcal,
                         dosecondary, sumifs, clversion, snversion, inbeamnames, 
                         targetnames, haveungated, ungateduvdata, dualphscal_setup, tabledir,
                         inbeamfilenums):
        """
        Functionality
        -------------
        Apply CALIB (self-calibration) solutions from in-beam calibrators.
        
        New features
        ------------
        1. Applying inbeamcalib solution to phsref is a new feature, which allows estimation
            of uncertainties of the inbeam absolute position.
        2. Can recognize and apply existent dual-phscal solutions.
        3. Allow targetonly==True for inverse referencing. That is to say, 
            targetonly means only applying the inbeam solutions to the inbeam-calibrator 'target'
            while inverse-referencing is requested. This feature is necessitated by the 
            PSR J1939+2134 work, where dual-phscal and inverse referencing are both requested.
            In such a case, dualphscal solutions are applied to all inbeamcals (which is not
            perfect).
            Note that for dual-phscal inverse referencing to more than one inbeam calibrators, 
            one needs to change yaml files to switch an inbeamcal to be the 'target', as it is
            very rare.

        Return parameters
        -----------------
        sncount : int
           Number of inbeam calibrators on which respective solutions are applied. 
        """
        phsrefnames = []
        for targetname in targetnames:
            [phscalname, junk] = self.expconfig2cals(expconfig)
            phsrefnames.append(phscalname)
        numinbeams = len(inbeamuvdatas)
        #for i in range(numtargets):
        sncount = 0
        ###########################################################################################
        ## delete old tables and load new table(s)
        ###########################################################################################
        for (inbeamsrc, targetfilenum) in zip(tocalnames, inbeamfilenums): ## already work for inverse-referencing
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
            print(inbeamsrc)
            if "CONCAT" in inbeamsrc:
                calibtablepath = "%sCONCAT%d.icalib.%s.sn" % (tabledir, sncount, calibstring)
            else:
                calibtablepath = "%s%s.icalib.%s.sn" % (tabledir, inbeamsrc, calibstring)
            ## for the first imbeamselfcalp1 or inbeamselfcalpn apply
            if not os.path.exists(calibtablepath):
                calibtablepath = calibtablepath.replace('.dualphscal', '')

            print("Applying table " + calibtablepath)
            if targetonly and (not os.path.exists(calibtablepath)):
                print("For target-only, the SN file must exist already - aborting!")
                print(calibtablepath)
                sys.exit(1)
            
            if targetfilenum >= 0:
                if not targetonly:
                    for i in range(numinbeams):
                        for j in range(20):
                            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+j)
                        vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath, snversion+sncount)
                if not calonly:
                    for j in range(20):
                        vlbatasks.deletetable(gateduvdata, 'SN', snversion+j)
                    vlbatasks.loadtable(gateduvdata, calibtablepath, snversion+sncount)
                    if haveungated:
                        for j in range(20):
                            vlbatasks.deletetable(ungateduvdata, 'SN', snversion+j)
                        vlbatasks.loadtable(ungateduvdata, calibtablepath, snversion+sncount)
            else: ## inverse referencing
                if not targetonly:
                    for j in range(20):
                        vlbatasks.deletetable(gateduvdata, 'SN', snversion+j)
                    vlbatasks.loadtable(gateduvdata, calibtablepath, snversion+sncount)
                    if haveungated:
                        for j in range(20):
                            vlbatasks.deletetable(ungateduvdata, 'SN', snversion+j)
                        vlbatasks.loadtable(ungateduvdata, calibtablepath, snversion+sncount)
                if not calonly:
                    for i in range(numinbeams): ## load table and apply to all inbeamcals
                        for j in range(20):
                            vlbatasks.deletetable(inbeamuvdatas[i], 'SN', snversion+j)
                        vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath, snversion+sncount)
            sncount += 1
        ######################################################################################
        ## merge the new tables (for diffferent target group)
        ######################################################################################
        print("Merging SN tables between " + str(snversion) + " and " + str(snversion + sncount -1))
        if any(targetfilenum >= 0 for targetfilenum in inbeamfilenums):
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.mergesntables(inbeamuvdatas[i], snversion, sncount, expconfig['refant'])
            if not calonly:
                vlbatasks.mergesntables(gateduvdata, snversion, sncount, expconfig['refant'])
                if haveungated:
                    vlbatasks.mergesntables(ungateduvdata, snversion, sncount, expconfig['refant'])
        elif any(targetfilenum < 0 for targetfilenum in inbeamfilenums):
            if not targetonly:
                vlbatasks.mergesntables(gateduvdata, snversion, sncount, expconfig['refant'])
                if haveungated:
                    vlbatasks.mergesntables(ungateduvdata, snversion, sncount, expconfig['refant'])
            if not calonly:
                for i in range(numinbeams):
                    vlbatasks.mergesntables(inbeamuvdatas[i], snversion, sncount, expconfig['refant'])
        else:
            print('You have one target group in the inverse referencing setup, while another in the \
                normal referencing setup. It is beyond the capability of this package. Aborting now.')
            sys.exit(1)
        ######################################################################################
        ## apply table(s)
        ######################################################################################
        if any(targetfilenum >= 0 for targetfilenum in inbeamfilenums):
            if not targetonly:
                for i in range(numinbeams):
                    for j in range(10):
                        vlbatasks.deletetable(inbeamuvdatas[i], 'CL', clversion+j+1)
                    sourcelist = []
                    for j in tocalindices:
                        if i < len(inbeamnames[j]):
                            sourcelist.append(inbeamnames[j][i]) 
                    print("Applying inbeamsn for ", sourcelist, " to file ", i)
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
                print(sourcelist)
                vlbatasks.applysntable(gateduvdata, snversion+sncount, '2PT', clversion, 
                                       expconfig['refant'], sourcelist, 'CALP')
                if haveungated:
                    for j in range(10):
                        vlbatasks.deletetable(ungateduvdata, 'CL', clversion+j+1)
                    vlbatasks.applysntable(ungateduvdata, snversion+sncount, '2PT', 
                                           clversion, expconfig['refant'], sourcelist, 'CALP')
        elif any(targetfilenum < 0 for targetfilenum in inbeamfilenums):
            if not targetonly:
                for j in range(10):
                    vlbatasks.deletetable(gateduvdata, 'CL', clversion+j+1)
                sourcelist = []
                for i in tocalindices:
                    sourcelist.append(targetnames[i])
                print(sourcelist)
                vlbatasks.applysntable(gateduvdata, snversion+sncount, '2PT', clversion, 
                                       expconfig['refant'], sourcelist, 'CALP')
                if haveungated:
                    for j in range(10):
                        vlbatasks.deletetable(ungateduvdata, 'CL', clversion+j+1)
                    vlbatasks.applysntable(ungateduvdata, snversion+sncount, '2PT', 
                                           clversion, expconfig['refant'], sourcelist, 'CALP')
            if not calonly:
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
        else:
            print('This should not happen. The code is broken.')
            sys.exit(1)
        
        return sncount + 1

    def __________applyinbeamcalib(self, tocalnames, tocalindices, inbeamuvdatas, gateduvdata, 
                         expconfig, targetconfigs, targetonly, calonly, doampcal,
                         dosecondary, sumifs, clversion, snversion, inbeamnames, 
                         targetnames, haveungated, ungateduvdata, dualphscal_setup, tabledir):
        """
        apply CALIB (self-calibration) solutions from in-beam calibrators

        Return parameters
        -----------------
        sncount : int
           Number of inbeam calibrators on which respective solutions are applied. 
        """
        phsrefnames = []
        for targetname in targetnames:
            [phscalname, junk] = self.expconfig2cals(expconfig)
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
            print("Merging SN tables between " + str(snversion) + " and " + str(snversion + sncount -1))
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
