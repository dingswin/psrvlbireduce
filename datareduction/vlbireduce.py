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
from support_vlbireduce import support_vlbireduce
from time import gmtime, strftime
from optparse import OptionParser
warnings.defaultaction = "always"

class vlbireduce(support_vlbireduce):
    """
    This class is the main part of the vlbireduce module. It includes all calibration steps.
    The functions are largely organized in the order of time.
    See more docstrings in the 'support_vlbireduce' class.
    """
    def __init__(self, runfromlevel, runtolevel):
        super(vlbireduce, self).__init__(runfromlevel, runtolevel) #python2 way to use super
    
    def load_uv_data(self, directory, experiment, expconfig, numinbeams, 
            inbeamfiles, inbeamuvdatas, calonly, targetonly, gateduvdata, ungateduvdata, 
            targetnames, haveungated, ungateduvfile, gateduvfile):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Loading UVDATA from disk")
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
            print("Skipping UVDATA load")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    def increase_the_number_of_IFs_if_requested(self, targetonly, numinbeams, experiment, uvsequence, inbeamuvdatas,
            klass, calonly, gateduvdata, haveungated, ungateduvdata, expconfig):
        try:
            moriffactor = expconfig['moriffactor']
        except KeyError:
            moriffactor = 0
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and moriffactor > 1:
            print("Runlevel " + str(self.runlevel) + ": Splitting each IF " + str(moriffactor) + " ways")
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
            print("Not running MORIF")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def unflag_Pie_Town_if_requested(self, expconfig, targetonly, numinbeams, inbeamuvdatas, calonly,
            ungateduvdata, gateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and expconfig['unflagpietown']:
            print("Runlevel " + str(self.runlevel) + ": Unflagging Pie Town")
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.unflagantenna(inbeamuvdatas[i], 1, 1, 9)
            if not calonly:
                if haveungated:
                    vlbatasks.unflagantenna(ungateduvdata, 1, 1, 9)
                vlbatasks.unflagantenna(gateduvdata, 1, 1, 9)
        else:
            print("Not unflagging Pie Town")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def load_the_user_flags_if_any(self, tabledir, targetonly, numinbeams, inbeamuvdatas, numtargets, 
            inbeamnames, calonly, ungateduvdata, gateduvdata, haveungated, targetnames):
        """
        flagging can be appointed to specific inbeamcal or target, but not for phscal; for phscal, just use additionaledit.flag
        the splitted uvfits file for phscal would also flag the additionaledit.inbeamname.flag
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Loading user flags")
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
                print("No user flag file - skipping")
            for userflagfile in userflagfiles:
                if "target" in userflagfile:
                    print("Warning: old-style additionaledit.target.flag file no longer supported!")
        else:
            print("Skipping flagging from user flag file")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def run_the_flagging_for_zero_fringe_rate_effects(self, targetconfigs, targetonly, inbeamuvdatas,
            haveungated, ungateduvdata, gateduvdata, calonly):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Running UVFLG for zero fringe rates")
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
            print("Skipping flagging times of zero fringe rate")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def run_the_autoflagger_if_requested(self, expconfig, tabledir, targetonly, numinbeams, inbeamuvdatas,
            calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and expconfig['autoflag']:
            print("Runlevel " + str(self.runlevel) + ": Running autoflagger")
            autoflagfile = tabledir + "auto.rawuvdata.flag"
            if targetonly and (not os.path.exists(autoflagfile)):
                print("For target-only, autoflag file %s must already exist!" % autoflagfile)
                print("Aborting!")
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
            print("Skipping autoflagging")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def duplicate_the_initial_CL_table(self, targetonly, numinbeams, inbeamuvdatas, gateduvdata,
            ungateduvdata, calonly, haveungated):
        """
        In cases where runlevel starts from 2, CL>=2 and SN>=2 are deleted here.
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Duplicating CL table 1")
            if not targetonly:
                for i in range(numinbeams):
                    for j in range(20):
                        vlbatasks.deletetable(inbeamuvdatas[i], 'SN', j+1+self.snversion) #in cases runlevel starts from 2.
                        vlbatasks.deletetable(inbeamuvdatas[i], 'CL', j+1+self.snversion)
                    vlbatasks.tacop(inbeamuvdatas[i], 'CL', 1, inbeamuvdatas[i], 2)
            if not calonly:
                for j in range(20):
                    vlbatasks.deletetable(gateduvdata, 'SN', j+1+self.snversion) #in cases runlevel starts from 2.
                    vlbatasks.deletetable(gateduvdata, 'CL', j+1+self.snversion)
                vlbatasks.tacop(gateduvdata, 'CL', 1, gateduvdata, 2)
                if haveungated:
                    for j in range(20):
                        vlbatasks.deletetable(ungateduvdata, 'SN', j+1+self.snversion) #in cases runlevel starts from 2.
                        vlbatasks.deletetable(ungateduvdata, 'CL', j+1+self.snversion)
                    vlbatasks.tacop(ungateduvdata, 'CL', 1, ungateduvdata, 2)
        else:
            print("Skipping duplication of initial CL table")
        self.runlevel += 1
        self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def correct_positions(self, alluvdatas, targetonly, gateduvdata, haveungated, 
            ungateduvdata, calonly, inbeamuvdatas, expconfig):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Running CLCOR to shift sources")
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
                print("No sources to shift")
            for s in shifts:
                if s == "": continue
                srcname  = s.split(',')[0]
                rashift  = float(s.split(',')[1])
                decshift = float(s.split(',')[2])
                print(("Shifting " + srcname + " by " + str(rashift) + "," + str(decshift)))
                foundone = False
                for dataset in uvdatas:
                    present = False
                    for row in dataset.table('SU', 1):
                        print(row.source.strip())
                        if row.source.strip() == srcname:
                            present = True
                            foundone = True
                    if present:
                        print("About to operate on ", dataset, srcname, rashift, decshift, self.clversion)
                        vlbatasks.shift_source(dataset, srcname, rashift, decshift,
                                               self.clversion)
                if not foundone:
                    print("Didn't find source " + srcname + " in any datasets!")
                    sys.exit()
        else:
            print("Skipping source shifting")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_TECOR(self, expconfig, targetonly, numinbeams, inbeamuvdatas, logdir,
            calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not expconfig['skiptecor']:
            print("Runlevel " + str(self.runlevel) + ": Running TECOR to correct ionosphere")
            try:
                follow = expconfig['tecorfollow']
            except KeyError:
                print("Follow not specified in expconfig file, defaulting to 0.2")
                follow = 0.2
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.correct_iono(inbeamuvdatas[i], logdir, self.clversion, follow)
            if not calonly:
                vlbatasks.correct_iono(gateduvdata, logdir, self.clversion, follow)
                if haveungated:
                    vlbatasks.correct_iono(ungateduvdata, logdir, self.clversion, follow)
        else:
            print("Skipping ionospheric corrections")
        
        self.runlevel += 1
        if not expconfig['skiptecor']:
            self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def run_CLCOR_to_correct_EOPs(self, targetonly, numinbeams, inbeamuvdatas, logdir,
            calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Running CLCOR to correct for EOPs")
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.correct_eops(inbeamuvdatas[i], logdir, self.clversion)
            if not calonly:
                vlbatasks.correct_eops(gateduvdata, logdir, self.clversion)
                if haveungated:
                    vlbatasks.correct_eops(ungateduvdata, logdir, self.clversion)
        else:
            print("Skipping EOP corrections")
        self.runlevel += 1
        self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def prepare_leakage_related_variables(self, expconfig, ampcalsrc):
        self.leakagedopol = 0
        try:
            self.xpolscan = expconfig['xpolscan']
        except KeyError:
            self.xpolscan = -1
        try:
            self.xpolsource = expconfig['xpolsource']
        except KeyError:
            self.xpolsource = ampcalsrc
        try:
            self.leakagescan = expconfig['leakagescan']
            self.leakagedopol = 2
        except KeyError:
            self.leakagescan = -1
        try:
            self.leakagesource = expconfig['leakagesource']
        except KeyError:
            self.leakagesource = ampcalsrc
        try:
            self.leakageuvrange = expconfig['leakageuvrange']
        except KeyError:
            self.leakageuvrange = [0,0]
        try:
            self.leakageweightit = expconfig['leakageweightit']
        except KeyError:
            self.leakageweightit = 0

    def run_CLCOR_to_correct_PANG(self, targetonly, numinbeams,
            inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and self.leakagescan > 0:
            print("Runlevel " + str(self.runlevel) + ": Running CLCOR to correct PANG")
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.clcor_pang(inbeamuvdatas[i], self.clversion)
            if not calonly:
                vlbatasks.clcor_pang(gateduvdata, self.clversion)
                if haveungated:
                    vlbatasks.clcor_pang(ungateduvdata, self.clversion)
        else:
            print("Skipping PANG corrections")
        self.runlevel += 1
        if self.leakagescan > 0:
            self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def correct_for_a_priori_delays_if_available_using_CLCOR(self, tabledir, alluvdatas,
            targetonly, gateduvdata, ungateduvdata, haveungated, calonly, inbeamuvdatas):
        adelayfile = tabledir + 'aprioridelays.txt'
        gatedadelayfile = tabledir + 'aprioridelays.gated.txt'
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            if os.path.exists(adelayfile):
                print("Runlevel " + str(self.runlevel) + ": Using CLCOR to correct measured delays")
            else:
                print("Runlevel " + str(self.runlevel) + ": Not correcting a priori delays, just duplicating CL table")
            uvdatas = alluvdatas
            if targetonly:
                uvdatas = [gateduvdata]
                if haveungated:
                    uvdatas = [gateduvdata, ungateduvdata]
            elif calonly:
                uvdatas = inbeamuvdatas
            for uvdata in uvdatas:
                vlbatasks.tacop(uvdata, 'CL', self.clversion, uvdata, self.clversion+1)
                thisfile = adelayfile
                if uvdata == gateduvdata and os.path.exists(gatedadelayfile):
                    useroverride = interaction.yesno("There is an override file for a priori delays for the gated dataset! Do you wish to use this file? ONLY IF THERE WAS A MISTAKE WITH THE CLOCKS AT CORRELATION TIME!!! : ")
                    if useroverride:
                        thisfile = gatedadelayfile
                if os.path.exists(thisfile):
                    vlbatasks.clcordelaysfromfile(uvdata, thisfile, self.clversion+1)
        else:
            print("Skipping correction of a priori clock delays")
        self.runlevel += 1
        self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def do_the_amplitude_calibration_with_existing_table_or_in_an_interative_way__and_inspect(self, expconfig, targetonly,
            tabledir, alluvdatas, gateduvdata, ungateduvdata, haveungated, calonly, inbeamuvdatas, logdir, experiment):
        try:
            skipapcalaccor = expconfig['skipapcalaccor']
        except KeyError:
            skipapcalaccor = False
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not skipapcalaccor:
            print("Runlevel " + str(self.runlevel) + ": Doing amplitude calibration")
            if targetonly and (not os.path.exists(tabledir + 'accor.sn') or \
                               not os.path.exists(tabledir + 'apcal.sn')):
                print("For target-only, the SN files must exist already - aborting!")
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
                    vlbatasks.deletetable(inbeamuvdatas[0], 'SN', i+self.snversion)
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
                    print("Can't do tsys fix and copytsys!")
                    sys.exit()
                vlbatasks.accor(inbeamuvdatas[0], accorinterval)
                if accorampsmo > 0:
                    vlbatasks.snsmo(inbeamuvdatas[0], 'BOTH', 20, accorampsmo, 0.0, 0.0, self.snversion, expconfig['refant'], True)
                else:
                    vlbatasks.tacop(inbeamuvdatas[0], 'SN', self.snversion, inbeamuvdatas[0], self.snversion+1)
                ### load .antab file if available
                antabfile = logdir + '/' + experiment + '.antab'
                if os.path.exists(antabfile):
                    for uvdata in alluvdatas:
                        vlbatasks.antab(uvdata, antabfile, 1, 1)
                vlbatasks.apcal(inbeamuvdatas[0], self.snversion+2)
                if apcalampsmo > 0:
                    vlbatasks.snsmo(inbeamuvdatas[0], 'BOTH', 20, apcalampsmo, 0.0, 0.0, self.snversion+2, expconfig['refant'], True)
                else:
                    vlbatasks.tacop(inbeamuvdatas[0], 'SN', self.snversion+2, inbeamuvdatas[0], self.snversion+3)
                if dotsysfix:
                    print("Doing tsys fix")
                    vlbatasks.fixtsys(inbeamuvdatas[0], self.snversion+3)
                elif copytsys != "":
                    groups = copytsys.split(':')
                    for g in groups:
                        copytsyssplit = g.split(',')
                        if not len(copytsyssplit) == 3:
                            print("Copytsys was " + copytsys + ", must be --copytsys=AN1,AN2,scalefactor")
                            sys.exit()
                        if copytsyssplit[0] == copytsyssplit[1]:
                            ifsplit = copytsyssplit[2].split('@')
                            if len(ifsplit) != 2 or not ifsplit[0][:2] == "if" or not ifsplit[1][:2] == "if":
                                print("If copying between IFs of the same antenna, must use")
                                print("--copytsys=AN,AN,ifX@ifY where X and Y are input and output IFs respectively (-1 for all)")
                            vlbatasks.copytsys(inbeamuvdatas[0], self.snversion+3, copytsyssplit[0],
                                               copytsyssplit[1], 1.0, int(ifsplit[0][2:]), int(ifsplit[1][2:]))
                        else:
                            vlbatasks.copytsys(inbeamuvdatas[0], self.snversion+3, copytsyssplit[0],
                                               copytsyssplit[1], float(copytsyssplit[2]))
                        if not g==groups[-1]:
                            vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+3)
                            vlbatasks.tacop(inbeamuvdatas[0], 'SN', self.snversion+4, inbeamuvdatas[0], self.snversion+3)
                            vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+4)
                else:
                    vlbatasks.tacop(inbeamuvdatas[0], 'SN', self.snversion+3, inbeamuvdatas[0], self.snversion+4)
                snoutver1 = self.snversion + 1
                snoutver2 = self.snversion + 4
                if not expconfig['skipsnedt']:
                    vlbatasks.snedt(inbeamuvdatas[0], self.snversion+1)
                    vlbatasks.snedt(inbeamuvdatas[0], self.snversion+4)
                    snoutver1 = self.snversion+5
                    snoutver2 = self.snversion+6
                if os.path.exists(tabledir + 'accor.sn'):
                    os.remove(tabledir + 'accor.sn')
                if os.path.exists(tabledir + 'apcal.sn'):
                    os.remove(tabledir + 'apcal.sn')
                vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver1, tabledir + 'accor.sn')
                vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver2, tabledir + 'apcal.sn')
                vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver1, 'AMP', 0, 2, 4, tabledir + 'accor.ps')
                vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver2, 'AMP', 0, 2, 4, tabledir + 'apcal.ps')
                for i in range(7):
                    vlbatasks.deletetable(inbeamuvdatas[0], 'SN', i+self.snversion)
            for uvdata in uvdatas:
                vlbatasks.loadtable(uvdata, tabledir + 'accor.sn', self.snversion)
                vlbatasks.loadtable(uvdata, tabledir + 'apcal.sn', self.snversion+1)
                vlbatasks.applysntable(uvdata, self.snversion, accorinterpol, self.clversion, expconfig['refant'])
                vlbatasks.applysntable(uvdata, self.snversion+1, '2PT', self.clversion+1, expconfig['refant'])
        else:
            print("Skipping amplitude calibration")
        if not skipapcalaccor:
            self.clversion += 2
            self.snversion += 2
        
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def correct_for_primary_beam_attenuation(self, expconfig, directory, experiment,
            ampcalsrc, targetonly, numinbeams, numtargets, phscalnames, inbeamnames, tabledir, inbeamuvdatas, calonly,
            targetnames, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not expconfig['skippbcor']:
            print("Runlevel " + str(self.runlevel) + ": Doing primary beam correction")
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
                    vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion)
                    vlbatasks.correct_primarybeam(inbeamuvdatas[i], self.snversion-1, i, scanlist, fieldsourcenames, False, False)
                    if os.path.exists(pbsntable):
                        os.remove(pbsntable)
                    vlbatasks.writetable(inbeamuvdatas[i], 'SN', self.snversion, pbsntable)
                    vlbatasks.plottops(inbeamuvdatas[i], 'SN', self.snversion, 'AMP', 0, 2, 4, tabledir + 'pbcor.cal' + str(i) + '.ps')
                    vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion)
                    vlbatasks.loadtable(inbeamuvdatas[i], pbsntable, self.snversion)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion, '2PT', self.clversion, expconfig['refant'])
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
                    vlbatasks.deletetable(uvdata, 'SN', self.snversion)
                    vlbatasks.correct_primarybeam(uvdata, self.snversion-1, 0, scanlist, fieldsourcenames, False, False)
                    if os.path.exists(pbsntable):
                        os.remove(pbsntable)
                    vlbatasks.writetable(uvdata, 'SN', self.snversion, pbsntable)
                    vlbatasks.deletetable(uvdata, 'SN', self.snversion)
                    vlbatasks.loadtable(uvdata, pbsntable, self.snversion)
                    vlbatasks.applysntable(uvdata, self.snversion, 'SELN', self.clversion, expconfig['refant'])
                vlbatasks.plottops(gateduvdata, 'SN', self.snversion, 'AMP', 0, 2, 4, tabledir + 'pbcor.target.ps')
        else:
            print("Skipping primary beam correction")
        self.runlevel  = self.runlevel + 1
        if not expconfig['skippbcor']:
            self.clversion = self.clversion + 1
            self.snversion = self.snversion + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def do_PCAL_correction_and_inspect(self, expconfig, targetonly, tabledir,
            inbeamuvdatas, ampcalsrc, calonly, gateduvdata, ungateduvdata, haveungated, numinbeams):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            not expconfig['skippcal'] and expconfig['ampcalscan'] > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing pulse cal calibration")
            if targetonly and (not os.path.exists(tabledir + 'pccor.sn')):
                print("For target-only, the PC file must exist already - aborting!")
                sys.exit(1)
            if  not (targetonly or (os.path.exists(tabledir + 'pccor.sn') and \
                                    interaction.yesno("Do you wish to used saved SN table " + \
                                          "for pulse cal?"))):
                try:
                    inbeamuvdatas[0].table('SN', self.snversion).zap()
                except IOError:
                    print("No need to delete old SN table")
                vlbatasks.pccor(inbeamuvdatas[0], ampcalsrc, self.snversion, expconfig['ampcalscan'], expconfig['refant'])
                if os.path.exists(tabledir + 'pccor.sn'):
                    os.remove(tabledir + 'pccor.sn')
                vlbatasks.writetable(inbeamuvdatas[0], 'SN', self.snversion, tabledir + 'pccor.sn')
                vlbatasks.plottops(inbeamuvdatas[0], 'SN', self.snversion, 'PHAS', 0, 2, 4, tabledir + 'pccor.ps')
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion)
            if not calonly:
                vlbatasks.loadtable(gateduvdata, tabledir + 'pccor.sn', self.snversion)
                vlbatasks.applysntable(gateduvdata, self.snversion, '2PT', self.clversion, expconfig['refant'])
                if haveungated:
                    vlbatasks.loadtable(ungateduvdata, tabledir + 'pccor.sn', self.snversion)
                    vlbatasks.applysntable(ungateduvdata, self.snversion, '2PT', self.clversion, expconfig['refant'])
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'pccor.sn', self.snversion)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion, '2PT', self.clversion, expconfig['refant'])
        else:
            print("Skipping pulse cal corrections")
        self.runlevel  = self.runlevel + 1
        if not expconfig['skippcal'] and expconfig['ampcalscan'] > 0:
            self.snversion = self.snversion + 1
            self.clversion = self.clversion + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_FRING_with_fringe_finder(self, targetonly, expconfig, tabledir,
            modeldir, ampcalsrc, modeltype, inbeamuvdatas): 
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not targetonly and \
            expconfig['ampcalscan'] > 0:
            if not os.path.exists(tabledir + 'ampcalfring.sn') or \
                   not interaction.yesno("Do you wish to used saved SN table for ampcal FRING?"):
                print("Runlevel " + str(self.runlevel) + ": FRING'ing amplitude calibrator")
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
                    print("Using " + modeltype + " model for " + ampcalsrc)
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
                vlbatasks.fring(inbeamuvdatas[0], self.snversion, self.clversion, 
                                expconfig['bandpassfringmins'], inttimesecs, ampcalsrc, expconfig['refant'],
                                False, expconfig['bandpassfringsnr'], sumifs, ampcalmodeldata,
                                sumrrll, bandpassuvrange, zerorates, delaywin, ratewin,
                                doexhaustive, halfbandwidth, dispersivefit)
                print("Still at runlevel " + str(self.runlevel) + ": editing FRING for ampcal calibrator")
                vlbatasks.snsmo(inbeamuvdatas[0], 'DELA', 20, 0.0, 0.0, 10.0, self.snversion,
                                expconfig['refant'])
                snoutver = self.snversion+1
                if not expconfig['skipsnedt']:
                    vlbatasks.snedt(inbeamuvdatas[0], self.snversion+1)
                    snoutver = self.snversion+2
                if os.path.exists(tabledir + 'ampcalfring.sn'):
                    os.remove(tabledir + 'ampcalfring.sn')
                vlbatasks.writetable(inbeamuvdatas[0], 'SN', snoutver, tabledir + \
                                     'ampcalfring.sn')
                vlbatasks.plottops(inbeamuvdatas[0], 'SN', snoutver, 'DELA', 0, 2, 4,
                                   tabledir + 'ampcalfring.ps')
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion)
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+1)
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+2)
        else:
            print("Skipping FRING and SNSMO/SNEDT of ampcal FRING results")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def copy_the_FRING_SN_table_around_and_apply_it(self, expconfig, targetonly,
            tabledir, numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and expconfig['ampcalscan'] > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading ampcal FRING SN table & calibrating")
            if targetonly and (not os.path.exists(tabledir + 'ampcalfring.sn')):
                print("For target-only, the SN file must exist already - aborting!")
                sys.exit(1)
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'ampcalfring.sn', self.snversion)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion, '2PT', self.clversion,
                                           expconfig['refant'], [], 'CALP')
            if not calonly:
                vlbatasks.loadtable(gateduvdata, tabledir + 'ampcalfring.sn', self.snversion)
                vlbatasks.applysntable(gateduvdata, self.snversion, '2PT', self.clversion, 
                                       expconfig['refant'], [], 'CALP')
                if haveungated:
                    vlbatasks.loadtable(ungateduvdata, tabledir + 'ampcalfring.sn', self.snversion)
                    vlbatasks.applysntable(ungateduvdata, self.snversion, '2PT', self.clversion, 
                                           expconfig['refant'], [], 'CALP')
        else:
            print("Skipping calibration of ampcal FRING results")
        self.runlevel  = self.runlevel + 1
        if expconfig['ampcalscan'] > 0:
            self.snversion = self.snversion + 1
            self.clversion = self.clversion + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_xpoldelaycal(self, tabledir, targetonly, expconfig, modeldir, inbeamuvdatas):
        xpolsnfilename = tabledir + '/xpolfring.sn'
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not targetonly and \
            self.xpolscan > 0:
            if not os.path.exists(xpolsnfilename) or \
               not interaction.yesno("Do you wish to used saved SN table for Xpol cal?"):
                print("Runlevel " + str(self.runlevel) + ": Calculating XPOL delays")
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
                xpolmodelfile = modeldir + '/' + self.xpolsource + ".clean.fits"
                if not os.path.exists(xpolmodelfile):
                    print(xpolmodelfile, "doesn't exit: must have model for xpol calibration")
                    sys.exit()
                vlbatasks.fitld_image(xpolmodelfile, xpolmodel)
                print("int time in seconds is ", inttimesecs)
                vlbatasks.xpoldelaycal(inbeamuvdatas[0], self.clversion, expconfig['refant'], 
                                       self.xpolsource, self.xpolscan, xpolmodel, xpolsolintmins, 
                                       inttimesecs, xpolsnfilename, delaywin, ratewin)
        else:
            print("Skipping determining xpoldelays")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        return xpolsnfilename
    
    def load_and_apply_xpoldelaycal(self, targetonly, numinbeams,
            inbeamuvdatas, xpolsnfilename, expconfig, gateduvdata, ungateduvdata, haveungated, calonly):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and self.xpolscan > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading XPOL CAL FRING SN table & calibrating")
            if targetonly and (not os.path.exists(xpolsnfilename)):
                print("For target-only, the SN file must exist already - aborting!")
                sys.exit(1)
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.loadtable(inbeamuvdatas[i], xpolsnfilename, self.snversion)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion, '2PT', self.clversion,
                                           expconfig['refant'])
            if not calonly:
                vlbatasks.loadtable(gateduvdata, xpolsnfilename, self.snversion)
                vlbatasks.applysntable(gateduvdata, self.snversion, '2PT', self.clversion,
                                       expconfig['refant'])
                if haveungated:
                    vlbatasks.loadtable(ungateduvdata, xpolsnfilename, self.snversion)
                    vlbatasks.applysntable(ungateduvdata, self.snversion, '2PT', self.clversion,
                                           expconfig['refant'])
        else:
            print("Skipping calibration of xpolcal cross-pol FRING results")
        self.runlevel  = self.runlevel + 1
        if self.xpolscan > 0:
            self.snversion += 1 
            self.clversion += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_leakagecal_on_the_leakage_calibrator_and_save_AN_file(self, targetonly,
            tabledir, modeldir, inbeamuvdatas, expconfig, directory, experiment):
        leakagefilename = tabledir + "/leakage.an"
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not targetonly and \
            self.leakagescan > 0:
            if not os.path.exists(leakagefilename) or \
                   not interaction.yesno("Do you wish to used saved AN table for leakage cal?"):
                print("Runlevel " + str(self.runlevel) + ": Running leakage cal")
                if os.path.exists(leakagefilename):
                    os.remove(leakagefilename)
                leakagemodel = AIPSImage("LEAKSRC","CLNMOD",1,1)
                if leakagemodel.exists():
                    leakagemodel.zap()
                leakagemodelfile = modeldir + '/' + self.leakagesource + ".clean.fits"
                leakageoutputfile = directory + '/' + experiment + "_" + self.leakagesource + "_leakagecal_uv.fits"
                if not os.path.exists(leakagemodelfile):
                    print(leakagemodelfile, "doesn't exit: must have model for leakage calibration")
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
                vlbatasks.leakagecalc(inbeamuvdatas[0], self.leakagesource, leakagemodel, leakagefilename, 
                            expconfig['refant'], leakageacalmins, leakagepcalmins, self.leakagescan, self.clversion,
                            hasbptable, leakageoutputfile, self.leakageuvrange, self.leakageweightit)
        else:
            print("Skipping calibration of leakage")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        return leakagefilename
    
    def load_up_the_leakage_calibrated_AN_table(self, targetonly, xpolsnfilename, numinbeams,
            inbeamuvdatas, leakagefilename, calonly, gateduvdata, ungateduvdata, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and self.leakagescan > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading leakage-calculated AN table")
            if targetonly and (not os.path.exists(xpolsnfilename)):
                print("For target-only, the SN file must exist already - aborting!")
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
            print("Skipping loading the leakage calibrated AN table")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_BPASS(self, targetonly, expconfig, tabledir, modeldir, ampcalsrc, 
            modeltype, inbeamuvdatas):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not targetonly and \
            expconfig['ampcalscan'] > 0:
            print("Runlevel " + str(self.runlevel) + ": Generating bandpass corrections")
            if not os.path.exists(tabledir + 'bpass.bp') or not \
               interaction.yesno("Do you wish to used saved BP table for bandpass?"):
                ampcalmodeldata = None
                ampcalmodelfile = modeldir + ampcalsrc + '.clean.fits'
                if os.path.exists(ampcalmodelfile):
                    aipscalname = ampcalsrc
                    if len(ampcalsrc) > 12:
                        aipscalname = ampcalsrc[:12]
                    print("Using " + modeltype + " model for " + ampcalsrc)
                    ampcalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                    if ampcalmodeldata.exists():
                        ampcalmodeldata.zap()
                    vlbatasks.fitld_image(ampcalmodelfile, ampcalmodeldata)
                vlbatasks.bpass(inbeamuvdatas[0], ampcalsrc, self.clversion, 
                                expconfig['ampcalscan'], ampcalmodeldata, self.leakagedopol)
                if os.path.exists(tabledir + 'bpass.bp'):
                    os.remove(tabledir + 'bpass.bp')
                vlbatasks.writetable(inbeamuvdatas[0], 'BP', 1, tabledir + 'bpass.bp')
                vlbatasks.deletetable(inbeamuvdatas[0], 'BP', 1)
        else:
            print("Skipping bandpass corrections")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
    
    def load_BPASS_solutions(self, expconfig, tabledir, calonly, 
            gateduvdata, ungateduvdata, targetonly, numinbeams, inbeamuvdatas, haveungated):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            expconfig['ampcalscan'] > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading bandpass corrections")
            if not os.path.exists(tabledir + 'bpass.bp'):
                print("Error - bandpass table " + tabledir + 'bpass.bp does not exist')
                sys.exit(1)
            if not calonly:
                vlbatasks.loadtable(gateduvdata, tabledir + 'bpass.bp', 1)
                if haveungated:
                    vlbatasks.loadtable(ungateduvdata, tabledir + 'bpass.bp', 1)
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'bpass.bp', 1)
        else:
            print("Skipping loading of bandpass corrections")
        bandpassclversion = self.clversion
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        return bandpassclversion 
    
    def plot_bandpass(self, plotbandpass, directory, inbeamuvdatas, ifs, chans):
        if plotbandpass:
            bandpassplot = directory + '/' + 'bandpass.ps' 
            #vlbatasks.plotbandpass(inbeamuvdatas[0], 1, True, 9, bandpassplot, self.clversion, ifs, 0, chans, "I")
            vlbatasks.plotbandpass(inbeamuvdatas[0], 1, True, 9, bandpassplot, self.clversion, ifs, 0, chans, "I")
        else:
            pass

    def prepare_variables_for_calibration_on_phase_calibrators(self, targetconfigs, phscalnames):
        """
        prepare variables for FRING (phase reference calibrators)
        """
        self.maxphsreffringmins = -1
        self.maxphsrefcalibapnmins = -1
        self.maxphsrefcalibpnmins = -1
        for config in targetconfigs:
            if config['phsreffringmins'] > self.maxphsreffringmins:
                self.maxphsreffringmins = config['phsreffringmins']
            if config['phsrefcalibapnmins'] > self.maxphsrefcalibapnmins:
                self.maxphsrefcalibapnmins = config['phsrefcalibapnmins']
            try:
                if config['phsrefcalibpnmins'] > self.maxphsrefcalibpnmins:
                    self.maxphsrefcalibpnmins = config['phsrefcalibpnmins']
            except KeyError:
                pass
        self.donephscalnames = []
        self.doneconfigs = []
        for i in range(len(phscalnames)):
            phscal = phscalnames[i]
            if phscal in self.donephscalnames:
                continue
            self.donephscalnames.append(phscal)
            self.doneconfigs.append(targetconfigs[i])

    def run_FRING_with_phase_reference_calibrators(self, targetonly,
            tabledir, modeldir, modeltype, expconfig, experiment, inbeamuvdatas):
        dophscaldump = False
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not targetonly and \
           self.maxphsreffringmins > 0:
            if not os.path.exists(tabledir + 'phsreffring.sn') or \
                   not interaction.yesno("Do you wish to used saved SN table for phsref FRING?"):
                print("Runlevel " + str(self.runlevel) + ": FRING'ing phase calibrator")
                for phscal, config in zip(self.donephscalnames, self.doneconfigs):
                    try:
                        phscalseparateifmodel = config['phscalseparateifmodel']
                    except KeyError:
                        phscalseparateifmodel = False
                    phscalmodeldata = None
                    phscalmodelfile = modeldir + phscal + self.cmband + '.clean.fits'
                    if not phscalseparateifmodel:
                        if os.path.exists(phscalmodelfile):
                            aipscalname = phscal
                            if len(phscal) > 12:
                                aipscalname = phscal[:12]
                            print("Using " + modeltype + " model for " + phscal)
                            phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                            if phscalmodeldata.exists():
                                phscalmodeldata.zap()
                            vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
                        else:
                            print("Currently no " + modeltype + " model for " + phscal)
                            if modeltype == "preliminary":
                                try:
                                    allowphscalpointmodel = expconfig['allowphscalpointmodel']
                                    if not allowphscalpointmodel:
                                        print("Using a point source, will dump output after FRING")
                                        dophscaldump = True
                                except KeyError:
                                    print("Using a point source, will dump output after FRING")
                                    dophscaldump = True
                            else:
                                print("Aborting!")
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
                            print("Overriding sumifs because no amp cal fring was done")
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
                    print("Delay and rate windows: ", delaywin, ratewin)
                    print("UV range: ", phsrefuvrange)
                    if phscalseparateifmodel:
                        tofringuvdata = self.normalise_UVData_with_separate_IF_model_and_concatenate(phscal, config, expconfig, inbeamuvdatas[0],\
                            modeldir, self.clversion)
                        #phscalmodeldata = None
                        #doband = False
                        vlbatasks.fring(tofringuvdata, 1, -1,  
                                        config['phsreffringmins'], inttimesecs, phscal, 
                                        expconfig['refant'], False, 
                                        config['phsreffringsnr'], sumifs, None,
                                        sumrrll, phsrefuvrange,False,delaywin,ratewin,
                                        doexhaustive, halfbandwidth, dispersivefit, self.leakagedopol)
                        vlbatasks.writetable(tofringuvdata, 'SN', 1, tabledir + \
                                             'multi_IF_raw_phsreffring.sn')
                        vlbatasks.loadtable(inbeamuvdatas[0], tabledir + 'multi_IF_raw_phsreffring.sn', self.snversion)
                    else:
                        vlbatasks.fring(inbeamuvdatas[0], self.snversion, self.clversion,  
                                        config['phsreffringmins'], inttimesecs, phscal, 
                                        expconfig['refant'], doband, 
                                        config['phsreffringsnr'], sumifs, phscalmodeldata,
                                        sumrrll, phsrefuvrange,False,delaywin,ratewin,
                                        doexhaustive, halfbandwidth, dispersivefit, self.leakagedopol)
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
                snoutver = self.snversion
                if dosnsmo:
                    #vlbatasks.snsmo(inbeamuvdatas[0], 'DELA', 20, 0.0, 0.0, 10.0, self.snversion,
                    #                smoothedrefant, True)
                    if sumifs:
                        doratesmoothing = False
                        ratesmoothingmins = 0.0
                    else:
                        doratesmoothing = True
                        ratesmoothingmins = 3.0
                    vlbatasks.snmwfclip(inbeamuvdatas[0], mwfminutes, 0.0, 0.0, delaymwf, ratemwf, self.snversion,
                                        smoothedrefant, doratesmoothing, ratesmoothingmins)
                    snoutver = self.snversion+1
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
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion)
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+1)
                vlbatasks.deletetable(inbeamuvdatas[0], 'SN', self.snversion+2)
        else:
            print("Skipping FRING and SNSMO/SNEDT of FRING results")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        return dophscaldump

    def copy_the_phsref_FRING_SN_table_around_and_apply_it(self, targetonly, tabledir, expconfig,
            numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated, directory, experiment, dophscaldump):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and self.maxphsreffringmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading phsref FRING SN table & calibrating")
            if targetonly and (not os.path.exists(tabledir + 'phsreffring.sn')):
                print("For target-only, the SN file must exist already - aborting!")
                sys.exit(1)
            try:
                interptype = expconfig['phsreffringinterp']
            except KeyError:
                interptype = 'AMBG'
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.loadtable(inbeamuvdatas[i], tabledir + 'phsreffring.sn', self.snversion)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion, interptype, self.clversion, 
                                           expconfig['refant'])
            if not calonly:
                vlbatasks.loadtable(gateduvdata, tabledir + 'phsreffring.sn', self.snversion)
                vlbatasks.applysntable(gateduvdata, self.snversion, interptype, self.clversion, 
                                       expconfig['refant'])
                if haveungated:
                    vlbatasks.loadtable(ungateduvdata, tabledir + 'phsreffring.sn', self.snversion)
                    vlbatasks.applysntable(ungateduvdata, self.snversion, interptype, self.clversion, 
                                           expconfig['refant'])
        else:
            print("Skipping calibration of FRING results")

        if self.maxphsreffringmins > 0:
            self.snversion = self.snversion + 1
            self.clversion = self.clversion + 1
        print(self.donephscalnames)
        
        if dophscaldump: # Need to dump out the phs cal sources so we can make models of them
            for phscal in self.donephscalnames:
                for i in range(20): #Clear any old CALIB split catalog entries
                    phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, i)
                    if phscal_uv_data.exists():
                        phscal_uv_data.zap()
                phscal_uv_data = AIPSUVData(phscal[:12], 'CALIB', 1, 1)
                rawuvoutputfile = directory + experiment.upper() + '_' + \
                                           phscal + self.cmband + ".formodeling.uv.fits"
                doband = False
                domulti = False
                if expconfig['ampcalscan'] > 0:
                    doband = True
                combineifs = False
                beginif = -1
                endif = -1
                vlbatasks.splittoseq(inbeamuvdatas[0], self.clversion, 'CALIB', phscal, 1, domulti,
                                     doband, beginif, endif, combineifs, self.leakagedopol)
                vlbatasks.writedata(phscal_uv_data, rawuvoutputfile, True)
            print("UV datasets of the phase reference sources have been written out to model")
            sys.exit()
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_phase_CALIB_on_the_phase_reference_sources(self,
            tabledir, targetonly, inbeamuvdatas, expconfig, modeltype, modeldir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxphsrefcalibpnmins > 0:
            haveall = True
            for phscal in self.donephscalnames:
                calibtablepath = tabledir + phscal + '.calibpn.sn'
                if not os.path.exists(calibtablepath):
                    haveall = False
            if not targetonly and (not haveall or \
                   not interaction.yesno("Do you wish to use saved SN table for phsref phase CALIB?")):
                print("Runlevel " + str(self.runlevel) + ": Running CALIB on phsref sources")
                for phscal, config in zip(self.donephscalnames, self.doneconfigs):
                    try:
                        phscalseparateifmodel = config['phscalseparateifmodel']
                    except KeyError:
                        phscalseparateifmodel = False
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
                    vlbatasks.splittoseq(inbeamuvdatas[0], self.clversion, 'CALIB', phscal, 1, domulti,
                                         doband, beginif, endif, combineifs, self.leakagedopol)
                    phscal_uv_data.table('NX', 1).zap()
                    phscalmodeldata = None
                    phscalmodelfile = modeldir + phscal + self.cmband + '.clean.fits'
                    try:
                        phsrefuvrange = config['phsrefuvrange']
                    except KeyError:
                        phsrefuvrange = [0, 0]
                    try:
                        phsrefweightit = expconfig['phsrefcalibpnweightit']
                    except KeyError:
                        phsrefweightit = 0
                    if os.path.exists(phscalmodelfile):
                        print("Using " + modeltype + " model for " + phscal)
                        aipscalname = phscal
                        if len(phscal) > 12:
                            aipscalname = phscal[:12]
                        phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                        if phscalmodeldata.exists():
                            phscalmodeldata.zap()
                        vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
                    else:
                        print("Currently no " + modeltype + " model for " + phscal + '; aborting!')
                        sys.exit()
                    calibsoltype = config['phsrefcalibpntype']
                    calibtablepath = tabledir + phscal + '.calibpn.sn'
                    try:
                        flagwheremodelbelow = expconfig['phsrefminmodelflux'] #Jy
                    except KeyError:
                        flagwheremodelbelow = -1
                    if phscalseparateifmodel:
                        phscal_uv_data = self.normalise_UVData_with_separate_IF_model_and_concatenate(phscal, config, expconfig,\
                            inbeamuvdatas[0], modeldir, self.clversion)
                        vlbatasks.singlesource_calib(phscal_uv_data, None,
                                                     1, expconfig['refant'], False, 
                                                     config['phsrefcalibpnmins'],
                                                     False, calibsoltype, 
                                                     config['phsrefcalibpnsnr'], False,
                                                     phsrefuvrange, phsrefweightit,
                                                     flagwheremodelbelow)
                    else:
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
            print("Skipping phase cal CALIB on the phase reference sources")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def load_all_the_phase_CALIB_solutions_obtained_with_phase_calibrators(self,
            phscalnames, targetonly, numinbeams, inbeamuvdatas, calonly, gateduvdata, haveungated, ungateduvdata,
            tabledir, expconfig):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxphsrefcalibpnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading phs ref PN CALIB " + \
                  "solutions and applying")
            sncount = 0
            for phscal in phscalnames:
                if not targetonly:
                    for i in range(numinbeams):
                        vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion+sncount)
                if not calonly:
                    vlbatasks.deletetable(gateduvdata, 'SN', self.snversion+sncount)
                    if haveungated:
                        vlbatasks.deletetable(ungateduvdata, 'SN', self.snversion+sncount)
                sncount += 1
            sncount = 0
            for phscal in self.donephscalnames:
                calibtablepath = tabledir + phscal + '.calibpn.sn'
                if targetonly and (not os.path.exists(calibtablepath)):
                    print("For target-only, the SN file must exist already - aborting!")
                    sys.exit(1)
                if not targetonly:
                    for i in range(numinbeams):
                        vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath,
                                            self.snversion+sncount)
                if not calonly:
                    vlbatasks.loadtable(gateduvdata, calibtablepath,
                                        self.snversion+sncount)
                    if haveungated:
                        vlbatasks.loadtable(ungateduvdata, calibtablepath,
                                            self.snversion+sncount)
                sncount += 1
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion+sncount)
                    vlbatasks.mergesntables(inbeamuvdatas[i], self.snversion, sncount, 
                                            expconfig['refant'])
            if not calonly:
                vlbatasks.deletetable(gateduvdata, 'SN', self.snversion+sncount)
                vlbatasks.mergesntables(gateduvdata, self.snversion, sncount, 
                                        expconfig['refant'])
                if haveungated:
                    vlbatasks.deletetable(ungateduvdata, 'SN', self.snversion+sncount)
                    vlbatasks.mergesntables(ungateduvdata, self.snversion, sncount, 
                                            expconfig['refant'])
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.deletetable(inbeamuvdatas[i], 'CL', self.clversion+1)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion+sncount, 'SELF', 
                                           self.clversion, expconfig['refant'])
            if not calonly:
                vlbatasks.deletetable(gateduvdata, 'CL', self.clversion+1)
                vlbatasks.applysntable(gateduvdata, self.snversion+sncount, 'SELF', 
                                       self.clversion, expconfig['refant'])
                if haveungated:
                    vlbatasks.deletetable(ungateduvdata, 'CL', self.clversion+1)
                    vlbatasks.applysntable(ungateduvdata, self.snversion+sncount, 'SELF', 
                                           self.clversion, expconfig['refant'])
        else:
            print("Skipping loading/application of phase reference source phase CALIB solutions")
            sncount = 0
            if self.maxphsrefcalibpnmins > 0:
                sncount += len(self.donephscalnames)

        if self.maxphsrefcalibpnmins > 0:
            self.snversion = self.snversion + sncount + 1
            self.clversion = self.clversion + 1
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def run_amp_CALIB_on_the_phase_reference_sources(self, targetconfigs,
            tabledir, targetonly, expconfig, inbeamuvdatas, modeldir, modeltype, phscalnames):
        #prepare two variables for amp CALIB on the phase calibrators
        self.donephscalnames = []
        self.doneconfigs = []
        for i in range(len(phscalnames)):
            phscal = phscalnames[i]
            try:
                if phscal in self.donephscalnames or targetconfigs[i]['phsrefcalibapnmins'] < 0:
                    continue
                self.donephscalnames.append(phscal)
                self.doneconfigs.append(targetconfigs[i])
            except KeyError:
                continue

        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxphsrefcalibapnmins > 0:
            haveall = True
            for phscal in self.donephscalnames:
                calibtablepath = tabledir + phscal + '.calibapn.sn'
                if not os.path.exists(calibtablepath):
                    haveall = False
            if not targetonly and (not haveall or \
                   not interaction.yesno("Do you wish to used saved SN table for phsref amp CALIB?")):
                print("Runlevel " + str(self.runlevel) + ": Running CALIB on phsref sources")
                for phscal, config in zip(self.donephscalnames, self.doneconfigs):
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
                    vlbatasks.splittoseq(inbeamuvdatas[0], self.clversion, 'CALIB', phscal, 1, domulti,
                                         doband, beginif, endif, combineifs, self.leakagedopol)
                    phscal_uv_data.table('NX', 1).zap()
                    phscalmodeldata = None
                    phscalmodelfile = modeldir + phscal + self.cmband + '.clean.fits'
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
                        print("Using " + modeltype + " model for " + phscal)
                        aipscalname = phscal
                        if len(phscal) > 12:
                            aipscalname = phscal[:12]
                        phscalmodeldata = AIPSImage(aipscalname, 'CLNMOD', 1, 1)
                        if phscalmodeldata.exists():
                            phscalmodeldata.zap()
                        vlbatasks.fitld_image(phscalmodelfile, phscalmodeldata)
                    else:
                        print("Currently no " + modeltype + " model for " + phscal + '; aborting!')
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
                           print("Overriding default value of phsrefcalibapnclip (", phsrefcalibapnclip, ") for this source with", tmpphsrefcalibapnclip)
                        phsrefcalibapnclip = tmpphsrefcalibapnclip
                    except KeyError:
                        pass
                    try:
                        phscalseparateifmodel = config['phscalseparateifmodel']
                    except KeyError:
                        phscalseparateifmodel = False
                    if phscalseparateifmodel:
                        phscal_uv_data = self.normalise_UVData_with_separate_IF_model_and_concatenate(phscal, config, expconfig,\
                            inbeamuvdatas[0], modeldir, self.clversion)
                        vlbatasks.singlesource_calib(phscal_uv_data, None,
                                                     1, expconfig['refant'], True,
                                                     config['phsrefcalibapnmins'],
                                                     False, calibsoltype,
                                                     config['phsrefcalibapnsnr'], False,
                                                     phsrefuvrange, phsrefweightit,
                                                     flagwheremodelbelow, normalise)
                    else:
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
            print("Skipping amp cal CALIB on the phase reference sources")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        #return self.donephscalnames

    def load_all_the_amp_CALIB_solutions_obtained_with_phase_calibrators(self,
            phscalnames, targetonly, numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated,
            tabledir, expconfig):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxphsrefcalibapnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Loading phs ref APN CALIB " + \
                  "solutions and applying")
            sncount = 0
            for phscal in phscalnames:
                if not targetonly:
                    for i in range(numinbeams):
                        vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion+sncount)
                if not calonly:
                    vlbatasks.deletetable(gateduvdata, 'SN', self.snversion+sncount)
                    if haveungated:
                        vlbatasks.deletetable(ungateduvdata, 'SN', self.snversion+sncount)
                sncount += 1
            sncount = 0
            for phscal in self.donephscalnames:
                calibtablepath = tabledir + phscal + '.calibapn.sn'
                if targetonly and (not os.path.exists(calibtablepath)):
                    print("For target-only, the SN file must exist already - aborting!")
                    sys.exit(1)
                if not targetonly:
                    for i in range(numinbeams):
                        vlbatasks.loadtable(inbeamuvdatas[i], calibtablepath,
                                            self.snversion+sncount)
                if not calonly:
                    vlbatasks.loadtable(gateduvdata, calibtablepath,
                                        self.snversion+sncount)
                    if haveungated:
                        vlbatasks.loadtable(ungateduvdata, calibtablepath,
                                            self.snversion+sncount)
                sncount += 1
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.deletetable(inbeamuvdatas[i], 'SN', self.snversion+sncount)
                    vlbatasks.mergesntables(inbeamuvdatas[i], self.snversion, sncount,
                                            expconfig['refant'])
            if not calonly:
                vlbatasks.deletetable(gateduvdata, 'SN', self.snversion+sncount)
                vlbatasks.mergesntables(gateduvdata, self.snversion, sncount,
                                        expconfig['refant'])
                if haveungated:
                    vlbatasks.deletetable(ungateduvdata, 'SN', self.snversion+sncount)
                    vlbatasks.mergesntables(ungateduvdata, self.snversion, sncount,
                                            expconfig['refant'])
            if not targetonly:
                for i in range(numinbeams):
                    vlbatasks.deletetable(inbeamuvdatas[i], 'CL', self.clversion+1)
                    vlbatasks.applysntable(inbeamuvdatas[i], self.snversion+sncount, 'SELF',
                                           self.clversion, expconfig['refant'])
            if not calonly:
                vlbatasks.deletetable(gateduvdata, 'CL', self.clversion+1)
                vlbatasks.applysntable(gateduvdata, self.snversion+sncount, 'SELF',
                                       self.clversion, expconfig['refant'])
                if haveungated:
                    vlbatasks.deletetable(ungateduvdata, 'CL', self.clversion+1)
                    vlbatasks.applysntable(ungateduvdata, self.snversion+sncount, 'SELF',
                                           self.clversion, expconfig['refant'])
        else:
            print("Skipping loading/application of phase reference source amp CALIB solutions")
            sncount = 0
            if self.maxphsrefcalibapnmins > 0:
                sncount += len(self.donephscalnames)

        if self.maxphsrefcalibapnmins > 0:
            self.snversion = self.snversion + sncount + 1
            self.clversion = self.clversion + 1
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def prepare_variables_for_inbeamselfcal(self, inbeamuvdatas, targetconfigs, numtargets, inbeamnames, targetnames):
        """
        Functionality
        -------------
        prepare some variables for operations on inbeam calibrators.
        
        Note
        ----
        1. must be run before any inbeamselfcal step.
        2. inverse referencing is made possible here.
        3. self.doneinbeams include inbeam calibrator for each target to do self-calibration.
        """
        print('\nNow preparing variables for inbeam self-calibrations...')
        self.doneinbeams = []
        self.inbeamfilenums = []
        self.secondaryinbeams = []
        self.secondaryfilenums = []
        self.numifs = 4
        self.maxinbeamcalibp1mins = -1 # Phase-only, primary inbeam, summed IFs
        self.maxinbeamcalibpnmins = -1 # Phase-only, primary inbeam, separate IFs
        self.maxinbeamcalibsp1mins = -1 # Phase-only, secondary inbeam (rarely used)
        self.maxinbeamcalibap1mins = -1 # Amplitude and phase, primary inbeam, combined IFs
        self.maxinbeamcalibapnmins = -1 # Amplitude and phase, primary inbeam, separate IFs
        if inbeamuvdatas[0].exists():
            self.numifs = vlbatasks.getNumIFs(inbeamuvdatas[0])
        self.beginif = 1
        self.endif = self.numifs
        for config in targetconfigs:
            if config['inbeamcalibp1mins'] > self.maxinbeamcalibp1mins:
                self.maxinbeamcalibp1mins = config['inbeamcalibp1mins']
            if config['inbeamcalibap1mins'] > self.maxinbeamcalibap1mins:
                self.maxinbeamcalibap1mins = config['inbeamcalibap1mins']
            try:
                if config['inbeamcalibsp1mins'] > self.maxinbeamcalibsp1mins:
                    self.maxinbeamcalibsp1mins = config['inbeamcalibsp1mins']
            except KeyError:
                print("No secondary inbeam self-calibration requested...")
            try:
                if config['inbeamcalibpnmins'] > self.maxinbeamcalibpnmins:
                    self.maxinbeamcalibpnmins = config['inbeamcalibpnmins']
            except KeyError:
                print("No separate IF inbeam selfcal...")
            try:
                if config['inbeamcalibapnmins'] > self.maxinbeamcalibapnmins:
                    self.maxinbeamcalibapnmins = config['inbeamcalibapnmins']
            except KeyError:
                print("No separate IF inbeam amp selfcal...")
        print("maxinbeamcalibp1mins", self.maxinbeamcalibp1mins)
        print("maxinbeamcalibpnmins", self.maxinbeamcalibpnmins)
        print("maxinbeamcalibap1mins", self.maxinbeamcalibap1mins)
        print("maxinbeamcalibapnmins", self.maxinbeamcalibapnmins)
        print("maxinbeamcalibsp1mins", self.maxinbeamcalibsp1mins)
        for i in range(numtargets):
            primaryinbeams = targetconfigs[i]['primaryinbeam'].split(',')
            try:
                secondaryinbeam = targetconfigs[i]["secondaryinbeam"]
                secondaryinbeam = secondaryinbeam.split(',')
            except KeyError:
                secondaryinbeam = "XXXXX"
            for j in range(len(inbeamnames[i])):
                if self.maxinbeamcalibsp1mins > 0:
                    for inbeamsrc in secondaryinbeam:
                        if not secondaryinbeam == None and (inbeamsrc.strip() == inbeamnames[i][j].strip()):
                            self.secondaryinbeams.append(inbeamsrc.strip())
                            self.secondaryfilenums.append(j)
                    if len(self.secondaryinbeams) == 0:
                        print("\nSecondary inbeam self-calibration is requested.\nHowever, either secondary inbeam calibrator(s) is not specified, or it does not exist. Aborted.")
                        sys.exit(1)
                for primaryinbeam in primaryinbeams:
                    if primaryinbeam.strip() == inbeamnames[i][j].strip():
                        if not primaryinbeam.strip() in self.doneinbeams:
                            self.doneinbeams.append(primaryinbeam.strip())
                            self.inbeamfilenums.append(j)
                    ## here is the first ingredient of inverse referencing
                    if primaryinbeam.strip() == targetnames[i].strip():
                        if not targetnames[i] in self.doneinbeams: 
                            self.doneinbeams.append(targetnames[i]) 
                            self.inbeamfilenums.append(-1) ## < 0 --> inverse referencing, see inbeamselfcal()
                    ##################
        for primaryinbeam in primaryinbeams:
            if not primaryinbeam.strip() in self.doneinbeams:
                print("Didn't find primary inbeam " + primaryinbeam.strip() + " amongst the data! Check the inbeam name in your config file!")
                print(self.doneinbeams)
                sys.exit()

    def generate_a_raw_inbeam_dataset_without_phase_selfcal_on_itself_if_requested(self, expconfig,
            targetonly, numtargets, targetconfigs, inbeamnames, directory, experiment, inbeamuvdatas):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
           expconfig['writerawinbeam']:
            print("Runlevel " + str(self.runlevel) + ": Writing raw inbeam outputs")
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
                        rawuvoutputfile =  directory + experiment.lower() + '_' + \
                                           inbeamsrc + self.cmband + ".formodeling.uv.fits"
                        splitdata = AIPSUVData(aipssrcname, 'NOIB', 1, 1)
                        if splitdata.exists():
                            splitdata.zap()
                        uvdata = inbeamuvdatas[count]
                        doband = False
                        domulti = False
                        if expconfig['ampcalscan'] > 0:
                            doband = True
                        combineifs = False
                        vlbatasks.splittoseq(uvdata, self.clversion, 'NOIB', aipssrcname, 1,
                                             domulti, doband, self.beginif, self.endif-subtractif, 
                                             combineifs, self.leakagedopol)
                        vlbatasks.writedata(splitdata, rawuvoutputfile, True)
                        count += 1
        else:
            print("Skipping dump of raw inbeam data (pre-selfcal)")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])

    def do_a_combined__IF_and_pol__phase_selfcal_on_the_inbeams_if_requested(self,
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype, targetonly,
            calonly, targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibp1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing phase-only inbeam selfcal (combined IFs)")
            tocalnames, tocalindices = self.inbeamselfcal(self.doneinbeams, self.inbeamfilenums, inbeamuvdatas, gateduvdata, 
                                       expconfig, targetconfigs, modeldir, modeltype, targetonly, 
                                       calonly, self.beginif, self.endif, False, False, True, self.clversion, targetnames, numtargets, 
                                       inbeamnames, directory, tabledir, alwayssaved, self.leakagedopol)
        else:
            print("Skipping inbeam phase-only selfcal (combined IFs)")
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
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion, inbeamuvdatas[0])
        return tocalnames, tocalindices

    def load_inbeam_CALIB_solutions_obtained_with__IF_and_pol__combined(self, tocalnames, 
            tocalindices, inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, 
            haveungated, ungateduvdata, tabledir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibp1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Applying inbeam CALIB p1 sols")
            sncount = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, 
                                       targetconfigs, targetonly, calonly, False, False, True,
                                       self.clversion, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                       ['-1','0'], tabledir, self.inbeamfilenums)
            ## sncount is used to point at SN table in post-phscal stage
        else:
            print("Skipping application of inbeam phase-only selfcal (combined IFs)")
            if self.maxinbeamcalibp1mins > 0:
                sncount = len(tocalnames) + 1
            else:
                sncount = 0
        if self.maxinbeamcalibp1mins > 0:
            self.snversion = self.snversion + sncount
            self.targetcl += 1
        self.runlevel = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def do_a_separate_IF_phase_selfcal_on_the_inbeams_if_requested(self, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype, targetonly,
            calonly, targetnames, numtargets, directory, tabledir, alwayssaved, inbeamnames):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibpnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing phase-only inbeam selfcal (separate IFs)")
            tocalnames, tocalindices = self.inbeamselfcal(self.doneinbeams, self.inbeamfilenums, inbeamuvdatas, gateduvdata,
                                       expconfig, targetconfigs, modeldir, modeltype, targetonly,
                                       calonly, self.beginif, self.endif, False, False, False, self.clversion+self.targetcl, 
                                       targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved, self.leakagedopol)
        else:
            print("Skipping inbeam phase-only selfcal (separate IFs)")
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
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])
        return tocalnames, tocalindices
    

    def do_dual_phscal_calibration_correcting_the_CALIB_solutions_on_inbeams_with__IF_and_pol__combined(self, dualphscal_setup, 
            directory, tabledir, inbeamuvdatas, gateduvdata, ungateduvdata, targetonly, calonly, haveungated, tocalnames, tocalindices, 
            expconfig, targetconfigs, inbeamnames, targetnames):
        """
        Functionality
        -------------
        Correct the inbeamcalibp1 solutions, then apply to the target only.
        
        Note
        ----
        So far dualphscal_setup applies for all target groups. But it can adapt easily if necessary.
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
           int(dualphscal_setup[0].strip()) > 0:
            print("Adopt dual-phscal mode now...")
            inbeamselfcal_phase_time_folder = directory+'/inbeamselfcal_phase_time_evolution'
            if not os.path.exists(inbeamselfcal_phase_time_folder):
                os.system('mkdir %s' % inbeamselfcal_phase_time_folder)
            
            for (tocalname, tocalindex) in zip(tocalnames, tocalindices):
                ## >>> note that this function does not work for the scenario where two inbeams are provided as primaryinbeam!
                inbeamselfcalp1sntable = tabledir + '/' + tocalname + '.icalib.p1.sn'
                ## <<<
                
                ## >>>> use the real inbeamcal data (target can work the same) as a host to produce dualphscal solutions
                filenum = self.inbeamfilenums[tocalindex]
                if filenum >= 0:
                    uvdata = inbeamuvdatas[filenum] ## uvdata here refers to the uvdata to do self-calibration
                else:
                    uvdata = gateduvdata
                
                for i in range(20):
                    vlbatasks.deletetable(uvdata, 'SN', self.snversion+i) ## supposed to be clear beforehand

                vlbatasks.loadtable(uvdata, inbeamselfcalp1sntable, self.snversion)
                dualphscalp1 = vlbatasks.calibrate_target_phase_with_two_colinear_phscals(uvdata)
                dualphscalp1.compile_into_table()
                #originalinbeamselfcalp1sntable = dualphscal.copy_inbeamselfcal_sntable(inbeamselfcalp1sntable) 
                
                final_inbeamselfcal_phase_edit = inbeamselfcal_phase_time_folder + '/.corrected_phases_inbeam_selfcal.final'
                dualphscal_edit = tabledir + '/dualphscal.edit'
                if (not os.path.exists(dualphscal_edit)) and (not os.path.exists(final_inbeamselfcal_phase_edit)): ## the old interactive approach
                    print("the final saved_inbeamselfcal_phase_edit not found, now heading to interactive phase correction. When you finalize the edit, make a copy of the output file, rename it to .corrected_phases_inbeam_selfcal.final and rerun the pipeline.")
                    dualphscalp1.interactively_solve_phase_ambiguity(inbeamselfcal_phase_time_folder)
                    sys.exit(0)
                
                if os.path.exists(dualphscal_edit):
                    dualphscalp1.correcting_inbeamcalib_phase_with_dualphscal_edit(tabledir, inbeamselfcal_phase_time_folder) ## this will make/overwrite final_inbeamselfcal_phase_edit
                
                phase_correction_factor = float(dualphscal_setup[1].strip())
                dualphscalp1.load_final_inbeamselfcal_phase_edit_and_prepare_for_edit_in_AIPS(final_inbeamselfcal_phase_edit,
                                                                                            phase_correction_factor)
                dualphscaloutputsn = inbeamselfcalp1sntable.replace('.sn', '.dualphscal.sn')
                dualphscalp1.edit_AIPS_sntable_and_write_out(self.snversion, dualphscaloutputsn)
                
                vlbatasks.deletetable(uvdata, 'SN', self.snversion) ## clean up after producing dualphscaloutputsn
                ## <<<<
                    
            ## >>> apply the corrected p1.sn only to the de-facto target.\
            ## the original clversion+1 CL table will be deleted and replaced in applyinbeamcalib() !!
            junk = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, 
                                   targetconfigs, True, calonly, False, False, True,
                                   self.clversion+self.targetcl-1, self.snversion, inbeamnames, targetnames, haveungated, 
                                   ungateduvdata, dualphscal_setup, tabledir, self.inbeamfilenums)
            ## <<< 
                    
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0]) ## do\
        ## not trust this printTableAndRunlevel result if you are requesting iverse referencing!
    
    def load_inbeam_CALIB_solutions_on_separate_IFs(self, tocalnames, tocalindices, inbeamuvdatas, 
            gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata,
            tabledir):
        """
        Functionality
        -------------
        If requested, apply *pn.sn to phscal(s), inbeamcals and target(s).
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibpnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Applying inbeam CALIB pn sols")
            ## >>> dualphscal_setup==['-1','0']; sumifs==False
            sncount = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                                       targetconfigs, targetonly, calonly, False, False, False,
                                       self.clversion+self.targetcl, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                       ['-1','0'], tabledir, self.inbeamfilenums)
            ## <<<
        else:
            print("Skipping application of inbeam phase-only selfcal (separate IFs)")
            if self.maxinbeamcalibpnmins > 0:
                sncount = len(tocalnames) + 1
            else:
                sncount = 0

        if self.maxinbeamcalibpnmins > 0:
            self.snversion = self.snversion + sncount
            self.targetcl += 1
        self.runlevel = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def do_dual_phscal_calibration_correcting_the_CALIB_solutions_on_inbeams_on_separate_IFs(self, dualphscal_setup, 
            tabledir, inbeamuvdatas, gateduvdata, tocalnames, tocalindices, expconfig, targetconfigs, calonly,
            inbeamnames, targetnames, haveungated, ungateduvdata):
        """
        Functionality
        -------------
        Correct the inbeamcalibpn solutions.

        Note
        ----
        So far dualphscal_setup applies for all target groups. But it can adapt easily if necessary.
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
           int(dualphscal_setup[0].strip()) > 0 and self.maxinbeamcalibpnmins > 0:
            phase_correction_factor = float(dualphscal_setup[1].strip())
            for (tocalname, tocalindex) in zip(tocalnames, tocalindices):
                ## >>> note that this function does not work for the scenario where two inbeams are provided as primaryinbeam!
                inbeamselfcalpnsntable = tabledir + '/' + tocalname + '.icalib.pn.sn'
                ## <<<
                
                ## >>>> use the real inbeamcal data (target can work the same) as a host to produce dualphscal solutions
                filenum = self.inbeamfilenums[tocalindex]
                if filenum >= 0:
                    uvdata = inbeamuvdatas[filenum] ## uvdata here refers to the uvdata to do self-calibration
                else:
                    uvdata = gateduvdata
                for i in range(20):
                    vlbatasks.deletetable(uvdata, 'SN', self.snversion+i) ## supposed to be clear beforehand
                vlbatasks.loadtable(uvdata, inbeamselfcalpnsntable, self.snversion)
                dualphscalpn = vlbatasks.calibrate_target_phase_with_two_colinear_phscals(uvdata)
                dualphscalpn.read_inbeamselfcalpn_solutions()
                pndualphscaloutputsn = inbeamselfcalpnsntable.replace('.sn', '.dualphscal.sn')
                dualphscalpn.edit_inbeamselfcalpn_in_AIPS_and_write_out(phase_correction_factor, self.snversion, pndualphscaloutputsn)
                vlbatasks.deletetable(uvdata, 'SN', self.snversion) ## clean up after producing pndualphscaloutputsn
                ## <<<<
            
            ## >>> apply the corrected pn.sn only to the de-facto target.\
            ## the original clversion+2 CL table will be deleted and replaced in applyinbeamcalib() !!
            junk = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, 
                                   targetconfigs, True, calonly, False, False, False,
                                   self.clversion+self.targetcl-1, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                   dualphscal_setup, tabledir, self.inbeamfilenums)
            ## <<< 
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0]) ## do\
        ## not trust this printTableAndRunlevel result if you are requesting iverse referencing!
    


    def do_a_combined_IF__amp_and_phase__self_calibration_on_the_inbeams_if_requested(self,
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir,
            modeltype, targetonly, calonly, targetnames, numtargets, directory, tabledir, alwayssaved):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibap1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing amp+phase inbeam selfcal")
            tocalnames, tocalindices = self.inbeamselfcal(self.doneinbeams, self.inbeamfilenums, inbeamuvdatas, gateduvdata,
                                       expconfig, targetconfigs, modeldir, modeltype, targetonly,
                                       calonly, self.beginif, self.endif, True, False, True, self.clversion+self.targetcl, 
                                       targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved, self.leakagedopol)
        else:
            print("Skipping inbeam amp+phase selfcal (combined IFs)")
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
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])
        return tocalnames, tocalindices

    def load_inbeam_CALIB_on__amp_and_phase__with__IFs_and_pols__combined(self, tocalnames, tocalindices, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata, 
            dualphscal_setup, tabledir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibap1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Applying inbeam CALIB ap1 sols")
            sncount = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                                       targetconfigs, targetonly, calonly, True, False, True,
                                       self.clversion+self.targetcl, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                       dualphscal_setup, tabledir, self.inbeamfilenums)
        else:
            print("Skipping application of inbeam amp+phase selfcal (combined IFs)")
            if self.maxinbeamcalibap1mins > 0:
                sncount = len(tocalnames) + 1
            else:
                sncount = 0

        if self.maxinbeamcalibap1mins > 0:
            self.snversion = self.snversion + sncount
            self.targetcl += 1
        self.runlevel = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def do_a_separate_IF__amp_plus_phase__self_calibration_on_inbeams_if_requested(self, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype, targetonly, calonly, 
            targetnames, numtargets, directory, tabledir, alwayssaved):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibapnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing amp+phase inbeam selfcal")
            tocalnames, tocalindices = self.inbeamselfcal(self.doneinbeams, self.inbeamfilenums, inbeamuvdatas, gateduvdata,
                                       expconfig, targetconfigs, modeldir, modeltype, targetonly,
                                       calonly, self.beginif, self.endif, True, False, False, self.clversion+self.targetcl,
                                       targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved, self.leakagedopol)
        else:
            print("Skipping inbeam amp+phase selfcal (separate IFs)")
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
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])
        return tocalnames, tocalindices

    def load_inbeam_CALIB_solutions_on__amp_plus_phase__on_separate_IFs(self, tocalnames, tocalindices, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata, 
            dualphscal_setup, tabledir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibapnmins > 0:
            print("Runlevel " + str(self.runlevel) + ": Applying inbeam CALIB apn sols")
            sncount = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                                       targetconfigs, targetonly, calonly, True, False, False,
                                       self.clversion+self.targetcl, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                       dualphscal_setup, tabledir, self.inbeamfilenums)
        else:
            print("Skipping application of inbeam amp+phase selfcal (separate IFs)")
            if self.maxinbeamcalibapnmins > 0:
                sncount = len(tocalnames) + 1
            else:
                sncount = 0

        if self.maxinbeamcalibapnmins > 0:
            self.snversion = self.snversion + sncount
            self.targetcl += 1
        self.runlevel = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def do_a_secondary_phase_selfcal_on_inbeam_with__IFs_and_pols__combined_if_requested(self,
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype,
            targetonly, calonly, targetnames, numtargets, directory, tabledir, alwayssaved, inbeamnames):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibsp1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Doing phase-only inbeam selfcal on secondary inbeam")
            tocalnames, tocalindices = self.inbeamselfcal(self.secondaryinbeams, self.secondaryfilenums, 
                                       inbeamuvdatas, gateduvdata, expconfig, targetconfigs, 
                                       modeldir, modeltype, targetonly,
                                       calonly, self.beginif, self.endif, False, True, True, 
                                       self.clversion+self.targetcl, targetnames, numtargets, inbeamnames, directory, 
                                       tabledir, alwayssaved, self.leakagedopol)
        else:
            print("Skipping secondary inbeam phase-only selfcal (combined IFs)")
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
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])
        return tocalnames, tocalindices

    def load_secondaryinbeam_CALIB_solutions_with__IFs_and_pols__combined(self,
            tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly,
            calonly, inbeamnames, targetnames, haveungated, ungateduvdata, dualphscal_setup, tabledir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and \
            self.maxinbeamcalibsp1mins > 0:
            print("Runlevel " + str(self.runlevel) + ": Applying secondary inbeam CALIB sp1 sols")
            sncount = self.applyinbeamcalib(tocalnames, tocalindices, inbeamuvdatas, gateduvdata, expconfig,
                                       targetconfigs, targetonly, calonly, False, True, True,
                                       self.clversion+self.targetcl, self.snversion, inbeamnames, targetnames, haveungated, ungateduvdata, 
                                       dualphscal_setup, tabledir, self.inbeamfilenums)
        else:
            print("Skipping application of secondary inbeam phase-only selfcal (combined IFs)")
            if self.maxinbeamcalibsp1mins > 0:
                sncount = len(tocalnames) + 1
            else:
                sncount = 0

        if self.maxinbeamcalibsp1mins > 0:
            self.snversion = self.snversion + sncount
            self.targetcl += 1
        self.runlevel = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def calculate_the_scintillation_correction(self, numtargets, targetconfigs, tabledir, targetnames, expconfig, gateduvdata, inbeamuvdatas):
        beginif = 1
        endif = self.numifs
        scinttablepaths = []
        maxscintcorrmins = -1
        for i in range(numtargets):
            scinttablepaths.append(tabledir + targetnames[i] + '.scintcorrect.sn')
            if targetconfigs[i]['scintcorrmins'] > maxscintcorrmins:
                maxscintcorrmins = targetconfigs[i]['scintcorrmins']

        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and maxscintcorrmins > 0:
            print("Doing scintillation correction")
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
                    vlbatasks.splittoseq(gateduvdata, self.clversion+self.targetcl, 'SCINT', 
                                         targetnames[i], 1, domulti, doband, 
                                         beginif, endif, combineifs, self.leakagedopol)
                    vlbatasks.wizCorrectScint(gateduvdata, 1, self.snversion, splituvdata,
                                              targetconfigs[i]['scintcorrmins'], 
                                              scinttablepaths[i])
                    vlbatasks.plottops(gateduvdata, 'SN', self.snversion, 'AMP', self.numifs, 1, 
                                       self.numifs, pspath)
                    vlbatasks.deletetable(gateduvdata, 'SN', self.snversion)
                else:
                    os.system("rm -f " + scinttablepaths[i])
                    os.system("rm -f " + pspath)

        else:
            print("Skipping scintillation correction calculation")
        self.runlevel += 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])
        return scinttablepaths, beginif, endif

    def prepare_variables_for_final_split(self, numtargets, inbeamnames, targetconfigs, expconfig, phscalnames, targetnames,
            directory, experiment):
        """
        Note
        ----
        1. must be run before the final split and centroid-location fitting.
        """
        for i in range(20): #Clear any old FINAL split catalog entries
            for j in range(numtargets):
                for inbeamsrc in inbeamnames[j]:
                    aipssrcname = inbeamsrc
                    if len(inbeamsrc) > 12:
                        aipssrcname = inbeamsrc[:12]
                    splitdata = AIPSUVData(aipssrcname, 'FINAL', 1, 1)
                    if splitdata.exists():
                        splitdata.zap()
                    if inbeamsrc in targetconfigs[j]['primaryinbeam']:
                        self.splitdata_PS = AIPSUVData(aipssrcname, 'PRESEL', 1, 1) #pre-selfcalibration
                        if self.splitdata_PS.exists():
                            self.splitdata_PS.zap()
                ## >>> for inverse referencing where target is used as primaryinbeam
                for targetname in targetnames: 
                    if targetname.strip() == targetconfigs[j]['primaryinbeam'].strip():
                        aipstargetname = targetname
                        if len(targetname) > 12:
                            aipstargetname = targetname[:12]
                        self.splitdata_PS = AIPSUVData(aipstargetname, 'PRESEL', 1, 1) #pre-selfcalibration
                        if self.splitdata_PS.exists():
                            self.splitdata_PS.zap()
                ## <<<
                splitdata = AIPSUVData(phscalnames[j], 'FINAL', 1, 1)
                if splitdata.exists():
                    splitdata.zap()
                splitdata = AIPSUVData(targetnames[j], 'UFINL', 1, 1)
                if splitdata.exists():
                    splitdata.zap()
                splitdata = AIPSUVData(targetnames[j], 'GFINL', 1, 1)
                if splitdata.exists():
                    splitdata.zap()
        self.phscaluvfiles = []
        self.inbeamuvfiles = []
        self.gateduvfiles = []
        self.dividedinbeamuvfiles = []
        self.ungateduvfiles = []
        self.ungatedpresent = []
        self.dividesourcelist = []
        self.ibshiftdivphscaluvfiles = []
        self.inbeampreselfcaluvfiles = []
        for config in targetconfigs:
            try:
                if not config['dividesources'] == None:
                    dividesources = [x.strip(' ') for x in config['dividesources'].split(',')]
                    self.dividesourcelist.extend(dividesources)
            except KeyError:
                pass
        for phscalsrc in phscalnames:
            if len(phscalsrc) > 12:
                phscalsrc = phscalsrc[:12]
            self.phscaluvfiles.append(directory + '/' + experiment + "_" + phscalsrc + \
                                 "_pipeline_uv.fits")
            self.ibshiftdivphscaluvfiles.append(directory + '/' + experiment + '_' + phscalsrc + '_ibshiftdiv_uv.fits')
            self.inbeamuvfiles.append([])
            self.inbeampreselfcaluvfiles.append('')
            self.dividedinbeamuvfiles.append([])
        for i in range(numtargets):
            targetoutname = targetnames[i]
            gfile = directory + '/' + experiment + "_" + targetoutname + \
                    "_pipeline_uv.gated.fits"
            ufile = directory + '/' + experiment + "_" + targetoutname + \
                    "_pipeline_uv.ungated.fits"
            if targetoutname in targetconfigs[i]['primaryinbeam']:
                ifile_preselfcal = directory + '/' + experiment + "_preselfcal_" + \
                    targetoutname + "_pipeline_uv.fits"
            if expconfig['dodefaultnames']:
                gfile = directory + '/' + experiment + "_pulsar_pipeline_uv.gated.fits"
                ufile = directory + '/' + experiment + "_pulsar_pipeline_uv.ungated.fits"
                if targetoutname in targetconfigs[i]['primaryinbeam']:
                    ifile_preselfcal = directory + '/' + experiment + "_preselfcal" + \
                        "_pulsar_pipeline_uv.ungated.fits"
            self.gateduvfiles.append(gfile)
            self.ungateduvfiles.append(ufile)
            if targetoutname in targetconfigs[i]['primaryinbeam']:
                self.inbeampreselfcaluvfiles[i] = ifile_preselfcal

        for i in range(numtargets):
            icount = 0
            for inbeamoutname in inbeamnames[i]:
                ifile = directory + '/' + experiment + "_" + \
                        inbeamoutname + "_pipeline_uv.fits"
                dfile = directory + '/' + experiment + "_" + \
                        inbeamoutname + "_pipeline_divided_uv.fits"
                if inbeamoutname in targetconfigs[i]['primaryinbeam']:
                    ifile_preselfcal = directory + '/' + experiment + "_preselfcal_" + \
                        inbeamoutname + "_pipeline_uv.fits"
                if expconfig['dodefaultnames']:
                    ifile = directory + '/' + experiment + "_inbeam-" + str(i) + \
                            "_" + str(icount) + "_pipeline_uv.fits"
                    if inbeamoutname in targetconfigs[i]['primaryinbeam']:
                        ifile_preselfcal = directory + '/' + experiment + "_preselfcal_inbeam-" + str(i) + \
                                "_" + str(icount) + "_pipeline_uv.fits"
                    dfile = directory + '/' + experiment + "_inbeam-" + str(i) + \
                            "_" + str(icount) + "_pipeline_divided_uv.fits"
                self.inbeamuvfiles[i].append(ifile)
                self.dividedinbeamuvfiles[i].append(dfile)
                if inbeamoutname in targetconfigs[i]['primaryinbeam']:
                    self.inbeampreselfcaluvfiles[i] = ifile_preselfcal
                icount += 1
    
    def split__normalize_on_individual_basis__then_write_out_uvdata_for_all_sources(self,
            numtargets, targetconfigs, expconfig, targetonly, imageoutofbeam, inbeamuvdatas, bandpassclversion, ampcalsrc,
            directory, experiment, phscalnames, modeldir, inbeamnames,
            calonly, targetnames, haveungated, ungateduvdata, scinttablepaths, gateduvdata):
        """
        Note
        ----
        1. An extra model-divided fits file will be made for IBCs added to 'dividesources' of the target yaml file; a statsfile
            will be acquired for that extra fitsfile using vlbatasks.jmfit()
        """
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Runlevel " + str(self.runlevel) + ": Splitting and writing final images")
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
                        ################################################
                        ## First the ampcal source, only do it once
                        ################################################
                        if i==0:
                            splitdata = AIPSUVData(ampcalsrc, 'FINAL', 1, 1)
                            if splitdata.exists():
                                splitdata.zap()
                            if expconfig['ampcalscan'] > 0:
                                vlbatasks.splittoseq(inbeamuvdatas[0], bandpassclversion, 'FINAL', 
                                                     ampcalsrc, splitseqno, splitmulti, splitband, splitbeginif, 
                                                     splitendif, combineifs, self.leakagedopol)
                                ampcaluvfile = directory + '/' + experiment + "_" + ampcalsrc + "_pipeline_uv.fits"
                                vlbatasks.writedata(splitdata, ampcaluvfile, True)
                            #else:
                            #    vlbatasks.splittoseqnoband(inbeamuvdatas[0], bandpassclversion, 'FINAL', 
                            #                               ampcalsrc, 1)
                            #ampcaluvfile = directory + '/' + experiment + "_" + ampcalsrc + "_pipeline_uv.fits"
                            #vlbatasks.writedata(splitdata, ampcaluvfile, True)
                        
                        ################################################
                        ## Then the phase reference source
                        ################################################
                        aipssrcname = phscalnames[i]
                        if len(aipssrcname) > 12:
                            aipssrcname = aipssrcname[:12]
                        for split_phscal_option in ['FINAL','IBS']: #ibshift
                            splitdata = AIPSUVData(aipssrcname, split_phscal_option, 1, 1)
                            if splitdata.exists():
                                splitdata.zap()
                            if split_phscal_option == 'FINAL':
                                vlbatasks.splittoseq(inbeamuvdatas[0], self.clversion, split_phscal_option, aipssrcname, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                                vlbatasks.writedata(splitdata, self.phscaluvfiles[i], True)
                            else:
                                try:
                                    phscalseparateifmodel = config['phscalseparateifmodel']
                                except KeyError:
                                    phscalseparateifmodel = False
                                if phscalseparateifmodel:
                                    divideddata = self.normalise_UVData_with_separate_IF_model_and_concatenate(phscalnames[i], config, expconfig,\
                                        inbeamuvdatas[0], modeldir, self.clversion+self.targetcl)
                                else:
                                    vlbatasks.splittoseq(inbeamuvdatas[0], self.clversion+self.targetcl, split_phscal_option, aipssrcname,\
                                        splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                                    phscal_image_file = modeldir + aipssrcname + self.cmband + ".clean.fits"
                                    if not os.path.exists(phscal_image_file):
                                        print("Need a model for " + aipssrcname + " since --divideinbeammodel=True")
                                        print("But " + phscal_image_file + " was not found.")
                                        sys.exit()
                                    modeldata = AIPSImage(aipssrcname, "CLEAN", 1, 1)
                                    if modeldata.exists():
                                        modeldata.zap()
                                    vlbatasks.fitld_image(phscal_image_file, modeldata)
                                    divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                                    if divideddata.exists():
                                        divideddata.zap()
                                    vlbatasks.normaliseUVData(splitdata, modeldata,  divideddata)
                                vlbatasks.writedata(divideddata, self.ibshiftdivphscaluvfiles[i], True)
                                #plotfile = directory + '/' + experiment + '_' + aipssrcname + '.clean.ps'
                                #if not skipplots:
                                #    vlbatasks.image(splitdata, 0.5, 512, 75, 0.5, phscalnames[i], plotfile, False,
                                #                    fullauto, stokesi)
                    ################################################
                    ## Then the inbeams
                    ################################################
                    for inbeamsrc in inbeamnames[i]:
                        aipssrcname = inbeamsrc
                        if len(inbeamsrc) > 12:
                            aipssrcname = inbeamsrc[:12]
                        splitdata = AIPSUVData(aipssrcname, 'FINAL', 1, 1)
                        if splitdata.exists():
                            splitdata.zap()
                        vlbatasks.splittoseq(inbeamuvdatas[count], self.clversion+self.targetcl, 'FINAL', inbeamsrc, 
                                             splitseqno, splitmulti, splitband, splitbeginif, splitendif, 
                                             combineifs, self.leakagedopol)
                        #plotfile = directory + '/' + experiment + '_' + aipssrcname + '.clean.ps'
                        #if not skipplots:
                        #    vlbatasks.image(splitdata, 0.5, 512, 75, 0.5, inbeamsrc, plotfile, False,
                        #                    fullauto, stokesi)
                        vlbatasks.writedata(splitdata, self.inbeamuvfiles[i][count], True)
                        #split primary in-beam calibrator referenced to phscal (pre-inbeamselfcal)
                        if inbeamsrc in config['primaryinbeam']:
                            vlbatasks.splittoseq(inbeamuvdatas[count], self.clversion, 'PRESEL', inbeamsrc, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                            #write out pre-selfcal in-beam cal data
                            tempfile = directory + '/temp.fits'
                            os.system("rm -f " + tempfile)
                            #if os.path.exists(temp_inbeampreselfcaluvfile):
                            #    os.remove(temp_inbeampreselfcaluvfile)
                            #print(inbeamuvdatas[count], inbeamsrc, splitseqno, splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                            vlbatasks.writedata(self.splitdata_PS, tempfile, True)
                            os.system("mv -f " + tempfile + " " + self.inbeampreselfcaluvfiles[i])
                            #os.system("mv -f " + tempdivfile + " " + self.inbeampreselfcaluvfiles[i][count])
                        #divide inbeamuvfiles by model
                        tempdivfile = directory + '/temp.fits'
                        if inbeamsrc in self.dividesourcelist:
                            if config['separateifmodel'] and inbeamsrc in config['primaryinbeam']:
                                ifdivideddatas = []
                                for j in range(beginif,endif+1):
                                    inbeam_image_file = modeldir + inbeamsrc + ".IF" + str(j) + self.cmband + ".clean.fits"
                                    if not os.path.exists(inbeam_image_file):
                                        print("Need a frequency-dependent model for " + inbeamsrc + " since --divideinbeammodel=True and separateifmodel=True")
                                        print("But " + inbeam_image_file + " was not found.")
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
                                inbeam_image_file = modeldir + inbeamsrc + self.cmband + ".clean.fits"
                                if not os.path.exists(inbeam_image_file):
                                    print("Need a model for " + inbeamsrc + " since --divideinbeammodel=True")
                                    print("But " + inbeam_image_file + " was not found.")
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
                            os.system("mv -f " + tempdivfile + " " + self.dividedinbeamuvfiles[i][count])
                        #divide pre-selfcal primary inbeam uvfile by model
                        """
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
                            vlbatasks.normaliseUVData(self.splitdata_PS, modeldata,  divideddata_PS)
                            #os.system("rm -f " + tempdivfile)
                            vlbatasks.writedata(divideddata_PS, tempdivfile, True)
                            os.system("mv -f " + tempdivfile + " " + self.inbeampreselfcaluvfiles[i][count])
                        """
                        count += 1
                ################################################
                ## finally gateduvdata and ungateduvdata
                ################################################
                if not calonly:
                    splitdata1 = AIPSUVData(targetnames[i], 'UFINL', 1, 1)
                    if splitdata1.exists():
                        splitdata1.zap()
                    if haveungated:
                        try:
                            vlbatasks.splittoseq(ungateduvdata, self.clversion+self.targetcl, 'UFINL', targetnames[i], 
                                                 splitseqno, splitmulti, splitband, splitbeginif, splitendif, 
                                                 combineifs, self.leakagedopol)
                            if os.path.exists(scinttablepaths[i]) and config['scintcorrmins'] > 0:
                                vlbatasks.loadtable(splitdata1, scinttablepaths[i], 1)
                                splitdataS = AIPSUVData(targetnames[i], 'UFINS', 1, 1)
                                if splitdataS.exists():
                                    splitdataS.zap()
                                vlbatasks.splat(splitdata1, 1, [0,0,0,0,0,0,0,0], splitdataS, 1)
                                splitdata1.zap()
                                splitdata1 = splitdataS
                            self.ungatedpresent.append(True)
                        except RuntimeError:
                            print("Guess there was no ungated for " + targetnames[i])
                            self.ungatedpresent.append(False)
                    splitdata2 = AIPSUVData(targetnames[i], 'GFINL', 1, 1)
                    if splitdata2.exists():
                        splitdata2.zap()
                    vlbatasks.splittoseq(gateduvdata, self.clversion+self.targetcl, 'GFINL', targetnames[i], splitseqno,
                                         splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                    ## >>> inverse referencing
                    if targetnames[i].strip() == config['primaryinbeam'].strip():
                        #splitdata_PS = AIPSUVData(targetnames[i], 'PRESEL', 1, 1) #pre-selfcalibration
                        if self.splitdata_PS.exists():
                            self.splitdata_PS.zap()
                        vlbatasks.splittoseq(gateduvdata, self.clversion, 'PRESEL', targetnames[i], splitseqno,
                             splitmulti, splitband, splitbeginif, splitendif, combineifs, self.leakagedopol)
                        vlbatasks.writedata(self.splitdata_PS, self.inbeampreselfcaluvfiles[i], True)
                    ## <<<
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
                    if haveungated and self.ungatedpresent[i]:
                        vlbatasks.writedata(splitdata1, self.ungateduvfiles[i], True)
                    vlbatasks.writedata(splitdata2, self.gateduvfiles[i], True)
        else:
            print("Skipping splitting/writing of final images")
            for i in range(numtargets):
                if os.path.exists(self.ungateduvfiles[i]):
                    self.ungatedpresent.append(True)
                else:
                    self.ungatedpresent.append(False)
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def image_targets_using_DIFMAP_and_fit_for_position(self, calonly, numtargets, targetconfigs, expconfig,
            directory, experiment, targetnames, beginif, endif, haveungated, phscalnames, inbeamnames, inbeamuvdatas,
            uvtaperstring, difmaptargetuvaverstring):
        #fullauto  = True
        stokesi   = True
        #gaussiantarget = False
        #gaussianinbeam = True
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel and not calonly:
            print("Runlevel " + str(self.runlevel) + ": Fitting target positions using difmap")
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
                ## >>> gated data first
                targetimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                                  '_difmap.gated.fits'
                jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                            '.gated.difmap.jmfit'
                if expconfig['dodefaultnames']:
                    targetimagefile = directory + '/' + experiment + '_pulsar' + \
                                      '_difmap.gated.fits'
                    jmfitfile = directory + '/' + experiment + '_pulsar' + \
                                    '.gated.difmap.jmfit'
                vlbatasks.difmap_maptarget(self.gateduvfiles[i], targetimagefile, fullauto, stokesi,
                                           config['difmappixelmas'], config['difmapnpixels'],
                                           config['difmapweightstring'], difmaptargetuvaverstring, 
                                           uvtaperstring, config['usegaussiantarget'],
                                           beginif, endif-subtractif)
                vlbatasks.jmfit(targetimagefile, jmfitfile, targetnames[i], stokesi, endif-subtractif)
                ## <<< gated data first

                ## >>> preselfcal gated in the case of inverse referencing
                if targetnames[i].strip() == config['primaryinbeam'].strip():
                    targetimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                        '_preselfcal.difmap.gated.fits'
                    jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                        '_preselfcal.gated.difmap.jmfit'
                    if expconfig['dodefaultnames']:
                        targetimagefile = directory + '/' + experiment + '_pulsar' + \
                            '_preselfcal.difmap.gated.fits'
                        jmfitfile = directory + '/' + experiment + '_pulsar_preselfcal' + \
                            '.gated.difmap.jmfit'
                    vlbatasks.difmap_maptarget(self.inbeampreselfcaluvfiles[i], targetimagefile, fullauto, stokesi,
                                               config['difmappixelmas'], config['difmapnpixels'],
                                               config['difmapweightstring'], difmaptargetuvaverstring,
                                               uvtaperstring, config['usegaussiantarget'],
                                               beginif, endif-subtractif)
                    vlbatasks.jmfit(targetimagefile, jmfitfile, targetnames[i], stokesi, endif-subtractif)
                ## <<< preselfcal gated in the case of inverse referencing

                ## >>> ungated data
                targetimagefile = directory + '/' + experiment + '_' + targetnames[i] + \
                                  '_difmap.ungated.fits'
                jmfitfile = directory + '/' + experiment + '_' + targetnames[i] + \
                            '.ungated.difmap.jmfit'
                if expconfig['dodefaultnames']:
                    targetimagefile = directory + '/' + experiment + '_pulsar' + \
                                      '_difmap.ungated.fits'
                    jmfitfile = directory + '/' + experiment + '_pulsar' + \
                                    '.ungated.difmap.jmfit'
                if haveungated and self.ungatedpresent[i]:
                    vlbatasks.difmap_maptarget(self.ungateduvfiles[i], targetimagefile, fullauto, 
                                               stokesi, config['difmappixelmas'], config['difmapnpixels'], 
                                               config['difmapweightstring'], difmaptargetuvaverstring, 
                                               uvtaperstring, config['usegaussiantarget'], 
                                               beginif, endif-subtractif)
                    vlbatasks.jmfit(targetimagefile, jmfitfile, targetnames[i], stokesi, endif-subtractif)
                ## <<< ungated data

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
                #vlbatasks.difmap_maptarget(self.ibshiftdivphscaluvfiles[i], targetimagefile, fullauto, stokesi, config['difmappixelmas'], config['difmapnpixels'], config['difmapweightstring'], config['usegaussiantarget'], beginif, endif-subtractif)
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
                    vlbatasks.difmap_maptarget(self.inbeamuvfiles[i][j], targetimagefile, 
                                               fullauto, stokesi, config['difmappixelmas'],
                                               config['difmapnpixels'], config['difmapweightstring'], '20, False',
                                               uvtaperstring, config['usegaussianinbeam'], 
                                               beginif, endif-subtractif)
                    vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                                    stokesi, endif-subtractif)
                    if inbeamnames[i][j] in self.dividesourcelist:
                        targetimagefile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                                          '_divided_difmap.fits'
                        jmfitfile = directory + '/' + experiment + '_' + inbeamnames[i][j] + \
                                    '_divided.difmap.jmfit'
                        if expconfig['dodefaultnames']:
                            targetimagefile = directory + '/' + experiment + '_inbeam' + str(i) + \
                                              "_" + str(icount) + '_pipeline_divided_uv.fits'
                            jmfitfile = directory + '/' + experiment + '_inbeam' + str(i) + \
                                        "_" + str(icount) + '_divided.difmap.jmfit'
                        vlbatasks.difmap_maptarget(self.dividedinbeamuvfiles[i][j], targetimagefile, 
                                                   fullauto, stokesi, config['difmappixelmas'], 
                                                   config['difmapnpixels'], config['difmapweightstring'], '20, False', 
                                                   uvtaperstring, config['usegaussianinbeam'],
                                                   beginif, endif-subtractif)
                        vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                                        stokesi, endif-subtractif)
                    if inbeamnames[i][j] in config['primaryinbeam']:
                        #targetimagefile = directory + '/' + experiment + '_' + inbeamnames[i][j] + '_preselfcal.divided_uv.fits'
                        #vlbatasks.difmap_maptarget(self.inbeampreselfcaluvfiles[i][j], targetimagefile, 
                        #                           fullauto, stokesi, config['difmappixelmas'],
                        #                           config['difmapnpixels'], config['difmapweightstring'],
                        #                           config['usegaussianinbeam'], beginif,
                        #                           endif-subtractif)
                        #vlbatasks.jmfit(targetimagefile, jmfitfile, inbeamnames[i][j], 
                        #                stokesi, endif-subtractif)
                        aipssrcname = inbeamnames[i][j]
                        if len(aipssrcname) > 12:
                            aipssrcname = aipssrcname[:12]
                        jmfitfile = directory + '/' + experiment + '_' + aipssrcname + '_preselfcal.difmap.jmfit' 
                        preselfcalinbeamimage = AIPSImage(aipssrcname, "ICL001", 1, 1)
                        if preselfcalinbeamimage.exists():
                            preselfcalinbeamimage.zap()
                        #divideddata = AIPSUVData(aipssrcname, 'DIV', 1, 1)
                        vlbatasks.widefieldimage(self.splitdata_PS, aipssrcname, 256, 0.75, True, 0.050,
                                   0, 0, 0, 100, 20)
                        vlbatasks.nonpulsarjmfit("", jmfitfile, aipssrcname, -1, -1, True,
                                           False, preselfcalinbeamimage, 48) ## imagedata == loadedfile == preselfcalinbeamimage

                    icount += 1
        else:
            print("Skipping imaging/position fitting of target")
        self.runlevel  = self.runlevel + 1
        self.printTableAndRunlevel(self.runlevel, self.snversion, self.clversion+self.targetcl, inbeamuvdatas[0])

    def make_diagnostic_plots(self, directory, codedir):
        if self.runfromlevel <= self.runlevel and self.runtolevel >= self.runlevel:
            print("Making final diagnostic plots")
            os.chdir(directory)
            os.system("%s/make_final_diagnostic.py" % codedir)
        else:
            print("Skipping making of diagnostic plots")
