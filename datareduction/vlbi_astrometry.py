#!/usr/bin/env ParselTongue
__doc__ = """
Functionality
-------------
A ParselTongue script for phase referencing VLBA data and getting positions 
of pulsars and inbeam calibrators.

Code structure
--------------
vlbi_astrometry.py :
    the user-level code, which contains the main() data reduction procedure provoking
    functions in vlbireduce.py.
vlbireduce.py :
    to be described
support_vlbireduce.py :
    to be described
vlbatasks.py :
    repetitively used basic ParselTongue scripts.

Versions
--------
initial version :
    finished at 08 Jan 2013 by Adam Deller, fully functional for the data
    reduction of the PSRPI project.
new versions : 
    have been maintained by Hao Ding since Nov. 2018. The changes include         
    re-organization of the code into object-baesd form and added new features
    (e.g. the dual-phscal function) since Nov 2018.

Commenting syntax
-----------------
Apart from the docstrings provided below the class/function titles, rich comments
are offered to help users understand and develop upon this code. All comments
will be gradually replaced by three major formats.
1. same-line comments :
    comments given right after a line of code, starting with '## '.
2. multi-line comments :
    comments given in two lines singling out a paragraph of code to be commented
    on. The top line of such comments begins with '## >>> ', while the bottom one
    stays as '## <<< ' to mark the end of the paragraph. If the singled-out
    paragraph is long, the comments in the top line will be repeated in the bottom
    comment. The indentation shall stay the same as the paragraph to be commented on.
    Furthermore, the '## >>>' needs to completed with '## <<<'. '## >>>>' is expected
    to be completed with '## <<<<'.
3. segment comments :
    the longer part of code is divided according its function into several segments.
    Such segments start with the following comments.
    #######################################
    ## SEGMENT COMMENTS
    #######################################
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
from vlbireduce import vlbireduce
from time import gmtime, strftime
from optparse import OptionParser
warnings.defaultaction = "always"

################################################################################
## Little logger class to put print statements to the log file
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
    def flush(self):
        pass

def main():
    """
    main program
    """
    ################################################################################
    ## Option parsing and defaulted global variables
    ################################################################################
    try:
        aipsver = os.environ['PSRVLBAIPSVER']
    except KeyError:
        try:
            aipsver = os.environ['AIPS_VERSION'].split('/')[-1]
        except KeyError:
            aipsver = '31DEC20'
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
    parser.add_option("-k", "--skipdiagnosticplots", dest="skipdiagnosticplots", default=False,
                      action="store_true", help="Do not make diagnostic plots")
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
    skipdiagnosticplots = options.skipdiagnosticplots
    runfromlevel    = int(runsplit[0])
    runtolevel      = 999
    if len(runsplit) > 1:
        runtolevel = int(runsplit[1])


    ################################################################################
    ## Check validity of inputs, load up the config files
    ################################################################################
    if experiment == "":
        parser.error("You must supply an experiment name")
    if runfromlevel > runtolevel:
        parser.error("Runtolevel (" + str(runtolevel) + ") must be greater " \
                     "than or equal to runfromlevel (" + str(runfromlevel) + ")")
    try:
        auxdir = os.environ['PSRVLBAUXDIR']
    except KeyError:
        print("PSRVLBAUXDIR is not defined - aborting!")
        sys.exit(1)
    try:
        codedir = os.environ['PSRVLBICODEDIR']
    except KeyError:
        print("PSRVLBICODEDIR is not defined - aborting!")
        sys.exit(1)
    configdir     = auxdir + '/configs/'
    finalmodeldir = auxdir + '/sourcemodels/final/'
    prelimmodeldir= auxdir + '/sourcemodels/preliminary/'
    expconfigfile = configdir + experiment + '.yaml'
    if not os.path.exists(expconfigfile):
        parser.error("Experiment config file %s does not exist!" % expconfigfile)
    expconfig     = yaml.safe_load(open(expconfigfile))
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
    #clversion     = 1
    #snversion     = 1
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
        #targetconfigfile = configdir + expconfig['rootdir'].split('/')[-1] + '.yaml'
        if not os.path.exists(targetconfigfile):
            parser.error("Target config file %s does not exist!" % targetconfigfile)
        targetconfigs.append(yaml.safe_load(open(targetconfigfile)))
    skiplastif = targetconfigs[0]['skiplastif']
    if numtargets > 1:
        for config in targetconfigs[1:]:
            if skiplastif != targetconfigs[i]['skiplastif']:
                print("All targets must have the same value for skiplastif! Aborting")
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
        plotbandpass_setup = expconfig['plotbandpass']
        plotbandpass = True
        plotbandpass_IFs = plotbandpass_setup[0:2]
        plotbandpass_chans = plotbandpass_setup[2:4]
    except KeyError:
        plotbandpass = False
        plotbandpass_IFs = None
        plotbandpass_chans = None
    try:
        uvtaperstring = expconfig['difmaptargetuvtaperstring']
    except KeyError:
        uvtaperstring = '0.99,999'
    try:
        dualphscal_setup = targetconfigs[0]['dualphscal'].split(',')
    except KeyError:
        dualphscal_setup = ['-1','0']
    try:
        secondary_dualphscal_setup = targetconfigs[0]['secondarydualphscal'].split(',')
    except KeyError:
        secondary_dualphscal_setup = ['-1','0']

    try:
        triphscal_setup = targetconfigs[0]['triphscal'].split(';')
        dotriphscal = True
    except KeyError:
        dotriphscal = False

    try:
        difmaptargetuvaverstring = expconfig['difmaptargetuvaverstring']
    except KeyError:
        difmaptargetuvaverstring = '20,false'

    if float(dualphscal_setup[0]) > 0 and dotriphscal:
        print("Can't do both dualphscal and triphscal; aborting...")
        sys.exit()


    ################################################################################
    ## get an instance of the vlbireduce class
    ################################################################################
    reducevlbi = vlbireduce(runfromlevel, runtolevel)


    ################################################################################
    ## Parse the source file and set some more variables
    ################################################################################
    try:
        gateduvfile   = expconfig['gateduvfile']
        if gateduvfile[0] != '/':
            gateduvfile = directory + gateduvfile
        ungateduvfile = expconfig['ungateduvfile']
        if ungateduvfile[0] != '/':
            ungateduvfile = directory + ungateduvfile
        inbeamfiles   = expconfig['inbeamfiles']
        if type(inbeamfiles) != list:
            inbeamfiles = [inbeamfiles]
        for inbeamfile in inbeamfiles:
            if inbeamfile[0] != '/':
                inbeamfile = directory + inbeamfile
        numinbeams    = len(inbeamfiles)
        targetnames   = expconfig['targets']
        if type(targetnames) != list:
            targetnames = [targetnames]
        inbeamnames   = expconfig['inbeamnames']
        if type(inbeamnames) != list:
            inbeamnames = [[inbeamnames]]
        if type(inbeamnames[0]) != list:
            inbeamnames = [inbeamnames]
        phscalnames   = expconfig['phscalnames']
        if type(phscalnames) != list:
            phscalnames = [phscalnames]
        ampcalsrc     = expconfig['ampcalsrc'] ## here we should consider to generalize to ampcalsrcs later
    except KeyError:
    ## the deprecated source file setup that would still be in use until old data are published
        print("some of the expconfig variables are missing, try to find it in source file as an alternative")
        sourcefile = directory + experiment.lower() + ".source"
        if not os.path.exists(sourcefile):
            print(("%s does not exist; aborting" % sourcefile))
            sys.exit()
        else:
            gateduvfile, ungateduvfile, numinbeams, inbeamfiles, inbeamuvdatas, \
            targetnames, inbeamnames, phscalnames, ampcalsrc = reducevlbi.parsesourcefile(sourcefile, experiment, klass, uvsequence)

    ## the deprecated part ends here
    reducevlbi.cmband = ""
    try:
        reducevlbi.cmband = "." + str(expconfig['cmband']).strip() + "cm"
    except KeyError:
        pass
    
    haveungated = True
    gateduvdata    = AIPSUVData(experiment.upper() + "_G", klass, 1, uvsequence)
    ungateduvdata  = AIPSUVData(experiment.upper() + "_U", klass, 1, uvsequence)
    alluvdatas     = [gateduvdata, ungateduvdata]
    inbeamuvdatas  = []
    for i in range(numinbeams):
        inbeamuvdatas.append(AIPSUVData(experiment.upper() + "_I" + str(i+1), klass, 1, uvsequence))
    if ungateduvfile == "":
        haveungated = False
        alluvdatas = [gateduvdata]
    if calonly:
        alluvdatas = []
    for uvdata in inbeamuvdatas:
        alluvdatas.append(uvdata)

    ################################################################################
    ## Save a record of the time and the arguments of this run, and start the log
    ## Echo inputs from yaml files to the log
    ################################################################################
    AIPS.log = open(logfile, logmode)
    sys.stdout = Logger(AIPS.log)
    logout = open(directory + 'runlog.txt', 'a')
    logout.write(100*'=' + '\n' + 100*'=' + '\n')
    logout.write(strftime("%a, %d %b %Y %H:%M:%S +0000\n", gmtime()))
    for a in sys.argv:
        logout.write(a + ' ')
    logout.write('\n\n****** %s.yaml inputs: ******\n' % experiment)
    readfile = open(expconfigfile)
    lines = readfile.readlines()
    readfile.close()
    logout.writelines(lines)
    for targetname in targetnames:
        logout.write('\n\n****** %s.yaml inputs: ******\n' % targetname)
        readfile = open(configdir + '/' + targetname + '.yaml')
        lines = readfile.readlines()
        readfile.close()
        logout.writelines(lines)
    logout.write('\n\n')
    logout.close()

    ################################################################################
    ## Zap the existing cal tables if requested
    ################################################################################
    if zapallcaltables and runfromlevel > 1:
        print("Zapping all SN, BP, and CL tables (except CL 1)!")
        for uvdata in alluvdatas:
            tables = uvdata.tables
            for table in tables:
                if table[1][-2:] == "SN":
                    print("Zapping SN " + str(table[0]) + " from " + uvdata.name)
                    uvdata.table('SN', table[0]).zap()
                elif table[1][-2:] == "BP":
                    print("Zapping BP " + str(table[0]) + " from " + uvdata.name)
                    uvdata.table('BP', table[0]).zap()
                elif table[1][-2:] == "CL" and table[0] > 1:
                    print("Zapping CL " + str(table[0]) + " from " + uvdata.name)
                    uvdata.table('CL', table[0]).zap()

    ################################################################################
    ## Begin the actual data processing
    ################################################################################

    ## Load the uv data ############################################################
    reducevlbi.load_uv_data(directory, experiment, expconfig, numinbeams, inbeamfiles, 
            inbeamuvdatas, calonly, targetonly, gateduvdata, ungateduvdata, targetnames, 
            haveungated, ungateduvfile, gateduvfile)
    ## Increase the number of IFs if requested #####################################
    reducevlbi.increase_the_number_of_IFs_if_requested(targetonly, numinbeams, experiment, uvsequence, inbeamuvdatas,
            klass, calonly, gateduvdata, haveungated, ungateduvdata, expconfig)
    ## Unflag Pie Town if requested ################################################
    reducevlbi.unflag_Pie_Town_if_requested(expconfig, targetonly, numinbeams, inbeamuvdatas, calonly,
                ungateduvdata, gateduvdata, haveungated)
    ## Load the user flags (if any) ################################################
    reducevlbi.load_the_user_flags_if_any(tabledir, targetonly, numinbeams, inbeamuvdatas, numtargets, 
            inbeamnames, calonly, ungateduvdata, gateduvdata, haveungated, targetnames)
    ## Run the flagging for zero-fringe rate effects ###############################
    reducevlbi.run_the_flagging_for_zero_fringe_rate_effects(targetconfigs, targetonly, inbeamuvdatas,
            haveungated, ungateduvdata, gateduvdata, calonly)
    ## Run the autoflagger (if requested) ##########################################
    reducevlbi.run_the_autoflagger_if_requested(expconfig, tabledir, targetonly, numinbeams,
            inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated)
    ## Duplicate the initial CL table ##############################################
    reducevlbi.duplicate_the_initial_CL_table(targetonly, numinbeams, inbeamuvdatas, gateduvdata,
                ungateduvdata, calonly, haveungated)
    ## Correct positions ###########################################################
    reducevlbi.correct_positions(alluvdatas, targetonly, gateduvdata, haveungated, 
                ungateduvdata, calonly, inbeamuvdatas, expconfig)
    ## Run TECOR ###################################################################
    reducevlbi.run_TECOR(expconfig, targetonly, numinbeams, inbeamuvdatas, logdir,
                calonly, gateduvdata, ungateduvdata, haveungated)
    ## Run CLCOR to correct EOPs ###################################################
    reducevlbi.run_CLCOR_to_correct_EOPs(targetonly, numinbeams, inbeamuvdatas, logdir,
                calonly, gateduvdata, ungateduvdata, haveungated)
    ## prepare leakage-related variables ############################################
    reducevlbi.prepare_leakage_related_variables(expconfig, ampcalsrc)
    ## Run CLCOR to correct PANG ###################################################
    reducevlbi.run_CLCOR_to_correct_PANG(targetonly, numinbeams,
                inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated)
    ## Correct for a priori delays if available using CLCOR ########################
    reducevlbi.correct_for_a_priori_delays_if_available_using_CLCOR(tabledir, alluvdatas,
                targetonly, gateduvdata, ungateduvdata, haveungated, calonly, inbeamuvdatas)
    ## Do the amplitude calibration (either load existing table or do and then edit) and inspect
    reducevlbi.do_the_amplitude_calibration_with_existing_table_or_in_an_interative_way__and_inspect(expconfig,
            targetonly, tabledir, alluvdatas, gateduvdata, ungateduvdata, haveungated, calonly, inbeamuvdatas, logdir, experiment)
    ## Correct for primary beam attenuation ########################################
    reducevlbi.correct_for_primary_beam_attenuation(expconfig, directory, experiment,
                ampcalsrc, targetonly, numinbeams, numtargets, phscalnames, inbeamnames, tabledir, inbeamuvdatas, calonly,
                targetnames, gateduvdata, ungateduvdata, haveungated)
    ## Do PCAL correction and inspect ##############################################
    reducevlbi.do_PCAL_correction_and_inspect(expconfig, targetonly, tabledir,
                inbeamuvdatas, ampcalsrc, calonly, gateduvdata, ungateduvdata, haveungated, numinbeams)
    ## Run FRING (bandpass calibrator) #############################################
    reducevlbi.run_FRING_with_fringe_finder(targetonly, expconfig, tabledir,
                modeldir, ampcalsrc, modeltype, inbeamuvdatas)
    ## Copy the FRING SN table around and apply it #################################
    reducevlbi.copy_the_FRING_SN_table_around_and_apply_it(expconfig, targetonly,
                tabledir, numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated)
    ## Run xpoldelaycal ############################################################
    xpolsnfilename = reducevlbi.run_xpoldelaycal(tabledir, targetonly, expconfig, modeldir, inbeamuvdatas)
    ## Load and apply xpoldelaycal #################################################
    reducevlbi.load_and_apply_xpoldelaycal(targetonly, numinbeams,
                inbeamuvdatas, xpolsnfilename, expconfig, gateduvdata, ungateduvdata, haveungated, calonly)
    ## Run leakagecal on the leakage calibrator and save AN file ###################
    leakagefilename = reducevlbi.run_leakagecal_on_the_leakage_calibrator_and_save_AN_file(targetonly,
            tabledir, modeldir, inbeamuvdatas, expconfig, directory, experiment)
    ## Load up the leakage-calibrated AN table #####################################
    reducevlbi.load_up_the_leakage_calibrated_AN_table(targetonly, xpolsnfilename, numinbeams,
                inbeamuvdatas, leakagefilename, calonly, gateduvdata, ungateduvdata, haveungated)
        # Remember to set dopol=2 on everything hereon!!!!

    ## Run BPASS ###################################################################
    reducevlbi.run_BPASS(targetonly, expconfig, tabledir, modeldir, 
            ampcalsrc, modeltype, inbeamuvdatas)
    ## Load BPASS ##################################################################
    bandpassclversion = reducevlbi.load_BPASS_solutions(expconfig, tabledir, 
            calonly, gateduvdata, ungateduvdata, targetonly, numinbeams, inbeamuvdatas, haveungated)
    ## Plot bandpass ###############################################################
    reducevlbi.plot_bandpass(plotbandpass, directory, 
            inbeamuvdatas, plotbandpass_IFs, plotbandpass_chans)
    ## prepare variables for FRING (phase reference calibrators) ####################
    reducevlbi.prepare_variables_for_calibration_on_phase_calibrators(targetconfigs, phscalnames)

    ## Run FRING (phase reference calibrators) #####################################
    dophscaldump = reducevlbi.run_FRING_with_phase_reference_calibrators(targetonly,
            tabledir, modeldir, modeltype, expconfig, experiment, inbeamuvdatas)
    ## Copy the phsref FRING SN table around and apply it ##########################
    reducevlbi.copy_the_phsref_FRING_SN_table_around_and_apply_it(targetonly,
            tabledir, expconfig, numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated,
            directory, experiment, dophscaldump)
    ## Run phase CALIB on the phase reference sources ##############################
    reducevlbi.run_phase_CALIB_on_the_phase_reference_sources(tabledir, targetonly, inbeamuvdatas, expconfig, modeltype, modeldir)
    ## Load all the phase CALIB solutions #########################################
    reducevlbi.load_all_the_phase_CALIB_solutions_obtained_with_phase_calibrators(phscalnames, 
            targetonly, numinbeams, inbeamuvdatas, calonly, gateduvdata, haveungated, ungateduvdata, tabledir, expconfig)
    ## Run amp CALIB on the phase reference sources ################################
    reducevlbi.run_amp_CALIB_on_the_phase_reference_sources(targetconfigs,
            tabledir, targetonly, expconfig, inbeamuvdatas, modeldir, modeltype, phscalnames)
    ## Load all the amp CALIB solutions ###########################################
    reducevlbi.load_all_the_amp_CALIB_solutions_obtained_with_phase_calibrators(
            phscalnames, targetonly, numinbeams, inbeamuvdatas, calonly, gateduvdata, ungateduvdata, haveungated,
            tabledir, expconfig)
    ## prepare some variables for operations on inbeam calibrators ###########################    
    reducevlbi.prepare_variables_for_inbeamselfcal(inbeamuvdatas, targetconfigs, numtargets, inbeamnames, targetnames)
    ## Generate a raw (no phase selfcal) inbeam dataset if requested ##############
    reducevlbi.generate_a_raw_inbeam_dataset_without_phase_selfcal_on_itself_if_requested(expconfig,
                targetonly, numtargets, targetconfigs, inbeamnames, directory, experiment, inbeamuvdatas)
    ## Do a combined IF (and pol) phase selfcal on the inbeams if requested #################
    [tocalnames, tocalindices] = reducevlbi.do_a_combined__IF_and_pol__phase_selfcal_on_the_inbeams_if_requested(
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs,
            modeldir, modeltype, targetonly, calonly, targetnames, numtargets, inbeamnames, directory, tabledir, alwayssaved)
    ## Load the p1 inbeam CALIB solutions ########################################
    reducevlbi.load_inbeam_CALIB_solutions_obtained_with__IF_and_pol__combined(tocalnames,
            tocalindices, inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, 
            ungateduvdata, tabledir)
    ## Do a separate IF phase selfcal on the inbeams if requested #################
    [tocalnames, tocalindices] = reducevlbi.do_a_separate_IF_phase_selfcal_on_the_inbeams_if_requested(
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype,
            targetonly, calonly, targetnames, numtargets, directory, tabledir, alwayssaved, inbeamnames)
    ## Do dual-phscal calibration if requested: 1). correct INBEAM.icalib.p1.sn ###################################
    reducevlbi.do_dual_phscal_calibration_correcting_the_CALIB_solutions_on_inbeams_with__IF_and_pol__combined(dualphscal_setup, directory,
            tabledir, inbeamuvdatas, gateduvdata, ungateduvdata, targetonly, calonly, haveungated, tocalnames, tocalindices, expconfig, 
            targetconfigs, inbeamnames, targetnames)
    ## Load the pn inbeam CALIB pn solutions ########################################################################
    reducevlbi.load_inbeam_CALIB_solutions_on_separate_IFs(tocalnames, tocalindices, inbeamuvdatas, 
            gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata, 
            tabledir)
    ## Do dual-phscal calibration if requested: 2). correct INBEAM.icalib.pn.sn ###################################
    reducevlbi.do_dual_phscal_calibration_correcting_the_CALIB_solutions_on_inbeams_on_separate_IFs(dualphscal_setup, 
            tabledir, inbeamuvdatas, gateduvdata, tocalnames, tocalindices, expconfig, targetconfigs, calonly,
            inbeamnames, targetnames, haveungated, ungateduvdata)
    ## Do a combined IF amp + phase selfcal on the inbeams if requested ###########################################
    [tocalnames, tocalindices] = reducevlbi.do_a_combined_IF__amp_and_phase__self_calibration_on_the_inbeams_if_requested(
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype, targetonly, calonly, 
            targetnames, numtargets, directory, tabledir, alwayssaved)
    ## Load amp+p1 inbeam CALIB solutions ########################################################################
    reducevlbi.load_inbeam_CALIB_on__amp_and_phase__with__IFs_and_pols__combined(tocalnames, 
            tocalindices, inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated,
            ungateduvdata, tabledir)
    ## Do a separate IF amp + phase selfcal on the inbeams if requested ###########################################
    [tocalnames, tocalindices] = reducevlbi.do_a_separate_IF__amp_plus_phase__self_calibration_on_inbeams_if_requested(
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype, targetonly, calonly, 
            targetnames, numtargets, directory, tabledir, alwayssaved)
    ## Load amp+pn inbeam CALIB solutions ########################################################################
    reducevlbi.load_inbeam_CALIB_solutions_on__amp_plus_phase__on_separate_IFs(tocalnames, tocalindices, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata, tabledir)
    ## Do a secondary phase selfcal on the inbeam(s) if requested #################################################
    [tocalnames, tocalindices] = reducevlbi.do_a_secondary_phase_selfcal_on_inbeam_with__IFs_and_pols__combined_if_requested(
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, modeldir, modeltype,
            targetonly, calonly, targetnames, numtargets, directory, tabledir, alwayssaved, inbeamnames)
    ## Load the secondary inbeam CALIB solutions ########################################################################
    reducevlbi.load_secondaryinbeam_CALIB_solutions_with__IFs_and_pols__combined(tocalnames, tocalindices, 
            inbeamuvdatas, gateduvdata, expconfig, targetconfigs, targetonly, calonly, inbeamnames, targetnames, haveungated, ungateduvdata, tabledir)
    ## Do dual-phscal calibration on the prIBC-secIBC line if requested: correct secIBC.icalib.sp1.sn #############
    reducevlbi.do_dual_phscal_calibration_correcting_the_CALIB_solutions_on_inbeams_with__IF_and_pol__combined(secondary_dualphscal_setup, directory,
            tabledir, inbeamuvdatas, gateduvdata, ungateduvdata, targetonly, calonly, haveungated, tocalnames, tocalindices, expconfig, 
            targetconfigs, inbeamnames, targetnames)

    ## Calculate the scintillation correction #####################################################################
    [scinttablepaths, beginif, endif] = reducevlbi.calculate_the_scintillation_correction(numtargets, targetconfigs,
                tabledir, targetnames, expconfig, gateduvdata, inbeamuvdatas)
        # Scintillation correction is applied later at the stage of the final split !!! 

    ## Prepare variables for final split ##########################################################################
    reducevlbi.prepare_variables_for_final_split(numtargets, inbeamnames, targetconfigs, expconfig, phscalnames, targetnames,
                directory, experiment)
    ## Split, image and write all three #############################################################################
    reducevlbi.split__normalize_on_individual_basis__then_write_out_uvdata_for_all_sources(numtargets, targetconfigs, expconfig, 
            targetonly, imageoutofbeam, inbeamuvdatas, bandpassclversion, ampcalsrc, directory, experiment, phscalnames, 
            modeldir, inbeamnames, calonly, targetnames, haveungated, ungateduvdata, scinttablepaths, gateduvdata)
    ## Image targets using Difmap and fit for position ###############################################################
    reducevlbi.image_targets_using_DIFMAP_and_fit_for_position(calonly, numtargets, targetconfigs, expconfig,
                directory, experiment, targetnames, beginif, endif, haveungated,
                phscalnames, inbeamnames, inbeamuvdatas, uvtaperstring, difmaptargetuvaverstring)
    ## Make some nice diagnostic plots #############################################
    if not skipdiagnosticplots:
        reducevlbi.make_diagnostic_plots(directory, codedir)

if __name__ == "__main__":
    main()
