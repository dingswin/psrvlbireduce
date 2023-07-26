# User Guide for PSRVLBIREDUCE

Overview of this repository
This project provides a pipeline for VLBI data reduction, with special focus on VLBI astrometry. The main data reduction code, written in python-based parseltongue (Kettenis et al. 2006), is vlbi_astrometry.py which calls classes and functions from vlbireduce.py and support_vlbireduce.py. The latter two use functions provided in vlbatasks.py, where some new features have also been added. The pipeline was originally developed and used by Dr. Deller for the PSRPI project (Deller et al. 2019) and other projects. Since 2018, it has undergone large upgrades (including systematic restructuring) made by Hao Ding. Currently, the vlbireduce_v3 branch is already compatible with python3, and will be merged to the main branch after being tested properly.

New Feature of this branch: This branch is now compatible with python3, and will be merged to master after being tested.
   
Should you have any inquiry about the pipeline or comments for improvements, feel free to contact Hao (hdingastro@hotmail.com) or Adam (adeller@astro.swin.edu.au).

This README.md was generated by Alice Curtin (alice.curtin@mail.mcgill.ca)

## Required Packages and Software for PSRVLBIReduce: AIPS and PARSELTONGUE
Below details the required packages and software that needs to be installed prior to downloading and running PSRVLBIREDUCE.

The first required package is, of course, AIPs. AIPs can be downloaded by following the instructions located [here](http://www.aips.nrao.edu/index.shtml) by clicking on the date for which you want to download AIPs. Note that as of recently, AIPs is now compatible with M1 macs.

Once you have downloaded AIPs, you should ensure it works as expected. Download any VLBA or other VLBI instrument file, and then start AIPs from that location by running either `AIPS` or `aips TV=local:0.0`. Note that sometimes the latter is necessary to get the TV to work properly. Then, confirm that you can load in a file use fitld, and check that the header of the file looks okay. If you aren’t familiar with basic AIPs tasks, then this would be the time to become familiar. There are three really nice AIPs tutorials: [A small guide to AIPS](http://www3.mpifr-bonn.mpg.de/staff/hrk/AIPS_TUTORIAL/HRK_AIPS_1.html), [simple VLBA project including self-calibration](https://casaguides.nrao.edu/index.php/AIPS-Simple-Self-Cal) and [spectral-line VLBA project plus astrometry](https://casaguides.nrao.edu/index.php/AIPS-Spectral_Lines_and_Astrometry) that are available to introduce you to the basics of AIPs. 

The second required package is ParselTongue. ParselTongue can be downloaded by following the instructions located [here](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue). Note that ParselTongue does work on the M1 macs, and the easiest way to download it is through using brew. 
Note that ParselTongue not only requires AIPS, but also requires OBIT. It is easiest to download and install Obit using the instructions on the ParselTongue website, as it can otherwise be a bit tricky to install.
Additionally, depending on your operating machine, ParselTongue requires a few additional dependencies such as gcc, glib, etc. Instructions on what is necessary can be found [here](https://www.jive.eu/parseltongue/codex.pdf).

Before moving on to the next step of actually setting up `PSRVLBIREDUCE`, you should confirm that ParselTongue is running nominally on your machine. To test whether or not ParselTongue is working nominally, start up ParselTongue from the command line. This should bring up the ParselTongue server. If it does not, something has either gone wrong with your installation or the insertion of the proper baths into your bashrc or bash_login file. 

Next, make sure you can load all of the AIPs software okay from within the ParselTongue window e.g., 

 ```
 from AIPS import AIPS
 from AIPSTask import AIPSTask, AIPSList
 from AIPSData import AIPSUVData, AIPSImage
 ```

Next, try loading in your data and confirming it exists a-okay e.g.. 

```
fitld = AIPSTask('FITLD') 
fitld.datain = "<path to data>/gated_test.idifits" # note even if in same directory, need to pass the full path length 
fitld.outname = "test_file" 
fitld.outclass=”UVDATA”
fitld.inp() # checks the input for fitld to make sure it looks okay
fitld.go()
```

This should load your data. Next, let’s confirm it exists a-okay.

```
uvdata = AIPSUVData("test_gated", "UVDATA", 1, 1)
uvdata.exists()
```

This should return the value True, and your level 0 testing is now complete. You are ready for PSRVLBIREDUCE installation 

## Installation
Below, we detail the instructions for installing the PSRVLBIREDUCE pipeline.

Download the `PSRPIVLBIREDUCE` github repo from https://github.com/dingswin/psrvlbireduce. You may need to request permission to view this github repository and/or edit it. 
Once you have downloaded the repository, there are a few paths that need to be set to get it running smoothly. The easiest way to do this is to create a file called `source_file_psrvlbireduce.sh` in the same directory as this `README.md`. Your source file needs to contain the following:

```
. <path to AIPs>/AIPs/LOGIN.sh
export PATH=$PATH:<path to psrvlbireduce>/psrvlbireduce
export PSRVLBAAUXDIR=<path to psrvlbireduce>/psrvlbireduce/examples
export PSRVLBICODEDIR=<path to psrvlbireduce>
export PSRVLBICODEDIR=<path to psrvlbireduce>/psrvlbireduce/datareduction
export PSRVLBAIPSVER=<version of AIPS you are running>
```

An example .sh file is listed below:

```
. /Users/Alice/software/AIPs/LOGIN.SH
export PATH=$PATH:/Users/Alice/Documents/PSR/psrvlbireduce
export PSRVLBAUXDIR=/Users/Alice/Documents/PSR/psrvlbireduce/examples
export PSRVLBICODEDIR=/Users/Alice/Documents/PSR/
export PSRVLBICODEDIR=/Users/Alice/Documents/PSR/psrvlbireduce/datareduction
export PSRVLBAIPSVER=31DEC23
export PATH=$PATH:/Users/Alice/Documents/PSR/psrvlbireduce/datareduction
```

You will also need to add the location of the psrvlbireduce code to your bash_login e.g., `export PYTHONPATH=$PYTHONPATH:<path to psrvlbireduce>/psrvlbireduce/`

Ideally, this should be everything you need to now start playing around with the code. 

## Example
For an example of how to run this pipeline on some previously acquired and publicly accessible data, please see the `README.md` in the `/examples` folder.

## Detailed Overview of the PSRVLBIREDUCE Pipeline
Below, we detail the main functions performed by the PSRVLBIREDUCE pipeline. Note that this was compiled on July 26, 2023, and hence might not apply for a much later and updated version of `psrvlbireduce`. 

### _**LOAD DATA AND PERFORM BASIC CORRECTIONS**_

1. Load the data
2. [OPTIONAL]: increase the number of IFs
    * Relevant parameters in .yaml file:
        * moriffactor
3. [OPTIONAL]: Unflag pitown
    * Relevant parameters in .yaml file:
        * unflagpietown
4. [OPTIONAL]: Load user flags for RFI
    * No specific parameters required in .yaml file. Instead, it will scrape the /tables director and look for any relevant files.
5. Run the flagging for zero-fringe rate effects
    * Flagging baselines when the natural fringe rate is low
        * These times are susceptible to corruption by radio frequency interference (RFI) and instrumental effects
    * Relevant parameters in yaml file:
        * fringerateflagsuppressionfactor
        * If this parameter doesn’t exist, suppression factor set to be 7
6. [OPTIONAL]: Run the auto flagger
    * Flags data per source, per baseline, and per polarization 
    * Finds the median of the data, and then flags data that is some sigma away e.g., 7*MAD, away from this median
    * Relevant parameters in yaml file:
        * autoflag
7. [OPTIONAL]: Correct positions
    * Shifts UV-data by amounts specified in the yaml file if the positions were initially off
    * Can be useful if your initial pulsar position was very incorrect, or if the in-beam position wasn’t well known (and can be better determined later)
    * Relevant parameters in yaml file:
        * shift 
        * e.g., ["J201810+2833, 60, 820"] where first you list the name, then the ra offset in thousands of mas and then the dec offset in thousands of mas 
8. Run TECOR to correct for ionosphere
    * Corrects for the ionosphere in the gated, ungated, and in-beam calibrator datasets using previously downloaded ionospheric model files
    * does not count for short time interval ionospheric changes, but larger, global changes
9. Run CLCOR to correct EOPs
    * Corrects for the EOPs using the most up-to-date usno_finals.erp file
10. Prepare leakage-related variables
    * Relevant .yaml parameters:
        * xpolscan, xpolsource, leakagescan, leakagesource, leakageuvrange, leakageweightit
        * If nothing in the yaml file, the leakage source is the amplitude calibrator, the leakage weight is 0, the leakage uvrange is [0,0]
11. Run CLCOR to correct PANG
    * Make parallactic angle corrections using CLCOR
12. Correct for a-priori delays if available using CLCOR
    * Loads and applies known apriori clock delays which must be supplied
    * Relevant files:
        * aprioridelays.txt located in /tables
13. Correct for primary beam attenuation
    * As neither the in-beam calibrator nor the pulsar will lie at the centre of the field of view. Hence, there is attenuation at both of their positions that is corrected by using a Bessel approximate to estimate the response of the VLBA beam
    * Relevant parameter in .yaml files:
        * skippbcor — if True, does not perform primary beam corrections
14. Do PCAL correction and inspect
    * Runs additional calibration using a pulse calibrator source
    * NOTE: The pulsar calibrator file must already exist, and hence must have been requested when submitting your observations 

### _**FRINGE FINDER CORRECTIONS AND LEAKAGE CORRECTIONS**_

15. Run FRINGE on fringe calibrator and then apply it
    * Runs fringe on the fringe finder source (called ampcalsrc). Can use a saved .sn table if it exists otherwise creates and then saves a .sn file to them be loaded and applied. 
16. [OPTIONAL]: Run xpoldelaycal, and then apply it
    * Performs cross-polar delay calibration
    * Requires the  variable Xpolscan needed to have been set in the yaml file, otherwise the earlier leakage stage set xpolscan to -1 which causes this part to not fun
    * Relevant parameter in .yaml files:
        * inttimesecs, xpolsolintmins, delaywindows, ratewindowmhz
    * Requires a model for xpolarization calibrator
17. Run leakage cal on leakage calibrator
    * Requires leakgescan variable in .yaml file (comes from previous preparing leakage-related variables step) 
    * Performs leakage corrections using leakage source through computing instrumental polarization corrections
18. Run BPASS, load it, and plot it
    * Corrects for bandpass using the fringe finder source (ampcalsrc)
19. [OPTIONAL]: Plot the bandpass solutions
    * Required variable in .yaml config:
        * plotbandpass 
        
### _**PHASE CALIBRATOR OPERATIONS**_

20. Run FRINGE on phase reference calibrator
    * Runs fringe on the phase reference calibrator
    * If no model for phase calibrator, program will break here after creating a .fits file for the phase calibrator which will need to then be imaged for the program to keep running
    * Optional variables in .yaml file:
        * ampcalscan, phsrefuvrange, delay windows, ratewindowmhz, phsreffringsumifs, phsreffringratemwfhz, phsreffringedelaymwfns, phsreffringsumpols, resreffringhalfbandwidth, phsreffringdispersivefit, inttimesecs, exhaustivefring
    * After running bring, we copy the SN table and apply it to the data
21. [OPTIONAL]: Run phase CALIB on phase reference calibrators
    * Runs phase self-calibration on the phase reference source
    * Relevant .yaml variables:
        * maxphsrefcalibpnmins
        * If this variable < 1, then self-calibration is not performed
22. [OPTIONAL]: Run amplitude CALIB on phase reference calibrators
    * Runs amplitude self-calibration on the phase reference source
    * Note you must already have a model for your in-beam calibrator 
    * Relevant .yaml variables:
        * maxphsrefcalibapnmins — if this variable < 1, then self-calibration is not performed
IN-BEAM CALIBRATOR OPERATIONS
23. [OPTIONAL]: generate a raw in-beam dataset pre-any self calibration
    * Prior to running any self-calibration on the in-beam source, writes to disk the inbeam
    * Relevant .yaml variables:
        * writerawinbeam 
24. [OPTIONAL]: Do a combined IF and polarization phase self-calibration on the in-beam calibrator
    * Perform one round of phase self-calibration on the in-beam calibrator with the IFs and polarizations combined
    * Relevant .yaml variables:
        * maxinbeamcalibp1mins — this must be above must be > 0 for this to run
        * inbeamuvrange, inbeamcalibp1mins, inbeamcalibp1snr, inbeamcalibp1type, inbeamminmodelflux
25. [OPTIONAL]: Do a separate IF phase self-calibration on the in-beam calibrators
    * Same as above except do self-calibration per IF rather than combining all together
26. [OPTIONAL]: Do dual-phasecal calibration with IFs and pol combined
    * Applies dual-phase cal mode on the primary in-beam calibrator 
    * Dual-phase calibrator mode refers to use two co-linear phase calibrators
    * Seems to be fairly rarely used 
27. [OPTIONAL]: Do dual-phasecal calibration with separate IFs
    * Same as above but now per IF rather than combing all IFs
28. [OPTIONAL]: Combined IF amplitude and phase self-calibration for in-beam calibrator
    * Amplitude and phase self-calibration for the in-beam calibrator with IFs combined
    * Relevant .yaml parameters:
        * maxinbeamcalibap1mins — Must be >0 for this to run
        * inbeamcalibap1mins, inbeamcalibap1snr, inbeamcalibap1type, inbeamcalibap1weightit
29. [OPTIONAL]: Separate IF amplitude and phase self-calibration for in-beam calibrator 
    * Same as above except now for separate IFs
    * All variables are the same except ap1 is now apn
    * Relevant .yaml parameters:
        * maxinbeamcalibapnmins — Must be >0 for this to run
        * inbeamcalibapnmins, inbeamcalibapnsnr, inbeamcalibapnweightit, inbeamcalibapnmins, inbeamcalibapnstokesi
30. [OPTIONAL]: Phase self-calibration on second in-beam calibrator
    * Rarely used
    * Relevant variables are the same as that for the first round except p1 is replaced with sp1 e.g., maxinbeamcalibsp1mins
31. [OPTIONAL]: Do dual-phscal calibration on the prIBC-secIBC line if requested
    * Applies dual-phase cal mode on the secondary in-beam calibrators
    * Dual-phase calibrator mode refers to use two co-linear phase calibrators
    * Seems to be fairly rarely used 
32. [OPTIONAL]: Calculate scintillation corrections
    * Correcting for scintillation that can cause a scattering disk in the final image
    * Relevant .yaml file parameters:
        * noscintcorr -> needs to be set to True
        * scintcorrmins — must be greater than zero for this to run
        * 
    * Dual-phase calibrator mode refers to use two co-linear phase calibrators
    * Seems to be fairly rarely used 

### _**IMAGE AND WRITE OUT DATA**_

33. Split, image, and write out all sources
34. Image the target sources using DIFMAP
    * Images the targets (pulsars, in-beam cals, and phase calibrators). Imaging technique depends on which source is being imaged e.g., phase calibrator uses a wide-field image while the other two don’t
    * Fits for positions using JMFIT
35. [OPTIONAL]: Produce some nice diagnostics
    * Puts together a Diagnostic.html file with all of the final images and a bunch of diagnostic plots
    * VERY helpful!

