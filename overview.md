## Detailed Overview of the PSRVLBIREDUCE Pipeline
Below, we detail the main functions performed by the PSRVLBIREDUCE pipeline. We note what functions automatically will run, and which are optional and hence must be indicated in the .yaml configuration files.

_Note that this was compiled on July 26, 2023, and hence might not apply for a much later and updated version of PSRVLBIREDUCE._

### _**LOAD DATA AND PERFORM BASIC CORRECTIONS**_

1. Load the data
2. [OPTIONAL]: increase the number of IFs
    * Relevant parameters in .yaml file:
        * `moriffactor`
3. [OPTIONAL]: Unflag pitown
    * Relevant parameters in .yaml file:
        * `unflagpietown`
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
        * `autoflag`
7. [OPTIONAL]: Correct positions
    * Shifts UV-data by amounts specified in the yaml file if the positions were initially off
    * Can be useful if your initial pulsar position was very incorrect, or if the in-beam position wasn’t well known (and can be better determined later)
    * Relevant parameters in yaml file:
        * shift 
        * e.g., `["J201810+2833, 60, 820"]` where first you list the name, then the ra offset in thousands of mas and then the dec offset in thousands of mas 
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
        * `aprioridelays.txt` located in `/tables`
13. Correct for primary beam attenuation
    * As neither the in-beam calibrator nor the pulsar will lie at the centre of the field of view. Hence, there is attenuation at both of their positions that is corrected by using a Bessel approximate to estimate the response of the VLBA beam
    * Relevant parameter in .yaml files:
        * `skippbcor` — if True, does not perform primary beam corrections
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
        * `inttimesecs`, `xpolsolintmins`, `delaywindows`, `ratewindowmhz`
    * Requires a model for xpolarization calibrator
17. Run leakage cal on leakage calibrator
    * Requires leakgescan variable in .yaml file (comes from previous preparing leakage-related variables step) 
    * Performs leakage corrections using leakage source through computing instrumental polarization corrections
18. Run BPASS, load it, and plot it
    * Corrects for bandpass using the fringe finder source (ampcalsrc)
19. [OPTIONAL]: Plot the bandpass solutions
    * Required variable in .yaml config:
        * `plotbandpass `
        
### _**PHASE CALIBRATOR OPERATIONS**_

20. Run FRINGE on phase reference calibrator
    * Runs fringe on the phase reference calibrator
    * If no model for phase calibrator, program will break here after creating a .fits file for the phase calibrator which will need to then be imaged for the program to keep running
    * Optional variables in .yaml file:
        * `ampcalscan`, `phsrefuvrange`, `delay windows`, `ratewindowmhz`, `phsreffringsumifs`, `phsreffringratemwfhz`, `phsreffringedelaymwfns`, `phsreffringsumpols`, `resreffringhalfbandwidth`, `phsreffringdispersivefit`, `inttimesecs`, `exhaustivefring`
    * After running bring, we copy the SN table and apply it to the data
21. [OPTIONAL]: Run phase CALIB on phase reference calibrators
    * Runs phase self-calibration on the phase reference source
    * Relevant .yaml variables:
        * `maxphsrefcalibpnmins`
        * If this variable < 1, then self-calibration is not performed
22. [OPTIONAL]: Run amplitude CALIB on phase reference calibrators
    * Runs amplitude self-calibration on the phase reference source
    * Note you must already have a model for your in-beam calibrator 
    * Relevant .yaml variables:
        * `maxphsrefcalibapnmins` — if this variable < 1, then self-calibration is not performed
IN-BEAM CALIBRATOR OPERATIONS
23. [OPTIONAL]: generate a raw in-beam dataset pre-any self calibration
    * Prior to running any self-calibration on the in-beam source, writes to disk the inbeam
    * Relevant .yaml variables:
        * `writerawinbeam`
24. [OPTIONAL]: Do a combined IF and polarization phase self-calibration on the in-beam calibrator
    * Perform one round of phase self-calibration on the in-beam calibrator with the IFs and polarizations combined
    * Relevant .yaml variables:
        * `maxinbeamcalibp1mins` — this must be above must be > 0 for this to run
        * `inbeamuvrange`, `inbeamcalibp1mins`, `inbeamcalibp1snr`, `inbeamcalibp1type`, `inbeamminmodelflux`
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
        * `maxinbeamcalibap1mins` — Must be >0 for this to run
        * `inbeamcalibap1mins`, `inbeamcalibap1snr`, `inbeamcalibap1type`, `inbeamcalibap1weightit`
29. [OPTIONAL]: Separate IF amplitude and phase self-calibration for in-beam calibrator 
    * Same as above except now for separate IFs
    * All variables are the same except ap1 is now apn
    * Relevant .yaml parameters:
        * `maxinbeamcalibapnmins` — Must be >0 for this to run
        * `inbeamcalibapnmins`, `inbeamcalibapnsnr`, `inbeamcalibapnweightit`, `inbeamcalibapnmins`, `inbeamcalibapnstokesi`
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
        * `noscintcorr` -> needs to be set to True
        * `scintcorrmins` — must be greater than zero for this to run
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

