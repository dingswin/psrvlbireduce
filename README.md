# User Guide for PSRVLBIREDUCE

Overview of this repository
This project provides a pipeline for VLBI data reduction, with special focus on VLBI astrometry. The main data reduction code, written in python-based parseltongue (Kettenis et al. 2006), is vlbi_astrometry.py which calls classes and functions from vlbireduce.py and support_vlbireduce.py. The latter two use functions provided in vlbatasks.py, where some new features have also been added. The pipeline was originally developed and used by Dr. Deller for the PSRPI project (Deller et al. 2019) and other projects. Since 2018, it has undergone large upgrades (including systematic restructuring) made by Hao Ding. Currently, the vlbireduce_v3 branch is already compatible with python3, and will be merged to the main branch after being tested properly.

New Feature of this branch:
1. This branch is now compatible with python3, and will be merged to master after being tested.
   
Should you have any inquiry about the pipeline or comments for improvements, feel free to contact Hao (hdingastro@hotmail.com) or Adam (adeller@astro.swin.edu.au).

## Required Packages and Software for PSRVLBIReduce: AIPS and PARSELTONGUE
Below details the required packages and software that needs to be installed prior to downloading and running PSRVLBIREDUCE.

The first required package is, of course, AIPs. AIPs can be downloaded by following the instructions located here by clicking on the date for which you want to download AIPs. Note that as of recently, AIPs is now compatible with M1 macs.

Once you have downloaded AIPs, you should ensure it works as expected. Download any VLBA or other VLBI instrument file, and then start AIPs from that location by running either AIPS or aips TV=local:0.0. Note that sometimes the latter is necessary to get the TV to work properly. Then, confirm that you can load in a file use fitld, and check that the header of the file looks okay. If you aren’t familiar with basic AIPs tasks, then this would be the time to become familiar. There are three really nice AIPs tutorials: A small guide to AIPS, simple VLBA project including self-calibration and spectral-line VLBA project plus astrometry that are available to introduce you to the basics of AIPs. 

The second required package is ParselTongue. ParselTongue can be downloaded by following the instructions located here. Note that ParselTongue does work on the M1 macs, and the easiest way to download it is through using brew. 
Note that ParselTongue not only requires AIPS, but also requires OBIT. It is easiest to download and install Obit using the instructions on the ParselTongue website, as it can otherwise be a bit tricky to install.
Additionally, depending on your operating machine, ParselTongue requires a few additional dependencies such as gcc, glib, etc. Instructions on what is necessary can be found here.

Before moving on to the next step of actually setting up PSRVLBIREDUCE, you should confirm that ParselTongue is running nominally on your machine. To test whether or not ParselTongue is working nominally, start up ParselTongue from the command line. This should bring up the ParselTongue server. If it does not, something has either gone wrong with your installation or the insertion of the proper baths into your bashrc or bash_login file. 

Next, make sure you can load all of the AIPs software okay from within the ParselTongue window e.g., 

 from AIPS import AIPS
 from AIPSTask import AIPSTask, AIPSList
 from AIPSData import AIPSUVData, AIPSImage

Next, try loading in your data and confirming it exists a-okay e.g.. 

fitld = AIPSTask('FITLD') 
fitld.datain = "<path to data>/gated_test.idifits" # note even if in same directory, need to pass the full path length 
fitld.outname = "test_file" 
fitld.outclass=”UVDATA”
fitld.inp() # checks the input for fitld to make sure it looks okay
fitld.go()

This should load your data. Next, let’s confirm it exists a-okay.

uvdata = AIPSUVData("test_gated", "UVDATA", 1, 1)
uvdata.exists() 

This should return the value True, and your level 0 testing is now complete. You are ready for PSRVLBIREDUCE installation 

## Installation
Below, we detail the instructions for installing the PSRVLBIREDUCE pipeline.

Download the PSRPIVLBIREDUCE github repo from https://github.com/dingswin/psrvlbireduce. You may need to request permission to view this github repository and/or edit it. 
Once you have downloaded the repository, there are a few paths that need to be set to get it running smoothly. The easiest way to do this is to create a file called source_file_psrvlbireduce.sh in the same directory as this README.md. Your source file needs to contain the following:

. <path to AIPs>/AIPs/LOGIN.sh
export PATH=$PATH:<path to psrvlbireduce>/psrvlbireduce
export PSRVLBAAUXDIR=<path to psrvlbireduce>/psrvlbireduce/examples
export PSRVLBICODEDIR=<path to psrvlbireduce>
export PSRVLBICODEDIR=<path to psrvlbireduce>/psrvlbireduce/datareduction
export PSRVLBAIPSVER=<version of AIPS you are running>

An example .sh file is listed below:

. /Users/Alice/software/AIPs/LOGIN.SH
export PATH=$PATH:/Users/Alice/Documents/PSR/psrvlbireduce
export PSRVLBAUXDIR=/Users/Alice/Documents/PSR/psrvlbireduce/examples
export PSRVLBICODEDIR=/Users/Alice/Documents/PSR/
export PSRVLBICODEDIR=/Users/Alice/Documents/PSR/psrvlbireduce/datareduction
export PSRVLBAIPSVER=31DEC23
export PATH=$PATH:/Users/Alice/Documents/PSR/psrvlbireduce/datareduction


You will also need to add the location of the psrvlbireduce code to your bash_login e.g., export PYTHONPATH=$PYTHONPATH:<path to psrvlbireduce>/psrvlbireduce/

Ideally, this should be everything you need to now start playing around with the code. 
Example Module #1
First, before getting started, make sure to run your source_file_psrvlbireduce.sh file. Otherwise, everything will break.

Downloading the data
Next, the example module uses data from the experiment bd179i0 for PSR J1738+0333. This data will need to be downloaded from the NRAO data archive site. One at the data archive site, search for bd179i0 and then download the files:

VLBA_BD179I0_ungatedi0_BIN0_SRC0_0_150824T164359.idifits
VLBA_BD179I0_gatedi0_BIN0_SRC0_0_150824T163822.idifits
VLBA_BD179I0_inbeam1i0_BIN0_SRC0_0_150824T164856
VLBA_BD179I0_inbeam2i0_BIN0_SRC0_0_150824T165207.idifits
VLBA_BD179I0_inbeam3i0_BIN0_SRC0_0_150824T165321.idifits

Note that in total, these files are ~12.7 Gb, so you will need adequate storage on your machine. You should place these files under /Users/Bob/PSR/examples/J1738+0333/bd179i0. 

Downloading auxiliary data files
Next, you need to acquire the EOP (earth-orientation parameter) file and the ionospheric files for that given day. This can easily be done by running the program prepare_astrometric_epoch.py under /datareduction. To run this file successfully, you will either need to supply it with the .vex file from your observations, or run it under the directory in which your data lives so it can grab the name of the .idifits file and figure out the date for which the files need to be grabbed. For example, 

In /Users/Alice/PSR/examples/J1738+0333/bd179i0 you can either run prepare_astrometric_epoch.py bd179i0.vex or you can just run prepare_astrometric_epoch.py. This should download all of the necessary files for you to be on your way.

If you would prefer to download the above files manually, you can do so. The files you will need are:
Usno_finals.erp
codg<>.15i.Z
esag<>.15.Z
igsg<>.15i
jplg<>.15i.Z
The usno_finals.erp file can be fftp’ed from ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp. The ionospheric files can be found here.

You will need to place these files under a new directory /Users/Alice/PSR/examples/J1738+0333/bd179i0/logs. Additionally, you will need to make the directories /Users/Alice/PSR/examples/J1738+0333/bd179i0/images and /Users/Alice/PSR/examples/J1738+0333/bd179i0/tables. 

## Running the program
You should not be all set to get to work running the program. Go to /Users/Alice/PSR/psrvlbireduce/datareduction/ and run ParselTongue vlbi_astrometry.py -e bd179i0

You will be prompted with many questions along the way, as the program runs. The initial set of questions will be related to the loading in of your data. You want to always answer yes to these questions, otherwise the program will abort (seeing as it has no data to use). This is the computationally longest part of the pipeline.

Many of the questions after this will be related to using saved .SN and .CL tables. If it is your first time running the program, you won’t have any saved products, so you won’t have any option but to create these tables. However, if you are running the program for the nth time, and only want to adjust something much further downstream, you should feel free to use your saved .SN and .CL tables.

After you get to the phase calibrator step, the program will write out a .fits image (in this case it will be called J1740_formodeling.fits of your phase calibrator, and then exit. This is because it requires a model for the phase calibrator. Hence, you will need to load the phase calibrator file (in this case the phase calibrator is J1740_0311) into a normal AIPs environment, clean it, and then save the image as J1740+0311.clean.fits under /Users/Alice/PSR/psrvlbireduce/examples/sourcemodels/preliminary/.

At some point later on, you will be asked to do the same for your primary in-beam calibrator. After creating the clean images, proceed to re-run the program from the start (here would be a good time to use those saved .SN and .CL tables!). 

The program should (ideally) run all the way to the end. Given all your paths are properly set-up, it should end by running make_diagnostics.py which will produce a file Diagnostic.html which will contain all of the relevant plots from your work. Additionally, all of the logs outputted on the command line will be stored under bd179i0.datacheck.log and the pipeline summary will be saved under bd179i0.pipelinesummary. 

Interpreting the output data

If everything has gone right, you should have a file Diagnostic.html with many different plots. Below we detail how to interpret those plots. 

The document should start with three different image plots. They should show the pulsar and two of the in-beam calibrators along with their positions. 

Taking it one step further and generating a file for PMPAR

Coming soon…. 

## Detailed Overview of the PSRVLBIREDUCE Pipeline

Coming soon…



Should you have any inquiry about the pipeline or comments for improvements, feel free to contact Hao (hdingastro@hotmail.com) or Adam (adeller@astro.swin.edu.au).
