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


