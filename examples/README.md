## Example Module #1
Below, we will detail an example of how to run this pipeline on data that was acquired for Deller et al. (2019). Remember, before getting started, make sure to run your `source_file_psrvlbireduce.sh file`. Otherwise, everything will break.

##### Downloading the data
Next, the example module uses data from the experiment bd179i0 for PSR J1738+0333. This data will need to be downloaded from the [NRAO data archive site](https://data.nrao.edu/portal/). One at the data archive site, search for bd179i0 and then download the files:

```
VLBA_BD179I0_ungatedi0_BIN0_SRC0_0_150824T164359.idifits
VLBA_BD179I0_gatedi0_BIN0_SRC0_0_150824T163822.idifits
VLBA_BD179I0_inbeam1i0_BIN0_SRC0_0_150824T164856
VLBA_BD179I0_inbeam2i0_BIN0_SRC0_0_150824T165207.idifits
VLBA_BD179I0_inbeam3i0_BIN0_SRC0_0_150824T165321.idifits
```
Note that in total, these files are ~12.7 Gb, so you will need adequate storage on your machine. You should place these files under `/Users/Bob/PSR/examples/J1738+0333/bd179i0`. 

#### Downloading auxiliary data files
Next, you need to acquire the EOP (earth-orientation parameter) file and the ionospheric files for that given day. This can easily be done by running the program prepare_astrometric_epoch.py under /datareduction. To run this file successfully, you will either need to supply it with the .vex file from your observations, or run it under the directory in which your data lives so it can grab the name of the .idifits file and figure out the date for which the files need to be grabbed. For example, 

In `/Users/Alice/PSR/examples/J1738+0333/bd179i0` you can either run `prepare_astrometric_epoch.py bd179i0.vex` or you can just run `prepare_astrometric_epoch.py`. This should download all of the necessary files for you to be on your way.

If you would prefer to download the above files manually, you can do so. The files you will need are:
```
Usno_finals.erp
codg<>.15i.Z
esag<>.15.Z
igsg<>.15i
jplg<>.15i.Z
```
The usno_finals.erp file can be fftp’ed from ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp. The ionospheric files can be found [here](https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html).

You will need to place these files under a new directory `/Users/Alice/PSR/examples/J1738+0333/bd179i0/logs`. Additionally, you will need to make the directories `/Users/Alice/PSR/examples/J1738+0333/bd179i0/images` and `/Users/Alice/PSR/examples/J1738+0333/bd179i0/tables`. 

#### Running the program
You should not be all set to get to work running the program. Go to `/Users/Alice/PSR/psrvlbireduce/datareduction/` and run `ParselTongue vlbi_astrometry.py -e bd179i0`

You will be prompted with many questions along the way, as the program runs. The initial set of questions will be related to the loading in of your data. You want to always answer yes to these questions, otherwise the program will abort (seeing as it has no data to use). This is the computationally longest part of the pipeline.

Many of the questions after this will be related to using saved .SN and .CL tables. If it is your first time running the program, you won’t have any saved products, so you won’t have any option but to create these tables. However, if you are running the program for the nth time, and only want to adjust something much further downstream, you should feel free to use your saved .SN and .CL tables.

After you get to the phase calibrator step, the program will write out a .fits image (in this case it will be called `J1740_formodeling.fits` of your phase calibrator, and then exit. This is because it requires a model for the phase calibrator. Hence, you will need to load the phase calibrator file (in this case the phase calibrator is J1740_0311) into a normal AIPs environment, clean it, and then save the image as `J1740+0311.clean.fits` under `/Users/Alice/PSR/psrvlbireduce/examples/sourcemodels/preliminary/`. Note that for this work to work, you need to specify to use a preliminary model in the config file using `useprelimmodels`. This sets the model type to "preliminary" rather than the default of "final."

At some point later on, you will be asked to do the same for your primary in-beam calibrator. After creating the clean images, proceed to re-run the program from the start (here would be a good time to use those saved .SN and .CL tables!). 

The program should (ideally) run all the way to the end. Given all your paths are properly set-up, it should end by running `make_diagnostics.py` which will produce a file `Diagnostic.html` which will contain all of the relevant plots from your work. Additionally, all of the logs outputted on the command line will be stored under `bd179i0.datacheck.log` and the pipeline summary will be saved under `bd179i0.pipelinesummary`. 

#### Interpreting the output data

If everything has gone right, you should have a file `Diagnostic.html` with many different plots. Below we detail how to interpret those plots. 

The document should start with three different image plots. They should show the pulsar and two of the in-beam calibrators along with their positions. 

#### Taking it one step further and generating a file for PMPAR

The necessary input for pmpar is a file whose format is:

`<date> <ra> <ra error in hms> <dec> <dec error in mas>`

Thus, you can either use the .stats file from your gated or ungated pulsar dataset to retrieve all of these quantities, or you can use the program `jmfit2pmpar.py`. At this time, however, detailed notes on how to use `jmfit2pmpar.py` do not exist.
