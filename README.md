About psrvlbireduce:

This project provides pipeline for VLBI data reduction, with special focus on VLBI astrometry.
The main data reduction code, written in python-based parseltongue (Kettenis et al. 2006), is vlbi_astrometry_reduce.py calling classes and functions from vlbireduce.py and support_vlbireduce.py. The latter two uses functions provided in vlbatasks.py, where some new features have also been added. 
The pipeline was originally developed and used by Dr. Deller for the PSRPI project (Deller et al. 2019) and other projects. Since 2018, it has gone through large upgrade (including systematic restructuring) made by Hao Ding. Currently, the master branch is already outdated. Please use the latest branch vlbireduce_v2, which is under continuous maintainance and development. vlbireduce_v2 will be merged to the master branch by mid-2022.

Should you have any inquiry about the pipeline or comments for improvements, feel free to contact Hao (haoding@swin.edu.au) or Adam (adeller@gmail.com).
