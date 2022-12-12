About psrvlbireduce:

This project provides pipeline for VLBI data reduction, with special focus on VLBI astrometry.
The main data reduction code, written in python-based parseltongue (Kettenis et
al. 2006), is vlbi_astrometry.py calling classes and functions from
vlbireduce.py and support_vlbireduce.py. The latter two use functions provided
in vlbatasks.py, where some new features have also been added. The pipeline was
originally developed and used by Dr. Deller for the PSRPI project (Deller et
al. 2019) and other projects. Since 2018, it has undergone large upgrade
(including systematic restructuring) made by Hao Ding. 

New Feature of this branch:
1. This branch is now compatible with python3, and will be merged to master after being tested.

Should you have any inquiry about the pipeline or comments for improvements, feel free to contact Hao (hdingastro@hotmail.com) or Adam (adeller@astro.swin.edu.au).
