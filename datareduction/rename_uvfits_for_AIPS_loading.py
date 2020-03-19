#!/usr/bin/env python
import os, glob
uvfitsfiles = glob.glob(r'*.uvfits')
uvfitsfiles.sort()
i = 1
for uvfitsfile in uvfitsfiles:
    os.system('cp %s file%d' % (uvfitsfile, i))
    i += 1
