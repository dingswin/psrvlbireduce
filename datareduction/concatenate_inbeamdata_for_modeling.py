#!/usr/bin/env python
"""
Usage
-----
concatenate_*py inbeamname

Input
-----
inbeamname, input as argv[1].

Note
----
Run under the targetdir.
"""
import glob,os,sys
[junk, inbeamname] = sys.argv
vexfiles = glob.glob(r'*/*.vex')
if len(vexfiles) == 0:
    print('No vex file is found; aborting')
    sys.exit()
if not os.path.exists(inbeamname):
    os.mkdir(inbeamname)
for vexfile in vexfiles:
    expno = vexfile.split('/')[0].strip()
    inbeamuvdatas = glob.glob(r'%s/%s_%s_pipeline_uv.fits' % (expno, expno, inbeamname))
    if len(inbeamuvdatas) != 1:
        print('No inbeamuvdata found or multiple found; aborting')
        sys.exit()
    inbeamuvdata = inbeamuvdatas[0]
    os.system('cp %s %s' % (inbeamuvdata, inbeamname))

os.chdir(inbeamname)
print('Combining single-epoch uvfitsfiles')
os.system('dbcon.py')
os.system('mv dbcon_uv.fits %s_formodeling_uv.fits' % inbeamname)
#print('Deleting single-epoch uvfitsfiles...')
#fitsfile2delete = glob.glob(r'*_pipeline_uv.fits')
#for fitsfile in fitsfile2delete:
#    os.remove(fitsfile)
