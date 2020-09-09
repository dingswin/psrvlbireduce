#!/usr/bin/env python
import os, sys
"""
copy the archive info to a downloads.txt at the webpage step to select 'Retrieve over internet',
then you have to copy and paste webpage link into the downloads.txt after a '#' symbol, and copy user and password in the same line separated with '#',
then run download*py downloads.txt
"""
usage = "download_listed_NRAO_archives.py downloads.txt\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()
[junk, downloadlistfile] = sys.argv
lines = open(downloadlistfile).readlines()
for line in lines:
    if line.startswith('#'):
        line = line.split('#')
        link = line[1].strip()
        user = line[2].strip()
        passwd = line[3].strip()
    #prompt = 'Public File available : '
    #if prompt in line:
    else:
        line = line.replace('  ', ' ').replace('  ', ' ')
        fitsfile = line.strip().split(' ')[0].strip()
        fitsfile = link + '/' + fitsfile
        os.system('wget -t 0 --user %s --password %s %s' % (user, passwd, fitsfile))
print "Downloads complete"
