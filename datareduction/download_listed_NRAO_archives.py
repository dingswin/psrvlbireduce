#!/usr/bin/env python
import os, sys
usage = "download_listed_NRAO_archives.py downloads.txt\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()
[junk, downloadlistfile] = sys.argv
lines = open(downloadlistfile).readlines()
for line in lines:
    prompt = 'Public File available : '
    if prompt in line:
        link = line.split(prompt)[-1].strip()
        os.system('wget -t 0 %s' % link)
print "Downloads complete"
