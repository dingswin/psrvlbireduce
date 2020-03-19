#!/usr/bin/env python
import sys, os
target = sys.argv[1]
nepoch = int(sys.argv[2])
for i in range(nepoch):
    print "remove the epoch %d and see the difference" % i
    os.system('generatepmparin.py -t %s -bp -d %d' % (target, i))
print "\ndiagnose finished\n"

