#!/usr/bin/env ParselTongue

import os, sys, math
import vlbatasks
import numpy
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData

usage  = "uvsubScaled.py <uvinputfile1> <uvinputfile2> <scalefactor> <uvoutputfile>\n"
usage += "  (subtracts uvinputfile2 from uvinputfile1, after scaling uvinputfile2 by scalefactor)"
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    aipsver = '31DEC23'
AIPS.userno = 101
uvdata1     = AIPSUVData("JUNK1", "JUNK1", 1, 1)
uvdata2     = AIPSUVData("JUNK2", "JUNK2", 1, 1)
uvfile1     = sys.argv[1]
uvfile2     = sys.argv[2]
outuvfile   = sys.argv[4]
scalefactor = float(sys.argv[3])
if uvdata1.exists():
    uvdata1.zap()
if uvdata2.exists():
    uvdata2.zap()

if uvfile1[0] != '/':
    uvfile1 = os.getcwd() + '/' + uvfile1
if uvfile2[0] != '/':
    uvfile2 = os.getcwd() + '/' + uvfile2
if outuvfile[0] != '/':
    outuvfile = os.getcwd() + '/' + outuvfile
tempoutfile = os.getcwd() + '/junkout.fits'
vlbatasks.fitld_uvfits(uvfile1, uvdata1, [])
vlbatasks.fitld_uvfits(uvfile2, uvdata2, [])

wizuvdata1 = WizAIPSUVData(uvdata1)
wizuvdata2 = WizAIPSUVData(uvdata2)

numif = 0
fqtable = uvdata1.table('FQ', 1)
for row in fqtable:
    try:
        for iffreq in row.if_freq:
            numif += 1
    except (AttributeError, TypeError):
        numif += 1
try:
    numchan = int(uvdata1.table('FQ', 1)[0].total_bandwidth[0]/\
              uvdata1.table('FQ', 1)[0].ch_width[0] + 0.5)
except TypeError:
    numchan = int(uvdata1.table('FQ', 1)[0].total_bandwidth/\
              uvdata1.table('FQ', 1)[0].ch_width + 0.5)
numstokes = len(uvdata1.stokes)
print("%d IFs, %d chans, %d stokes" % (numif, numchan, numstokes))

times = []
baselines = []
rowcount = 0
output2 = open("junk2.txt", "w")
for row in wizuvdata2:
    rowcount += 1
    times.append(row.time)
    baselines.append([row.baseline[0],row.baseline[1]])
    output2.write("%.6f %d-%d\n" % (row.time, row.baseline[0], row.baseline[1]))
output2.close()

print("Finished reading file 1...")

visibilities = numpy.zeros(rowcount*numif*numchan*numstokes*3).reshape(rowcount, numif, numchan, numstokes, 3)

print("Allocated space for visibilities")
rowcount = 0
for row in wizuvdata2:
    for i in range(numif):
        for j in range(numchan):
            for k in range(numstokes):
                for l in range(3):
                    visibilities[rowcount][i][j][k][l] = row.visibility[i][j][k][l]
    rowcount += 1

print("Now we've read in file 1's visibilities, time to skip through file 2")
numvisibilities = rowcount

times1 = []
output = open("junk.txt", "w")
for row in wizuvdata1:
    times1.append(row.time)
    output.write("%.6f %d-%d\n" % (row.time, row.baseline[0], row.baseline[1]))
output.close()

print("Loaded all of the data for subtraction...")
print("Master data has " + str(len(times1)) + " visibilities")
print("Data to subtract has " + str(numvisibilities) + " visibilities")

atindex = 0
numskipped = 0
numskipped2 = 0
numvis = 0
rowcount = 0
TINY = 0.000001
for row in wizuvdata1:
    if atindex == len(times):
        print("Ran out of times")
        break
    numvis += 1
    while atindex < len(times) and row.time > times[atindex]+TINY:
        print("Skipping a visibility due to time mismatch! %.10f %.10f" % (row.time, times[atindex]))
        atindex += 1
        numskipped += 1
#    while atindex < len(times) and (row.baseline[0] > baselines[atindex][0] or ((row.baseline[0] == baselines[atindex][0]) and (row.baseline[1] > baselines[atindex][1]))):
    #while atindex < len(times) and row.time < times[atindex]-TINY and (row.baseline[0] > baselines[atindex][0] or ((row.baseline[0] == baselines[atindex][0]) and (row.baseline[1] > baselines[atindex][1]))):
    ## WARNING - CHANGING THIS FROM ABOVE, SINCE THE TIME CHECK SEEMED INCORRECT
    while ( atindex < len(times) # Not running off the end of the times array
             and 
              row.time < times[atindex]+TINY # Current row is earlier than or equal to the time we are currently looking at - ensures we don't race off into the future
             and  
              (row.baseline[0] > baselines[atindex][0] # The row's baseline is a higher number than the one we are currently looking at
               or 
                ((row.baseline[0] == baselines[atindex][0]) and (row.baseline[1] > baselines[atindex][1])) # The row's baseline is a higher number than the one we're currently looking at
               or 
                (row.baseline[0] != baselines[atindex][0] and row.baseline[1] != baselines[atindex][1] and baselines[atindex][0] == baselines[atindex][1]) # The baseline doesn't match and we're currently looking at an autocorrelation, which might be present in the to-subtract dataset and not the primary dataset
              )
           ):
        print("Skipping a visibility due to baseline mismatch!")
        print(row.baseline, baselines[atindex], atindex, len(baselines), len(times))
        atindex += 1
        numskipped += 1
    if atindex == len(times):
        print("Ran out of times")
        break
    rowcount += 1
    if row.time < times[atindex]-TINY or row.baseline[0] != baselines[atindex][0] or row.baseline[1] != baselines[atindex][1]:
        print("Skipping a visibility because times no longer match after baseline mismatch")
        print(row.baseline, baselines[atindex], row.time, times[atindex])
        for i in range(numif):
            for j in range(numchan):
                for k in range(numstokes):
                    row.visibility[i][j][k][2] = 0.0
        row.update()
        continue
    #print("\n\n\n" + str(baselines[atindex]) + " " + str(row.baseline) + ", Before: " + str(row.visibility))
    if row.baseline[0] != baselines[atindex][0] or row.baseline[1] != baselines[atindex][1] or row.time != times[atindex]:
        print("ARGGH!")
    #if row.visibility[0][0][0][2] > 0:
    #    weightratio = row.visibility[0][0][0][2]/visibilities[atindex][0][0][0][2]
    #    if weightratio > 1.1 or weightratio < 0.9:
    #        for i in range(numif):
    #            for j in range(numchan):
    #                for k in range(numstokes):
    #                    row.visibility[i][j][k][2] = 0.0
    #        numskipped2 += 1
    #        atindex += 1
    #        continue
    for i in range(numif):
        for j in range(numchan):
            for k in range(numstokes):
                row.visibility[i][j][k][0] -= scalefactor*visibilities[atindex][i][j][k][0]
                row.visibility[i][j][k][1] -= scalefactor*visibilities[atindex][i][j][k][1]
    row.update()
    #if math.sqrt(row.visibility[0][0][0][0]**2 + row.visibility[0][0][0][1]**2) > 0.035 and math.sqrt(row.uvw[0]**2+row.uvw[1]**2+row.uvw[2]**2) > 50000000:
    #    print(row.visibility)
    #print("After: " + str(row.visibility))
    atindex += 1

wizuvdata1.update()

print("Skipped " + str(numskipped) + " out of " + str(numvis) + " visibilities, plus another " + str(numskipped2) + " because of weight mismatch")
vlbatasks.writedata(uvdata1, tempoutfile, True)
os.system("mv %s %s" % (tempoutfile,outuvfile))
