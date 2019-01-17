#!/usr/bin/env python
import os, sys

def parse_timerange(response, formatstring, dayno):
    splitresponse = response.split(',')
    if len(splitresponse) != 2:
        print "Malformed timerange string " + response + " - aborting"
        sys.exit()
    splittimes = [splitresponse[0].split(':'), splitresponse[1].split(':')]
    if len(splittimes[0]) != 3 or len(splittimes[1]) != 3:
        print "Malformed timerange string " + response + " - aborting"
        sys.exit()
    entries = []
    for s in splittimes:
        entries.append(dayno + int(s[0])/24)
        entries.append(int(s[0])%24)
        entries.append(int(s[1]))
        entries.append(int(s[2]))
    return formatstring % (entries[0], entries[1], entries[2], entries[3], 
                           entries[4], entries[5], entries[6], entries[7])

wd = os.getcwd()
splitwd = wd.split('/')
exp = splitwd[-1]
tabledir = wd + '/tables/'
vexfile = wd + '/' + exp + '.vex'

if not os.path.exists(tabledir):
    print tabledir + " does not exist - aborting!"
    sys.exit()

if not os.path.exists(vexfile):
    print vexfile + " does not exist - aborting!"
    sys.exit()

flagsource = raw_input('Flag source: ')
if flagsource=="":
    flagfile = tabledir + 'additionaledit.flag'
else:
    flagfile = tabledir + 'additionaledit.' + flagsource + '.flag'

dayno = -1
vexlines = open(vexfile).readlines()
for line in vexlines:
    if 'doy' in line:
        dayno = int(line.split()[-1])
        break
if dayno < 0:
    print "Couldn't fine day number in vex file! Aborting"
    sys.exit()

print "Going to add a flag for the observation " + exp
print "At any time, enter q to indicate that no more variables need to be set"
if not os.path.exists(flagfile):
    print "Creating additionaledit flag file " + flagfile
    output = open(flagfile, "w")
    output.write("opcode = 'FLAG'\n")
    output.write("dtimrang = 1  timeoff = 0\n")
else:
    print "Appending to " + flagfile
    output = open(flagfile, "a")

flagline = ""
ant = ""
bas = ""
time1 = ""
time2 = ""
if1 = ""
if2 = ""
reason = "manual/additional"

antennas = ["BR","FD","HN","KP","LA","MK","NL","OV","PT","SC"]
queries = []
entries = []
queries.append("Flag antenna")
entries.append("ant_name='%s' ")
queries.append("Flag baselines from %s to")
entries.append("bas_name='%s' ")
queries.append("Flag timerange (hh:mm:ss,hh:mm:ss)")
entries.append("timerang=%03d,%02d,%02d,%02d, %03d,%02d,%02d,%02d ")
#queries.append("Flag sources")
#entries.append("sources='%s' ")
queries.append("Flag IFs (begin,end)")
entries.append("bif=%d eif=%d ")
queries.append("Flag stokes")
entries.append("stokes='%s' ")

answer = ""
index = 0
flagants = ""
while index < len(queries):
    query = queries[index]
    if "Flag baselines from " in query:
        if flagants == "":
            index += 1
            continue
        query = query % flagants
    answer = raw_input(query + ': ')
    if answer == "q": break
    if query == "Flag antenna":
       flagants = answer
    splitanswer = answer.split(',')
    if answer != "":
        if "timerange" in query:
            flagline += parse_timerange(answer, entries[index], dayno)
            #flagline += parse_timerange(answer, entries[index], 0)
        elif "IFs" in query:
            if len(splitanswer) == 2:
                flagline += (entries[index] % (int(splitanswer[0]), int(splitanswer[1])))
            else:
                print "Malformed IF pair " + answer
                sys.exit()
        else:
            flagline += (entries[index] % answer)
    index += 1
#flagline += ("sources='%s' " % flagsource)

if flagline == "":
    print "You didn't enter anything to flag..?   I'm not going to do anything."
    sys.exit()

if not "timerang" in flagline:
    flagline += parse_timerange("0:0:0,47:00:00","timerang=%03d,%02d,%02d,%02d, %03d,%02d,%02d,%02d ",dayno) 
    #flagline += parse_timerange("0:0:0,47:00:00","timerang=%03d,%02d,%02d,%02d, %03d,%02d,%02d,%02d ",0)
if not "bas_name" in flagline:
    flagline += " bas_name=''"
if not "sources" in flagline:
    flagline += " sources=''"
if not "bif" in flagline:
    flagline += " bif=1 eif=0 "
if not "stokes" in flagline:
    flagline += " stokes=''"

if "ant_name" in flagline:
    output.write(flagline + " Reason='manual/additional' / \n")
else:
    for a in antennas:
        output.write(flagline + " ant_name='" + a + "' Reason='manual/additional' / \n")
output.close()
