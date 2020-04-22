#!/usr/bin/env python
###########################################################################
###-- for astrometric runs of bin-gated data
###-- tailored mainly for bd223 data
###-- produce source file automatically and copy yaml file
###-- running prepare_astrometric_epoch.py unless -n
###-- NOTICE: RUN UNDER EXPDIR!!! (e.g. bd223d)
############################################################################
import os,sys,glob,yaml
from optparse import OptionParser

usage = "usage: run under EXPDIR\n%prog []\n-n --noprepare\n-b --binno\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-n", "--noprepare", dest="prepare", default="True",
                  action="store_false",help="NOT running prepare_astrometric_epoch.py beforehand")
parser.add_option("-b", "--binno", dest="binno", default=-1,
                  help="choose only one bin to reduce, input bin number")
(options, junk) = parser.parse_args()
prepare         = options.prepare
binno           = int(options.binno)

expdir = os.getcwd()
experiment = expdir.split('/')[-1]
gatedfitsfiles = glob.glob(r'./*_gated*idifits')
if len(gatedfitsfiles) < 2:
    gatedfitsfiles = glob.glob(r'./*_binbd*idifits')    
gatedfitsfiles.sort()
if binno != -1:
    gatedfitsfiles = glob.glob(r'./*_gated*BIN%d*idifits' % binno)
    if len(gatedfiltsfiles) == 0:
        gatedfitsfiles = glob.glob(r'./*_gated*BIN0%d*idifits' % binno)
ungatedfitsfile = glob.glob(r'./*_ungated*idifits')
if len(ungatedfitsfile) > 1:
    print "There are more than 1 ungatedfitsfile; aborting\n"
    sys.exit()
ungatedfitsfile = ungatedfitsfile[0].split('/')[-1].strip()
ungatedfitsfile = expdir + '/' + ungatedfitsfile

if prepare:
    print("run prepare_astrometric_epoch.py...\n")
    os.system('prepare_astrometric_epoch.py %s.vex' % experiment)

auxdir    = os.environ['PSRVLBAUXDIR']
configdir = auxdir + '/configs'
yamlfile  = configdir + '/' + experiment + '.yaml'
if not os.path.exists(yamlfile):
    print("%s doesn't exist\n" % yamlfile)
    yamlmodels = glob.glob(r'%s/%s*.yaml' % (configdir, experiment[:5]))
    yamlmodels.sort()
    if len(yamlmodels) == 0:
        print("no %s*.yaml model is found" % experiment[:5])
        alternative_model_expno = raw_input("Input the alternative expno for model yaml file: ")
        yamlmodels = [configdir + '/' + alternative_model_expno + '.yaml']
        if not os.path.exists(yamlmodels[-1]):
            print('model for yaml file does not exist; aborting')
            sys.exit()
    print("copy %s to %s...\n" % (yamlmodels[-1], yamlfile))
    os.system('cp %s %s' % (yamlmodels[-1], yamlfile))

sourcefile = expdir + '/' + experiment + '.source'
print sourcefile
bindir = expdir + '/bin_astrometric_run_results'
if not os.path.exists(bindir):
    os.system('mkdir %s' % bindir)

for gatedfitsfile in gatedfitsfiles:
    gatedfitsfile = gatedfitsfile.split('/')[-1].strip()
    binumber = gatedfitsfile.split('_')[3].lower()
    gatedfitsfile = expdir + '/' + gatedfitsfile
    print gatedfitsfile
    print("writing %s.source using %s as gatedfitsfile...\n" % (experiment, gatedfitsfile))
    writefile = open(sourcefile,'w')
    writefile.write('GATED FITSFILE: %s\n' % gatedfitsfile)
    writefile.write('UNGATED FITSFILE: %s\n' % ungatedfitsfile)
    writefile.write('INBEAM 0 FITSFILE: %s\n' % ungatedfitsfile)
    writefile.write('BANDPASS CAL: J1733-1304\n')
    writefile.write('NUMBER OF TARGETS: 1\n')
    writefile.write('TARGET 0 NAME: XTEJ1810-197\n')
    writefile.write('TARGET 0 PHSREF: J1753-1843\n')
    writefile.write('TARGET 0 INBEAM 0 NAME: J1819-2036\n')
    writefile.close()
    print("run final_astrometric_reduce.py using %s as gatedfitsfile\n" % gatedfitsfile)
    os.system('rerunfinalastrometry.py -e %s' % experiment)
    resultdir = bindir + '/' + binumber
    if not os.path.exists(resultdir):
        os.system('mkdir %s' % resultdir)
    print("copy stats and fits file to %s" % resultdir)
    os.system('cp {*.gated*stats,*pipeline*.gated*fits,} %s' % resultdir)
print "successfully fly over"
