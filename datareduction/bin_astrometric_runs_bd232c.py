#!/usr/bin/env python
###########################################################################
###-- for astrometric runs of bin-gated data
###-- tailored for bd232c data
###-- produce source file automatically and copy yaml file
###-- running prepare_astrometric_epoch.py unless -n
###-- NOTICE: RUN UNDER EXPDIR!!! (e.g. bd232c)
############################################################################
import os,sys,glob,yaml
from optparse import OptionParser

usage = "usage: run under EXPDIR\n%prog []\n-n --noprepare\n-b --binno\n-f --freq\n-h or --help for more"
parser = OptionParser(usage)
parser.add_option("-f", "--freq", dest="freq", default="2", help="2 and 4 represent 2cm and 4cm respectively")
parser.add_option("-n", "--noprepare", dest="prepare", default="True",
                  action="store_false",help="NOT running prepare_astrometric_epoch.py beforehand")
parser.add_option("-b", "--binno", dest="binno", default=-1,
                  help="choose only one bin to reduce, input bin number")
(options, junk) = parser.parse_args()
prepare         = options.prepare
binno           = int(options.binno)
freq            = options.freq

expdir = os.getcwd()
experiment = expdir.split('/')[-1]
gatedfitsfiles = glob.glob(r'./*_gated*%scm*idifits' % freq)
if len(gatedfitsfiles) == 0:
    gatedfitsfiles = glob.glob(r'./*_binbd*%scm*idifits' % freq)    
gatedfitsfiles.sort()
#print(gatedfitsfiles)
#sys.exit()
if binno != -1:
    gatedfitsfiles = glob.glob(r'./*_gated*%scm*BIN%d*idifits' % (freq, binno))
    if len(gatedfitsfiles) != 1:
        gatedfitsfiles = glob.glob(r'./*_gated*%scm*BIN0%d*idifits' % (freq, binno))
ungatedfitsfiles = glob.glob(r'./*_ungated*%scm*idifits' % freq)
if len(ungatedfitsfiles) != 1:
    print "There are other than 1 ungatedfitsfile; aborting\n"
    sys.exit()
ungatedfitsfile = ungatedfitsfiles[0]
#for fitsfile in ungatedfitsfiles:
#    if 'J1818' in fitsfile.split('/')[-1]:
#        ungatedfitsfile = fitsfile.split('/')[-1].strip()
#        ungatedfitsfile = expdir + '/' + ungatedfitsfile
#    else:
#        IGR_fitsfile = fitsfile.split('/')[-1].strip()
#        IGR_fitsfile = expdir + '/' + IGR_fitsfile

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
bindir = expdir + '/bin_astrometric_run_results_' + freq + 'cm'
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
    #writefile.write('INBEAM 1 FITSFILE: %s\n' % IGR_fitsfile)
    writefile.write('BANDPASS CAL: J1549+0237\n')
    writefile.write('NUMBER OF TARGETS: 1\n')
    writefile.write('TARGET 0 NAME: J1818-1607\n')
    writefile.write('TARGET 0 PHSREF: J1825-1718\n')
    writefile.write('TARGET 0 INBEAM 0 NAME: J1805-1408\n')
    #writefile.write('TARGET 0 INBEAM 1 NAME: J1817-1621\n')
    writefile.close()
    print("run final_astrometric_reduce.py using %s as gatedfitsfile\n" % gatedfitsfile)
    os.system('rerunfinalastrometry.py -e %s' % experiment)
    resultdir = bindir + '/' + binumber
    if not os.path.exists(resultdir):
        os.system('mkdir %s' % resultdir)
    print("copy stats and fits file to %s" % resultdir)
    os.system('cp {*.gated*stats,*pipeline*.gated*fits,} %s' % resultdir)
print "successfully fly over"
