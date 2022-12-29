#!/usr/bin/env python
#
# Script to make a diagnostic Web document.
#
from __future__ import print_function ## this import has to be placed on the top; it makes python3 print works for python2
import os, glob, re, pwd
import subprocess
import datetime
from MarkupPy import markup

months = { 1:'January',  2:'February',  3:'March',
           4:'April',    5:'May',       6:'June',
           7:'July',     8:'August',    9:'September',
          10:'October', 11:'November', 12:'December'}

#convert  = '/usr/bin/convert'
#psselect = '/usr/bin/psselect'
#difmap   = 'difmap'
#mv       = '/bin/mv'

##############################

def plt_uvdata(code):
    """
    Plot the u-v data in a directory, using difmap.
    """
    filelist_a = glob.glob(code.lower()+'*_pipeline_uv.fits')
    filelist_b = glob.glob(code.lower()+'*_pipeline_uv.gated.fits')
    filelist = filelist_a + filelist_b

    plotfile = []

    for file in filelist:
        plotfile.append(file.rstrip('fits') + 'png')
        plotfile.append(file.rstrip('fits') + 'avg.png')

        dfile = open('diagnostic_tmp.par', 'w')
        dfile.write('!\n')
        dfile.write('observe %s\n'%file)
        dfile.write('select I\n')
        pgotfile = plotfile[-2] + '/PNG'
        dfile.write('device %s\n'%pgotfile)
        dfile.write('radplot\n')
        dfile.write('uvaver 20\n')
        pgotfile = plotfile[-1] + '/PNG'
        dfile.write('device %s\n'%pgotfile)
        dfile.write('radplot\n')
        dfile.write('quit\n')
        dfile.close()

        dfile = open('diagnostic_tmp.par', 'r')
        status = subprocess.call('difmap',
                                 stdin=dfile, stdout=subprocess.PIPE)
        dfile.close()

        status = subprocess.call(['mv', plotfile[-2], 'images'])
        status = subprocess.call(['mv', plotfile[-1], 'images'])
        os.unlink('diagnostic_tmp.par')
        os.unlink('difmap.log')

    return(plotfile)
#
def plt_auto():
    """
    Make images from automagically generated plots.

    Gets complicated because convert does not handle multi-page 
    PostScript files well.
    """
    filelist = glob.glob('tables/*.ps')

    for fl in filelist:
        tmp = fl.split('/')
        basefile = tmp[1].rstrip('.ps')

        f = open(fl, 'r')
        il = 0
        for l in f:
            if l.startswith('%%Pages:'):
                il = il + 1
                if il == 2:
                    pages = l.split(':')
        f.close()
        
        for ipage in range(int(pages[1])):
            status = subprocess.call(['psselect', '-p%d'%(ipage+1), fl, 'junk.ps']) ## call works for both python2/3
            outfile = 'images/%s-%02d.png'%(basefile, ipage+1)
            status = subprocess.call(['convert', 'junk.ps', '-colorspace',
                                     'RGB', outfile]) 
    
    os.unlink('junk.ps')

    filelist = glob.glob('images/*-*.png')
    return(filelist)
#
def plt_images(code):
    """
    Make figures of images.
    """
    filelist_a = glob.glob(code+'_*_difmap*fits*ii.a.clean')
    filelist_b = glob.glob(code+'_*_difmap.*gated.fits.ii.a.')
    filelist = filelist_a + filelist_b
    plotfile = []

    for f in filelist:
        plotfile.append('images/' + f + '.png')
        infile = 'ps:' + f
        status = subprocess.call([convert, '-rotate', '90', infile, plotfile[-1]])

    return(plotfile)

#
def parse_log(code):
    """
    Parse the operator log, looking for stations at which data 
    were lost or potentially affected.
    """
    dataloss = []

    databad = re.compile('\*[BFHKLMNOPS][A-Z]')
    datamay = re.compile('\%[BFHKLMNOPS][A-Z]')

    infile = os.environ['PSRVLBAUXDIR'] + '/runlogs/' + code.lower() + 'log.vlba'
    if os.path.exists(infile):
        f = open(infile, 'r')
        for l in f:
            if databad.search(l):
                dataloss.append(l)
            if datamay.search(l):
                dataloss.append(l)
    else:
        dataloss.append("No log found")

    return(dataloss)

##############################
#
# Main routine
#
if __name__ == "__main__":

    print("Making diagnostic plots ...")
    #
    # Basic stuff
    #
    code = os.getcwd().split('/')[-1]
    code = code.upper()

    #tmp = glob.glob('VLBA_'+code+'_gated*.idifits')[0]
    #tmp = tmp.split('_')[6]
    #tmp = tmp.split('.')[0]
    #year  = int(tmp[:2])
    #month = int(tmp[2:4])
    #day   = int(tmp[4:6])
    tmp = glob.glob('*.vex')[0]
    for line in open(tmp):
        if 'date     :' in line:
            datestring = line.split()
            year  = int(datestring[-1])
            month = datestring[-2]
            day   = int(datestring[-3])
            break

    # Zap any preexisting images
    os.system("rm -f tables/*.png")
    os.system("rm -f images/*png")

    #
    # Make Web document.
    #
    # Start with basic Web doc stuff.

    diag = markup.page()

    title = "Diagnostic output for %s, %d %s %d"%(code, year, month, day)
    diag.init( title=title )

    #
    # operator log parsing
    #
    logcode = code
    if len(logcode) == 8:
        logcode = logcode[:7]
    log_info = parse_log(logcode)

    diag.hr()
    diag.h1('Operator Log Comments')
    if len(log_info) == 0:
        diag.p('No stations reported with potential data issues.')
    else:
        for l in log_info:
            diag.pre(l)
    
    #
    # Make image plots.
    #
    images_images = plt_images(code)

    diag.hr()
    diag.h1('images')
    for image in images_images:
        diag.object('%s'%image,
                    data='%s'%image, type='image/png', width=512, style='border-width: medium ')
    #
    # Make u-v plots.
    #
    images_uv = plt_uvdata(code)

    diag.hr()
    diag.h1('u-v data')
    for image in images_uv:
        if not 'avg.png' in image:
            diag.object('%s'%image, data='images/%s'%image, type='image/png', width=512)
    diag.hr()
    diag.h1('u-v data (averaged to 20 seconds)')
    for image in images_uv:
        if 'avg.png' in image:
            diag.object('%s'%image, data='images/%s'%image, type='image/png', width=512)
    #
    #
    #
    images_auto = plt_auto()
    images_auto.sort()

    diag.hr()
    diag.h1('ACCOR plots')
    for image in images_auto:
        if image.find('accor') > 0:
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('APCAL plots')
    for image in images_auto:
        if image.find('apcal') > 0:
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('fringe finder amplitude calibration plots')
    for image in images_auto:
        if image.find('ampcalfring') > 0:
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('phase reference calibration plots')
    for image in images_auto:
        if image.find('phsreffring') > 0:
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('phase reference amplitude selfcal plots (amplitude)')
    lastsrc = ''
    for image in images_auto:
        if (image.find('calibapn.am')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('phase reference amplitude selfcal plots (phase)')
    lastsrc = ''
    for image in images_auto:
        if (image.find('calibapn.ph')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('primary in-beam phase calibration plots')
    lastsrc = ''
    for image in images_auto:
        if (image.find('icalib.p1')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('primary in-beam phase calibration plots (separate IFs)')
    lastsrc = ''
    for image in images_auto:
        if (image.find('icalib.pn')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('primary in-beam amp+phase calibration plots')
    lastsrc = ''
    for image in images_auto:
        if (image.find('icalib.ap1')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('secondary (daisy-chained) in-beam amp+phase calibration plots')
    lastsrc = ''
    for image in images_auto:
        if (image.find('icalib.sp1')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    diag.h1('Scintillation correction plots')
    lastsrc = ''
    for image in images_auto:
        if (image.find('scintcorr')) > 0:
            srcname = image.split('/')[-1].split('.')[0]
            if not srcname in lastsrc:
                lastsrc = srcname
                diag.p((srcname + ":"))
            diag.object('%s'%image, data='%s'%image, type='image/png', width=384)

    #
    # Get relevant info from the log
    # 
    os.system("make_pipeline_summary.py %s %s/../" % (code.lower(), os.getcwd()))
    diag.hr()
    diag.h1('Important messages from the pipeline log')
    summarylines = open(code.lower() + '.pipelinesummary').readlines()
    oneline = ""
    for line in summarylines:
        oneline += line
    diag.pre(oneline)

    #
    # Footer
    #
    diag.hr()
    now = datetime.date.today()
    myname = pwd.getpwuid(os.getuid()).pw_name
    timestamp = 'make_diagnostic: automagically generated %d %s %02d by %s running in %s'%( now.year, months[now.month], now.day, myname, os.getcwd())
    diag.p(timestamp)
    #        
    # Print out HTML.
    #
    print(diag, file=open("Diagnostic.html", 'w')) ## python3 way
    #print >> open("Diagnostic.html", 'w'), diag ## python2 way
