#!/usr/bin/env python
import sys,os,glob
class making_frequency_dependent_calibrator_model_using_CASA:
    """
    load the *pipeline_uv.fits (uvfits files should have been copied to the folder);
    TCLEAN jointly
    usage:
    start from table_number=0, then incrementally moving further, then TCLEAN(N, 100, 'data') in the end to make the image 
    """
    def __init__(s):
        s.uvfitsfiles = glob.glob(r'*pipeline_uv.fits')
        s.uvfitsfiles.sort()
        print((s.uvfitsfiles))
        s.MS0s = [uvfitsfile.replace('fits', 'ms') for uvfitsfile in s.uvfitsfiles]
        print(s.MS0s)
        i = 20
        s.MSs = []
        while s.MSs == [] or len(s.MSs) != len(s.MS0s):
            s.MSs = glob.glob(r'*pipeline_uv.v%d.ms' % i)
            i -= 1
            if i < 0:
                s.MSs = s.MS0s
                break
        print(s.MSs)
    def source_name(s):
        path = os.getcwd()
        sourcename = path.split('/')[-1]
        return sourcename
    def load_uvfitsfiles_into_CASA(s):
        script_write = open('.junkscript.py', 'w')
        for uvfitsfile in s.uvfitsfiles:
            MS = uvfitsfile.replace('fits', 'ms')
            print(MS)
            os.system("rm -rf %s" % MS)
            script_write.write("importuvfits(fitsfile='%s', vis='%s')\n" % (uvfitsfile, MS))
        script_write.close()
        os.system('casa -c .junkscript.py')
    def TCLEAN_jointly(s, image_number=1, niter_value=100, savemodel_value='modelcolumn', overwrite=False):
        """
        image_number always 0 unless wanting to start from scratch
        """
        imagefoldername = s.source_name()+'v'+str(image_number)
        print(imagefoldername)
        if overwrite:
            os.system("rm -rf %s*" % imagefoldername)
        script_write = open('.junkscript.py', 'w')
        #if table_number > 0:
        #    datacolumn_value = 'data'
        #    savemodel_value = 'modelcolumn'
        #else:
        #    datacolumn_value = 'corrected'
        #    savemodel_value = 'none'
        script_write.write("tclean(vis=%s, imagename='%s', datacolumn='data', imsize=[256,256], cell=[0.0001,0.0001], deconvolver='mtmfs',\
            nterms=2, threshold='0mJy', interactive=True, niter=%d, gain=0.1, savemodel='%s', restart=True)\n" % (s.MSs, imagefoldername, niter_value, savemodel_value))
        script_write.close()
        os.system('casa -c .junkscript.py')
    def GAINCAL_and_APPLYCAL_and_SPLIT(s, table_number, calmode_value='p', solint_value='30s'):
        """
        table_number=0,1,2,...
        """
        MS0s = s.MS0s
        script_write = open('.junkscript.py', 'w')
        if table_number > 0:
            suffix = '.v' + str(table_number) + '.ms'
            old_MSs = [MS0.replace('.ms', suffix) for MS0 in MS0s]
        else:
            suffix = '.ms'
            old_MSs = MS0s
        print(old_MSs)
        s.MSs = []
        for MS in old_MSs:
            epoch = MS.split('_')[0]
            caltable_value = epoch + '.' + calmode_value + str(table_number) + '.table'
            print(caltable_value)
            os.system("rm -rf %s" % caltable_value)
            if calmode_value == 'p':
                script_write.write("gaincal(vis='%s', calmode='%s', solint='%s', gaintype='G', caltable='%s', minsnr=0)\n" % (MS,\
                    calmode_value, solint_value, caltable_value))
            elif calmode_value == 'ap':
                script_write.write("gaincal(vis='%s', calmode='%s', solint='%s', gaintype='T', caltable='%s', minsnr=0, solnorm=True, normtype='median')\n" % (MS,\
                    calmode_value, solint_value, caltable_value))
            script_write.write("applycal(vis='%s', gaintable='%s')\n" % (MS, caltable_value))
            next_MS = MS.replace(suffix, '.v'+str(table_number+1)+'.ms')
            s.MSs.append(next_MS)
            print(next_MS)
            os.system("rm -rf %s " % next_MS)
            script_write.write("split(vis='%s', outputvis='%s', datacolumn='corrected')\n" % (MS, next_MS))
        script_write.close()
        os.system('casa -c .junkscript.py')
    def plotms(s, vis_value, xaxis_value, yaxis_value):
        script_write = open('.junkscript.py', 'w')
        script_write.write("plotms(vis='%s', xaxis='%s', yaxis='%s')\n" % (vis_value, xaxis_value, yaxis_value))
        script_write.close()
        os.system('casa -c .junkscript.py')

        
        

