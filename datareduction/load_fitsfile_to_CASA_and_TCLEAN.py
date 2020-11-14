#!/usr/bin/env python
import sys,os,glob
class making_frequency_dependent_calibrator_model_using_CASA:
    """
    load the *pipeline_uv.fits (uvfits files should have been copied to the folder);
    TCLEAN jointly 
    """
    def __init__(s):
        s.uvfitsfiles = glob.glob(r'*pipeline_uv.fits')
        s.uvfitsfiles.sort()
        print(s.uvfitsfiles)
        s.MSs = [uvfitsfile.replace('fits', 'ms') for uvfitsfile in s.uvfitsfiles]
        print s.MSs
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
    def TCLEAN_jointly(s, table_number=1, niter_value=100):
        imagefoldername = s.source_name()+'v'+str(table_number)
        print(imagefoldername)
        os.system("rm -rf %s" % imagefoldername)
        script_write = open('.junkscript.py', 'w')
        script_write.write("tclean(vis=%s, imagename='%s', imsize=[512,512], cell=[0.0005,0.0005], deconvolver='mtmfs',\
            nterms=2, threshold='0mJy', interactive=True, niter=%d, gain=0.1, savemodel='modelcolumn')\n" % (s.MSs, imagefoldername, niter_value))
        script_write.close()
        os.system('casa -c .junkscript.py')
    def GAINCAL_and_APPLYCAL_and_SPLIT(s, table_number=1, calmode_value='p', solint_value='6s'):
        script_write = open('.junkscript.py', 'w')
        for MS in s.MSs:
            epoch = MS.split('_')[0]
            caltable_value = epoch + '.' + calmode_value + str(table_number) + '.table'
            print(caltable_value)
            os.system("rm -rf %s" % caltable_value)
            script_write.write("gaincal(vis='%s', calmode='%s', solint='%s', gaintype='G', caltable='%s', minsnr=5)\n" % (MS,\
                calmode_value, solint_value, caltable_value))
            script_write.write("applycal(vis='%s', gaintable='%s')\n" % (MS, caltable_value))
            splitted_MS = MS.replace('.ms', '.'+calmode_value+str(table_number)+'.ms')
            print(splitted_MS)
            os.system("rm -rf %s " % splitted_MS)
            script_write.write("split(vis='%s', outputvis='%s', datacolumn='corrected')\n" % (MS, splitted_MS))
        script_write.close()
        os.system('casa -c .junkscript.py')
    def plotms(s, vis_value, xaxis_value, yaxis_value):
        script_write = open('.junkscript.py', 'w')
        script_write.write("plotms(vis='%s', xaxis='%s', yaxis='%s')\n" % (vis_value, xaxis_value, yaxis_value))
        script_write.close()
        os.system('casa -c .junkscript.py')

        
        

