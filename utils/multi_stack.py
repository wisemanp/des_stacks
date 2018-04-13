import pathos.pools as pp
import multiprocessing
from multiprocessing import Process
import os
import subprocess
from des_stacks.utils.stack_tools import make_swarp_cmd
from des_stacks.utils.sex_tools import sex_for_psfex, psfex, sex_for_cat
import time
import numpy as np
import os
from itertools import repeat
from functools import partial
import logging
from des_stacks.analysis.astro import astrometry,init_phot
import pandas as pd
def stack_worker(arg_pair):
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    '''if zp_cut>0:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:'''
    logger.setLevel(logging.INFO)
    ch.setLevel(logging.INFO)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    chip,args =arg_pair[1],arg_pair[0]
    s,y,field,band,cuts,final,logger2= [args[i]for i in range(len(args))]
    started = float(time.time())
    #logger.info('Stacking chip %s; starting by creating mini-stacks to save time'%chip)
    cmd_list = make_swarp_cmd(s,y,field,chip,band,s.logger,cuts,final)
    #logger.info("Pulled commands list")
    staged_imgs = []
    for key,value in cmd_list.items():

        cmd,outname = value
        if outname != False:

            staged_imgs.append(outname)
        if cmd == False:
            print ('Already stacked this chip with these cuts, going straight to the next chip')
            #logger.info("Already stacked this chip with these cuts, going straight to astrometry")
            pass
        else:
            #s.logger.info('Stacking... please be patient.'.format(cmd))
            os.chdir(s.temp_dir)
            try:

                print ('Stacking chip %s, part %s'%(chip,key))
                starttime=float(time.time())
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                outs,errs = p.communicate()
                endtime=float(time.time())

            except (OSError, IOError):
                #s.logger.warn("Swarp failed.", exc_info=1)
                print ('Swarp failed for some reason in chip %s'%chip)


            print('Finished doing chip %s, part %s. Took %.3f seconds' %(chip,key,endtime-starttime))
        print('Added %s to list of images to make final stack' %outname)
    print ("Finished first SWarp")
    print('Now combining mini-stacks into final science frame')
    staged_list = np.array(staged_imgs)
    print('Combining these frames:')
    print(staged_list)
    staged_listname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_staged.lst'%(y,field,band,chip,s.cutstring))
    np.savetxt(staged_listname,staged_list,fmt='%s')
    resamp_cmd =['swarp','@%s'%staged_listname,'-COMBINE','N','-RESAMPLE','Y','-c','default.swarp']
    os.chdir(s.band_dir)
    #s.logger.info('Resampling and weighting the intermediate images:\n %s'%resamp_cmd)
    print ('Resampling and weighting intermediate images')
    res_start = float(time.time())
    rf = subprocess.Popen(resamp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    r_out,r_errs = rf.communicate()
    res_end = float(time.time())
    print ('Done resampling intermediate images on chip %s in %.3f seconds'%(chip,(res_end-res_start)))
    #s.logger.info('Done resampling intermediate stacks, took %.3f seconds'%(res_end-res_start))
    resamplist = []
    weightlist = []
    for img in staged_list:
        imgname = os.path.split(img)[-1]
        imgnameroot = imgname[:-5]
        resamplist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.fits'))
        weightlist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.weight.fits'))
    final_resampname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_final.lst'%(y,field,band,chip,s.cutstring))
    final_weightname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_final.wgt.lst'%(y,field,band,chip,s.cutstring))
    np.savetxt(final_resampname,resamplist,fmt='%s')
    np.savetxt(final_weightname,weightlist,fmt='%s')
    imgout_name = staged_list[0][:-7]+'_sci.fits'
    weightout_name = staged_list[0][:-7]+'_wgt.fits'
    final_cmd = ['swarp','@%s'%final_resampname,'-IMAGEOUT_NAME',imgout_name,'-c','default.swarp','-WEIGHTOUT_NAME',weightout_name,'-COMBINE_TYPE','WEIGHTED','-WEIGHT_IMAGE','@%s'%final_weightname]
    #s.logger.info('Doing this command to do the final stack:\n %s'%final_cmd)
    print ('Doing final stack!')
    final_start = float(time.time())
    pf = subprocess.Popen(final_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    f_out,f_errs = pf.communicate()
    final_end = float(time.time())

    print("Done combining mini-stacks, took %.3f seconds"%(final_end -final_start))
    print("Saved final science frame at %s"%imgout_name)
    print("And final weightmap at %s"%weightout_name)
    t_tot = float(time.time()) - started
    return t_tot

def sex_worker(arg_pair):
    print ('Attempting to run a sex worker')
    chip,args =arg_pair[1],arg_pair[0]
    s,y,field,band,cuts,final,logger2= [args[i]for i in range(len(args))]
    started = float(time.time())

    sex_for_psfex(s,chip,cuts)

    model_fwhm = psfex(s,chip,retval='FWHM',cuts=cuts)
    os.chdir(s.ana_dir)

    sexcat = sex_for_cat(s,chip,cuts)

    zp,sex_fwhm = astrometry(s,chip,sexcat)

    zp_psf,psf_fwhm = astrometry(s,chip,sexcat,phot_type = 'PSF')

    qual = open(os.path.join(s.ana_dir,'%s_ana.qual'%s.cutstring),'w')
    print ('# Quality parameters for %s %s %s %s' %(s.my,s.field,s.band,chip),file =qual)
    print ('# Parameters:',file=qual)
    print ('# Zeropoint from sextractor',file=qual)
    print ('# Zeropoint from PSF matches', file = qual)
    print ('# FWHM from PSFex',file = qual)
    print ('# FWHM from SExtractor using PSF model',file = qual)
    print ('%s %s %s %s'%(zp,zp_psf,model_fwhm,sex_fwhm),file=qual)
    qual.close()
    print("Written quality factors to %s_ana.qual" %s.cutstring)
    qual_dict = {'zp':zp,'fwhm_psfex':model_fwhm,'fwhm_sex':sex_fwhm}
    qual_df = pd.DataFrame([qual_dict],index = [chip])
    print("Quality of stack:\n %s" %qual_df)
    print("********** Done measuring quality of stack! **********")
    print("******************************************************")
    return sexcat
def multitask(s,y,field,band,cuts,final,w='stack'):
    n_chips = len(s.chips)

    args = [s,y,field,band,cuts,final]

    pool_size = multiprocessing.cpu_count()*2
    #s.logger.info("Starting %s processes"%pool_size)
    act = multiprocessing.active_children()
    pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
    pool._clear()
    pool._serve()

    chips = list(s.chips)
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    args.append(logger)
    all_args = []
    for c in chips:
        all_args.append([args,c])
        #p = Process(target=worker,args=(args,c))
        #p.start()
        #p.join()
    if w =='stack':
        results = pool.map(stack_worker,all_args)
    elif w=='sex':
        results = pool.map(sex_worker,all_args)
    pool.close()
    pool.join()
    return results
