# -*- coding: utf-8 -*-

import pathos.pools as pp
import multiprocessing
from multiprocessing import Process
import os
import subprocess
import time
import numpy as np
import os
import logging
import pandas as pd
import astropy.io.fits as fits

from des_stacks.utils.stack_tools import make_swarp_cmds, combine_mask_weight
from des_stacks.utils.source_tools import source_for_psfex, psfex, source_for_cat
from des_stacks.analysis.astro import init_phot, init_calib

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
    s,logger2= [args[i]for i in range(len(args))]
    field,band,y,cuts,final = s.field,s.band,s.my,s.cuts,s.final
    started = float(time.time())
    #logger.info('Stacking chip %s; starting by creating mini-stacks to save time'%chip)
    cmd_list = make_swarp_cmds(s,chip,s.logger,cuts,final)
    #logger.info("Pulled commands list")
    staged_imgs = []
    for key,value in cmd_list.items():

        clip_cmd,wgt_cmd,outname = value
        print (clip_cmd,wgt_cmd,outname)
        if outname:

            staged_imgs.append(outname.replace('clipped','weighted'))
        if clip_cmd == False:
            print ('Already stacked chip %s, part %s with these cuts, going straight to the next chip'%(chip,key))
            #logger.info("Already stacked this chip with these cuts, going straight to astrometry")
            pass
        else:
            #s.logger.info('Stacking... please be patient.'.format(cmd))
            os.chdir(s.temp_dir)
            try:

                print ('Stacking chip %s, part %s, clipped'%(chip,key))
                print ('Current dir: %s'%os.curdir)
                print ('This command: %s'%clip_cmd)
                starttime=float(time.time())
                pcli = subprocess.Popen(clip_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                outs,errs = pcli.communicate()
                endtime=float(time.time())
                print('Finished doing chip %s, part %s, clipped. Took %.3f seconds' %(chip,key,endtime-starttime))
            except (OSError, IOError):
                #s.logger.warn("Swarp failed.", exc_info=1)
                print ('Swarp failed for some reason in chip %s'%chip)
        if outname==False:
            print ('No outname - assume you are skipping this for a reason')
        elif wgt_cmd==False:
            print ('Already done the weighted stack of chip %s, part %s, going to the next, or the combination'%(chip,key))
        else:
            ## Start by making a mask
            mm_conf_name = os.path.join(s.temp_dir,'cliptabs','%s_%s_%s_%s_%s_%s_mm.config'%(y,field,band,chip,s.cutstring,key))
            mm_conf = open(mm_conf_name, 'w')
            stackhead = fits.getheader(outname)
            stackhead_name = outname.replace('.fits','.head')
            try:
                stackhead.totextfile(stackhead_name)
            except:
                pass
            outhead = fits.getheader(outname)
            headlistlst = np.genfromtxt(os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.head.lst'%(y,field,band,chip,s.cutstring,key)),dtype=str,delimiter='\n')
            widths, lengths = [],[]
            for headfname in headlistlst:
                inth = fits.Header.fromtextfile(headfname)
                widths.append(inth['NAXIS1'])
                lengths.append(inth['NAXIS2'])
            xsize,ysize = max(widths),max(lengths)
            params = {
            'outliers': os.path.join(s.temp_dir,'cliptabs','%s_%s_%s_%s_%s_%s_clipped.tab'%(y,field,band,chip,s.cutstring,key)),
            'stackhead': stackhead_name,
            'headlist': os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.head.lst'%(y,field,band,chip,s.cutstring,key)),
            'mask':os.path.join(s.temp_dir,'mask.conf'),
            'masksuffix':'.mask.fits',
            'xsize':'%s'%int(xsize),
            'ysize':'%s'%int(ysize)
            }
            for p in params.keys():
                mm_conf.write('%s  = %s     ;\n'%(p,params[p]))
            mm_conf.close()
            print ('Wrote mm_conf to: %s'%mm_conf_name)
            maskmap_cmd = ['/home/wiseman/software/cliputils/MaskMap']
            try:
                print ('Making mask for chip %s, part %s'%(chip,key))
                print ('Current dir: %s'%os.curdir)
                print ('This command: %s'%maskmap_cmd)
                config_file = open(mm_conf_name)
                print ('Stdin: %s'%mm_conf_name)
                starttime=float(time.time())
                maskp = subprocess.Popen(maskmap_cmd,stdin=config_file,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                outs,errs = maskp.communicate()
                endtime=float(time.time())

                print('Finished masking chip %s, part %s. Took %.3f seconds' %(chip,key,endtime-starttime))
            except (OSError, IOError):
                #s.logger.warn("Swarp failed.", exc_info=1)
                print ('MaskMap failed for some reason in chip %s'%chip)
            ## Combine the masks with the weightmaps
            print ('Combining masks with weightmaps for chip %s, part %s'%(chip,key))
            starttime=float(time.time())
            combine_mask_weight(s,chip,key)
            endtime=float(time.time())
            print('Finished creating masked weightmaps for chip %s, part %s. Took %.3f seconds' %(chip,key,endtime-starttime))
            ##  And here the weighted stack
            try:
                print ('Stacking chip %s, part %s, weighted'%(chip,key))
                print ('Current dir: %s'%os.curdir)
                print ('This command: %s'%wgt_cmd)
                starttime=float(time.time())
                pwgt = subprocess.Popen(wgt_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                outs,errs = pwgt.communicate()
                endtime=float(time.time())
                print('Finished doing chip %s, part %s, weighted. Took %.3f seconds' %(chip,key,endtime-starttime))
            except (OSError, IOError):
                #s.logger.warn("Swarp failed.", exc_info=1)
                print ('Swarp failed for some reason in chip %s'%chip)

        if outname!=False:
            print('Added %s to list of images to make final stack' %outname.replace('clipped','weighted'))

    print('Now combining mini-stacks into final science frame for chip %s'%chip)
    staged_list = np.array(staged_imgs)
    print('Combining these frames for chip %s:'%chip)
    print(staged_list)
    staged_listname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_staged.lst'%(y,field,band,chip,s.cutstring))
    np.savetxt(staged_listname,staged_list,fmt='%s')
    imgout_name = staged_list[0][:-16]+'_clipweighted_sci.fits'
    resamp_cmd =['swarp',
    '@%s'%staged_listname,
    '-COMBINE','N',
    '-RESAMPLE','Y']
    os.chdir(s.band_dir)
    #s.logger.info('Resampling and weighting the intermediate images:\n %s'%resamp_cmd)
    print ('Making weights for intermediate images with: %s'%resamp_cmd)
    res_start = float(time.time())
    rf = subprocess.Popen(resamp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    r_out,r_errs = rf.communicate()
    res_end = float(time.time())
    print ('Done weighting intermediate images on chip %s in %.3f seconds'%(chip,(res_end-res_start)))
    #s.logger.info('Done resampling intermediate stacks, took %.3f seconds'%(res_end-res_start))
    resamplist = []
    weightlist = []
    for img in staged_list:
        imgname = os.path.split(img)[-1]
        imgnameroot = imgname[:-5]
        resamplist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.fits'))
        weightlist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.weight.fits'))
    final_resampname = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_final.lst'%(y,field,band,chip,s.cutstring))
    final_weightname = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_final.wgt.lst'%(y,field,band,chip,s.cutstring))
    np.savetxt(final_resampname,resamplist,fmt='%s')
    np.savetxt(final_weightname,weightlist,fmt='%s')

    weightout_name = imgout_name.replace('_sci.fits','_wgt.fits')
    final_cmd = ['swarp','@%s'%final_resampname,'-IMAGEOUT_NAME',imgout_name,
    '-WEIGHTOUT_NAME',weightout_name,
    '-COMBINE_TYPE','WEIGHTED',
    '-WEIGHT_TYPE','MAP_WEIGHT',
    '-WEIGHT_IMAGE','@%s'%final_weightname]
    #s.logger.info('Doing this command to do the final stack:\n %s'%final_cmd)
    print ('Doing final stack with command: %s'%final_cmd)
    final_start = float(time.time())
    pf = subprocess.Popen(final_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    f_out,f_errs = pf.communicate()
    final_end = float(time.time())

    print("Done combining intermediate stacks, took %.3f seconds"%(final_end -final_start))
    print("Saved final science frame at %s"%imgout_name)
    #print("And final weightmap at %s"%weightout_name)
    t_tot = float(time.time()) - started
    return t_tot

def source_worker(arg_pair):
    print ('Attempting to run a source worker')
    chip,args =arg_pair[1],arg_pair[0]
    s,band,logger2= [args[i]for i in range(len(args))]
    field,band,y,cuts,final = s.field,s.band,s.my,s.cuts,s.final
    started = float(time.time())

    source_for_psfex(s,chip,cuts)

    model_fwhm = psfex(s,chip,retval='FWHM',cuts=cuts)
    os.chdir(os.path.join(s.band_dir,str(chip)))

    sourcecat = source_for_cat(s,chip,cuts)

    zp,zp_sig,source_fwhm,source_fwhm_sig = init_calib(s,chip,sourcecat)

    qual = np.array([zp,zp_sig,source_fwhm,source_fwhm_sig])
    qual_fn = os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring)
    np.savetxt(qual_fn,qual)

    print("Written quality factors to %s" %qual_fn)
    qual_dict = {'zp':zp,'fwhm_psfex':model_fwhm,'fwhm_source':source_fwhm}
    qual_df = pd.DataFrame([qual_dict],index = [chip])
    print("Quality of stack:\n %s" %qual_df)
    print("********** Done measuring quality of stack! **********")
    print("******************************************************")
    return sourcecat
def multitask(s,w='stack'):
    args = [s]
    pool_size = multiprocessing.cpu_count()*2
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
    elif w=='source':
        results = pool.map(source_worker,all_args)
    pool.close()
    pool.join()
    return results
