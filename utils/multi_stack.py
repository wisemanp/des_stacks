import multiprocessing
import os
import subprocess
from des_stacks.utils.stack_tools import make_swarp_cmd
import time
import numpy as np
import os

class chip(object):
    def __init__(self, chip):
        self.chip = chip

    def __call__(self):
        return 'Done Stacking chip %s' % (self.chip)
    def __str__(self):
        return 'Stacking chip %s' % (self.chip)

def worker(args,chip):
    started = float(time.time)
    y,field,band,logger,cuts,final,chip=[args[i]for i in range(len(args))]
    s.logger.info('Stacking CCD {0}; starting by creating mini-stacks to save time'.format(chip))
    cmd_list = make_swarp_cmd(s,y,field,chip,band,s.logger,cuts,final)
    staged_imgs = []
    for key,value in cmd_list.items():

        cmd,outname = value
        staged_imgs.append(outname)
        if cmd == False:
            s.logger.info("Already stacked this chip with these cuts, going straight to astrometry")
        else:
            s.logger.info('Stacking... please be patient.'.format(cmd))
            os.chdir(s.temp_dir)
            try:
                starttime=float(time.time())
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                outs,errs = p.communicate()
                endtime=float(time.time())

            except (OSError, IOError):
                s.logger.warn("Swarp failed.", exc_info=1)
            s.logger.info('Finish stacking chip {0}'.format(chip))
            s.logger.info('Took %.3f seconds' % (endtime-starttime))
        s.logger.info('Added %s to list of images to make final stack' %outname)
    s.logger.info('Now combining mini-stacks into final science frame')
    staged_list = np.array(staged_imgs)
    s.logger.info('Combining these frames:')
    s.logger.info(staged_list)
    staged_listname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_staged.lst'%(y,s.field,s.band,chip,s.cutstring))
    np.savetxt(staged_listname,staged_list,fmt='%s')
    resamp_cmd =['swarp','@%s'%staged_listname,'-COMBINE','N','-RESAMPLE','Y','-c','default.swarp']
    os.chdir(s.band_dir)
    s.logger.info('Resampling and weighting the intermediate images:\n %s'%resamp_cmd)
    res_start = float(time.time())
    rf = subprocess.Popen(resamp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    r_out,r_errs = rf.communicate()
    res_end = float(time.time())
    s.logger.info('Done resampling intermediate stacks, took %.3f seconds'%(res_end-res_start))
    resamplist = []
    weightlist = []
    for img in staged_list:
        imgname = os.path.split(img)[-1]
        imgnameroot = imgname[:-5]
        resamplist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.fits'))
        weightlist.append(os.path.join(s.band_dir,imgnameroot+'.resamp.weight.fits'))
    final_resampname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_final.lst'%(y,s.field,s.band,chip,s.cutstring))
    final_weightname = os.path.join(s.temp_dir,'%s_%s_%s_%s_%s_final.wgt.lst'%(y,s.field,s.band,chip,s.cutstring))
    np.savetxt(final_resampname,resamplist,fmt='%s')
    np.savetxt(final_weightname,weightlist,fmt='%s')
    imgout_name = staged_list[0][:-7]+'_sci.fits'
    weightout_name = staged_list[0][:-7]+'_wgt.fits'
    final_cmd = ['swarp','@%s'%final_resampname,'-IMAGEOUT_NAME',imgout_name,'-c','default.swarp','-WEIGHTOUT_NAME',weightout_name,'-COMBINE_TYPE','WEIGHTED','-WEIGHT_IMAGE','@%s'%final_weightname]
    s.logger.info('Doing this command to do the final stack:\n %s'%final_cmd)
    final_start = float(time.time())
    pf = subprocess.Popen(final_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    f_out,f_errs = pf.communicate()
    final_end = float(time.time())
    s.logger.info("Done combining mini-stacks, took %.3f seconds"%(final_end -final_start))
    s.logger.info("Saved final science frame at %s"%imgout_name)
    s.logger.info("And final weightmap at %s"%weightout_name)
    t_tot = float(time.time()) - started
    return t_tot

def start_swarp(s):
    s.logger.info("Starting", multiprocessing.current_process().name)

def multitask(s,y,field,band,logger,cuts,final):

    n_chips = len(s.chips)
    pool_size = multiprocessing.cpu_count()*2
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer = start_swarp(s),
                                maxtasksperchild=2,
                                )
    args = (s,y,field,band,logger,cuts,final)

    results = [pool.apply(worker,args=(args,chip,)) for chip in s.chips]
    pool.close()
    pool.join()
    return results
                            
