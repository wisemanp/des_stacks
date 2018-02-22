# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import seaborn as sns
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import datetime
import configparser
import os
import glob
import logging
import time
import subprocess
def reduce_info(info,**kwargs):
    pass
    #Get a smaller info dataframe based on kwargs

def make_good_frame_list(stack,field,band,zp_cut = -0.15, psf_cut = 2.5):
    """Returns a list of images for a certain chip that are of quality better than zp_cut and psf_cut.
    Arguments:
    field (str): (e.g. 'X2')
    band (str): (e.g. 'g')

    Keyword arguments:
    sig (float): the residual from mean ZP to clip at (float, default = -0.15)

    Returns:
    good_exps (DataFrame): A reduced DataFrame of good exposures for each (field,band)
    """
    ## 1
    ## Get median ZP for each exposure
    ## And calculate residual for each chip compared to that median
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    if zp_cut>0:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info('Initiating make_good_frame_list.py')
    logger.info('Finding good frames for {0}, in {1}, clipping ZP at {2}'.format(field,band,zp_cut))

    logger.info('Getting median zeropoint for each exposure, calculating residual for each image')

    info = stack.info_df
    import math
    info = info[info['FIELD']==field]
    logger.debug('These are the bands available for field {0}'.format(field))
    logger.debug(info.BAND.unique())
    info = info[info['BAND']==band]
    info['ZP_EXPRES']=''
    for counter,exp in enumerate(info.EXPNUM.unique()):
        this_exp = info[info['EXPNUM']==exp]
        exp_idx = this_exp.index
        med_zp = this_exp['CHIP_ZERO_POINT'].median()
        info.loc[exp_idx,'ZP_EXPRES'] = this_exp['CHIP_ZERO_POINT']-med_zp

    logger.debug('Here is the column of initial residuals for each exposure')
    logger.debug(info['ZP_EXPRES'])
    ######################################################
    ## 2
    ## Calculate the average residual from part 1 over all exposures for that chip (field,band)
    logger.info('Getting the median zeropoint residual for all exposures for each chip')
    info['ZP_ADJ1'] = ''
    info['ZP_SIG_ADJ1'] = ''
    for counter, chip in enumerate(info.CCDNUM.unique()):

        this_chip = info[info['CCDNUM']==chip]
        chip_idx = this_chip.index
        med_zp = this_chip['ZP_EXPRES'].median()
        sig_chip_zps = this_chip['ZP_EXPRES'].std()
        ## 3
        ## Subtract that average residual from the given ZP to give an adjusted ZP
        info.loc[chip_idx,'ZP_ADJ1']=this_chip['CHIP_ZERO_POINT']-med_zp
        info.loc[chip_idx,'ZP_SIG_ADJ1']=sig_chip_zps
    logger.debug('Here is the column of adjusted zeropoints')
    logger.debug(info['ZP_ADJ1'])
    #####################################################
    ## 4
    ## Subtract the average ZP for each year off the individual zps for that year

    years =['Y1','Y2','Y3','Y4']
    info['ZP_RES']=''
    for year in years:
        logger.info('Subtracting the median zeropoint for {0} from all exposures that year'.format(year))
        this_year = info[info['YEAR']==year]
        year_idx = this_year.index
        year_med = this_year['ZP_ADJ1'].median()
        year_sig = this_year['ZP_ADJ1'].std()
        final_resid = this_year['ZP_ADJ1']-year_med
        info.loc[year_idx,'ZP_RES']=final_resid

    #####################################################################
    ## 5
    ## Now cut exposures (field,band,chip) based on whether they make the cut and return them
    logger.info('Getting rid of exposures whose ZP residual is below {0}'.format(zp_cut))
    exps = info.EXPNUM.unique()
    zp_cut     = float(zp_cut)
    psf_cut = float(psf_cut)
    nbad = 15
    good_exps = []
    good_frame = pd.DataFrame()
    for exp in exps:
        this_exp = info[info['EXPNUM']==exp]
        logger.debug('Cutting in exposure {0}'.format(exp))
        resids = this_exp['ZP_RES']
        resids = resids.as_matrix()

        bad_resids = 0
        reformatted_resids = []
        for i in range(len(resids)):
            res = resids[i]
            res = float(res)
            reformatted_resids.append(res)
            if float(res)<zp_cut:
                bad_resids +=1
        this_exp['ZP_RES']=np.array(reformatted_resids)
        logger.debug('Number of frames in exposure {0} that fail the ZP cut: {1}'.format(exp,bad_resids))
        bads =bad_resids+len(this_exp[this_exp['PSF_NEA']>psf_cut])
        logger.info("Number of frames failing the combined ZP and PSF cut: %s" %bads)
        if bads <nbad:
            #is a good frame
            logger.debug('...is a good frame!')
            good_exps.append(exp)

            good_frame = good_frame.append(this_exp)
    logger.info("%s exposures were rejected!" %(len(exps)-len(good_exps)))
    ## Save results
    np.savetxt(os.path.join(stack.list_dir,'good_exps_%s_%s.txt'%(field,band)),good_exps,fmt='%s')
    try:
        good_table = Table.from_pandas(good_frame.drop(['ZP_RES','ZP_EXPRES','ZP_ADJ1','ZP_SIG_ADJ1'],axis=1))
    except ValueError:
        logger.exception("Not enough good frames!!")
    logger.debug('Here is the good_table, to write to fits format')
    logger.debug(good_table)
    good_fn = os.path.join(stack.list_dir,'good_exps_%s_%s_%s_%s.fits'%(field,band,zp_cut,psf_cut))
    logger.info('Writing out good exposure list to {0}'.format(good_fn))
    try:
        os.remove(good_fn)
    except OSError:
        pass
    good_table.write(good_fn)
    return good_frame

def make_swarp_cmd(stack,MY,field,chip,band,logger = None,zp_cut = -0.15,psf_cut =2.5,final=True):
    if not os.path.isdir(os.path.join(stack.out_dir,'MY%s'%MY,field,band)):
        os.mkdir(os.path.join(stack.out_dir,'MY%s'%MY,field,band))
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    """function to make swarp command to stack Nminus1_year, field chip, band"""
    logger.info('Preparing the frames to be stacked')
    #band = band + '    '
    ## Get the list of exposures for MY, field, chip, band.
    good = stack.good_frame
    good_band = good[good['BAND']==band]
    if MY == 'none':
        good_band_my =good_band
    else:

        good_band_my = good_band[good_band['YEAR']!='Y{0}'.format(MY)]
    good_my_exps = good_band_my['EXPNUM'].unique()
    #for each good exposure, find it's file
    stack_fns = {}
    logger.info('Adding files to the stack')
    good_band_my.sort_values('CHIP_ZERO_POINT',ascending=False,inplace=True)
    n,l=0,0
    stack_fns[0]=[]
    for counter,exp in enumerate(good_my_exps):

        this_exp = good_band_my[good_band_my['EXPNUM']==exp]
        first = this_exp.iloc[0]
        night = first['NITE']
        #chip = first['CCDNUM']
        this_exp_fn = get_dessn_obs(stack,field,band,night,exp,chip,logger)
        #logger.info("Adding file from %s" %night)
        if this_exp_fn:
            for fn in this_exp_fn:
                n+=1
                stack_fns[l].append(fn)
                if n%100.0 == 0:
                    l+=1
                    stack_fns[l]=[]
    logger.info('Added {} files'.format(counter))
    cmd_list = {}
    for j in range(0,l+1):
        fns = np.array(stack_fns[j])
        fn_list = os.path.join(stack.temp_dir,\
        'stack_fns_MY%s_%s_%s_%s_%.3f_%s_%s.lst' %(MY,field,band,chip,zp_cut,psf_cut,j))
        logger.info('Saving list of files to stack at {0}'.format(fn_list))
        np.savetxt(fn_list,fns,fmt='%s')

        if final == True:
            fn_out = os.path.join(stack.out_dir,'MY%s'%MY,field,band)+\
            '/ccd_%s_%s_%.3f_%s_%s.fits'%(chip,band,zp_cut,psf_cut,j)
        else:
            fn_out = os.path.join(stack.out_dir,'MY%s'%MY,field,band)+\
            '/ccd_%s_%s_%.3f_%s_%s_temp.fits'%(chip,band,zp_cut,psf_cut,j)
        
        weightlist_name = os.path.join(stack.list_dir,'%s_%s_%s_%s_%.3f_%s_%s.wgt.lst'%(MY,stack.field,stack.band,chip,zp_cut,psf_cut,j))
        resamplist_name = os.path.join(stack.list_dir,'%s_%s_%s_%s_%.3f_%s_%s.resamp.lst'%(MY,stack.field,stack.band,chip,zp_cut,psf_cut,j))
        if not os.path.isfile(weightlist_name):
            weightlist_name,resamplist_name = make_weightmap(stack,fn_list,MY,chip,zp_cut,psf_cut,j,logger)
        if os.path.isfile(fn_out):
            cmd_list[j] = (False,fn_out)
        else:
            cmd_list[j]=(['swarp','-IMAGEOUT_NAME','{0}'.format(fn_out),
         '@%s'%resamplist_name,'-c','default.swarp','-COMBINE_TYPE',
         'WEIGHTED','-WEIGHT_TYPE','MAP_WEIGHT',
         '-RESCALE_WEIGHTS','N','-WEIGHT_IMAGE','@%s'%weightlist_name],fn_out)

    logger.info(cmd_list)
    return cmd_list
#############################################
def get_des_obs_year(night,logger=None):
    if not logger:
        logger = logging.getLogger(__name__)
        logger.handlers =[]
        logger.setLevel(logging.DEBUG)
        formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    night = int(night)
    cp=configparser.ConfigParser()
    # read the .ini file
    cp.read('/media/data1/wiseman/des/coadding/config/snobs_params.ini')
    # Make a list of years
    years= ['Y1','Y2','Y3','Y4']
    year_night_lims = {}
    for y in years:
        year_night_lim = cp.get('year_night_lims',y)
        year_night_lims[y]=[int(lim.strip()) for lim in year_night_lim.split(',')]
    if ((night > year_night_lims['Y1'][0]) and
        (night < year_night_lims['Y1'][1])):
        year = 'Y1'
    elif ((night > year_night_lims['Y2'][0]) and
          (night < year_night_lims['Y2'][1])):
        year = 'Y2'
    elif ((night > year_night_lims['Y3'][0]) and
          (night < year_night_lims['Y3'][1])):
        year = 'Y3'
    elif ((night > year_night_lims['Y4'][0]) and
          (night < year_night_lims['Y4'][1])):
        year = 'Y4'
    else:
        raise ValueError
    return year
###############################################

def get_dessn_obs(stack, field, band, night, expnum, chipnum,logger=None):
    '''Function to get the filename for a DES image for a
       given field, band, night, chip, and expnum.
       Uses an object of the Stack class.
       Returns path and name of the file requested.'''
    if not logger:
        logger.handlers =[]
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    #------------------------------------
    chipnum = int(chipnum)
    # step 1 - get the year of the observation
    year = get_des_obs_year(night,logger)
    #------------------------------------
    # step 2 - find the directory for this night
    year_dir = stack.data_dirs[year]
    glob_str = year_dir+night+'-r????'
    subdir = glob.glob(glob_str)[-1]+'/'
    # then this field
    field_glob_str = subdir+'*%s_%s*' % (field, band)
    try:
        field_subdir = glob.glob(field_glob_str)[-1]+'/'
    except IndexError:
        #that subdir doesn't have obs for our field in it
        subdir=glob.glob(glob_str)[0]+'/'
        field_glob_str = subdir+'*%s_%s*' % (field, band)
        field_subdir = glob.glob(field_glob_str)[-1]+'/'
    field_subsubdir = field_subdir + os.listdir(field_subdir)[-1]+'/'
    # then the chip
    chip_subdir = '%sccd%02d' % (field_subsubdir, chipnum)
    #------------------------------------
    # step 3 - find the final fits file at the end of the dir structure
    subdir_list = [chip_subdir]
    rabbit_hole = True
    while rabbit_hole:
        curr_subdir = '/'.join(subdir_list)
        try:
            curr_dir_cont = os.listdir(curr_subdir)
        except:
            logger.info('CCD {0} not observed on that night in that band'.format(chipnum))
            return None
        next_subdir = '/'.join([
            curr_subdir,
            curr_dir_cont[0]])
        try:
            next_subdir_cont = os.listdir(next_subdir)
            subdir_list.append(curr_dir_cont[0])
        except:
            obs_dir = next_subdir
            rabbit_hole = False
    #------------------------------------
    # step 4 - CHECK THIS IS THE CORRECT EXPOSURE NUMBER!!!
    obs_fns = []
    for base_obs_fn in os.listdir(os.path.dirname(obs_dir)):
        try:
            obs_fn = os.path.join(os.path.dirname(obs_dir),os.path.basename(base_obs_fn))
            fits_expnum = fits.getheader(obs_fn)['EXPNUM']
            #logger.info('Attempted EXPNUM: {0}; Got: {1}'.format(expnum,fits_expnum))
        except:
            #logger.info(base_obs_fn)
            #logger.info(fits.getheader(obs_fn)['EXPNUM'])
            continue

        if year == 'Y4':
            obs_fn = obs_fn+'[0]'
        if base_obs_fn[:3]!='DSN':
            obs_fns.append(obs_fn)

    return obs_fns

def rn_stack(stack):
    pass
def chiplims(stack):
    import os
    import astropy.io.fits as fits
    from astropy.wcs import WCS

    root_dir = stack.out_dir+'/MYnone'
    fields = ['E1','E2','S1','S2','C1','C2','C3','X1','X2','X3']
    chips = range(1,63)
    chiplims={}
    for f in fields:
        chiplims[f]={}
        for chip in chips:

            if chip not in [2,31,61]:
                print ('Field:', f,'Chip',chip)
                try:
                    img = os.path.join(root_dir,'SN-%s'%f,'r','ccd_%s_r_-0.150_2.5.fits'%chip)
                    h =fits.getheader(img)
                except:
                    img = os.path.join(root_dir,'SN-%s'%f,'r','ccd_%s.fits'%chip)
                    h =fits.getheader(img)
                lenra = h['NAXIS1']
                lendec= h['NAXIS2']
                wcs = WCS(h)
                lims =wcs.all_pix2world([[0,0],[0,lendec],[lenra,0],[lenra,lendec]],1)
                chiplims[f][chip]=lims
def get_y3a1():
    import easyaccess as ea
    conn = ea.connect()
    fcents ={'C1':(54.2743,-27.1116),
    'C2':(54.2743,-29.0884),
    'C3':(52.6484,-28.1000),
    'E1':(7.8744,-43.0096),
    'E2':(9.5000,-43.9980),
    'S1':(42.8200,0.0000),
    'S2':(41.1944,-0.9884),
    'X1':(34.4757,-4.9295),
    'X2':(35.6645,-6.4121),
    'X3':(36.4500,-4.6000)}
    fields = ['E1','E2','S1','S2','C1','C2','C3','X1','X2','X3']
    for f in fields:
        ra,dec = fcents[f]
        RAMIN,RAMAX,DECMIN,DECMAX = ra-1.6,ra+1.6,dec-1.6,dec+1.6
        q = 'select Y3A1_COADD_OBJECT_SUMMARY.COADD_OBJECT_ID, RA, DEC, \
        MAGERR_AUTO_G, MAGERR_AUTO_R, MAGERR_AUTO_I, MAGERR_AUTO_Z,\
        MAGERR_DETMODEL_G, MAGERR_DETMODEL_R, MAGERR_DETMODEL_I, MAGERR_DETMODEL_Z,\
        MAG_AUTO_G, MAG_AUTO_R, MAG_AUTO_I, MAG_AUTO_Z, \
        MAG_DETMODEL_G, MAG_DETMODEL_R, MAG_DETMODEL_I, MAG_DETMODEL_Z, \
        SPREAD_MODEL_G, SPREAD_MODEL_R, SPREAD_MODEL_I, SPREAD_MODEL_Z, \
        SPREADERR_MODEL_G, SPREADERR_MODEL_R, SPREADERR_MODEL_I, SPREADERR_MODEL_Z, \
        FLAGS_G, FLAGS_I, FLAGS_R, FLAGS_Z \
        from Y3A1_COADD_OBJECT_SUMMARY where RA between %s and %s and DEC between %s AND %s and \
        flags_g = 0 and flags_r = 0 and flags_i = 0 and flags_z = 0 and \
        imaflags_iso_g = 0 and imaflags_iso_r = 0 and imaflags_iso_i = 0 and imaflags_iso_z = 0\
        order by DEC, RA'%(RAMIN,RAMAX,DECMIN,DECMAX)

        dat = conn.query_to_pandas(q)
        dat.to_csv('/home/wiseman/y3a1_%s_summary.csv'%f)

def make_weightmap(s,lst,y,chip,zp_cut,psf_cut,j,logger):
    img_list = np.loadtxt(lst,dtype='str')
    starttime=float(time.time())
    logger.info('Creating weightmaps for individual input exposures')
    weightlist = []
    resamplist = []
    for img in img_list:
        imgname = os.path.split(img)[-1]
        imgroot = imgname[:-5]
        if imgroot[-2:]=='fi':
            imgroot = imgroot[:-3]
        weightlist.append(os.path.join(s.temp_dir,imgroot +'.resamp.weight.fits'))
        resamplist.append(os.path.join(s.temp_dir,imgroot+'.resamp.fits'))
    os.chdir(s.temp_dir)
    swarp_cmd = ['swarp','@%s'%lst,'-COMBINE','N','-RESAMPLE','Y','-c','default.swarp']
    p = subprocess.Popen(swarp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    endtime=float(time.time())
    logger.info('Finished creating weightmaps, took %.3f seconds'%(endtime-starttime))
    weightlist = np.array(weightlist)
    weightlist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%.3f_%s_%s.wgt.lst'%(y,s.field,s.band,chip,zp_cut,psf_cut,j))
    np.savetxt(weightlist_name,weightlist,fmt='%s')
    resamplist = np.array(resamplist)
    resamplist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%.3f_%s_%s.resamp.lst'%(y,s.field,s.band,chip,zp_cut,psf_cut,j))
    np.savetxt(resamplist_name,resamplist,fmt='%s')
    return (weightlist_name,resamplist_name)
