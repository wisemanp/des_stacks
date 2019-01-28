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
import _pickle as cpickle
def reduce_info(info,**kwargs):
    pass
    #Get a smaller info dataframe based on kwargs

def make_good_frame_list(s,field,band,cuts={'teff':0.2, 'zp':None,'psf':None}):
    """Returns a list of images for a certain chip that are of quality better than a given cut.
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
    '''if zp_cut>0:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:'''
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info('Initiating make_good_frame_list.py')
    info = s.info_df
    import math
    info = info[info['FIELD']==field]
    logger.info('These are the bands available for field {0}'.format(field))
    logger.info(info.BAND.unique())
    info = info[info['BAND']==band]
    logger.info(cuts)

    if cuts['zp']!= 'None' and cuts['zp']!= None:
        logger.warning('Gone to do ZP residuals, not sure you want this.')
        logger.info(cuts['zp'])
        info['ZP_EXPRES']=''
        for counter,exp in enumerate(info.EXPNUM.unique()):
            this_exp = info[info['EXPNUM']==exp]
            exp_idx = this_exp.index
            med_zp = this_exp['CHIP_ZERO_POINT'].median()
            info.loc[exp_idx,'ZP_EXPRES'] = this_exp['CHIP_ZERO_POINT']-med_zp

        logger.info('Here is the column of initial residuals for each exposure')
        logger.info(info['ZP_EXPRES'])
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

        years =[1,2,3,4,5]
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
        if cuts['psf']:
            psf_cut = float(cuts['psf'])
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
        good_fn = os.path.join(s.list_dir,'good_exps_%s_%s_%s_%s.fits'%(field,band,zp_cut,psf_cut))
        logger.info("%s exposures were rejected!" %(len(exps)-len(good_exps)))
    ## Save results
    elif -1 < cuts['teff'] < 500:
        logger.info('Doing the cut based on T_eff > %s'%cuts['teff'])
        good_frame = pd.DataFrame()
        good_exps = []

        for counter,exp in enumerate(info.EXPNUM.unique()):

            this_exp = info[info['EXPNUM']==exp]
            exp_idx = this_exp.index
            try:
                t_eff = this_exp['T_EFF'].values[0]
            except:
                t_eff = 2
            if t_eff < -1:
                t_eff = 2
            elif not t_eff:
                t_eff = 2
            try:
                fwhm = this_exp['FWHM_ASEC'].values[0]
            except:
                fwhm = 2
            if fwhm < -1:
                fwhm = 5
            elif not fwhm:
                fwhm = 5
            if not cuts['psf']:
                psf_cut=5
            else:
                psf_cut = cuts['psf']

            if t_eff > cuts['teff'] and len(this_exp[this_exp['PSF_NEA']>psf_cut])<15 :
                this_exp['T_EFF']= t_eff
                good_frame = good_frame.append(this_exp)
                good_exps.append(exp)
        good_fn = os.path.join(s.list_dir,'good_exps_%s_%s_%s_%s.csv'%(field,band,cuts['teff'],cuts['psf']))
    txtname = good_fn[:-4]+'txt'
    np.savetxt(txtname,good_exps,fmt='%s')
    logger.info('Writing out good exposure list to {0}'.format(good_fn))

    good_frame.drop(['ZP_RES','ZP_EXPRES','ZP_ADJ1','ZP_SIG_ADJ1'],axis=1).to_csv(good_fn)
    return good_frame

def make_swarp_cmds(s,MY,field,chip,band,logger = None,cuts={'teff':0.2, 'zp':None,'psf':None},final=True):
    if not os.path.isdir(os.path.join(s.out_dir,'MY%s'%MY,field,band)):
        os.mkdir(os.path.join(s.out_dir,'MY%s'%MY,field,band))
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    """function to make swarp command to stack Nminus1_year, field chip, band"""
    logger.info('Preparing the frames for %s, %s band, chip %s'%(field,band, chip))
    #band = band + '    '
    ## Get the list of exposures for MY, field, chip, band.
    good = s.good_frame
    good_band = good[good['BAND']==band]
    if MY == 'none':
        good_band_my =good_band
    else:

        good_band_my = good_band[good_band['SEASON']!=int(MY)]
    good_my_exps = good_band_my['EXPNUM'].unique()
    #for each good exposure, find it's file
    stack_fns = {}
    logger.info('Adding files to the %s, %s band, chip %s stack'%(field, band,chip))
    good_band_my.sort_values('CHIP_ZERO_POINT',ascending=False,inplace=True)
    n,l=0,0
    stack_fns[0]=[]
    nights = []
    for counter,exp in enumerate(good_my_exps):

        this_exp = good_band_my[good_band_my['EXPNUM']==exp]
        first = this_exp.iloc[0]
        night = first['NITE']
        #chip = first['CCDNUM']
        this_exp_fn = get_dessn_obs(s,field,band,night,exp,chip,logger)
        #logger.info("Adding file from %s" %night)
        if night not in nights:
            if this_exp_fn:
                nights.append(night)
                for fn in this_exp_fn:
                    n+=1
                    stack_fns[l].append(fn)

                    if n%100.0 == 0:
                        l+=1
                        stack_fns[l]=[]
    cmd_list = {}
    for j in range(0,l+1):
        fns = np.array(stack_fns[j])
        fn_list = os.path.join(s.temp_dir,\
        'stack_fns_MY%s_%s_%s_%s_%s_%s.lst' %(MY,field,band,chip,s.cutstring,j))
        logger.info('Saving list of files to stack at {0}'.format(fn_list))
        np.savetxt(fn_list,fns,fmt='%s')

        if final == True:
            fn_out = os.path.join(s.out_dir,'MY%s'%MY,field,band)+\
            '/ccd_%s_%s_%s_%s_clipped.fits'%(chip,band,s.cutstring,j)
        else:
            fn_out = os.path.join(s.out_dir,'MY%s'%MY,field,band)+\
            '/ccd_%s_%s_%s_%s_temp_clipped.fits'%(chip,band,s.cutstring,j)

        weightlist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.wgt.lst'%(MY,s.field,s.band,chip,s.cutstring,j))
        resamplist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.resamp.lst'%(MY,s.field,s.band,chip,s.cutstring,j))
        headerlist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.head.lst'%(MY,s.field,s.band,chip,s.cutstring,j))

        weightout_name = fn_out[:-4]+'wgt.fits'
        nofiles = 0
        if not os.path.isfile(headerlist_name):
            logger.info("%s, %s band, chip %s: Going to do resampling!"%(field,band,chip))
            try:
                weightlist_name,resamplist_name,headerlist_name = resample(s,fn_list,MY,chip,cuts,j,logger)
            except TypeError:
                logger.info("No files in list: %s" %fn_list)
                nofiles = 1
        else:
            logger.info("Header list exists: %s, \n Going to make initial stack."%headerlist_name)
        if os.path.isfile(fn_out.replace('clipped','weighted')):
            cmd_list[j] = (False,False,fn_out)
        else:
            if nofiles ==0:
                cliptab_name = os.path.join(s.temp_dir,'cliptabs','%s_%s_%s_%s_%s_%s_clipped.tab'%(MY,s.field,s.band,chip,s.cutstring,j))
                clip_sigs = {
                'SN-X1':15,
                'SN-X2':15,
                'SN-C1':15,
                'SN-C2':15,
                'SN-S1':15,
                'SN-S2':15,
                'SN-E1':15,
                'SN-E2':15,
                'SN-X3':10,
                'SN-C3':10,
                }
                swarp_clip = [
                'swarp',
                '@%s'%resamplist_name,
                '-RESAMPLE','N',
                '-IMAGEOUT_NAME',fn_out,
                '-CLIP_LOGNAME',cliptab_name,
                '-CLIP_SIGMA',str(clip_sigs[s.field])
                ]
                maskweightlist_name = resamplist_name.replace('resamp','maskweight')


                swarp_weight = [
                'swarp',
                '@%s'%resamplist_name,
                '-RESAMPLE', 'N',
                '-COMBINE_TYPE', 'WEIGHTED',
                '-WEIGHT_TYPE', 'MAP_WEIGHT',
                '-WEIGHT_IMAGE', '@%s'%maskweightlist_name,
                '-IMAGEOUT_NAME', fn_out.replace('clipped.fits','weighted.fits')
                ]
                if os.path.isfile(fn_out):
                    cmd_list[j] = (False,swarp_weight,fn_out)
                else:

                    cmd_list[j]=(swarp_clip,swarp_weight,fn_out)
            else:
                cmd_list[j]=(False,False,False)

    #logger.info(cmd_list)
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
    cp.read('/media/data3/wiseman/des/coadding/config/snobs_params.ini')
    # Make a list of years
    years= ['Y1','Y2','Y3','Y4','Y5']
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
    elif ((night > year_night_lims['Y5'][0]) and
          (night < year_night_lims['Y5'][1])):
        year = 'Y5'
    else:
        raise ValueError
    return year
###############################################

def get_dessn_obs(s, field, band, night, expnum, chipnum,logger=None):
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
    year_dir = s.data_dirs[year]
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
    procs = []
    counts  = []
    if len(os.listdir(field_subdir))==0:
        return None

    for counter,ext in enumerate(os.listdir(field_subdir)):
        procs.append(int(ext[1:]))
        counts.append(counter)
    which = np.argmax(np.array(procs))
    try:
        field_subsubdir = field_subdir + os.listdir(field_subdir)[counts[which]]+'/'
    except:
        field_subsubdir = field_subdir + os.listdir(field_subdir)[-1]+'/'
    #field_subsubdir = field_subdir + os.listdir(field_subdir)[-1]+'/'
    # then the chip
    chip_subdir = '%sccd%02d' % (field_subsubdir, chipnum)
    #------------------------------------
    # step 3 - find the final fits file at the end of the dir structure
    curr_obs_dir = os.path.join(chip_subdir,'red')
    rabbit_hole=True
    while rabbit_hole:
        obs_subdir_list = os.listdir(curr_obs_dir)
        try:
            if os.path.isfile(os.path.join(curr_obs_dir,obs_subdir_list[0])):
                obs_dir = curr_obs_dir
                rabbit_hole=False
            else:
                curr_obs_dir = os.path.join(curr_obs_dir,obs_subdir_list[0])
        except:
            logger.warning('No files in %s'%curr_obs_dir)
            return None


    #------------------------------------
    # step 4 - CHECK THIS IS THE CORRECT EXPOSURE NUMBER!!!
    obs_fns = []
    if os.path.split(obs_dir)[-1] =='red':
        list_of_imgs = glob.glob(os.path.join(obs_dir,'*immasked*'))
        if len(list_of_imgs) == 0:
            list_of_imgs = glob.glob(os.path.join(obs_dir,'*'))
    else:
        list_of_imgs = glob.glob(os.path.join(obs_dir,'*'))
    for obs_fn in list_of_imgs:
        try:
            fits_expnum = fits.getheader(obs_fn)['EXPNUM']
        except:
            continue

        if obs_fn[-9:] =='sked.fits':
            obs_fn = obs_fn+'[0]'
        fn_ext = os.path.split(obs_fn)[-1]
        if fn_ext[:3]!='DSN':
            obs_fns.append(obs_fn)
    return obs_fns

def rn_stack(stack):
    pass
def chiplims(stack):
    import os
    import astropy.io.fits as fits
    from astropy.wcs import WCS

    root_dir = s.out_dir+'/MYnone'
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

def resample(s,lst,y,chip,cuts,j,logger,stamp_sizex=4200,stamp_sizey=2200):

    img_list = np.genfromtxt(lst,dtype='str',delimiter='\n')
    if len(img_list)==0:
        logger.info('Empty list: %s \n %s'%(lst,img_list))
        return False

    starttime=float(time.time())
    logger.info('Creating weightmaps for individual input exposures in %s, %s band, chip %s'%(s.field,s.band,chip))
    weightlist,resamplist,headerlist,masklist  = [],[],[],[]
    os.chdir(s.temp_dir)
    #try:
    #    pixel_scale = 3600.0*abs(fits.getheader(img_list[0])['CD1_1'])
    #except:
        #pixel_scale = 3600.0*abs(fits.getheader(img_list[0][:-3])['CD1_1'])
    ra_cent,dec_cent = get_chip_vals(s.field,int(chip))
    try:
        if len(img_list)>1:

            for img in img_list:
                imgname = os.path.split(img)[-1]
                imgroot = imgname[:-5]
                if imgroot[-2:]=='fi':
                    imgroot = imgroot[:-3]
                weightlist.append(os.path.join(s.temp_dir,imgroot +'.resamp.weight.fits'))
                resamplist.append(os.path.join(s.temp_dir,imgroot+'.resamp.fits'))
                headerlist.append(os.path.join(s.temp_dir,imgroot+'.resamp.head'))
                masklist.append(os.path.join(s.temp_dir,imgroot+'.resamp.mask.fits'))


            swarp_cmd = [
            'swarp',
            '@%s'%lst,
            '-COMBINE','N',
            '-RESAMPLE','Y',
            #'-IMAGE_SIZE','%s,%s'%(4113,2058),
            #'-CENTER_TYPE','MANUAL',
            #'-CENTER','%f,%f'%(ra_cent,dec_cent),
            #'-PIXELSCALE_TYPE','MANUAL',
            #'-PIXEL_SCALE','0.2363',#'%.03f'%pixel_scale,
            '-BACK_SIZE','256',
            ]

    except TypeError:
        img = str(img_list)
        swarp_cmd = [
        'swarp','%s'%img,
        '-COMBINE','N',
        '-RESAMPLE','Y',
        '-IMAGE_SIZE','%s,%s'%(stamp_sizex,stamp_sizey),
        '-CENTER_TYPE','MANUAL',
        '-CENTER','%f,%f'%(ra_cent,dec_cent),
        '-PIXELSCALE_TYPE','MANUAL',
        '-PIXEL_SCALE','0.2363',#'%.03f'%pixel_scale,
        '-BACK_SIZE','256',
        ]
        imgname = os.path.split(img)[-1]
        imgroot = imgname[:-5]
        if imgroot[-2:]=='fi':
            imgroot = imgroot[:-3]
        weightlist.append(os.path.join(s.temp_dir,imgroot +'.resamp.weight.fits'))
        resamplist.append(os.path.join(s.temp_dir,imgroot+'.resamp.fits'))
        headerlist.append(os.path.join(s.temp_dir,imgroot+'.resamp.head'))
        masklist.append(os.path.join(s.temp_dir,imgroot+'.resamp.mask.fits'))

    logger.info('Resampling with command: %s'%swarp_cmd)
    p = subprocess.Popen(swarp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    endtime=float(time.time())
    logger.info('Finished creating weightmaps for %s, %s band, chip %s; took %.3f seconds'%(s.field,s.band,chip,endtime-starttime))

    weightlist = np.array(weightlist)
    weightlist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.wgt.lst'%(y,s.field,s.band,chip,s.cutstring,j))
    np.savetxt(weightlist_name,weightlist,fmt='%s')
    resamplist = np.array(resamplist)
    resamplist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.resamp.lst'%(y,s.field,s.band,chip,s.cutstring,j))
    np.savetxt(resamplist_name,resamplist,fmt='%s')
    headerlist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.head.lst'%(y,s.field,s.band,chip,s.cutstring,j))
    np.savetxt(headerlist_name,headerlist,fmt='%s')

    try:
        for img in resamplist:
            imgname = os.path.split(img)[-1]
            imgroot = imgname[:-5]

            if imgroot[-2:]=='fi':
                imgroot = imgroot[:-3]
            try:
                imgn = img[:-3]
                h = fits.getheader(imgn)
            except:
                h = fits.getheader(img)
            header_name = os.path.join(s.temp_dir,imgroot+'.head')
            if os.path.isfile(header_name):
                os.remove(header_name)
            h.totextfile(header_name)

    except TypeError:
        img = str(img_list)
        h = fits.getheader(img)
        h.totextfile(os.path.join(s.temp_dir,imgroot+'.head'))
    return (weightlist_name,resamplist_name,headerlist)

def make_cap_stamps(sg,sr,si,sz,chip,sn_name,ra,dec,stamp_sizex=4100,stamp_sizey=2100):
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # resample to the same grid using SWarp
    # start by getting the science frames for each band

    sci_frames = []
    for s in [sg,sr,si,sz]:
        bd = s.band_dir
        # assume we don't have multiple versions of the science frame
        glob_string = os.path.join(bd,'ccd_%s_%s_*clipweighted_sci.fits'%(str(chip),s.band))
        logger.info("Looking for things that look like: '%s'"%glob_string)
        glob_list = glob.glob(glob_string)
        sci_frames.append(glob_list[0])
        logger.info("Found the correct coadd, exists at: '%s'"%glob_list[0])
    pixel_scale = 3600.0*abs(fits.getheader(sci_frames[0])['CD1_1'])
    sci_frame_str = sci_frames[0]+' '+sci_frames[1]+' '+sci_frames[2]+' '+sci_frames[3]

    # set up the directory if it doesn't already exist
    cap_dir = os.path.join(sg.out_dir,'CAP')
    if not os.path.isdir(cap_dir):
        os.mkdir(cap_dir)
    sn_dir = os.path.join(cap_dir,sn_name)
    if not os.path.isdir(sn_dir):
        os.mkdir(sn_dir)
    os.chdir(sn_dir)
    # make a white stamp as a det image
    logger.info('Resampling all bands in a stamp around %s'%sn_name)
    resamp_cmd = ['swarp',
    '-COMBINE','N',
    '-RESAMPLE','Y',
    '-IMAGE_SIZE','%s,%s'%(stamp_sizex,stamp_sizey),
    '-CENTER_TYPE','MANUAL',
    '-CENTER','%f,%f'%(ra,dec),
    '-PIXELSCALE_TYPE','MANUAL',
    '-PIXEL_SCALE','%.03f'%pixel_scale,
    '-BACK_SIZE','512',
    sci_frame_str]
    logger.info("Swarping with the command:\n '%s'"%resamp_cmd)
    starttime=float(time.time())
    p = subprocess.Popen(resamp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    endtime=float(time.time())
    logger.info('Done resampling a stamp around %s, took %.3f seconds'%(sn_name,endtime-starttime))

    # Now resample the  image
    resamp_frames = []
    for s in [sg,sr,si,sz]:
        glob_string = os.path.join(sn_dir,'ccd_%s_%s_*_sci.resamp.fits'%(str(chip),s.band))
        glob_list = glob.glob(glob_string)
        resamp_frames.append(glob_list[0])
        cmd = ['swarp',
        '-IMAGE_SIZE','%s,%s'%(stamp_sizex,stamp_sizey),
        '-CENTER_TYPE','MANUAL',
        '-CENTER','%f,%f'%(ra,dec),
        '-PIXELSCALE_TYPE','MANUAL',
        '-PIXEL_SCALE','%.03f'%pixel_scale,
        '-BACK_SIZE','512',
        '-IMAGEOUT_NAME',glob_list[0],
        glob_list[0]]
        pr = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = pr.communicate()
    resamp_frame_str = resamp_frames[0]+' '+resamp_frames[1]+' '+resamp_frames[2]+' '+resamp_frames[3]

    return

def get_chip_vals(f,chip,vals = 'center'):
    lim_file=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
    chiplims = cpickle.load(lim_file)
    short_field = f[-2:]
    field_lims = chiplims[short_field]
    this_chip_lims = field_lims[chip]
    if vals == 'center':
        ras = [this_chip_lims[i][0] for i in range(4)]
        decs = [this_chip_lims[i][1] for i in range(4)]
        return (np.mean(ras),np.mean(decs))
    elif vals == 'lims':
        return this_chip_lims

def resample_chip_for_cap(sg,sr,si,sz,chip,stamp_sizex=4300,stamp_sizey=2300,npix_off1 = 0,npix_off2 = 0):
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # resample to the same grid using SWarp
    # start by getting the science frames for each band

    sci_frames = []
    naxis1s,naxis2s = [],[]
    for s in [sg,sr,si,sz]:
        bd = s.band_dir
        # assume we don't have multiple versions of the science frame
        glob_string = os.path.join(bd,'ccd_%s_%s_*_clipweighted_sci.fits'%(str(chip),s.band))
        logger.info("Looking for things that look like: '%s'"%glob_string)
        glob_list = glob.glob(glob_string)
        sci_frames.append(glob_list[0])
        logger.info("Found the correct coadd, exists at: '%s'"%glob_list[0])
        naxis1s.append(fits.getheader(glob_list[0])['NAXIS1'])
        naxis2s.append(fits.getheader(glob_list[0])['NAXIS2'])
    ghead = fits.getheader(sci_frames[0])
    pixel_scale = 3600.0*abs(ghead['CD1_1'])
    sci_frame_str = sci_frames[0]+' ' +sci_frames[1]+' '+sci_frames[2]+' '+sci_frames[3]

    # set up the directory if it doesn't already exist
    cap_dir = os.path.join(sg.out_dir,'MY%s'%sg.my,sg.field,'CAP')
    if not os.path.isdir(cap_dir):
        os.mkdir(cap_dir)
    cap_chip_dir = os.path.join(cap_dir,str(chip))
    if not os.path.isdir(cap_chip_dir):
        os.mkdir(cap_chip_dir)
    os.chdir(cap_chip_dir)
    # find the center of the chip
    ra_cent,dec_cent = ghead['CRVAL1'],ghead['CRVAL2']
    smallest1,smallest2 = min(naxis1s),min(naxis2s)
    stamp_sizex,stamp_sizey = smallest1-80-npix_off1,smallest2-80-npix_off2
    # make a riz stamp as a det image
    logger.info('Resampling all bands in MY%s, %s, chip %s'%(sg.my,sg.field,chip))
    logger.info('Hopefully making images of size %s x %s'%(stamp_sizex,stamp_sizey))
    logger.info('Having taken %s off Axis 1 and % off Axis 2' %(npix_off1,npix_off2))
    resamp_cmd = ['swarp',
    '-COMBINE','N',
    '-RESAMPLE','Y',
    '-IMAGE_SIZE','%s,%s'%(stamp_sizex,stamp_sizey),
    '-CENTER_TYPE','MANUAL',
    '-CENTER','%f,%f'%(ra_cent,dec_cent),
    '-PIXELSCALE_TYPE','MANUAL',
    '-PIXEL_SCALE','%.03f'%pixel_scale,
    '-BACK_SIZE','512',
    sci_frame_str]
    logger.info("Swarping with the command:\n '%s'"%resamp_cmd)
    starttime=float(time.time())
    p = subprocess.Popen(resamp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    endtime=float(time.time())
    logger.info('Done resampling MY%s, %s, chip %s, took %.3f seconds'%(sg.my,sg.field,chip,endtime-starttime))
    # Now resample the  image
    resamp_frames = []
    for s in [sg,sr,si,sz]:
        glob_string = os.path.join(cap_chip_dir,'ccd_%s_%s_*_clipweighted*.resamp.fits'%(str(chip),s.band))
        glob_list = glob.glob(glob_string)
        resamp_frames.append(glob_list[0])

    resamp_frame_str = resamp_frames[1]+' '+resamp_frames[2]+' '+resamp_frames[3]
    logger.info('Creating riz with size %s x %s'%(stamp_sizex,stamp_sizey))
    riz_cmd = ['swarp',
    '-IMAGE_SIZE','%s,%s'%(stamp_sizex,stamp_sizey),
    '-CENTER_TYPE','MANUAL',
    '-CENTER','%f,%f'%(ra_cent,dec_cent),
    '-PIXELSCALE_TYPE','MANUAL',
    '-PIXEL_SCALE','%.03f'%pixel_scale,
    '-COMBINE','Y',
    '-RESAMPLE','N',
    '-BACK_SIZE','512',
    '-IMAGEOUT_NAME','%s_%s_%s_riz.fits'%(sg.my,sg.field,chip),
    resamp_frame_str]
    logger.info("Making a detection image for MY%s, %s, chip %s for CAP"%(sg.my,sg.field,chip))
    logger.info("Using this command: %s"%riz_cmd)
    starttime=float(time.time())
    p = subprocess.Popen(riz_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    endtime=float(time.time())
    logger.info('Done making detection image for for MY%s, %s, chip %s, took %.3f seconds'%(sg.my,sg.field,chip,endtime-starttime))
    logger.info('Checking they are the correct size')
    n_off1,n_off2 = check_resamps('%s_%s_%s_riz.fits'%(sg.my,sg.field,chip),resamp_frames)
    return '%s_%s_%s_riz.fits'%(sg.my,sg.field,chip),n_off1,n_off2

def check_resamps(riz_fn,resamp_frames):
    riz_h = fits.getheader(riz_fn)
    n1riz,n2riz = riz_h['NAXIS1'],riz_h['NAXIS2']
    #print ('Length of riz: %s x %s'%(n1riz,n2riz))
    n_off1 = 0
    n_off2 = 0
    for i in range(len(resamp_frames)):
        print (resamp_frames[i])
        h = fits.getheader(resamp_frames[i])
        n1,n2 = h['NAXIS1'],h['NAXIS2']
        n_diff1 = n1riz - n1
        n_diff2 = n2riz - n2
        print (n_diff1,n_diff2)
        if n_diff1>0:
            if n_diff1>n_off1:
                n_off1 = n_diff1
        if n_diff2>0:
            if n_diff2>n_off2:
                n_off2 = n_diff2
    print ('Returning',n_off1*2,n_off2*2)
    return (n_off1*2,n_off2*2)


def get_cuts(f,b):
    cp=configparser.ConfigParser()
    # read the .ini file
    cp.read('/media/data3/wiseman/des/coadding/config/snobs_params.ini')

    if f[-1]!='3':
        cuts = cp['%s_shallow'%b]
    else:
        cuts = cp['%s_deep'%b]
    cuts = dict(cuts)
    try:
        cuts['teff']=float(cuts['teff'])
    except:
        pass
    try:
        cuts['psf']=float(cuts['psf'])
    except:
        pass

    return cuts

def combine_mask_weight(s,chip,j):
    maskweightlist,masklist = [],[]
    resamplist_name = os.path.join(s.list_dir,'%s_%s_%s_%s_%s_%s.resamp.lst'%(s.my,s.field,s.band,chip,s.cutstring,j))
    for f in np.genfromtxt(resamplist_name,dtype='str',delimiter='\n'):
        if f[-3:]=='[0]':
            f = f[:-3]
        #try:
        wn = os.path.join(s.temp_dir,f.replace('fits','weight.fits'))
        mn = os.path.join(s.temp_dir,f.replace('fits','head.mask.fits'))
        masklist.append(mn)
        w = fits.open(wn)
        wd = w[0].data
        m = fits.getdata(mn)
        try:
            maskweight = np.multiply(wd,m)
        except ValueError:
            if wd.shape[0]<m.shape[0] or wd.shape[1]<m.shape[1]:
                new_mask = m[:wd.shape[0],:wd.shape[1]]
            elif wd.shape[0]>m.shape[0] or wd.shape[1]>m.shape[1]:
                new_mask = np.ones_like(wd)
                new_mask[:m.shape[0],:m.shape[1]]=m
            maskweight = np.multiply(wd,new_mask)
        w[0].data = maskweight
        #print ('Weight size',w[0].header['NAXIS1'],w[0].header['NAXIS2'])
        w.writeto(os.path.join(s.temp_dir,f.replace('fits','maskweight.fits')),overwrite=True)
        maskweightlist.append(os.path.join(s.temp_dir,f.replace('fits','maskweight.fits')))

            #except:
                #wn = os.path.join(s.temp_dir,f.replace('fits','weight.fits'))
                #w = fits.open(wn)
                #w.writeto(os.path.join(s.temp_dir,f.replace('fits','maskweight.fits')),overwrite=True)
                #maskweightlist.append(os.path.join(s.temp_dir,f.replace('fits','maskweight.fits')))

    np.savetxt(resamplist_name.replace('resamp','maskweight'),np.array(maskweightlist),fmt='%s')
    for mn in masklist:
        os.remove(mn)
    return maskweightlist
