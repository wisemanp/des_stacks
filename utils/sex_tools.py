# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from astropy.io import fits
import datetime
import os
import logging
from shutil import copyfile
import time
import subprocess
import _pickle as cpickle
import easyaccess as ea
import glob
from astropy.table import Table

def sex_for_psfex(s,chip,cuts=None):
    '''Runs SExtractor on a certain stacked frame, to send to PSFex'''
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
    logger.info('Starting sex_for_psfex')
    band_dir = os.path.join(s.out_dir,'MY%s'%s.my,s.field,s.band)
    sexcat = os.path.join(band_dir,chip,'%s_%s_temp.sexcat'%(chip,s.band))
    try:
        zp_cut,psf_cut,t_cut = cuts['zp'],cuts['psf'],cuts['teff']
    except:
        cuts = None
    if s.final == False:
        if not cuts:
            img = band_dir+'/ccd_%s_%s_temp.fits'%(chip,s.band)
        else:
            img = band_dir+'/ccd_%s_%s_%s_temp.fits'%(chip,s.band,s.cutstring)
    else:

        if not cuts:
            img = band_dir+'/ccd_%s.fits'%chip
        else:
            img = os.path.join(band_dir,'ccd_%s_%s_%s_clipweighted_sci.fits'%(chip,s.band,s.cutstring))
    os.chdir(os.path.join(band_dir,chip,'psf'))
    #run SExtractor
    logger.info('Got as far as starting SExtractor')
    sex_cmd = ['sex','-CATALOG_NAME',sexcat,'-CATALOG_TYPE','FITS_LDAC',img]
    start = float(time.time())
    try:
        logger.info("Running SExtractor on {0}".format(img))
        sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = sex_process.communicate()
        end = float(time.time())
        logger.info("Successfully ran SExtractor on {0}".format(img))
        logger.info("Took %.2f seconds" %(end - start))
        logger.info("Saved at {0}".format(sexcat))
        return sexcat
    except (OSError, IOError):
        logger.info("SExtractor failed...")
        return None

def psfex(s,chip,retval='FWHM',cuts=None):
    '''Runs PSFex on a SExtracted catalog to get the PSF Model'''
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
    band_dir = os.path.join(s.out_dir, 'MY%s' %s.my, s.field, s.band)
    init_cat = os.path.join(band_dir,chip,'%s_%s_temp.sexcat'%(chip,s.band))
    logger.info("Getting the PSF from the stack")
    logger.info("Running PSFex on %s" %init_cat)
    start =float(time.time())
    psfex_cmd = ['psfex', init_cat]
    os.chdir(os.path.join(band_dir,chip))
    psf_process = subprocess.Popen(psfex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,errs = psf_process.communicate()
    end = float(time.time())
    logger.info("Successfully made the PSF model")
    logger.info("Took %.2f secconds" % (end-start))
    copyfile(os.path.join(band_dir,chip,'%s_%s_temp.psf'%(chip,s.band)),os.path.join(band_dir,chip,'ana','default.psf'))
    if retval == 'FWHM':

        psf_out = os.path.join(band_dir,chip,'ana','default.psf')
        h = fits.getheader(psf_out,ext=1)
        fwhm = h['PSF_FWHM']*0.27 #convert from pixels to arcsec using DES chip pixel scale of 0.27 pix/arcsec
        return fwhm

def sex_for_cat(s,chip,cuts = None):
    '''Runs SExtractor on a certain stacked frame, given a PSF model from PSFex'''
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

    band_dir = os.path.join(s.out_dir, 'MY%s' %s.my, s.field, s.band)
    ana_dir = os.path.join(band_dir,chip,'ana')
    os.chdir(ana_dir)

    if not cuts:
        sexcat = os.path.join(ana_dir,'MY%s_%s_%s_%s.sexcat' %(s.my,s.field,s.band,chip))
    else:
        sexcat = os.path.join(ana_dir,'MY%s_%s_%s_%s_%s_clipweighted_sci.sexcat' %(s.my,s.field,s.band,chip,s.cutstring))

    try:
        zp_cut,psf_cut = cuts['zp'],cuts['psf']
    except:
        cuts = None
    if s.final == False:
        if not cuts:
            img = band_dir+'/ccd_%s_%s_temp.fits'%(chip,s.band)
        else:
            img = band_dir+'/ccd_%s_%s_%s_temp.fits'%(chip,s.band,s.cutstring)
    else:

        if not cuts:
            img = band_dir+'/ccd_%s.fits'%chip
        else:
            img = os.path.join(band_dir,'ccd_%s_%s_%s_clipweighted_sci.fits'%(chip,s.band,s.cutstring))
    logger.info("Starting source extraction using the modelled PSF")
    start = float(time.time())
    sex_cmd = ['sex','-CATALOG_NAME',sexcat,img]
    try:
        logger.info("Running SExtractor on {0}".format(img))
        sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = sex_process.communicate()
        logger.info("Successfully ran SExtractor on {0}".format(img))
        end = float(time.time())
        logger.info("Took %.2f seconds" % (end-start))
        logger.info("Saved at {0}".format(sexcat))
        return sexcat
    except (OSError, IOError):
        logger.info("SExtractor failed...")
        return None

def get_sn_dat(sn):
    f=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
    chiplims = cpickle.load(f)

    sncand = Table.read('/media/data3/wiseman/des/coadding/catalogs/sn_cand.fits').to_pandas()
    gap = ' '
    ngaps = (11-len(sn))*gap
    dat = sncand[sncand['TRANSIENT_NAME']==sn+ngaps]

    ra,dec =dat[['RA','DEC']].iloc[0].values
    y = dat['SEASON'].values[0]

    #################
    obj_field = dat['FIELD'].values[0].strip(' ')
    if len(obj_field)>2:
        obj_field = obj_field.split(',')[0]
    the_field = chiplims[obj_field]
    try:
        for ccd in the_field.keys():
            if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
                if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                    return (ra,dec,'SN-%s'%obj_field,y,ccd)

        if len(obj_field)>2:
            obj_field = obj_field.split(',')[1]
        the_field = chiplims[obj_field]
        for ccd in the_field.keys():
            if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
                if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                    return (ra,dec,'SN-%s'%obj_field,y,ccd)
    except KeyError:
        return False
def cap_sex_sn(sg,sr,si,sz,chip,sn_name,leave_if_done = False):
    '''Runs SExtractor in dual image mode to get common aperture photometry'''
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    logger.setLevel(logging.INFO)
    ch.setLevel(logging.INFO)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # get to the right directory
    sn_dir = os.path.join(sg.out_dir,'CAP',sn_name)
    os.chdir(sn_dir)
    white_name = sn_name+'_white_stamp.fits'
    # get the right config files in the directory
    for ext in ['sex','param','conv','nnw']:
        copyfile(os.path.join(sg.config_dir,'cap','default.%s'%ext),os.path.join(sn_dir,'default.%s'%ext))
    sexcats ={}
    for s in [sg,sr,si,sz]:
        sexcat = os.path.join(sn_dir,'%s_%s_cap_sci.sexcat'%(sn_name,s.band))
        if os.path.isfile(sexcat) and leave_if_done:
            logger.info("Already done the photometry in the %s band!"%s.band)
        else:
            glob_string = os.path.join(sn_dir,'ccd_%s_%s_*_sci.resamp.fits'%(str(chip),s.band))
            resamp_name = glob.glob(glob_string)[0]
            sex_cmd = ['sex','-CATALOG_NAME',sexcat,'%s,%s'%(white_name,resamp_name)]
            logger.info('Running SExtractor in dual image mode in order to get common aperture photometry in the %s band'%s.band)
            sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,errs = sex_process.communicate()
            logger.info('Dual image SExtractor complete in the %s band: you now have common aperture photometry on chip %s!'%(s.band,chip))
        sexcats[s.band]=sexcat
    return sexcats

def cap_sex_chip(sg,sr,si,sz,chip):
    '''Runs SExtractor in dual image mode to get common aperture
    photometry on a chip that has already had all bands resampled to the same pixels'''
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # get to the right directory
    cap_chip_dir = os.path.join(sg.out_dir,'MY%s'%sg.my,sg.field,'CAP',str(chip))
    if not os.path.isdir(cap_chip_dir):
        os.mkdir(cap_chip_dir)
    os.chdir(cap_chip_dir)
    # get the right config files in the directory
    for ext in ['sex','param','conv','nnw']:
        copyfile(os.path.join(sg.config_dir,'cap','default.%s'%ext),os.path.join(cap_chip_dir,'default.%s'%ext))
    sexcats ={}
    for s in [sg,sr,si,sz]:
        try:
            quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        except:
            s.run_stack_sex(cuts=s.cuts,final=True)
            quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        zp = float(quals[0])
        riz_name = '%s_%s_%s_riz.fits'%(s.my,s.field,chip)
        sexcat = os.path.join(cap_chip_dir,'%s_%s_%s_%s_cap_sci.sexcat'%(s.my,s.field,chip,s.band))
        redo = False
        if os.path.isfile(sexcat):
            logger.info("Already done the photometry in the %s band!"%s.band)
            test_cat = Table.read(sexcat).to_pandas()
            if 'FLUX_AUTO' in test_cat.columns:
                pass
            else:
                redo = True
                logger.info("There was a .sexcat file, but it didn't have the right parameters, so running SExtractor again")

        if os.path.isfile(sexcat) and redo ==False:
            pass
        else:
            glob_string = os.path.join(cap_chip_dir,'ccd_%s_%s_*_clipweighted*.resamp.fits'%(str(chip),s.band))
            resamp_name = glob.glob(glob_string)[0]
            check_name = os.path.join(cap_chip_dir,'%s_%s_%s_%s_check_aper.fits'%(s.my,s.field,chip,s.band))
            sex_cmd = [
            'sex',
            '-MAG_ZEROPOINT',str(zp),
            '-CATALOG_NAME',sexcat,
            '-CHECKIMAGE_TYPE','APERTURES',
            '-CHECKIMAGE_NAME',check_name,
            '%s,%s'%(riz_name,resamp_name)
            ]
            logger.info('Running SExtractor in dual image mode in order to get common aperture photometry in the %s band'%s.band)
            logger.debug(sex_cmd)
            sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,errs = sex_process.communicate()
            logger.info('Dual image SExtractor complete in the %s band: you now have common aperture photometry on chip %s!'%(s.band,chip))
        sexcats[s.band]=sexcat
        logger.info('Returning sexcats for %s,%s,%s,%s'%(sg.my,sg.field,chip,sg.band))
    return sexcats
