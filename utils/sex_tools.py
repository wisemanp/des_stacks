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
            img = os.path.join(band_dir,'ccd_%s_%s_%s_sci.fits'%(chip,s.band,s.cutstring))
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
        sexcat = os.path.join(ana_dir,'MY%s_%s_%s_%s_%s_sci.sexcat' %(s.my,s.field,s.band,chip,s.cutstring))

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
            img = os.path.join(band_dir,'ccd_%s_%s_%s_sci.fits'%(chip,s.band,s.cutstring))
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

def cap_sex(sg,sr,si,sz,y,f,chip,white_name):
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
    chip_cap_dir = os.path.join(sg.out_dir,'MY%s'%y,f,'cap',str(chip))
    os.chdir(chip_cap_dir)

    # get the right config files in the directory
    for ext in ['sex','param','conv','nnw']:
        copyfile(os.path.join(sg.config_dir,'cap','default.%s'%ext),os.path.join(chip_cap_dir,'default.%s'%ext))
    sexcats ={}
    for s in [sg,sr,si,sz]:
        sexcat = os.path.join(chip_cap_dir,'MY%s_%s_ccd_%s_%s_cap_sci.sexcat'%(y,f,chip,s.band))
        sexcats[s.band]=sexcat
        glob_string = os.path.join(cap_chip_dir,'ccd_%s_%s_*_sci.resamp.fits'%(str(chip),s.band))
        resamp_name = glob.glob(glob_string)[0]
        sex_cmd = ['sex','-CATALOG_NAME',sexcat,'%s,%s'%(white_name,resamp_name)]
        logger.info('Running SExtractor in dual image mode in order to get common aperture photometry in the %s band'%s.band)
        sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = sex_process.communicate()
        logger.info('Dual image SExtractor complete in the %s band: you now have common aperture photometry on chip %s!'%(s.band,chip))
    return sexcats
