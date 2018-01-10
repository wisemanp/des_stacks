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

def sex_for_psfex(stack,chip,cuts=None):
    '''Runs SExtractor on a certain stacked frame, to send to PSFex'''
    band_dir = os.path.join(stack.out_dir,'MY%s'%stack.my,stack.field,stack.band)
    sexcat = os.path.join(band_dir,chip,'%s_%s_temp.sexcat'%(chip,stack.band))
    try:
        zp_cut,psf_cut = cuts['zp'],cuts['psf']
    except:
        pass
    if stack.final == False:
        img = band_dir+'/ccd_%s_%s_%.3f_%s_temp.fits'%(chip,stack.band,zp_cut,psf_cut)
    else:
        img = os.path.join(band_dir,'ccd_%s_%s_%.3f_%s.fits'%(chip,stack.band,zp_cut,psf_cut))
    os.chdir(os.path.join(band_dir,chip,'psf'))
    #run SExtractor
    sex_cmd = ['sex','-CATALOG_NAME',sexcat,'-CATALOG_TYPE','FITS_LDAC',img]
    start = float(time.time())
    try:
        stack.logger.info("Running SExtractor on {0}".format(img))
        sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = sex_process.communicate()
        end = float(time.time())
        stack.logger.info("Successfully ran SExtractor on {0}".format(img))
        stack.logger.info("Took %10.2f seconds" %(end - start))
        stack.logger.info("Saved at {0}".format(sexcat))
        return sexcat
    except (OSError, IOError):
        stack.logger.info("SExtractor failed...")
        return None

def psfex(stack,chip,retval='FWHM',cuts=None):
    '''Runs PSFex on a SExtracted catalog to get the PSF Model'''

    band_dir = os.path.join(stack.out_dir, 'MY%s' %stack.my, stack.field, stack.band)
    init_cat = os.path.join(band_dir,chip,'%s_%s_temp.sexcat'%(chip,stack.band))
    stack.logger.info("Getting the PSF from the stack")
    stack.logger.info("Running PSFex on %s" %init_cat)
    start =float(time.time())
    psfex_cmd = ['psfex', init_cat]
    os.chdir(os.path.join(band_dir,chip))
    psf_process = subprocess.Popen(psfex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,errs = psf_process.communicate()
    end = float(time.time())
    stack.logger.info("Successfully made the PSF model")
    stack.logger.info("Took %10.2f secconds" % (end-start))
    copyfile(os.path.join(band_dir,chip,'%s_%s_temp.psf'%(chip,stack.band)),os.path.join(band_dir,chip,'ana','default.psf'))
    if retval == 'FWHM':

        psf_out = os.path.join(band_dir,chip,'ana','default.psf')
        h = fits.getheader(psf_out,ext=1)
        fwhm = h['PSF_FWHM']*0.27 #convert from pixels to arcsec using DES chip pixel scale of 0.27 pix/arcsec
        return fwhm

def sex_for_cat(stack,chip,cuts = None):
    '''Runs SExtractor on a certain stacked frame, given a PSF model from PSFex'''
    band_dir = os.path.join(stack.out_dir, 'MY%s' %stack.my, stack.field, stack.band)
    ana_dir = os.path.join(band_dir,chip,'ana')
    os.chdir(ana_dir)
    if not cuts:
        sexcat = os.path.join(ana_dir,'MY%s_%s_%s_%s.sexcat' %(stack.my,stack.field,stack.band,chip))
    else:
        sexcat = os.path.join(ana_dir,'MY%s_%s_%s_%s_%s_%s.sexcat' %(stack.my,stack.field,stack.band,chip,cuts['zp'],cuts['psf']))

    try:
        zp_cut,psf_cut = cuts['zp'],cuts['psf']
    except:
        pass
    if stack.final == False:
        img = band_dir+'/ccd_%s_%s_%.3f_%s_temp.fits'%(chip,stack.band,zp_cut,psf_cut)
    else:
        img = os.path.join(band_dir,'ccd_%s_%s_%.3f_%s.fits'%(chip,stack.band,zp_cut,psf_cut))
    stack.logger.info("Starting source extraction using the modelled PSF")
    start = float(time.time())
    sex_cmd = ['sex','-CATALOG_NAME',sexcat,img]
    try:
        stack.logger.info("Running SExtractor on {0}".format(img))
        sex_process = subprocess.Popen(sex_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,errs = sex_process.communicate()
        stack.logger.info("Successfully ran SExtractor on {0}".format(img))
        end = float(time.time())
        stack.logger.info("Took %10.2f seconds" % (end-start))
        stack.logger.info("Saved at {0}".format(sexcat))
        return sexcat
    except (OSError, IOError):
        stack.logger.info("SExtractor failed...")
        return None
