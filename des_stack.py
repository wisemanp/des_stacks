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
import subprocess
import logging
from shutil import copyfile

from des_stacks.utils.stack_tools import make_good_frame_list, make_swarp_cmd, get_dessn_obs, get_des_obs_year
from des_stacks.utils.sex_tools import sex_for_psfex, psfex, sex_for_cat
from des_stacks.analysis.astro import astrometry,init_phot

class Stack():
    def __init__(self, field, band, my, chips ,working_dir):
        self.field = field
        self.band = band
        self.my =my
        self.chips=chips
        self.coadding_dir =working_dir
        self._define_paths()
        self._init_log()
        self._get_info()
        self._get_configs()
        self.logger.info("Doing work in: %s as a root directory" % self.coadding_dir)

    ###############################################
    def _define_paths(self):
        '''Sets the base paths in which to do all of the work'''
        if self.coadding_dir == 'coadding':
            self.coadding_dir = '/media/data1/wiseman/des/coadding/'
        elif self.coadding_dir == 'current':
            self.coadding_dir = os.curdir
        self.config_dir = os.path.join(self.coadding_dir,'config')
        self.log_dir = os.path.join(self.coadding_dir,'log')
        self.cat_dir = os.path.join(self.coadding_dir,'catalogs')
        self.out_dir = os.path.join(self.coadding_dir,'stacks')
        self.temp_dir = os.path.join(self.coadding_dir,'temp')
        self.db_dir = os.path.join(self.coadding_dir,'db')
        self.list_dir = os.path.join(self.coadding_dir,'good_img_lists')
        self.band_dir = os.path.join(self.out_dir,'MY%s'%self.my,self.field,self.band)
    ########################################################
    def _init_log(self):
        '''Sets up the logger'''
        logger = logging.getLogger('des_stack.py')
        logger.handlers =[]
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        fh = logging.FileHandler(os.path.join(self.log_dir,'stack_%s%s%s%s.log'%(self.field,self.band,self.my,self.chips)))
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.addHandler(fh)
        self.logger = logger
    ###########################################################
    def _get_configs(self,dessn_ini = 'snobs_params.ini'):

        '''Get config files from the config directory
        parameters:
        dessn_ini: (str, optional, default = "snobs_params.ini"). name of
            the .ini file for the DES SN data, including paths for each
            year, plus year lims
        # sex_files: (string, or list of strings,optional, default =
            ["default.sex","params.sex","conv.sex"]). List of the default
            files for SExtractor'''
        if os.path.split(dessn_ini)[0]=='':
            #if only the filename for dessn_ini is given, assume it's in /coadding/config
            ini_fn = os.path.join(self.config_dir,dessn_ini)
        else:
            #assume I'm being given a totally new path
            ini_fn = dessn_ini
        cp=configparser.ConfigParser()
        # read the .ini file
        cp.read(ini_fn)
        # Make a list of years
        years = ['Y1','Y2','Y3','Y4']
        self.year_lims = {}
        self.data_dirs = {}
        for y in years:
            year_mjd_lim =cp.get('year_mjd_lims',y)
            year_night_lim = cp.get('year_night_lims',y)
            self.year_lims[y]={'mjd':year_mjd_lim,'night':year_night_lim}
            self.data_dirs[y]=cp.get('data_dirs',y)
        self.logger.info('Successfully pulled configuration from %s' %self.config_dir)
    ############################################################################
    def _get_info(self):
        '''Gets the current info file (originally from the DESDB)'''
        info_tab = Table.read(os.path.join(self.config_dir,'snobsinfo.fits'))
        self.info_df = info_tab.to_pandas()
    def do_my_stack(self, cuts={'zp':-0.15,'psf':2.5},final=True):
        '''Does a stack defined by the parameters from the Stack object it is passed.
        keyword arguments:
        zp_cut (float): the zeropoint cut to be used for the stack (default = -0.15)
        psf_cut (float): the seeing cut to be used for the stack (default = 2.5)
        final (Bool): whether the stack is final or intermediate (default = True)
        returns:
        none
        '''
        self.final=final
        zp_cut = cuts['zp']
        psf_cut = cuts['psf']
        self.zp_cut = zp_cut
        self.psf_cut = psf_cut
        field = self.field
        band = self.band
        my = self.my
        self.logger.info('******************************************************')
        self.logger.info('Initiating stack on {0} in {1} band'.format(field,band))
        self.logger.info('******************************************************')
        #does list of good frames exist?
        if not os.path.isfile(os.path.join(self.list_dir,'good_exps_%s_%s_%s_%s.fits'%(field,band,zp_cut,psf_cut))):
            #get the list of good frames
            self.logger.info('No good frame list yet, making a new one with ZP < %s and PSF < %s cut' %(zp_cut,psf_cut))
            self.good_frame = make_good_frame_list(self,field,band,zp_cut,psf_cut)

        else:
            good_fn = os.path.join(self.list_dir,'good_exps_%s_%s_%s_%s.fits'%(field,band,zp_cut,psf_cut))
            self.logger.info('Reading in list of good frames from {0}'.format(good_fn))
            good_table = Table.read(good_fn)
            self.good_frame = good_table.to_pandas()

        # how many chips are we stacking?
        chips = self.chips
        if chips == 'All':
            self.chips = self.info_df.CCDNUM.sort_values().unique()
        # get swarp commands
        log = open(os.path.join(self.log_dir,'swarp_%s_%s_%s_%s.log' %(field, band, my, self.chips)),'a')
        log.flush()
        for y in my:
            # catch silly string issue
            if y in 'none':
                y = 'none'
            self.logger.info('Stacking {0} in {1} band, skipping year {2}'.format(field,band,y))
            for chip in self.chips:

                self.logger.info('Stacking CCD {0}'.format(chip))
                cmd = make_swarp_cmd(self,y,field,chip,band,self.logger,zp_cut,psf_cut,final)
                if cmd == False:
                    self.logger.info("Already stacked this chip with these cuts, going straight to astrometry")

                else:

                    self.logger.info('Executing command: {0}'.format(cmd))
                    os.chdir(self.temp_dir)
                    try:
                        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        outs,errs = p.communicate()

                    except (OSError, IOError):
                        self.logger.warn("Swarp failed.", exc_info=1)
                    self.logger.info('Finish stacking chip {0}'.format(chip))

            self.logger.info('Finished stacking chips {0} for MY {1}'.format(self.chips,y))
            if y == 'none':
                break
        log.close()
        self.logger.info('Stacking of %s %s %s %s complete!' %(field, band, my, self.chips))
        self.logger.info('******************************************************')
    def run_stack_sex(self,cuts={'zp':-0.15,'psf':2.5},final=None):
        self.logger.info("******** Preparing to extract and match soures *******")
        if final ==None:
            #check to make sure we haven't just forgotten to say this is a temporary run of Sextractor
            try:
                final = self.final
            except AttributeError:
                final = True
        self.final = final
        zp_cut = cuts['zp']
        psf_cut = cuts['psf']
        self.zp_cut = zp_cut
        self.psf_cut = psf_cut
        qual_df = pd.DataFrame()
        self.sexcats=[]
        for chip in self.chips:
            # create file structure and copy defaults accross
            chip_dir = os.path.join(self.out_dir,'MY%s'%self.my,self.field,self.band,chip)
            if not os.path.isdir(chip_dir):
                os.mkdir(chip_dir)
            psf_dir = os.path.join(chip_dir,'psf')
            if not os.path.isdir(psf_dir):
                os.mkdir(psf_dir)
            copyfile(os.path.join(self.config_dir,'psf','default.sex'),os.path.join(psf_dir,'default.sex'))
            copyfile(os.path.join(self.config_dir,'psf','default.param'),os.path.join(psf_dir,'default.param'))
            copyfile(os.path.join(self.config_dir,'psf','default.conv'),os.path.join(psf_dir,'default.conv'))

            # do sex for psfex
            sex_for_psfex(self,chip,cuts)

            ana_dir = os.path.join(chip_dir,'ana')
            if not os.path.isdir(ana_dir):
                os.mkdir(ana_dir)
            self.chip_dir = chip_dir
            self.ana_dir = ana_dir
            copyfile(os.path.join(self.config_dir,'default.sex'),os.path.join(ana_dir,'default.sex'))
            copyfile(os.path.join(self.config_dir,'default.param'),os.path.join(ana_dir,'default.param'))
            copyfile(os.path.join(self.config_dir,'default.conv'),os.path.join(ana_dir,'default.conv'))

            # do psfex on sex, and get the fwhm from there
            model_fwhm = psfex(self,chip,retval='FWHM',cuts=cuts)
            # do sex on psf
            os.chdir(ana_dir)
            # Do SExtractor on the complete stacks
            sexcat = sex_for_cat(self,chip,cuts)
            self.sexcats.append(sexcat)
            # Compare new catalog to old one, get the ZP and FWHM out
            zp,sex_fwhm = astrometry(self,chip,sexcat)
            zp_psf,psf_fwhm = astrometry(self,chip,sexcat,phot_type = 'PSF')
            qual = open(os.path.join(ana_dir,'%s_%s_ana.qual'%(zp_cut,psf_cut)),'w')
            print ('# Quality parameters for %s %s %s %s' %(self.my,self.field,self.band,chip),file =qual)
            print ('# Parameters:',file=qual)
            print ('# Zeropoint from sextractor',file=qual)
            print ('# Zeropoint from PSF matches', file = qual)
            print ('# FWHM from PSFex',file = qual)
            print ('# FWHM from SExtractor using PSF model',file = qual)
            print ('%s %s %s %s'%(zp,zp_psf,model_fwhm,sex_fwhm),file=qual)
            qual.close()
            self.logger.info("Written quality factors to %s_%s_ana.qual" %(zp_cut,psf_cut))
            qual_dict = {'zp':zp,'fwhm_psfex':model_fwhm,'fwhm_sex':sex_fwhm}
            qual_df = qual_df.append(pd.DataFrame([qual_dict],index = [chip]))
            self.logger.info("And the DataFrame %s" %qual_df)
        self.logger.info("********** Done measuring quality of stack! **********")
        self.logger.info("******************************************************")
        return qual_df

    def init_phot(self):
        limmags = {}
        for counter,chip in enumerate(self.chips):
            sexcat = self.sexcats[counter]
            cat = Table.read(sexcat).to_pandas()
            limmags[chip]=init_phot(self,chip,cat)
        return limmags
