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
import multiprocessing
import logging
from shutil import copyfile
import time
from des_stacks.utils.stack_tools import make_good_frame_list, make_swarp_cmd, get_dessn_obs, get_des_obs_year,make_weightmap
from des_stacks.utils.sex_tools import sex_for_psfex, psfex, sex_for_cat
import des_stacks.utils.multi_stack as multi_stack
from des_stacks.analysis.astro import calib,init_phot

class Stack():
    def __init__(self, field, band, my, chips ,working_dir,cuts={'none':None}):
        self.field = field
        self.band = band
        self.my =my
        self.chips=chips
        self.coadding_dir =working_dir

        if cuts['none']:
            if self.band in ['g','r']:
                cuts ={'teff':0.15, 'zp':None,'psf':None}

            elif self.band in ['i','z']:
                cuts ={'teff':0.25, 'zp':None,'psf':None}

        if cuts['teff']:
            self.cutstring = '%s'%cuts['teff']
            
        elif cuts['zp'] and cuts['psf'] and not cuts['teff']:
            self.cutstring = '%.3f_%s'%(cuts['zp'],cuts['psf'])
        self._define_paths()
        self._init_log()
        self._get_info()
        self._get_configs()
        self.logger.info("Doing work in: %s as a root directory" % self.coadding_dir)

    ###############################################
    def _define_paths(self):
        '''Sets the base paths in which to do all of the work'''
        if self.coadding_dir == 'coadding':
            self.coadding_dir = '/media/data3/wiseman/des/coadding/'
        elif self.coadding_dir == 'current':
            self.coadding_dir = os.curdir
        self.config_dir = os.path.join(self.coadding_dir,'config')
        if not os.path.isdir(self.config_dir):
            os.mkdir(self.config_dir)
        self.log_dir = os.path.join(self.coadding_dir,'log')
        if not os.path.isdir(self.log_dir):
            os.mkdir(self.log_dir)
        self.cat_dir = os.path.join(self.coadding_dir,'catalogs')
        if not os.path.isdir(self.cat_dir):
            os.mkdir(self.cat_dir)
        self.out_dir = os.path.join(self.coadding_dir,'stacks')
        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)
        self.temp_dir = os.path.join(self.coadding_dir,'temp')
        if not os.path.isdir(self.temp_dir):
            os.mkdir(self.temp_dir)
        self.db_dir = os.path.join(self.coadding_dir,'db')
        if not os.path.isdir(self.db_dir):
            os.mkdir(self.db_dir)
        self.list_dir = os.path.join(self.coadding_dir,'good_img_lists')
        if not os.path.isdir(self.list_dir):
            os.mkdir(self.list_dir)
        self.weight_dir = os.path.join(self.coadding_dir,'weights')
        if not os.path.isdir(self.weight_dir):
            os.mkdir(self.weight_dir)
        self.band_dir = os.path.join(self.out_dir,'MY%s'%self.my,self.field,self.band)
        if not os.path.isdir(self.band_dir):
            if os.path.isdir(os.path.join(self.out_dir,'MY%s'%self.my,self.field)):
                os.mkdir(self.band_dir)
            else:
                if os.path.isdir(os.path.join(self.out_dir,'MY%s'%self.my)):
                    os.mkdir(os.path.join(self.out_dir,'MY%s'%self.my,self.field))
                else:
                    os.mkdir(os.path.join(self.out_dir,'MY%s'%self.my))
                    os.mkdir(os.path.join(self.out_dir,'MY%s'%self.my,self.field))
                    os.mkdir(self.band_dir)
    ########################################################
    def _init_log(self):
        '''Sets up the logger'''
        logger = logging.getLogger('des_stack.py')
        logger.handlers =[]
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        fh = logging.FileHandler(os.path.join(self.log_dir,'stack_%s%s%s.log'%(self.field,self.band,self.my)))
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
        self.qual_tab = pd.read_csv(os.path.join(self.config_dir,'SN_TABLE_EVAL.txt'),sep=' ')
    def do_my_stack(self, cuts={'teff':0.2, 'zp':None,'psf':None},final=True):
        '''Does a stack defined by the parameters from the Stack object it is passed.
        keyword arguments:
        zp_cut (float): the zeropoint cut to be used for the stack (default = -0.15)
        psf_cut (float): the seeing cut to be used for the stack (default = 2.5)
        final (Bool): whether the stack is final or intermediate (default = True)
        returns:
        none
        '''

        self.final=final
        self.logger.info("Cuts: %s"%cuts)
        self.zp_cut,self.psf_cut,self.t_cut = cuts['zp'], cuts['psf'],cuts['teff']
        if cuts['zp'] and cuts['psf'] and not cuts['teff']:
            self.cutstring = '%.3f_%s'%(cuts['zp'],cuts['psf'])
        elif cuts['teff']:
            self.cutstring = '%s'%cuts['teff']
        field = self.field
        band = self.band
        my = self.my
        self.logger.info('******************************************************')
        self.logger.info('Initiating stack on {0} in {1} band'.format(field,band))
        self.logger.info('******************************************************')
        #does list of good frames exist?
        if not os.path.isfile(os.path.join(self.list_dir,'good_exps_%s_%s_%s.fits'%(field,band,self.cutstring))):
            #get the list of good frames
            self.logger.info('No good frame list yet, making a new one with T_eff> %s, ZP < %s. and PSF < %s' %(self.t_cut,self.zp_cut,self.psf_cut))
            self.good_frame = make_good_frame_list(self,field,band,cuts)

        else:
            good_fn = os.path.join(self.list_dir,'good_exps_%s_%s_%s.fits'%(field,band,self.cutstring))
            self.logger.info('Reading in list of good frames from {0}'.format(good_fn))
            good_table = Table.read(good_fn)
            self.good_frame = good_table.to_pandas()

        # how many chips are we stacking?
        chips = self.chips
        if chips == 'All':
            self.chips =[]
            chips = self.info_df.CCDNUM.sort_values().unique()
            for c in chips:
                self.chips.append(str(c))
        # get swarp commands
        log = open(os.path.join(self.log_dir,'swarp_%s_%s_%s.log' %(field, band, my)),'a')
        log.flush()
        for y in my:
            # catch silly string issue
            if y in 'none':
                y = 'none'
            self.logger.info('Stacking {0} in {1} band, skipping year {2}'.format(field,band,y))
            multi_stack.multitask(self,y,field,band,cuts,final)
            self.logger.info('Finished stacking chips {0} for MY {1}'.format(self.chips,y))
            if y == 'none':
                break
        log.close()
        self.logger.info('Stacking of %s %s %s %s complete!' %(field, band, my, self.chips))
        self.logger.info('******************************************************')
    def run_stack_sex(self,cuts={'teff':0.2, 'zp':None,'psf':None},final=None):
        self.cuts=cuts
        self.logger.info("******** Preparing to extract and match soures *******")
        if not self.cutstring:
            if cuts['zp'] and cuts['psf'] and not cuts['teff']:
                self.cutstring = '%.3f_%s'%(cuts['zp'],cuts['psf'])
            elif cuts['teff']:
                self.cutstring = '%s'%cuts['teff']
        if final ==None:
            #check to make sure we haven't just forgotten to say this is a temporary run of Sextractor
            try:
                final = self.final
            except AttributeError:
                final = True
        self.final = final
        try:
            zp_cut = cuts['zp']
            psf_cut = cuts['psf']
            t_cut = cuts['teff']
        except TypeError:
            zp_cut,psf_cut,t_cut=None,None,None
        self.zp_cut,self.psf_cut,self.t_cut = zp_cut,psf_cut,t_cut
        qual_df = pd.DataFrame()
        self.sexcats=[]
        for chip in self.chips:
            # create file structure and copy defaults accross
            chip_dir = os.path.join(self.out_dir,'MY%s'%self.my,self.field,self.band,str(chip))
            if not os.path.isdir(chip_dir):
                os.mkdir(chip_dir)
            psf_dir = os.path.join(chip_dir,'psf')
            if not os.path.isdir(psf_dir):
                os.mkdir(psf_dir)
            copyfile(os.path.join(self.config_dir,'psf','default.sex'),os.path.join(psf_dir,'default.sex'))
            copyfile(os.path.join(self.config_dir,'psf','default.param'),os.path.join(psf_dir,'default.param'))
            copyfile(os.path.join(self.config_dir,'psf','default.conv'),os.path.join(psf_dir,'default.conv'))

            ana_dir = os.path.join(chip_dir,'ana')
            if not os.path.isdir(ana_dir):
                os.mkdir(ana_dir)
            self.chip_dir = chip_dir
            self.ana_dir = ana_dir
            copyfile(os.path.join(self.config_dir,'default.sex'),os.path.join(ana_dir,'default.sex'))
            copyfile(os.path.join(self.config_dir,'default.param'),os.path.join(ana_dir,'default.param'))
            copyfile(os.path.join(self.config_dir,'default.conv'),os.path.join(ana_dir,'default.conv'))
            copyfile(os.path.join(self.config_dir,'default.nnw'),os.path.join(ana_dir,'default.nnw'))

        self.sexcats = multi_stack.multitask(self,self.my,self.field,self.band,self.cuts,self.final,w='sex')

            # Compare new catalog to old one, get the ZP and FWHM out



    def init_phot(self,pl='n'):
        limmags = {}
        for counter,chip in enumerate(self.chips):
            chip_dir = os.path.join(self.out_dir,'MY%s'%self.my,self.field,self.band,chip)
            ana_dir = os.path.join(chip_dir,'ana')
            os.chdir(ana_dir)
            sexcat = self.sexcats[counter]
            cat = Table.read(sexcat).to_pandas()
            limmags[chip]=init_phot(self,chip,cat)
        return limmags
