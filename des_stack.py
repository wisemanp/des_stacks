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
import logging

from utils.stack_tools import make_good_frame_list, make_swarp_cmd, get_dessn_obs, get_des_obs_year

class Stack():
    def __init__(self, logger, field, band, my, chips ,working_dir):
        self.logger = logger
        self._define_paths()
        self._get_info()
        self._get_configs()
        self.logger.info('Initialising stack successful')
        self.field = field
        self.band = band
        self.my =my
        self.chips=chips
        self.coadding_dir =working_dir

    def _define_paths(self):
        '''Sets the base paths in which to do all of the work'''
        self.logger.info('Defining the paths so we can make our stack nice and tidily')
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
        print ('Writing some stuff')
    ########################################################
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
    def do_my_stack(self, cut=-0.15):
        '''Does a stack defined by the parameters from the Stack object it is passed.
        parameters:
        cut (float, optional): the zeropoint cut to be used for the stack
        returns:
        none
        '''
        field = self.field
        band = self.band
        my = self.my
        self.logger.info('Initiating stack on {0} in {1} band'.format(field,band))
        #does list of good frames exist?
        if not os.path.isfile(os.path.join(self.list_dir,'good_exps_%s_%s_%s.fits'%(field,band,cut))):
            #get the list of good frames
            self.good_frame = make_good_frame_list(self,field,band)
            self.logger.info('No good frame list yet, making a new one')
        else:
            good_fn = os.path.join(self.list_dir,'good_exps_%s_%s_%s.fits'%(field,band,cut))
            self.logger.info('Reading in list of good frames from {0}'.format(good_fn))
            good_table = Table.read(good_fn)
            self.good_frame = good_table.to_pandas()

        # how many chips are we stacking?
        chips = self.chips
        if chips == 'All':
            chips = self.info_df.CCDNUM.sort_values().unique()
        # get swarp commands
        for y in my:
            self.logger.info('Stacking {0} in {1} band, skipping year {2}'.format(field,band,my))
            for chip in chips:
                os.chdir(self.temp_dir)
                self.logger.info('Stacking CCD {0}'.format(chip))
                cmd = make_swarp_cmd(self,y,field,chip,band)
                self.logger.info('Executing command: {0}'.format(cmd))
                p = subprocess.Popen(cmd)
                p.wait()
                self.logger.info('Finish stacking chip {0}'.format(chip))
            self.logger.info('Finished stacking chips {0} for MY {1}'.format(chips,y))
        self.logger.info('Stacking complete!')
    def run_stack_sex(self):
        pass
        # Do SExtractor on the complete stacks

    def loop_after_sex(self):
        pass
        # Have a look at the sex files and compare them to something.
