#!/home/wiseman/anaconda3/bin/python
# -*- coding: utf-8 -*-
##############################################################################
#
# Copyright (c) 2019 University of Southampton
# All Rights Reserved.
# 12/05/2018
##############################################################################

__author__ = "Philip Wiseman <p.s.wiseman@soton.ac.uk>"
__version__ = "0.2"
__date__ = "26/11/19"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
import configparser
import os
import logging
import argparse
import glob

from time import gmtime, strftime
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time

from des_stacks import des_stack as stack
from des_stacks.utils.loop_stack import iterate_source_loop, init_source_loop
sns.set_color_codes(palette='colorblind')
# define some DES specific lists
all_years = ['none','1','2','3','4'] # add 5 when available
all_fields = ['SN-X1','SN-X2','SN-X3','SN-C1','SN-C2','SN-C3','SN-E1','SN-E2','SN-S1','SN-S2']
all_chips = np.arange(1,62)
all_bands = ['g','r','i','z']

class optimiser():

    def __init__(self):
        parsed = self._parser()

    def _parser(self):
        parser = argparse.ArgumentParser(description='Stack some DES SN images')
        parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
        parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
        parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
        parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
        parser.add_argument('-wd','--workdir', help = 'Working directory [coadding]', default = 'coadding')
        parser.add_argument('-l','--looptype', help ='Parameters to optimize (can be "psf", "depth", or a comma separated list of those")',required = False, default = 'depth')
        parser.add_argument('-pr','--psfrange',help = 'Range to optimize psf in (min,max): [1.5,3]',required=False,default= '1.5,3.0')
        parser.add_argument('-tr','--teffrange',help = 'Range to optimize teff in (min,max): [0,0.5]',required=False,default= '0.0,0.5')
        parser.add_argument('-st','--step',help = 'Size of step in the cut you want to optimize over (psf,teff): [0.25,0.01]',required = False, default = '0.25,0.01')
        parser.add_argument('-pl','--plot',help='Plot a heatmap of where the best cuts are?',required = False,action = 'store_true')
        parser.add_argument('-t','--tidy',help = 'Tidy up temporary files after?',action = 'store_true')
        args=parser.parse_args()
        parsed = {}
        try:
            fields = args.field.split(',')

        except:
            try:
                fields = args.field[0].split(' ')
            except:
                fields =args.field

        for i in range(len(fields)):
            try:
                field = fields[i]
                field = 'SN-'+field
                fields[i]=field
            except:
                fields = 'SN-'+fields[0]
        parsed['fields']=fields

        try:
            bands = args.band.split(',')
        except:
            try:
                bands = args.band[0].split(' ')
            except:
                bands = args.band
        parsed['bands']=bands
        try:
            mys = args.minusyears.split(',')
        except:
            try:
                mys = args.minusyears[0].split(' ')
            except:
                mys = args.minusyears
        parsed['mys']=mys

        if args.chips != 'All':
            try:
                chips = args.chips[0].split(',')
            except:
                if args.chips[0][0]== '[':
                    chip_bounds = args.chips[0][1:-1].split(',')
                    chips = np.arange(int(chip_bounds[0]), int(chip_bounds[-1]))
                else:
                    chips = args.chips[0].split(' ')

        else:
            chips = args.chips
        parsed['chips']=chips
        print ('Parsed chips as %s'%chips)
        if not args.workdir:
            workdir = 'current'
        else:
            workdir = args.workdir

        parsed['workdir']=workdir

        try:
            loop_types = args.looptype.split(',')
            parsed['looptype']=loop_types
        except:
            parsed['looptype']='depth'

        try:
            parsed['teffrange'] = args.teffrange.split(',')
        except:
            parsed['teffrange'] = [0.0,0.5]
        try:
            parsed['psfrange'] = args.psfrange.split(',')
        except:
            parsed['psfrange'] = [1.5,3.0]
        try:
            parsed['step'] = args.step.split(',')
        except:
            parsed['step'] = [0.25,0.01]

        parsed['tidy']=args.tidy
        self.parsed = parsed
        self.plot = args.plot

    def optimise(self,f,b,y,ch):
        # a function that iterates through stacks until the best one is reached
        t0,t1,ts = float(self.parsed['teffrange'][0]),float(self.parsed['teffrange'][1]),float(self.parsed['step'][1])
        p0,p1,ps = float(self.parsed['psfrange'][0]),float(self.parsed['psfrange'][1]),float(self.parsed['step'][0])
        wd,lt = self.parsed['workdir'],self.parsed['looptype'][0]
        print(t0,t1,ts)
        print(p0,p1,ps)
        print (lt)
        teff_range = np.arange(t0,t1,ts)
        psf_range = np.arange(p0,p1,ps)
        lim_df = pd.DataFrame(index = [str(r) for r in psf_range],columns=[str(r) for r in teff_range])
        psf_df = pd.DataFrame(index = [str(r) for r in psf_range],columns=[str(r) for r in teff_range])#create the DataFrame to put the quality measurements in
        lim_df.name = 'depth'
        psf_df.name = 'psf'

        for psf_cut in psf_range:
            for teff_cut in teff_range:
                lim,psf = self.do_stack(f,b,y,ch,wd,cuts = {'zp':None,'teff':teff_cut,'psf':psf_cut})
                lim_df.loc[str(psf_cut),str(teff_cut)] = lim
                psf_df.loc[str(psf_cut),str(teff_cut)] = psf
        best={'depth':None,'psf':None}

        '''smaller_teff_step = ts/5
        smaller_psf_step = ps/5
        if lt=='depth':
            teff_start = best['depth'][1]
            psf_start = best['depth'][0]
        elif lt=='psf':
            teff_start = best['psf'][1]
             psf_start = best['psf'][0]
        elif lt=='both':
            teff_start = np.mean(best['depth'][1],best['psf'][1])
            psf_start = np.mean(best['depth'][0],best['psf'][0])
        zoomed_teffrange = np.arange(teff_start-float(ts)*5,teff_start+float(ts)*5,smaller_teff_step)
        zoomed_psfrange = np.arange(psf_start-float(ps)*5,psf_start+float(ps)*5,smaller_psf_step)
        for newpsf in zoomed_psfrange:
            lim_df = lim_df.append(pd.DataFrame(index=[str(newpsf)],columns=lim_df.columns))
            psf_df = psf_df.append(pd.DataFrame(index=[str(newpsf)],columns=psf_df.columns))
            for newteff in zoomed_teffrange:
                lim_df[str(newteff)] = ''
                psf_df[str(newteff)] = ''
                lim,psf = do_stack(f,b,y,ch,wd,cuts = {'zp':None,'teff':newteff,'psf':newpsf})
                lim_df.loc[str(newpsf),str(newteff)] = lim
                psf_df.loc[str(newpsf),str(newteff)] = psf'''
        for df in [lim_df,psf_df]:
            best[df.name] = [np.float(np.argmax(df.max(axis=1))),np.float(np.argmax(df.max(axis=0)))]
            # ADD TO PLOT!
        if self.plot:

            f1,ax1 = plt.subplots()
            depthmin = np.min(lim_df.min().values)
            depthmax = np.max(lim_df.max().values)
            depthrang = depthmax-depthmin
            lim_df = lim_df.astype(float)
            psf_df = psf_df.astype(float)

            sns.heatmap(lim_df,ax=ax1,cmap='Oranges',cbar_kws={'label': 'Limiting Magnitude'})
            ax1.set_xlabel('$\\tau_{effective} cut$')
            ax1.set_ylabel('PSF cut')
            plt.savefig('/media/data3/wiseman/des/coadding/optimise/optimize_teff_%s_%s_%s_%s.pdf'%(f,b,y,ch[0]))
            plt.close()
            f2,ax2 = plt.subplots()

            sns.heatmap(psf_df,ax=ax2,cmap='Blues',cbar_kws={'label': 'Limiting Magnitude'})
            ax2.set_xlabel('$\\tau_{effective} cut$')
            ax2.set_ylabel('PSF cut')
            plt.savefig('/media/data3/wiseman/des/coadding/optimise/optimize_psf_%s_%s_%s_%s.pdf'%(f,b,y,ch[0]))

        return best
    def do_stack(self,f,b,y,ch,wd,cuts):
        #Performs the actual stack for a given set of cuts, and returns the limiting magnitudes and psf
        print ('Making stack of',f,b,y,ch,wd,cuts)
        s = stack.Stack(f,b,y,ch,wd,cuts,db=True)
        scifile = os.path.join(s.band_dir,'ccd_%s_%s_%.2f_%s_clipweighted_sci.fits'%(ch[0],b,cuts['teff'],cuts['psf']))
        if not os.path.isfile(scifile):
            print ('Did not find a file for these cuts; doing stack')
            s.do_my_stack(cuts=cuts,final=True)
        else:
            print ('Found a stacked file for these cuts; going to source')
        s.ana_dir = os.path.join(s.band_dir,ch[0],'ana')
        sourcename = os.path.join(s.ana_dir,'MY%s_%s_%s_%s_%.2f_%s_clipweighted_sci.sourcecat' %(y,f,b,ch[0],cuts['teff'],cuts['psf']))
        print ('Looking for file under the name: %s'%sourcename)
        if os.path.isfile(sourcename):
            print ('Found a sourcecat for these cuts at: %s'%sourcename)
            s.sourcecats = [sourcename]
            s.cuts=cuts
        else:
            print ('No sourcecat yet; running source extractor')
            print ('Sending %s to run_stack_source'%cuts)
            s.run_stack_source(cuts=cuts,final=True)
        s.cutstring = '%s_%s'%(cuts['teff'],cuts['psf'])
        #lim = np.median(s.init_phot()[ch[0]][-1])
        skylim = s.init_phot()[ch[0]][2]
        psf = np.loadtxt(os.path.join(s.band_dir,ch[0],'ana','%s_ana.qual'%s.cutstring))[2]
        psf_err = np.loadtxt(os.path.join(s.band_dir,ch[0],'ana','%s_ana.qual'%s.cutstring))[3]
        np.savetxt(os.path.join(s.ana_dir,'%s_limmags.txt'%s.cutstring),np.array([skylim,psf,psf_err]))
        return (skylim,psf)

def main():
    o = optimiser()
    parsed = o.parsed
    chips = [[str(chip)] for chip in parsed['chips'][0].split(',')]
    best_teff_df = pd.read_csv('/media/data3/wiseman/des/coadding/optimise/best_teff.csv',header=0)
    best_psf_df = pd.read_csv('/media/data3/wiseman/des/coadding/optimise/best_psf.csv',header=0)

    for y in parsed['mys']:
        for f in parsed['fields']:
            for b in parsed['bands']:

                for ch in chips:

                    print ('Sending chip %s to optimize'%ch)
                    best = o.optimise(f,b,y,ch)

                    best_teff_df = best_teff_df.append(pd.DataFrame([[f,b,ch,best['depth'][0],best['depth'][1]]],columns=best_teff_df.columns))
                    best_psf_df = best_psf_df.append(pd.DataFrame([[f,b,ch,best['psf'][0],best['psf'][1]]],columns=best_psf_df.columns))
    best_teff_df.to_csv('/media/data3/wiseman/des/coadding/optimise/best_teff.csv',index=False)
    best_psf_df.to_csv('/media/data3/wiseman/des/coadding/optimise/best_psf.csv',index=False)
if __name__=="__main__":
    main()
