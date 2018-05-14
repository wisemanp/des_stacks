# -*- coding: utf-8 -*-
##############################################################################
#
# Copyright (c) 2018 University of Southampton
# All Rights Reserved.
# 12/05/2018
##############################################################################

__author__ = "Philip Wiseman <p.s.wiseman@soton.ac.uk>"
__version__ = "0.1"
__date__ = "12/05/18"

import numpy as np
import pandas as pd
import seaborn as sns
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
from des_stacks.utils.loop_stack import iterate_sex_loop, init_sex_loop

# define some DES specific lists
all_years = ['none','1','2','3','4'] # add 5 when available
all_fields = ['SN-X1','SN-X2','SN-X3','SN-C1','SN-C2','SN-C3','SN-E1','SN-E2','SN-S1','SN-S2']
all_chips = np.arange(1,62)
all_bands = ['g','r','i','z']

def parser():
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

    return parsed

def optimize(f,b,y,ch,wd,t0,t1,ts,p0,p1,ps):
    # a function that iterates through stacks until the best one is reached
    print(t0,t1,ts)
    print(p0,p1,ps)
    teff_range = np.arange(t0,t1,ts)
    psf_range = np.arange(p0,p1,ps)
    lim_df = pd.DataFrame(index = [str(r) for r in psf_range],columns=[str(r) for r in teff_range])
    psf_df = pd.DataFrame(index = [str(r) for r in psf_range],columns=[str(r) for r in teff_range])#create the DataFrame to put the quality measurements in
    lim_df.name = 'depth'
    psf_df.name = 'psf'

    for psf_cut in psf_range:
        for teff_cut in teff_range:
            lim,psf = do_stack(f,b,y,ch,wd,cuts = {'zp':None,'teff':teff_cut,'psf':psf_cut})
            lim_df.loc[str(psf_cut),str(teff_cut)] = lim
            psf_df.loc[str(psf_cut),str(teff_cut)] = psf
    best={'depth':None,'psf':None}
    for df in [lim_df,psf_df]:
        best[df.name] = [np.float(np.argmax(df.max(axis=1))),np.float(np.argmax(df.max(axis=0)))]
        # TO BE ADDED: MAKE A PLOT
    smaller_teff_step = float(parsed['step'][1])/5
    smaller_psf_step = float(parsed['step'][0])/5
    if parsed['looptype']=='depth':
        teff_start = best['depth'][1]
        psf_start = best['depth'][0]
    elif parsed['looptype']=='psf':
        teff_start = best['psf'][1]
        psf_start = best['psf'][0]
    elif parsed['looptype']=='both':
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
            psf_df.loc[str(newpsf),str(newteff)] = psf
    for df in [lim_df,psf_df]:
        best[df.name] = [np.float(np.argmax(df.max(axis=1))),np.float(np.argmax(df.max(axis=0)))]
        # ADD TO PLOT!

    return best
def do_stack(f,b,y,ch,wd,cuts):
    #Performs the actual stack for a given set of cuts, and returns the limiting magnitudes and psf
    print ('Making stack of',f,b,y,ch,wd,cuts)
    s = stack.Stack(f,b,y,ch,wd,cuts)
    s.do_my_stack(cuts=cuts,final=False)
    s.run_stack_sex(cuts=cuts,final=False)
    lim = np.median(s.init_phot()[ch])
    psf = np.mean(np.loadtxt(os.path.join(s.band_dir,str(ch),'ana','%s_ana.qual'%s.cutstring))[1:])
    return (lim,psf)

def main():
    parsed = parser()
    for y in parsed['mys']:
        for f in parsed['fields']:
            for b in parsed['bands']:
                for ch in parsed['chips']:
                    t0,t1,ts = float(parsed['teffrange'][0]),float(parsed['teffrange'][1]),float(parsed['step'][1])
                    p0,p1,ps = float(parsed['psfrange'][0]),float(parsed['psfrange'][1]),float(parsed['step'][0])

                    best = optimize(f,b,y,ch,parsed['workdir'],t0,t1,ts,p0,p1,ps)
    print (best)

if __name__=="__main__":
    main()
