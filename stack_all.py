# -*- coding: utf-8 -*-
'''Stack Everything'''
# exectuable to run through the entire stack process
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
import argparse
from time import gmtime, strftime
from des_stacks import des_stack as stack
from des_stacks.utils.loop_stack import iterate_sex_loop, init_sex_loop

def parser():
    parser = argparse.ArgumentParser(description='Stack some DES SN images')
    parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
    parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
    parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
    parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
    parser.add_argument('-wd','--workdir', help = 'Working directory', default = './')
    parser.add_argument('-l','--looptype', help ='Type of loop (can be "psf", "zp", "b" or "both")',required = False, default = 'both')
    parser.add_argument('-ps','--psfstep', help ='The size of the cut step if using psf',required=False, default =0.25)
    parser.add_argument('-zs','--zpstep', help ='Size of the cut step for zeropoint cuts',required=False,default = 0.025)
    args=parser.parse_args()
    parsed = {}
    try:
        fields = args.field.split(',')

    except:
        try:
            fields = args.field[0].split(' ')
        except:
            fields =args.field

    for i in range(len(fields[0])):
        try:
            field = fields[i]
            field = 'SN-'+field
            fields[i]=field
        except:
            field = 'SN-'+fields[0]
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
    print ('Parsed chips as: %s' % chips)

    if not args.workdir:
        workdir = 'current'
    else:
        workdir = args.workdir

    parsed['workdir']=workdir

    try:
        loop_type = args.looptype
        parsed['looptype']=loop_type
    except:
        pass
    try:
        parsed['psf_step']=args.psfstep
    except:
        parsed['psf_step']=0.25

    try:
        parsed['zp_step']=args.zpstep
    except:
        parsed['zp_step'] = 0.05
    return parsed

def simple_stack(logger,parsed):
    '''code to run the stacks'''
    #read in parameters from the command line
    fields = parsed['fields']
    bands = parsed['bands']
    mys = parsed['mys']
    chips = parsed['chips']
    workdir = parsed['workdir']
    logger.info("Parsed command line and will work on:\n Fields: %s \n Bands: %s \n MYs: %s \n Chips: %s"%(fields,bands,mys,chips))
    for f in fields:
        for b in bands:
            for my in mys:
                s = stack.Stack(f,b,my,chips,workdir)
                s.do_my_stack()
                s.run_stack_sex(cut=-0.15)

def looped_stack(logger,parsed,init_cut,init_step):
    fields = parsed['fields']
    bands = parsed['bands']
    mys = parsed['mys']
    chips = parsed['chips']
    workdir = parsed['workdir']
    loop_type = parsed['looptype']
    psfstep = parsed['psf_step']
    for f in fields:
        for b in bands:
            for my in mys:
                logger.info("Doing the initial stack")
                norm,stringent,generous= init_sex_loop(logger,f,b,my,chips,loop_type,init_cut,init_step,workdir)
                for chip in chips:
                    logger.info("Now starting the loop to get the best stack for chip %s" %chip)
                    # find out which of the stacks is deeper
                    if generous.loc[chip,'zp']>norm.loc[chip,'zp']:
                        logger.info("Generous cut gives deeper stack than normal stack")
                        if generous.loc[chip,'zp']>stringent.loc[chip,'zp']:
                            logger.info("Generous cut gives deeper stack than stringent stack")
                            logger.info("Heading in the generous direction!")
                            # The generous stack is deepest, so go a bit further
                            cuts,qual = iterate_sex_loop(logger, f, b, my, chip, loop_type,init_cut,init_step,workdir,'ge',norm)
                            logger.info("Saved stack and exiting loop for chip %s" %chip)
                        elif generous.loc[chip,'zp']<=stringent.loc[chip,'zp']:
                            logger.info("Generous cut not as deep as stringent stack. Exploring options.")
                            # explore some options. Unlikely to get to this point.
                    elif generous.loc[chip,'zp']<norm.loc[chip,'zp']:
                        logger.info("Stack with a generous cut not as deep as the normal stack")
                        if norm.loc[chip,'zp']<stringent.loc[chip,'zp']:
                            # The stringent cut is deepest, so go more stringent
                            logger.info("Stack with a normal cut also not as deep as the stringent cut")
                            logger.info("Heading in the stringent direction!")
                            cuts,qual = iterate_sex_loop(logger, f, b, my, chip, loop_type,init_cut,init_step,workdir,'st',norm)
                            logger.info("Saved stack and exiting loop for chip %s" %chip)
                        elif norm.loc[chip,'zp']>stringent.loc[chip,'zp']:
                            logger.info("Stack with normal cut is deeper than the stringent cut. Exploring options")

if __name__=="__main__":
    logger = logging.getLogger('stack_all.py')
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("***********************************")
    logger.info("Initialising *** stack_all.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    logger.info("***********************************")
    parsed = parser()
    if 'looptype' in parsed.keys():
        if parsed['looptype']in ['b','both','zp','z','psf']:
            looped_stack(logger,parsed,[-0.15,2.5],[parsed['zp_step'],parsed['psf_step']])
        else:
            simple_stack(logger,parsed)
    else:
        simple_stack(logger,parsed)
