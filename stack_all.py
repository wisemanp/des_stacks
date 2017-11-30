# -*- coding: utf-8 -*-
'''Stack Everything'''

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


def parser():
    parser = argparse.ArgumentParser(description='Stack some DES SN images')
    parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
    parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
    parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
    parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
    parser.add_argument('-wd','--workdir', help = 'Working directory', default = './')
    parser.add_argument('-l','--looptype', help ='Type of loop (can be "psf", "zp", "b" or "both")',required = False, default = 'both')
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
        print ('Not looping')
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

def looped_stack(logger,parsed):
    fields = parsed['fields']
    bands = parsed['bands']
    mys = parsed['mys']
    chips = parsed['chips']
    workdir = parsed['workdir']
    looptype = parsed['looptype']
    for f in fields:
        for b in bands:
            for my in mys:
                #do initial stack
                s = stack.Stack(f,b,my,chips,workdir)
                cut = -0.14
                s.do_my_stack(cut)
                
                qual = s.run_stack_sex(cut)
                # do stack with more generous cut
                s.do_my_stack(cut = -0.2)
                qual_plus = s.run_stack_sex(cut= cut-0.1)
                # do stack with more stringent cut
                s.do_my_stack(cut = cut-0.1)
                qual_minus = s.run_stack_sex(cut = cut-0.1)
                #test which is better
                #try going further in that direction
                #until we get no further.
                logger.info("Quality of normal stack %s" %qual)
                logger.info("Quality of generous stack %s" %qual_plus)
                logger.info("Quality of stringent stack %s" %qual_minus)

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
        looped_stack(logger,parsed)
    else:
        do_stack(logger,parsed)
    
