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
all_logger = logging.getLogger('stack_all.py')
all_logger.setLevel(logging.DEBUG)
logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def parser():
    parser = argparse.ArgumentParser(description='Stack some DES SN images')
    parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
    parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
    parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
    parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
    parser.add_argument('-wd','--workdir', help = 'Working directory', default = './')
    args=parser.parse_args()
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

    try:
        bands = args.band.split(',')
    except:
        try:
            bands = args.band[0].split(' ')
        except:
            bands = args.band
    try:
        mys = args.minusyears.split(',')
    except:
        try:
            mys = args.minusyears[0].split(' ')
        except:
            mys = args.minusyears

    if args.chips != 'All':
        try:
            chips = args.chips.split(',')
        except:
            if args.chips[0][0]== '[':
                chip_bounds = args.chips[0][1:-1].split(',')
                print (chip_bounds)
                chips = np.arange(int(chip_bounds[0]), int(chip_bounds[-1]))
                print (chips)
            else:
                chips = args.chips[0].split(' ')

    else:
        chips = args.chips

    if not args.workdir:
        workdir = 'current'
    else:
        workdir = args.workdir
    return fields, bands, mys, chips, workdir

def do_stack():
    '''code to run the stacks'''
    #read in parameters from the command line
    fields,bands,mys,chips,workdir = parser()
    all_logger.info("Parsed command line and will work on:\n Fields: %s \n Bands: %s \n MYs: %s \n Chips: %s"%(fields,bands,mys,chips))
    for f in fields:
        for b in bands:
            for my in mys:
                s = stack.Stack(f,b,my,chips,workdir)
                s.do_my_stack()

if __name__=="__main__":
    all_logger.info("***********************************")
    all_logger.info("Initialising *** stack_all.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    all_logger.info("***********************************")
    do_stack()
