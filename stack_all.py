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


import .des_stack as stack

def parser(cmdline):
    parser = argparse.ArgumentParser(description='Stack some DES SN images')
    parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
    parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
    parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
    parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
    parser.add_argument('-wd','--workdir', help = 'Working directory', default = './')
    args=parser.parse_args(cmdline)
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
    #start the logger
    logging.basicConfig(level=logging.DEBUG)
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logger = logging.getLogger(__name__)
    
    #read in parameters from the command line
    fields,bands,mys,chips,workdir = parser()
    for f in fields:
        for b in bands:
            for my in mys:
                s = stack.Stack(logger,f,b,my,chips,workdir)
                s.do_my_stack()
