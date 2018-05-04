# -*- coding: utf-8 -*-
'''Wrapper to do common aperture photometry on DES SN hosts'''
# exectuable to run through the entire stack process
import numpy as np
import logging
import argparse
from time import gmtime, strftime
#Note: import this first else it crashes importing sub-modules
from des_stacks import des_stack as stack
from des_stacks.analysis.astro import cap_phot_sn

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--sn_name',help='Full DES supernova ID, e.g. DES16C2nm (case sensitive)',default=None)
    parser.add_argument('-l','--namelist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-wd','--workdir',help='Path to directory to work in',default = '/media/data3/wiseman/des/coadding')
    return parser.parse_args()

def cap(args,logger):
    if args.sn_name:
        logger.info("Doing common aperture photometry on %s"%args.sn_name)
        cap_phot_sn(args.sn_name,args.workdir)
    else:
        logger.info("Pulling list of SN on which to do common aperture photometry")
        for sn_name in np.loadtxt(args.namelist,dtype='str'):
            logger.info("Doing common aperture photometry on %s"%sn_name)
            cap_phot_sn(sn_name,args.workdir)

if __name__ == "__main__":
    logger = logging.getLogger('sn_cap.py')
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("***********************************")
    logger.info("Initialising *** sn_cap.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    logger.info("***********************************")
    args = parser()
    cap(args,logger)
