#!/home/wiseman/anaconda3/bin/python
# -*- coding: utf-8 -*-
'''Wrapper to do common aperture photometry on DES SN hosts'''
# exectuable to run through the entire stack process
import numpy as np
import logging
import argparse
import pandas as pd
from time import gmtime, strftime
import multiprocessing
from multiprocessing import Process
import pathos.pools as pp
#Note: import this first else it crashes importing sub-modules
from des_stacks import des_stack as stack
from des_stacks.analysis.astro import cap_phot_sn,cap_sn_lookup
import os

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--sn_name',help='Full DES supernova ID, e.g. DES16C2nm (case sensitive)',default=None)
    parser.add_argument('-l','--namelist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-av','--avoid',help='Avoid these SN, if it is in the list given by l. e.g. [DES16C2nm]',default=None)
    parser.add_argument('-wd','--workdir',help='Path to directory to work in',default = '/media/data3/wiseman/des/coadding')
    parser.add_argument('-sf','--savename',help='Filename to save results to',default=None)
    parser.add_argument('-ow','--overwrite',help='Overwrite existing results?',action = 'store_true')
    parser.add_argument('-th','--threshold',help='Distance threshold for host galaxy searching (arcsecs)',default=60)
    parser.add_argument('-v','--version',help='Way of getting photometry. 1 = Do own CAP; 2 = Go into chip results file',default=2)
    return parser.parse_args()


def cap_worker(arg_pair):
    args,sn = arg_pair[0],arg_pair[1]
    wd,dist_thresh,savename = [args[i]for i in range(len(args))]

    cap_sn_lookup(sn, wd = wd,savename = savename, dist_thresh = dist_thresh,autocuts=True)


def cap(args,logger):
    avoid_list = []
    args.version = int(args.version)
    if args.avoid:
        avoid_list = [i for i in args.avoid.split(',')]
    else:
       avoid_list = [None]
    if args.sn_name:
        logger.info('Have been given a name, doing CAP for just %s'%args.sn_name)
        if not args.savename:
            cap_sn_lookup(args.sn_name,wd=args.workdir,dist_thresh =  args.threshold,autocuts=True)
        else:
            logger.info('Been given a savename, so doing cap now and saving it to that')
            cap_sn_lookup(args.sn_name,wd=args.workdir,savename = args.savename,dist_thresh =  args.threshold,autocuts=True)
    else:
        logger.info("Pulling list of SN on which to do common aperture photometry")
        sn_list = np.genfromtxt(args.namelist,dtype=str,delimiter='\n')
        logger.info("Doing CAP on following input list")
        logger.info(sn_list)
        pool_size = 8#multiprocessing.cpu_count()*2
        act = multiprocessing.active_children()
        pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
        pool._clear()
        pool._serve()

        if not args.savename:

            try:
                done_sn = pd.read_csv('/media/data3/wiseman/des/coadding/results/all_sn_phot.csv',index_col=0)
            except:
                done_sn = pd.DataFrame(columns=['BAND', 'CLASS_STAR', 'ELONGATION', 'FWHM_WORLD', 'KRON_RADIUS', 'MAGERR_APER', 'MAGERR_AUTO', 'MAG_APER', 'MAG_AUTO', 'SN_NAME', 'X_WORLD', 'Y_WORLD','LIMMAG'])
        else:

            try:
                logger.info('Reading in results file to find out which ones I still need to do')
                logger.info(os.path.join('/media/data3/wiseman/des/coadding/results/',args.savename))
                done_sn = pd.read_csv(os.path.join('/media/data3/wiseman/des/coadding/results/',args.savename),index_col=0)
                logger.info('Read in %s'%os.path.join('/media/data3/wiseman/des/coadding/results',args.savename))
            except:
                done_sn = pd.DataFrame(columns=['BAND', 'CLASS_STAR', 'ELONGATION', 'FWHM_WORLD', 'KRON_RADIUS', 'MAGERR_APER', 'MAGERR_AUTO', 'MAG_APER', 'MAG_AUTO', 'SN_NAME', 'X_WORLD', 'Y_WORLD','LIMMAG'])

        all_args = []
        worker_args = [args.workdir,args.threshold,args.savename]
        for sn_name in sn_list :
            logger.info("Doing common aperture photometry on %s"%sn_name)
            logger.info(sn_name)

            if sn_name not in done_sn.SN_NAME.unique():


                if sn_name not in avoid_list:
                    all_args.append([worker_args,sn_name])

            elif args.overwrite == True:
                if sn_name not in avoid_list:
                    if not args.savename:

                        all_args.append([worker_args,sn_name])
                    else:
                        all_args.append([worker_args,sn_name])

            else:
                logger.info("Result for %s already in result file, and you told me not to overwrite it. Going to next one!"%sn_name)
        results = pool.map(cap_worker,all_args)

        pool.close()
        pool.join()
        return results
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
