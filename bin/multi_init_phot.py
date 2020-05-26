import numpy as np
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
from astropy.coordinates import SkyCoord
import logging
from astropy.table import Table
import astropy.io.fits as fits
import os
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM 
from astropy import wcs
from des_stacks import des_stack as stack
from des_stacks.bin import stack_all
from des_stacks.utils import stack_tools,source_tools,gen_tools
from des_stacks.analysis import astro
from des_stacks.utils.gen_tools import mc_robust_median as r_median

import time
import _pickle as cpickle
import itertools
import multiprocessing
from multiprocessing import Process
import pathos.pools as pp
bands = gen_tools.get_des_bands()


good_des_chips = []
for c in range(1,63):
    if c not in [2,31,61]:
        good_des_chips.append(c)
fields = ['E1']#,'E2','S1','S2','C1','C2','C3','X1','X2','X3']
bands = ['g','r']#,'i','z']

def init_phot_worker(arg_pair):
    args, chip = arg_pair[0],arg_pair[1]
    my,f,b,cuts = [args[i] for i in range(len(args))]
    ch = int(chip)
    bd = os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/MY%s/'%my,f,b)
    
    cat_fn = os.path.join(bd,str(chip),'ana',
                       'MY%s_%s_%s_%s_%.2f_%.1f_clipweighted_sci.sourcecat'%(my,f,b,
                                                                        str(ch),cuts['teff'],cuts['psf']))
    cat = Table.read(cat_fn).to_pandas()
    s = stack.Stack(f,b,my,ch,'coadding',cuts,db=False,new=True)
    s.cuts = cuts
    astro.init_phot(s,str(chip),cat)
    return

def multi_init_phot(my,f,b,chips):
    #cuts = {'psf':1.3,'teff':0.02}
    cuts =stack_tools.get_cuts(f,b) 
    args = [my,f,b,cuts]
    pool_size = multiprocessing.cpu_count()*2
    act = multiprocessing.active_children()
    pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
    pool._clear()
    pool._serve()

    chips = list(chips)

    all_args = []
    for c in chips:
        all_args.append([args,c])
        #p = Process(target=worker,args=(args,c))
        #p.start()
        #p.join()

    results = pool.map(init_phot_worker,all_args)

    pool.close()
    pool.join()
    return results

def main():
    for f in fields:
        f = 'SN-'+f
        for b in bands:
            cuts =stack_tools.get_cuts(f,b) 
            for y in ['1']:#,2,3,4,5]:
                #cuts = {'teff':0.02,'psf':1.3}
                
            
                multi_init_phot(y,f,b,good_des_chips)
if __name__=="__main__":
    main()
