import pandas as pd
import numpy as np
import os
import argparse
import glob
import subprocess
import pathos.pools as pp
import multiprocessing
from multiprocessing import Process

from des_stacks.utils.gen_tools import get_good_des_chips, get_des_bands

good_des_chips = get_good_des_chips()
bands = get_des_bands()

def worker(img):
    swarp_cmd = [
    'swarp',
    '%s'%img,
    '-COMBINE','N',
    '-RESAMPLE','Y',
    '-WEIGHTOUT_NAME','%s'%img.replace('.fits','wgt.fits')
    ]
    print ('Doing following: \n %s '%swarp_cmd)
    p = subprocess.Popen(swarp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    if os.path.isfile(img.replace('.fits','.resamp.fits')):
        os.remove(img.replace('.fits','.resamp.fits'))
    return
def multi_fn(lst):
    pool_size = multiprocessing.cpu_count()*2
    pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
    pool._clear()
    pool._serve()

    results = pool.map(worker,lst)

def main():
    fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
    mys = ['1','2','3','4','5','none']
    for my in mys:
        for f in fields:
            f = 'SN-'+f
            for b in bands:
                os.chdir(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/','MY%s'%my,f,b))
                list_of_scis = glob.glob(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/','MY%s'%my,f,b,'*clipweighted*'))
                multi_fn(list_of_scis)

if __name__=="__main__":
    main()
