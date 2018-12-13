import pandas as pd
import numpy as np
import os
import argparse
import glob
import subprocess
import pathos.pools as pp
import multiprocessing
import astropy.io.fits as fits
from multiprocessing import Process

from des_stacks.utils.gen_tools import get_good_des_chips, get_des_bands

good_des_chips = get_good_des_chips()
bands = get_des_bands()

def worker(img):
    if os.path.isfile(img.replace('.fits','.weight.fits')):
        os.remove(img.replace('.fits','.weight.fits'))
    swarp_cmd = [
    'swarp',
    '%s'%img,
    '-WEIGHTOUT_NAME','%s'%img.replace('.fits','.wgt.fits'),
    '-COMBINE','N',
    '-RESAMPLE','Y',
    '-DELETE_TMPFILES','N',
    '-BACK_SIZE','32'
    ]
    print ('Doing following: \n %s '%swarp_cmd)
    p = subprocess.Popen(swarp_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outs,errs = p.communicate()
    if os.path.isfile(img.replace('.fits','.resamp.fits')):
        os.remove(img.replace('.fits','.resamp.fits'))
    return
def worker2(img):
    fl = fits.open(img)
    wm = fits.open(img.replace('.fits','.resamp.weight.fits'))
    wgt = np.median(wm[0].data)
    wgtarr = np.ones_like(fl[0].data)*wgt
    fl[0].data = wgtarr
    fl[0].header['TYPE'] = 'WGTMAP  '
    os.remove(img.replace('.fits','.resamp.weight.fits'))
    fl.writeto(img.replace('.fits','.resamp.weight.fits'))
    fl.close()
    return
def multi_fn(lst):
    pool_size = multiprocessing.cpu_count()*2
    pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
    pool._clear()
    pool._serve()

    results = pool.map(worker2,lst)

def main():
    fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
    mys = ['1','2','3','4','5','none']
    for my in mys:
        for f in fields:
            f = 'SN-'+f
            for b in bands:
                os.chdir(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/','MY%s'%my,f,b))
                list_of_scis = glob.glob(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/','MY%s'%my,f,b,'*clipweighted_sci.fits'))
                multi_fn(list_of_scis)

if __name__=="__main__":
    main()
