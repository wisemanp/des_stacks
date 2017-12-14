
import numpy as np
import pandas as pd
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import datetime
import os
import logging
import time

def astrometry(stack,chip,sexcat,phot_type='AUTO'):
    '''Load in the existing DES and the newly SExtracted catalogs'''
    old_cat = os.path.join(stack.cat_dir,'%s_All_filters_3.csv'%(stack.field[3]))
    old = pd.DataFrame.from_csv(old_cat)
    sexdat = fits.getdata(sexcat,ext=1)
    new =pd.DataFrame(sexdat)
    new_obj = SkyCoord(ra=new['X_WORLD']*u.degree,dec =new['Y_WORLD']*u.degree)
    old_obj = SkyCoord(ra=old['RA_%s'%stack.band]*u.degree,dec =old['DEC_%s'%stack.band]*u.degree)
    # match the catalogs
    idx, d2d, d3d = new_obj.match_to_catalog_sky(old_obj)
    match_ids = idx
    match_dists = d2d.arcsec
    # get old cat mags of the matched objects
    init_match_cat_mag =old['CLIPPED_MEAN_%s'%stack.band].iloc[match_ids]
    init_match_cat_magerr =old['CLIPPED_SIGMA_%s'%stack.band].iloc[match_ids]
    # get indices of the objects that are within a specified distance of their matches
    dist_cut =2.0
    good_inds = np.nonzero(match_dists < dist_cut)[0]

    good_new_ra = new['X_WORLD'].iloc[good_inds]
    good_new_dec = new['Y_WORLD'].iloc[good_inds]
    reg = open(os.path.join(stack.out_dir,'MY%s'%stack.my,stack.field,stack.band,'ccd_%s'%chip+'_calib.reg'),'w')
    print ('global color=green',file=reg)
    for i in range(len(good_new_ra.values)):
        print ('fk5; circle(%s,%s,3p)'%(good_new_ra.values[i],good_new_dec.values[i]),file=reg)
    reg.close()
    # find the new mags that correspond to that
    good_new_mag = new['MAG_%s'%phot_type].iloc[good_inds]
    good_new_magerr = new['MAGERR_%s'%phot_type].iloc[good_inds]
    # and the old ones
    good_cat_mag = init_match_cat_mag.iloc[good_inds]
    good_cat_magerr = init_match_cat_magerr.iloc[good_inds]
    # subtract to get the frame ZP
    zp = np.median(good_cat_mag.values - good_new_mag.values)
    psf = np.median(new['FWHM_WORLD']*3600)
    return zp,psf
