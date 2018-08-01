import numpy as np
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.coordinates import SkyCoord
import logging
from astropy.table import Table
import astropy.io.fits as fits
import os
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy import wcs
from des_stacks import des_stack, stack_all
from des_stacks.utils import stack_tools,sex_tools
from des_stacks.analysis import astro
sns.set_color_codes(palette='colorblind')
import time
import _pickle as cpickle
import itertools
import sys
import os
import numpy as np
import math as m

#####################   DEFINE GLOBAL VARIABLES   #################################################

rad  = np.pi/180                   # convert deg to rad
pix_arcsec = 0.264                 # pixel scale (arcsec per pixel)
pix2_arcsec2 = 0.264**2            # pix^2 to arcsec^2 conversion factor
pix2_deg2 = pix2_arcsec2/(3600**2) # pix^2 to deg^2 conversion factor

def get_DLR_ABT(RA_SN, DEC_SN, RA, DEC, A_IMAGE, B_IMAGE, THETA_IMAGE, angsep):
    # inputs are arrays

    global numFailed
    rPHI = np.empty_like(angsep)
    d_DLR = np.empty_like(angsep)

    # convert from IMAGE units (pixels) to WORLD (arcsec^2)
    A_ARCSEC = A_IMAGE*pix_arcsec
    B_ARCSEC = B_IMAGE*pix_arcsec

    # angle between RA-axis and SN-host vector
    GAMMA = np.arctan((DEC_SN - DEC)/(np.cos(DEC_SN*rad)*(RA_SN - RA)))

    # angle between semi-major axis of host and SN-host vector
    PHI = np.radians(THETA_IMAGE) + GAMMA # angle between semi-major axis of host and SN-host vector

    rPHI = A_ARCSEC*B_ARCSEC/np.sqrt((A_ARCSEC*np.sin(PHI))**2 +
                                     (B_ARCSEC*np.cos(PHI))**2)

    # directional light radius
    #  where 2nd moments are bad, set d_DLR = 99.99
    d_DLR = angsep/rPHI

    return [d_DLR, A_ARCSEC, B_ARCSEC, rPHI]



sncand = Table.read('/media/data3/wiseman/des/coadding/catalogs/sn_cand.fits').to_pandas()
sngals = pd.read_csv('/media/data3/wiseman/des/coadding/catalogs/sngals2.csv',index_col=0)
sngals['DLR_reprod']=''
sngals['DLR_diff']=''
sngals['DLR_RANK_new']=''
n=0
galindex = []
for i in sngals.index:
    n+=1
    print ('Doing %s of %s'%(n,len(sngals)))
    snra,sndec,sn_name = sncand[['RA','DEC','TRANSIENT_NAME']].loc[i]
    sn_name = sn_name.strip(' ')
    match = sngals[sngals['TRANSIENT_NAME']==sn_name]
    sncoord = SkyCoord(ra=snra*u.degree,dec=sndec*u.degree)
    galcoords = SkyCoord(ra=match['RA'].values*u.deg,dec=match['DEC'].values*u.deg)
    dist = sncoord.separation(galcoords)
    dists = np.array([float(dist[i].to_string(unit=u.arcsec,decimal=True)) for i in range(len(dist))])
    new_dlr = get_DLR_ABT(snra,sndec,match['RA'].values,match['DEC'].values,
            match['A_IMAGE'].values,match['B_IMAGE'].values,
           match['THETA_IMAGE'].values,dists)[0]
    if len(new_dlr)>3:
        print (sn_name, len)


    for ind in match.index:
        galindex.append(ind)
    sngals['DLR_reprod'].loc[match.index]=new_dlr
    new_rank = sngals['DLR_reprod'].loc[match.index].rank()
    for r in range(len(new_rank)):
        if new_dlr[r]>4:
            new_rank.iloc[r]=new_rank.iloc[r]*-1
    sngals['DLR_RANK_new'].loc[match.index]=new_rank
    sngals['DLR_diff'].loc[match.index]=sngals['DLR'].loc[match.index]-new_dlr
sngals.to_csv('/media/data3/wiseman/coadding/catalogs/sngals_updated.csv')
