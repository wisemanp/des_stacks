import numpy as np
import pandas as pd
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import NDData
from astroimtools import nddata_stats
import datetime
import os
import sys
import logging
import time
import seaborn as sns
import matplotlib
import glob
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import UnivariateSpline as spln

resname = sys.argv[1]
def get_zs():
    survey_flags = {
    'DES_AAOmega':['1','2','3','4','6'],
    'ZFIRE_UDS':['3'],
    'NOAO_0522':['3','4','6'],
    'NOAO_0334':['3','4','6'],
    'N17B331':['4','6'],
    'MOSDEF':['Any'],
    'SpARCS':['1','2'],
    'PanSTARRS_AAOmega   ':['3','4','6'],
    'PanSTARRS_MMT': ['3','4','6'],
    'PRIMUS': ['3','4'],
    'NED': ['Any'],
    'UDS_FORS2':['A','B'],
    'UDS_VIMOS':['3','4'],
    'ACES': ['3','4'],
    'SDSS': ['0'],
    '6dF': ['4'],
    'ATLAS':['Any'],
    '2dFGRS':['3','4'],
    'GAMA':['4'],
    'SNLS_FORS           ':['1','2'],
    'CDB':['Any'],
    'VVDS_DEEP':['3','4','13','14','23','24','213','214'],
    'VVDS_CDFS':['3','4','13','14','23','24'],
    'MUSE':['3','2'],
    'SAGA':['4'],
    'SNLS_AAOmega':['3','4','6'],
    'VIPERS':['23.4', ' 2.4', ' 4.4', ' 3.5', ' 4.5', ' 2.2', ' 3.2', ' 4.2',
       ' 2.5', ' 9.5', ' 3.4', '19.5', '12.2', ' 9.4', ' 9.2', '13.2',
       '22.5', '24.2', '24.4', '14.2', '12.4', '24.5', '12.5', '22.2',
       '29.2', '23.5', '29.1', '22.1', '19.2', '13.5', '22.4', '29.5',
       '14.4', '23.2', '13.4', '14.5', '19.4', '23.1', '29.4', ' 2.1',
       '24.1', ' 4.1', ' 3.1', '219.', '13.1', '14.1', ' 9.1', '19.1',
       '12.1'],
    'DEEP2_DR4':['-1','3','4'],
    'VUDS_COSMOS':['3','4','13','14','23','24','43','44'],
    'VUDS_ECDFS':['3','4','13','14','23','24','43','44'],
    }
    grc = Table.read(os.path.join(s.cat_dir,'ozdes_grc.fits'))
    grc['ID'] = grc['ID'].astype(str)
    grc['flag'] = grc['flag'].astype(str)
    grc['source'] = grc['source'].astype(str)
    grc['comments'] = grc['comments'].astype(str)
    grc = grc.to_pandas()
    grc['flag'] = grc['flag'].str.strip(' ')
    good_redshifts = pd.DataFrame()
    for survey,flags in survey_flags.items():
        if flags !=['Any']:
            for flag in flags:
                good_redshifts = good_redshifts.append(grc[(grc['source']==survey)&(grc['flag']==flag)])
        else:
            good_redshifts = good_redshifts.append(grc[grc['source']==survey])

    z_gals = SkyCoord(ra=good_redshifts['RA'].values*u.degree,dec = good_redshifts['DEC'].values*u.degree)
    return good_redshifts,z_gals


def match_gals(catcoord,galscoord,cat,gals,dist_thresh):
    ordered_surveys = [
    'PRIMUS',
    'NED',
    'UDS_FORS2',
    'UDS_VIMOS',
    'ZFIRE_UDS',
    'ACES',
    'SDSS',
    '6dF',
    'ATLAS',
    '2dFGRS',
    'GAMA',
    'SNLS_FORS           ',
    'CDB',
    'VVDS_DEEP',
    'VVDS_CDFS',
    'MUSE',
    'SAGA',
    'DEEP2_DR4',
    'VUDS_COSMOS',
    'VUDS_ECDFS',
    'NOAO_0522',
    'NOAO_0334',
    'N17B331',
    'MOSDEF',
    'SpARCS',
    'VIPERS',
    'PanSTARRS_AAOmega   ',
    'PanSTARRS_MMT',
    'SNLS_AAOmega',
    'DES_AAOmega']

    inds,d2d,d3d = galscoord.match_to_catalog_sky(catcoord)
    init_matches = cat.iloc[inds]
    close_match_inds = d2d< dist_thresh*u.arcsec
    stack_gals_with_z = gals.iloc[close_match_inds]
    stack_gal_zs = init_matches[close_match_inds]

    logger.info('Matched %s galaxies with redshifts'%len(stack_gals_with_z))
    logger.debug('Going through them one-by-one to get the priority right')


    for c,i in enumerate(stack_gals_with_z.index):
        try:
            g = stack_gals_with_z.iloc[c:c+1]
        except:
            g = stack_gals_with_z.iloc[c:]
        gobj= SkyCoord(ra=g['X_WORLD'].values*u.deg,dec = g['Y_WORLD'].values*u.deg)
        idxc,idxgals,d2d,d3d = gobj.search_around_sky(catcoord,1*u.arcsec)
        hereitis=False
        grcres_full = cat.iloc[idxc]
        #print('So, around that galaxy, I found this: \n',grcres_full)
        '''if g['X_WORLD'].values<34.718 and g['X_WORLD'].values>34.716 and g['Y_WORLD'].values<-4.032 and g['Y_WORLD'].values>-4.034:
            hereitis=True
            logger.info(grcres)'''

        for survey in ordered_surveys:
            grcres = grcres_full[grcres_full['source']==survey]
            canskip = True
            for row in grcres[grcres['source']=='DES_AAOmega'].index:

                if grcres['ID'].loc[row][:10] =='SVA1_COADD':
                    ins = grcres[['z','z_Err','flag','source']].loc[row].values
                    stack_gals_with_z.loc[i,['SPECZ','SPECZ_ERR','SPECZ_FLAG','SPECZ_CATALOG']] = ins

                    if grcres['flag'].loc[row] in ['3','4']:
                        canskip=False
                    else:
                        canskip = True
                        #print('There is an ozdes source with name SVA1_COADD, but it has flag 1 or 2, so allowing further searching')
                else:
                    if canskip ==True:
                        if grcres['flag'].loc[row] in ['3','4']:
                            #print('I am going to insert an OzDES source that does not have name SVA1_COADD but does have a good flag')
                            ins = grcres[['z','z_Err','flag','source']].loc[row].values
                            stack_gals_with_z.loc[i,['SPECZ','SPECZ_ERR','SPECZ_FLAG','SPECZ_CATALOG']] = ins

            for row in grcres[grcres['source']!='DES_AAOmega'].index:
                bad_ozdes = 0
                for ozrow in grcres[grcres['source']=='DES_AAOmega'].index:
                    if grcres['flag'].loc[ozrow] in ['1','2']:
                        bad_ozdes =1
                if bad_ozdes ==1:
                    ins = grcres[['z','z_Err','flag','source']].loc[row].values
                    stack_gals_with_z.loc[i,['SPECZ','SPECZ_ERR','SPECZ_FLAG','SPECZ_CATALOG']] = ins
    gals.loc[stack_gals_with_z.index]=stack_gals_with_z

    return gals
if __name__ == "__main__":
    res = pd.read_csv(resname,index_col=0)
    noz_res = res[res['SPECZ']<0]
    gals_with_z,gals_with_z_coords = get_zs()
    nozrescoords = SkyCoord(ra=noz_res['RA'].values*u.deg,dec= noz_res['DEC'].values*u.deg)
    final_res = match_gals(gals_with_z_coords,nozrescoords,
                                    gals_with_z,noz_res,dist_thresh=1.5)
    final_res.to_csv(resname)
