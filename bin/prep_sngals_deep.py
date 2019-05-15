import sys
import pandas as pd
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
readfile = sys.argv[1]
version = sys.argv[2]
res = pd.read_csv(readfile,
                  names=['ANGSEP','A_IMAGE','B_IMAGE','CCDNUM','CLASS_STAR_g',
                         'CLASS_STAR_i','CLASS_STAR_r','CLASS_STAR_z','CXX_IMAGE',
                         'CXY_IMAGE','CYY_IMAGE','DLR','DLR_RANK','EDGE_FLAG',
                         'ELONGATION','FIELD','FLUXERR_APER_g','FLUXERR_APER_i',
                         'FLUXERR_APER_r','FLUXERR_APER_z','FLUXERR_AUTO_g','FLUXERR_AUTO_i',
                         'FLUXERR_AUTO_r','FLUXERR_AUTO_z','FLUX_APER_g','FLUX_APER_i',
                         'FLUX_APER_r','FLUX_APER_z','FLUX_AUTO_g','FLUX_AUTO_i','FLUX_AUTO_r',
                         'FLUX_AUTO_z','FLUX_RADIUS_g','FLUX_RADIUS_i','FLUX_RADIUS_r',
                         'FLUX_RADIUS_z','FWHM_WORLD_g','FWHM_WORLD_i','FWHM_WORLD_r',
                         'FWHM_WORLD_z','KRON_RADIUS','LIMFLUX_g','LIMFLUX_i','LIMFLUX_r',
                         'LIMFLUX_z','LIMMAG_g','LIMMAG_i','LIMMAG_r','LIMMAG_z','MAGERR_APER_g',
                         'MAGERR_APER_i','MAGERR_APER_r','MAGERR_APER_z','MAGERR_AUTO_g',
                         'MAGERR_AUTO_i','MAGERR_AUTO_r','MAGERR_AUTO_z','MAGERR_STATSYST_APER_g',
                         'MAGERR_STATSYST_APER_i','MAGERR_STATSYST_APER_r','MAGERR_STATSYST_APER_z',
                         'MAGERR_STATSYST_AUTO_g','MAGERR_STATSYST_AUTO_i','MAGERR_STATSYST_AUTO_r',
                         'MAGERR_STATSYST_AUTO_z','MAGERR_SYST_APER_g','MAGERR_SYST_APER_i',
                         'MAGERR_SYST_APER_r','MAGERR_SYST_APER_z','MAGERR_SYST_AUTO_g',
                         'MAGERR_SYST_AUTO_i','MAGERR_SYST_AUTO_r','MAGERR_SYST_AUTO_z','MAG_APER_g',
                         'MAG_APER_i','MAG_APER_r','MAG_APER_z','MAG_AUTO_g','MAG_AUTO_i','MAG_AUTO_r',
                         'MAG_AUTO_z','MAG_ZEROPOINT_ERR_g','MAG_ZEROPOINT_ERR_i','MAG_ZEROPOINT_ERR_r',
                         'MAG_ZEROPOINT_ERR_z','MAG_ZEROPOINT_g','MAG_ZEROPOINT_i','MAG_ZEROPOINT_r',
                         'MAG_ZEROPOINT_z','MY','PHOTOZ','PHOTOZ_ERR','SNID','THETA_IMAGE',
                         'X_IMAGE','X_WORLD','Y_IMAGE','Y_WORLD','flag','source','z','z_Err'])
res= res.reset_index(drop=True)

#res = res.drop([0,300150])

res = res.rename(index=str,columns={'X_WORLD':'RA','Y_WORLD':'DEC',
                                    'MY':'SEASON',
                                    'ANGSEP':'SEPARATION',
                                   'z':'SPECZ',
                                   'z_Err':'SPECZ_ERR',
                                   'source':'SPECZ_CATALOG',
                                   'flag':'SPECZ_FLAG',
                                   'MAG_APER_g':'MAG_APER_4_G',
                                   'MAG_APER_r':'MAG_APER_4_R',
                                   'MAG_APER_i':'MAG_APER_4_I',
                                   'MAG_APER_z':'MAG_APER_4_Z',
                                   'MAGERR_APER_g':'MAGERR_APER_4_G',
                                   'MAGERR_APER_r':'MAGERR_APER_4_R',
                                   'MAGERR_APER_i':'MAGERR_APER_4_I',
                                   'MAGERR_APER_z':'MAGERR_APER_4_Z',
                                   'MAGERR_SYST_APER_g':'MAGERR_SYST_APER_4_G',
                                   'MAGERR_SYST_APER_r':'MAGERR_SYST_APER_4_R',
                                   'MAGERR_SYST_APER_i':'MAGERR_SYST_APER_4_I',
                                   'MAGERR_SYST_APER_z':'MAGERR_SYST_APER_4_Z',
                                   'MAGERR_STATSYST_APER_g':'MAGERR_STATSYST_APER_4_G',
                                   'MAGERR_STATSYST_APER_r':'MAGERR_STATSYST_APER_4_R',
                                   'MAGERR_STATSYST_APER_i':'MAGERR_STATSYST_APER_4_I',
                                   'MAGERR_STATSYST_APER_z':'MAGERR_STATSYST_APER_4_Z',
                                   'MAG_AUTO_g':'MAG_AUTO_G',
                                   'MAG_AUTO_r':'MAG_AUTO_R',
                                   'MAG_AUTO_i':'MAG_AUTO_I',
                                   'MAG_AUTO_z':'MAG_AUTO_Z',
                                   'MAGERR_AUTO_g':'MAGERR_AUTO_G',
                                   'MAGERR_AUTO_r':'MAGERR_AUTO_R',
                                   'MAGERR_AUTO_i':'MAGERR_AUTO_I',
                                   'MAGERR_AUTO_z':'MAGERR_AUTO_Z',
                                   'MAGERR_SYST_AUTO_g':'MAGERR_SYST_AUTO_G',
                                   'MAGERR_SYST_AUTO_r':'MAGERR_SYST_AUTO_R',
                                   'MAGERR_SYST_AUTO_i':'MAGERR_SYST_AUTO_I',
                                   'MAGERR_SYST_AUTO_z':'MAGERR_SYST_AUTO_Z',
                                   'MAGERR_STATSYST_AUTO_g':'MAGERR_STATSYST_AUTO_G',
                                   'MAGERR_STATSYST_AUTO_r':'MAGERR_STATSYST_AUTO_R',
                                   'MAGERR_STATSYST_AUTO_i':'MAGERR_STATSYST_AUTO_I',
                                   'MAGERR_SSTATYST_AUTO_z':'MAGERR_STATSYST_AUTO_Z',
                                   'FLUX_AUTO_g':'FLUX_AUTO_G',
                                    'FLUX_AUTO_r':'FLUX_AUTO_R',
                                    'FLUX_AUTO_i':'FLUX_AUTO_I',
                                    'FLUX_AUTO_z':'FLUX_AUTO_Z',
                                   'FLUXERR_AUTO_g':'FLUXERR_AUTO_G',
                                    'FLUXERR_AUTO_r':'FLUXERR_AUTO_R',
                                    'FLUXERR_AUTO_i':'FLUXERR_AUTO_I',
                                    'FLUXERR_AUTO_z':'FLUXERR_AUTO_Z',
                                   'FLUX_APER_g':'FLUX_APER_4_G',
                                    'FLUX_APER_r':'FLUX_APER_4_R',
                                    'FLUX_APER_i':'FLUX_APER_4_I',
                                    'FLUX_APER_z':'FLUX_APER_4_Z',
                                    'FLUXERR_APER_g':'FLUXERR_APER_4_G',
                                    'FLUXERR_APER_r':'FLUXERR_APER_4_R',
                                    'FLUXERR_APER_i':'FLUXERR_APER_4_I',
                                    'FLUXERR_APER_z':'FLUXERR_APER_4_Z',
                                   'CLASS_STAR_g':'CLASS_STAR_G',
                                   'CLASS_STAR_r':'CLASS_STAR_R',
                                   'CLASS_STAR_i':'CLASS_STAR_I',
                                   'CLASS_STAR_z':'CLASS_STAR_Z',
                                   'MAG_ZEROPOINT_g':'MAG_ZEROPOINT_G',
                                   'MAG_ZEROPOINT_r':'MAG_ZEROPOINT_R',
                                   'MAG_ZEROPOINT_i':'MAG_ZEROPOINT_I',
                                   'MAG_ZEROPOINT_z':'MAG_ZEROPOINT_Z',
                                   'LIMMAG_g':'LIMMAG_G',
                                   'LIMMAG_r':'LIMMAG_R',
                                   'LIMMAG_i':'LIMMAG_I',
                                   'LIMMAG_z':'LIMMAG_Z',
                                   'LIMFLUX_g':'LIMFLUX_G',
                                   'LIMFLUX_r':'LIMFLUX_R',
                                   'LIMFLUX_i':'LIMFLUX_I',
                                   'LIMFLUX_z':'LIMFLUX_Z',
                                   'MAG_ZEROPOINT_ERR_g':'MAG_ZEROPOINT_ERR_G',
                                   'MAG_ZEROPOINT_ERR_r':'MAG_ZEROPOINT_ERR_R',
                                   'MAG_ZEROPOINT_ERR_i':'MAG_ZEROPOINT_ERR_I',
                                   'MAG_ZEROPOINT_ERR_z':'MAG_ZEROPOINT_ERR_Z'})

res = res.drop([
                'FLUX_RADIUS_g','FLUX_RADIUS_r','FLUX_RADIUS_i','FLUX_RADIUS_z',
               'FWHM_WORLD_g','FWHM_WORLD_r','FWHM_WORLD_i','FWHM_WORLD_z','EDGE_FLAG'],axis=1)

sngals = pd.read_csv('/media/data3/wiseman/des/coadding/catalogs/sngals_db.csv',index_col=0)
sngals.replace(-9999,np.NaN,inplace=True)
deep = res
deep['SNGALID'] = np.arange(len(deep))*-1
deep.dropna(subset=['SEPARATION'],inplace=True)
deep['DEC'] = deep['DEC'].astype(float)
deep.index = deep.index.astype(int)
deep[['CLASS_STAR_G','CLASS_STAR_R','CLASS_STAR_I','CLASS_STAR_Z']].replace(-9999,np.NaN,inplace=True)
sngals = sngals.dropna(axis=0,subset=['ra','dec'])
deep = deep.dropna(axis=0,subset=['RA','DEC'])

sngals_coords = SkyCoord(ra=sngals['ra'].values*u.deg,dec=sngals['dec'].values*u.deg)
deep_coords = SkyCoord(ra=deep['RA'].values*u.deg,dec=deep['DEC'].values*u.deg)

idx,d2d,d3d = deep_coords.match_to_catalog_sky(sngals_coords)

init_good_sngals = sngals.iloc[idx]
good_match_inds = np.nonzero(d2d.arcsec <2)[0]
good_deep  = deep.iloc[good_match_inds]
good_sngals = init_good_sngals.iloc[good_match_inds]
deep['SNGALID'].loc[good_deep.index] = good_sngals['sngalid'].values
deep['SNGALID'].loc[~deep.index.isin(good_deep.index)] = np.arange(len(deep.loc[~deep.index.isin(good_deep.index)]))+sngals['sngalid'].max()

duplicate_snglids = deep[deep.duplicated(['SNGALID']).values]
dupe_deep_coords = SkyCoord(ra=duplicate_snglids['RA'].values*u.deg,dec=duplicate_snglids['DEC'].values*u.deg)

idx,d2d,d3d = dupe_deep_coords.match_to_catalog_sky(sngals_coords,2)

init_good_sngals = sngals.iloc[idx]
good_match_inds = np.nonzero(d2d.arcsec <2)[0]
good_dupe_deep  = duplicate_snglids.iloc[good_match_inds]
good_sngals2 = init_good_sngals.iloc[good_match_inds]
deep['SNGALID'].loc[good_dupe_deep.index] = good_sngals2['sngalid'].values
deep['SNGALID'].loc[~duplicate_snglids.index.isin(good_dupe_deep.index)] = np.arange(len(duplicate_snglids.loc[~duplicate_snglids.index.isin(good_dupe_deep.index)]))+sngals['sngalid'].max()


deep['COADD_OBJECTS_ID']=np.NaN
deep['COADD_OBJECTS_ID'].loc[good_deep.index] = good_sngals['coadd_objects_id'].values
sncand = pd.read_csv('/media/data3/wiseman/des/coadding/catalogs/sncand_db.csv',index_col=0)
ids = sncand[sncand['snfake_id']==0][['transient_name','snid']]
ids=ids.rename(str,columns={'transient_name':'TRANSIENT_NAME','snid':'SNID'})
#ids = ids[ids['TRANSIENT_NAME']!='-9999']
deep = deep.merge(ids,on='SNID',how='inner')
deep['VERSION']=''
deep['VERSION']=str(version)
deep['GALFLAG']=1
deep['HOST']= 0
deep['HOST'].loc[deep[(deep['GALFLAG']==1) &(deep['DLR_RANK']==1)].index]=1
deep = deep.replace(np.NaN, -9.998)
deep.to_csv('/media/data3/wiseman/des/coadding/results/sngals_deep_v%s.csv'%version,index=False)
