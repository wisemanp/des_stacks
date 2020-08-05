import os
import sys
import pandas as pd
import easyaccess as ea
from tqdm import tqdm

deep = pd.read_hdf(sys.argv[1],key='main')
new_sncand = pd.read_csv('/media/data3/wiseman/des/coadding/catalogs/all_snids.csv')

#new_snids = new_sncand[new_sncand['SNFAKE_ID']==0]['SNID'].values
conn = ea.connect('desoper')
cursor = conn.cursor()
try:
    start = sys.argv[2]
except:
    start = 0
for i in tqdm(range(int(start),len(deep))): #len(deep)

    query =("INSERT INTO SNGALS_DEEP "
      "( A_IMAGE, B_IMAGE, CCDNUM, CLASS_STAR_G, CLASS_STAR_I,"
       "CLASS_STAR_R, CLASS_STAR_Z, CXX_IMAGE, CXY_IMAGE, CYY_IMAGE,"
       "DLR, DLR_RANK, ELONGATION, FIELD, FLUXERR_APER_4_G,"
       "FLUXERR_APER_4_I, FLUXERR_APER_4_R, FLUXERR_APER_4_Z,"
       "FLUXERR_AUTO_G, FLUXERR_AUTO_I, FLUXERR_AUTO_R, FLUXERR_AUTO_Z,"
       "FLUX_APER_4_G, FLUX_APER_4_I, FLUX_APER_4_R, FLUX_APER_4_Z,"
       "FLUX_AUTO_G, FLUX_AUTO_I, FLUX_AUTO_R, FLUX_AUTO_Z,"
       "KRON_RADIUS, LIMFLUX_G, LIMFLUX_I, LIMFLUX_R, LIMFLUX_Z,"
       "LIMMAG_G, LIMMAG_I, LIMMAG_R, LIMMAG_Z, MAGERR_APER_4_G,"
       "MAGERR_APER_4_I, MAGERR_APER_4_R, MAGERR_APER_4_Z,"
       "MAGERR_AUTO_G, MAGERR_AUTO_I, MAGERR_AUTO_R, MAGERR_AUTO_Z,"
       "MAGERR_STATSYST_APER_4_G, MAGERR_STATSYST_APER_4_I,"
       "MAGERR_STATSYST_APER_4_R, MAGERR_STATSYST_APER_4_Z,"
       "MAGERR_STATSYST_AUTO_G, MAGERR_STATSYST_AUTO_I,"
       "MAGERR_STATSYST_AUTO_R, MAGERR_STATSYST_AUTO_Z,"
       "MAGERR_SYST_APER_4_G, MAGERR_SYST_APER_4_I, MAGERR_SYST_APER_4_R,"
       "MAGERR_SYST_APER_4_Z, MAGERR_SYST_AUTO_G, MAGERR_SYST_AUTO_I,"
       "MAGERR_SYST_AUTO_R, MAGERR_SYST_AUTO_Z, MAG_APER_4_G,"
       "MAG_APER_4_I, MAG_APER_4_R, MAG_APER_4_Z, MAG_AUTO_G,"
       "MAG_AUTO_I, MAG_AUTO_R, MAG_AUTO_Z, MAG_ZEROPOINT_ERR_G,"
       "MAG_ZEROPOINT_ERR_I, MAG_ZEROPOINT_ERR_R, MAG_ZEROPOINT_ERR_Z,"
       "MAG_ZEROPOINT_G, MAG_ZEROPOINT_I, MAG_ZEROPOINT_R,"
       "MAG_ZEROPOINT_Z, SEASON, PHOTOZ, PHOTOZ_ERR, TRANSIENT_NAME,"
       "THETA_IMAGE, X_IMAGE, RA, Y_IMAGE, DEC, SPECZ_FLAG,"
       "SPECZ_CATALOG, SPECZ, SPECZ_ERR, COADD_OBJECTS_ID, SNGALID,"
       "VERSION, SNID, GALFLAG, HOST, SEPARATION,"
       "SFR, SPECSFR, SPECSFR_ERRMINUS ) "
       "VALUES ("
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %i, %f, '%-20s',"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f,"
       "%f, %f, %f, %f, %f, %5.0g, %f, %f, '%-11s',"
       "%f, %f, %9.6f, %f, %9.6f, '%-8s', '%-20s', %f,"
       "%f, %11.0f, %7.0f, '%-40s', %9.0f, %4.0g, %4.0g, %f, "
       "%f, %f, %f)"%(
       deep['A_IMAGE'].iloc[i], deep['B_IMAGE'].iloc[i], deep['CCDNUM'].iloc[i], deep['CLASS_STAR_G'].iloc[i], deep['CLASS_STAR_I'].iloc[i],
       deep['CLASS_STAR_R'].iloc[i], deep['CLASS_STAR_Z'].iloc[i], deep['CXX_IMAGE'].iloc[i], deep['CXY_IMAGE'].iloc[i], deep['CYY_IMAGE'].iloc[i],
       deep['DLR'].iloc[i], deep['DLR_RANK'].iloc[i], deep['ELONGATION'].iloc[i], deep['FIELD'].iloc[i], deep['FLUXERR_APER_4_G'].iloc[i],
       deep['FLUXERR_APER_4_I'].iloc[i], deep['FLUXERR_APER_4_R'].iloc[i], deep['FLUXERR_APER_4_Z'].iloc[i],
       deep['FLUXERR_AUTO_G'].iloc[i], deep['FLUXERR_AUTO_I'].iloc[i], deep['FLUXERR_AUTO_R'].iloc[i], deep['FLUXERR_AUTO_Z'].iloc[i],
       deep['FLUX_APER_4_G'].iloc[i], deep['FLUX_APER_4_I'].iloc[i], deep['FLUX_APER_4_R'].iloc[i], deep['FLUX_APER_4_Z'].iloc[i],
       deep['FLUX_AUTO_G'].iloc[i], deep['FLUX_AUTO_I'].iloc[i], deep['FLUX_AUTO_R'].iloc[i], deep['FLUX_AUTO_Z'].iloc[i],
       deep['KRON_RADIUS'].iloc[i], deep['LIMFLUX_G'].iloc[i], deep['LIMFLUX_I'].iloc[i], deep['LIMFLUX_R'].iloc[i], deep['LIMFLUX_Z'].iloc[i],
       deep['LIMMAG_G'].iloc[i], deep['LIMMAG_I'].iloc[i], deep['LIMMAG_R'].iloc[i], deep['LIMMAG_Z'].iloc[i], deep['MAGERR_APER_4_G'].iloc[i],
       deep['MAGERR_APER_4_I'].iloc[i], deep['MAGERR_APER_4_R'].iloc[i], deep['MAGERR_APER_4_Z'].iloc[i],
       deep['MAGERR_AUTO_G'].iloc[i], deep['MAGERR_AUTO_I'].iloc[i], deep['MAGERR_AUTO_R'].iloc[i], deep['MAGERR_AUTO_Z'].iloc[i],
       deep['MAGERR_STATSYST_APER_4_G'].iloc[i], deep['MAGERR_STATSYST_APER_4_I'].iloc[i],
       deep['MAGERR_STATSYST_APER_4_R'].iloc[i], deep['MAGERR_STATSYST_APER_4_Z'].iloc[i],
       deep['MAGERR_STATSYST_AUTO_G'].iloc[i], deep['MAGERR_STATSYST_AUTO_I'].iloc[i],
       deep['MAGERR_STATSYST_AUTO_R'].iloc[i], deep['MAGERR_STATSYST_AUTO_Z'].iloc[i],
       deep['MAGERR_SYST_APER_4_G'].iloc[i], deep['MAGERR_SYST_APER_4_I'].iloc[i], deep['MAGERR_SYST_APER_4_R'].iloc[i],
       deep['MAGERR_SYST_APER_4_Z'].iloc[i], deep['MAGERR_SYST_AUTO_G'].iloc[i], deep['MAGERR_SYST_AUTO_I'].iloc[i],
       deep['MAGERR_SYST_AUTO_R'].iloc[i], deep['MAGERR_SYST_AUTO_Z'].iloc[i], deep['MAG_APER_4_G'].iloc[i],
       deep['MAG_APER_4_I'].iloc[i], deep['MAG_APER_4_R'].iloc[i], deep['MAG_APER_4_Z'].iloc[i], deep['MAG_AUTO_G'].iloc[i],
       deep['MAG_AUTO_I'].iloc[i], deep['MAG_AUTO_R'].iloc[i], deep['MAG_AUTO_Z'].iloc[i], deep['MAG_ZEROPOINT_ERR_G'].iloc[i],
       deep['MAG_ZEROPOINT_ERR_I'].iloc[i], deep['MAG_ZEROPOINT_ERR_R'].iloc[i], deep['MAG_ZEROPOINT_ERR_Z'].iloc[i],
       deep['MAG_ZEROPOINT_G'].iloc[i], deep['MAG_ZEROPOINT_I'].iloc[i], deep['MAG_ZEROPOINT_R'].iloc[i],
       deep['MAG_ZEROPOINT_Z'].iloc[i], deep['SEASON'].iloc[i], deep['PHOTOZ'].iloc[i], deep['PHOTOZ_ERR'].iloc[i], deep['TRANSIENT_NAME'].iloc[i],
       deep['THETA_IMAGE'].iloc[i], deep['X_IMAGE'].iloc[i], deep['RA'].iloc[i], deep['Y_IMAGE'].iloc[i], deep['DEC'].iloc[i], deep['SPECZ_FLAG'].iloc[i],
       deep['SPECZ_CATALOG'].iloc[i], deep['SPECZ'].iloc[i], deep['SPECZ_ERR'].iloc[i], deep['COADD_OBJECTS_ID'].iloc[i], deep['SNGALID'].iloc[i],
       deep['VERSION'].iloc[i], deep['SNID'].iloc[i], deep['GALFLAG'].iloc[i], deep['HOST'].iloc[i],deep['SEPARATION'].iloc[i],
       deep['EDGE_FLAG'].iloc[i], deep['Z_RANK'].iloc[i],deep['Z_USE'].iloc[i]))
    #print (int(deep['SNID'].iloc[i]))
    try:
        cursor.execute(query)
    except:
        pass
        #print ('Successfully pushed row %s of %s to SNGALS_DEEP'%(i,len(deep)))
print ('DONE!')
