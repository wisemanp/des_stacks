#!/usr/bin/env python
"""get_sn_info.py: Python script to get information on hosts for single SN or a list of SN names."""
"""Simply type get_sn_info -n DESYYFFabc, where YY = year, FF = field, abc = name."""

__author__      = "Phil Wiseman"
__email__       = "P.S.Wiseman@soton.ac.uk"
__version__     = "0.9"
__date__        = "25/09/2018"

# CHANGELOG:
#   V 0.9: 2018-09-25, PW: initial version
#   V 0.91: 2018-11-05, PW: changed name from get_sn_info to sn_host_phot


import numpy as np
import pandas as pd
import _pickle as cpickle
import os
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--sn_name',help='Full DES supernova ID, e.g. DES16C2nm (case sensitive)',default=None)
    parser.add_argument('-l','--namelist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-av','--avoid',help='Avoid these SN, if it is in the list given by l. e.g. [DES16C2nm]',default=None)
    parser.add_argument('-sf','--savename',help='Filename to save results to',default=None)
    parser.add_argument('-ow','--overwrite',help='When given, I will overwrite existing results for SN that already have results in the named results file',action = 'store_true')
    parser.add_argument('-th','--threshold',help='Distance threshold for host galaxy searching (arcsecs)',default=15)
    return parser.parse_args()

def get_DLR_ABT(RA_SN, DEC_SN, RA, DEC, A_IMAGE, B_IMAGE, THETA_IMAGE, angsep):
    # inputs are arrays
    rad  = np.pi/180                   # convert deg to rad
    pix_arcsec = 0.264                 # pixel scale (arcsec per pixel)
    pix2_arcsec2 = 0.264**2            # pix^2 to arcsec^2 conversion factor
    pix2_deg2 = pix2_arcsec2/(3600**2) # pix^2 to deg^2 conversion factor
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

class sn():
    def __init__(self,sn_name,args):
        self.sn_name = sn_name
        self.args = args

    def get_sn_dat(self):

        f=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
        chiplims = cpickle.load(f)

        sncand = Table.read('/media/data3/wiseman/des/coadding/catalogs/sn_cand.fits').to_pandas()
        gap = ' '
        ngaps = (11-len(self.sn_name))*gap
        dat = sncand[sncand['TRANSIENT_NAME']==self.sn_name+ngaps]

        try:
            ra,dec =dat[['RA','DEC']].iloc[0].values
        except IndexError:
            print ('WARNING: it looks like this transient does not exist in SNCAND. Programme will crash')
        y = dat['SEASON'].values[0]
        self.ra,self.dec,self.y = ra,dec,y
        #################
        obj_field = self.sn_name[5:7]
        self.field=obj_field
        the_field = chiplims[obj_field]
        for ccd in the_field.keys():
            if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
                if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                    self.chip = ccd
                    return (ra,dec,'SN-%s'%obj_field,y,ccd)

    def check_res(self):
        sn_cap_dir = '/media/data3/wiseman/des/coadding/5yr_stacks/CAP/%s'%self.sn_name
        if not os.path.isdir(sn_cap_dir):
            print ('No directory for Common Aperture Photometry on %s yet'%self.sn_name)
            return None
        else:
            all_trans = pd.read_csv('/media/data3/wiseman/des/coadding/results/all_transients.result',index_col=0)
            sn_res_row = all_trans[all_trans['SN_NAME']==self.sn_name]
            if len(sn_res_row['SN_NAME'])==0:
                return None
            else:
                print ('Found the results in the main transient results file, /media/data3/wiseman/des/coadding/results/all_transients.result')

                return (sn_res_row)

    def check_chip_cap(self):
        chip_cap_dir = '/media/data3/wiseman/des/coadding/5yr_stacks/MY%s/SN-%s/CAP/%s'%(self.y,self.field,self.chip)
        chip_cap_res_fn = os.path.join(chip_cap_dir,'%s_SN-%s_%s_obj_deep.cat'%(self.y,self.field,self.chip))
        chip_cap_res = pd.read_csv(chip_cap_res_fn,index_col=0)
        cap_res_coords = SkyCoord(ra=chip_cap_res['X_WORLD'].values*u.deg,dec = chip_cap_res['Y_WORLD'].values*u.deg)
        sn_coords = SkyCoord(ra=self.ra*u.deg,dec =self.dec*u.deg)
        d2d= sn_coords.separation(cap_res_coords)
        close_inds = d2d <float(self.args.threshold)*u.arcsec
        dists = d2d[close_inds]
        match = chip_cap_res.iloc[close_inds]
        angsep = np.array([float(d2d[close_inds][j].to_string(unit=u.arcsec,decimal=True)) for j in range(len(d2d[close_inds]))])
        match['ANGSEP'] = angsep
        dlr = get_DLR_ABT(self.ra,self.dec, match.X_WORLD, match.Y_WORLD, match['A_IMAGE'], match['B_IMAGE'],  match['THETA_IMAGE'], angsep)[0]
        match['DLR'] = np.array(dlr)
        print ('Went into the Common Aperture Photometry and found galaxies with \n the following angular separations and DLRs')
        print (match[['ANGSEP','DLR']])
        return match

def get_results(sn_name,args):
    s = sn(sn_name,args)
    s.get_sn_dat()
    is_res = s.check_res()
    if not isinstance(is_res,pd.DataFrame):
        res = s.check_chip_cap()
    else:
        res = is_res
    return res

def main(args):
    avoid_list = []
    if args.avoid:
        avoid_list = [i for i in args.avoid.split(',')]
    else:
       avoid_list = [None]
    res = pd.DataFrame()
    if args.sn_name:
        print('Have been given a name, looking for photometry for %s only'%args.sn_name)
        res = res.append(get_results(args.sn_name,args))

    else:

        sn_list = np.genfromtxt(args.namelist,dtype=str,delimiter='\n')
        print("Looking for photometry on following input list")
        print(sn_list)
        if not args.savename:

            try:
                done_sn = pd.read_csv('/media/data3/wiseman/des/coadding/results/all_transients.result',index_col=0)
            except:
                done_sn = pd.DataFrame(columns=['BAND', 'CLASS_STAR', 'ELONGATION', 'FWHM_WORLD', 'KRON_RADIUS', 'MAGERR_APER', 'MAGERR_AUTO', 'MAG_APER', 'MAG_AUTO', 'SN_NAME', 'X_WORLD', 'Y_WORLD','LIMMAG'])
        else:

            try:
                print('Reading in results file to find out which ones I still need to do')
                print(os.path.join('/media/data3/wiseman/des/coadding/results/',args.savename))
                done_sn = pd.read_csv(os.path.join('/media/data3/wiseman/des/coadding/results/',args.savename),index_col=0)

                print('Read in %s'%os.path.join('/media/data3/wiseman/des/coadding/results',args.savename))
            except:
                done_sn = pd.DataFrame(columns=['BAND', 'CLASS_STAR', 'ELONGATION', 'FWHM_WORLD', 'KRON_RADIUS', 'MAGERR_APER', 'MAGERR_AUTO', 'MAG_APER', 'MAG_AUTO', 'SN_NAME', 'X_WORLD', 'Y_WORLD','LIMMAG'])

        for sn_name in sn_list :
            print("Getting photometry on %s"%sn_name)


            if sn_name not in done_sn.SN_NAME.unique():


                if sn_name not in avoid_list:
                    res = res.append(get_results(sn_name,args))
            elif args.overwrite == True:
                if sn_name not in avoid_list:
                    res = res.append(get_results(sn_name,args))
                else:
                    print("Result for %s already in result file, and you told me not to overwrite it. Going to next one!"%sn_name)
    if args.savename:
        res.to_csv(args.savename)
        print ('Saved to %s'%args.savename)
    else:
        print ('Here is the result. You did not specifcy a savename, so this wont be saved')
        print (res)

if __name__=="__main__":
    args = parser()
    main(args)
