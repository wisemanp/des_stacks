# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from pyraf import iraf
import datetime
import os
import logging
import time
import seaborn as sns
import matplotlib
import glob
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import UnivariateSpline as spln

from des_stacks import des_stack as stack
from des_stacks.utils.stack_tools import make_cap_stamps, resample_chip_for_cap, get_chip_vals, get_cuts
from des_stacks.utils.sex_tools import cap_sex_sn, cap_sex_chip, get_sn_dat

sns.set_palette('Dark2')
sns.set_color_codes(palette='colorblind')

def init_calib(s,chip,sexcat,phot_type='AUTO'):
    '''Load in the existing DES and the newly SExtracted catalogs'''
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    '''if zp_cut>0:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:'''
    logger.setLevel(logging.INFO)
    ch.setLevel(logging.INFO)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Reading in catalog in order to do photometry")
    cmap = {'PSF':'red','AUTO':'green','cat':'blue','APER':'purple'}
    y3a1_fn = os.path.join(s.cat_dir,'y3a1_%s_%s.csv'%(s.field[3],s.band))
    y3a1 = pd.DataFrame.from_csv(y3a1_fn)
    sexdat = fits.getdata(sexcat,ext=1)
    logger.info("Successfully read in catalog: %s" %y3a1_fn)
    Band = s.band.capitalize()
    star_inds = ((y3a1['MAG_AUTO_%s'%Band]<18) & (y3a1['CLASS_STAR_%s'%Band]>0.3) |((y3a1['SPREAD_MODEL_%s'%Band] + 3*y3a1['SPREADERR_MODEL_%s'%Band])<0.003)) & (y3a1['MAG_AUTO_%s'%Band]<22)

    y3a1_stars = y3a1[star_inds]

    logger.info("Matching objects...")
    new =pd.DataFrame(sexdat)
    new_obj = SkyCoord(ra=new['X_WORLD']*u.degree,dec =new['Y_WORLD']*u.degree)
    old_obj = SkyCoord(ra=y3a1_stars['RA']*u.degree,dec =y3a1_stars['DEC']*u.degree)
    # match the catalogs
    idx, d2d, d3d = new_obj.match_to_catalog_sky(old_obj)
    match_ids = idx
    match_dists = d2d.arcsec
    logger.info("Successfully matched %s objects!" %len(match_ids))
    # get old cat mags of the matched objects
    init_match_cat_mag =y3a1_stars['MAG_AUTO_%s'%Band].iloc[match_ids]
    init_match_cat_magerr =y3a1_stars['MAGERR_AUTO_%s'%Band].iloc[match_ids]
    # get indices of the objects that are within a specified distance of their matches
    dist_cut =2.0
    good_inds = np.nonzero(match_dists < dist_cut)[0]

    good_new_ra = new['X_WORLD'].iloc[good_inds]
    good_new_dec = new['Y_WORLD'].iloc[good_inds]
    logger.info("Using catalog magnitudes to calibrate photometry and get zeropoint")
    # find the new mags that correspond to that
    good_new_mag = new['MAG_AUTO'].iloc[good_inds]
    good_new_magerr = new['MAGERR_AUTO'].iloc[good_inds]
    # and the old ones
    good_cat_mag = init_match_cat_mag.iloc[good_inds]
    good_cat_magerr = init_match_cat_magerr.iloc[good_inds]
    # subtract to get the frame ZP
    zp = np.median(good_cat_mag.values - good_new_mag.values)
    psf = np.median(new['FWHM_WORLD']*3600)
    logger.info("Successfully calbirated this DES stack of: %s, MY %s, %s band, CCD %s" %(s.field,s.my,s.band,chip))
    return zp,psf

def init_phot(s,chip,cat,pl='n'):
    s.logger.info("Entered 'init_phot.py' to get Kron and PSF photometry and provide limiting magnitudes")

    ana_dir = os.path.join(s.band_dir,chip,'ana')
    try:
        final = s.final
    except AttributeError:
        final = True
    # first, get the raw magnitudes and add zero-points to make them proper magnitudes

    if not s.cuts:

        if final ==True:
            imgname = os.path.join(s.band_dir,'ccd_%s_sci.fits'%chip)
        else:
            imgname = s.band_dir+'/ccd_%s_temp.fits'%chip


    else:
        if final ==True:
            imgname = os.path.join(s.band_dir,'ccd_%s_%s_%s_sci.fits'%(chip,s.band,s.cutstring))
        else:
            imgname = s.band_dir+'/ccd_%s_%s_%s_temp.fits'%(chip,s.band,s.cutstring)
    cuts = imgname.split('_')
    quals= np.loadtxt(os.path.join(ana_dir,'%s_ana.qual'%s.cutstring))
    zp = float(quals[0])
    av_fwhm = float(quals[1])
    cat = cat.sort_values(by='X_WORLD')
    cat['MAG_AUTO']=cat['MAG_AUTO']+zp
    try:
        cat['MAG_APER']=cat['MAG_APER']+zp
    except:
        s.logger.info('Aperture photometry appears not to have been done yet; consider doing it')
    # get rid of clearly wrong values
    truth =cat['MAG_AUTO']<35
    cat = cat.iloc[truth.values]

    # make region files for ds9
    krreg = open(os.path.join(ana_dir,'%s_%s_%s_%s_auto.reg'%(s.my,s.field,s.band,chip)),'w')

    for i in range(len(cat['X_WORLD'].values)):
        print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f}'%(cat['X_WORLD'].iloc[i],cat['Y_WORLD'].iloc[i],cat['MAG_AUTO'].iloc[i],cat['MAGERR_AUTO'].iloc[i]),file=krreg)
    krreg.close()
    s.logger.info("Saved ds9 region files in /ana directory")
    sns.set_palette('Dark2')
    sns.set_color_codes(palette='colorblind')
    if pl == 'y':
        f,ax=plt.subplots()
        alp= 0.75
        cat.hist(column='MAG_AUTO',bins=150,normed=True,ax=ax,alpha=alp+0.25,label='Kron Magnitudes',color='r')

        ax.set_xlabel('Mag')
        ax.set_ylabel('Frequency Density')
        ax.set_title('Magnitude Distribution in MY %s, %s, CCD %s, %s' %(s.my,s.field,chip,s.band))
    #ax.set_yscale('log')
    hst,bin_edges = np.histogram(cat['MAG_AUTO'],bins=150,density=True)

    splkron = spln(bin_edges[1:],hst,s=0.02)

    x2 = np.linspace(bin_edges[0],bin_edges[-1],200)
    y2= splkron(x2)

    kr_lim = x2[np.argmax(y2)]

    limsig = 10
    errthresh = 2.5*np.log10(1+(1/limsig))
    if pl == 'y':
        ax.plot(x2,y2,c='r')
        ax.set_xlim(17,30)
        ax.vlines(kr_lim,0,1.1*np.max(y2),linestyle='--',label='Limiting Kron magnitude',color='r')
        ax.legend()
        f.savefig(os.path.join(ana_dir,'%s_%s_%s_%s_hist.jpg'%(s.my,s.field,s.band,chip)))
        f2,ax2 = plt.subplots()
        cat.plot.scatter('MAG_AUTO','MAGERR_AUTO',s=5,ax=ax2,label='Kron Magnitudes',color='r')
        ax2.set_xlabel('Magnitude')
        ax2.set_ylabel('Magnitude Error')
        ax2.hlines(errthresh,15,30,linestyle='--',color='#7570b3')
        ax2.set_xlim(17,30)
        ax2.set_ylim(-0.03,0.35)
        ax2.legend()
        f2.savefig(os.path.join(ana_dir,'%s_%s_%s_%s_mag_vs_err.jpg'%(s.my,s.field,s.band,chip)))
        plt.close('all')
    b_hi = errthresh +(errthresh/500)
    b_lo = errthresh -(errthresh/500)
    c2 = cat[cat['MAGERR_AUTO']<b_hi]
    c2 = c2[c2['MAGERR_AUTO']>b_lo]
    kr_lim2 = c2['MAG_AUTO'].median()

    nclip=50
    s.logger.info("Running iraf.imstat on %s in order to get sky noise" %imgname)
    out = iraf.imstat(imgname,fields='midpt,stddev',format=0,Stdout=1,nclip=nclip,usigma=2.8,lsigma=2.8)
    try:
        mean,skynoise = map(float, out[0].split())
    except ValueError:
        s.logger.error("iraf.imstat failed on %s with following output: %s" %(imgname,out))
        mean,skynoise = None,None
    h = fits.getheader(imgname)
    exptime= h['EXPTIME']
    pixscale=0.27

    thresh = 5
    skyflux = skynoise*np.sqrt(np.pi*(av_fwhm/pixscale)**2)
    skymag = 2.5*np.log10(thresh*skyflux)
    skylim = zp -skymag
    s.logger.info("Limiting Kron magnitude based on matched objects: %.3f\n"% kr_lim)
    s.logger.info("%s sigma limiting magnitude based on matched objects: %.3f\n"%(limsig,kr_lim2))
    s.logger.info("%s sigma limiting magnitude using zeropoint %.3f: %.3f\n "%(thresh,zp,skylim))

    resfile = open(os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(s.my,s.field,s.band,chip)),'w')
    cat['FWHM_WORLD'] = cat['FWHM_WORLD']*3600
    for i in range(len(cat['FWHM_WORLD'].values)):
        cat['FWHM_WORLD'].values[i] = float(cat['FWHM_WORLD'].values[i])
    radec=cat[['X_WORLD','Y_WORLD']].applymap("{0:7.5f}".format)
    try:
        rest = cat[['MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','MAG_APER','MAGERR_APER','FWHM_WORLD','ELONGATION']].applymap("{0:4.3f}".format)
    except:
        rest = cat[['MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','FWHM_WORLD','ELONGATION']].applymap("{0:4.3f}".format)

    rest[['X_WORLD','Y_WORLD']]=radec[['X_WORLD','Y_WORLD']]
    rest['CLASS_STAR']=cat['CLASS_STAR']
    rest['FLUX_RADIUS']=cat['FLUX_RADIUS']
    cols = rest.columns.tolist()
    rearranged = cols[-2:]+cols[:-2]
    re = rest[rearranged]
    re.to_csv(os.path.join(s.temp_dir,'temp_cat.csv'),index=False,sep=' ')
    stringthing = open(os.path.join(s.temp_dir,'temp_cat.csv'),'r')
    psfstring = stringthing.read()
    stringthing.close()
    reshead = '# Result file for a stack of Dark Energy Survey data taken by DECam\n'
    reshead +='# Field: %s\n'% s.field
    reshead +='# Minus year: %s\n'% s.my
    reshead +='# Band: %s\n' % s.band
    reshead +='# CCD Number: %s\n' % chip
    reshead +='# Total exposure time: %s s\n' %exptime
    reshead +='# Zeropoint based on AUTO photometry: %s \n'%zp
    reshead +='# Limiting Kron magnitude based on matched objects: %.3f\n'% kr_lim
    reshead +='# %s sigma limiting magnitude based on matched objects: %.3f\n'%(limsig,kr_lim2)
    reshead +='# %s sigma limiting magnitude using zeropoint %.3f: %.3f\n' %(thresh,zp,skylim)
    reshead +='# Columns:\n'
    reshead +='# Dec (J2000)\n'
    reshead +='# Kron Magnitude\n'
    reshead +='# Kron Magnitude error\n'
    reshead +='# PSF Magnitude\n'
    reshead +='# PSF Magnitude error\n'
    reshead +='# FWHM of the source (arcsec)\n'
    reshead +='# Elongation of source\n'
    reshead +='# Flux Radius\n'
    resfile.write(reshead)
    resfile.write(psfstring)
    savestring = os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(s.my,s.field,s.band,chip))
    s.logger.info("Saved result file to: %s"%savestring)
    return (kr_lim,kr_lim2,skylim,np.mean([kr_lim,kr_lim2,skylim]))

#####################################################################################################
<<<<<<< HEAD
def cap_phot_sn(sn_name,wd = 'coadding',savename = 'all_sn_phot.csv',dist_thresh = 15):
=======
def cap_phot_sn(sn_name,wd = 'coadding',savename = 'all_sn_phot.csv',dist_thresh = 5,autocuts=False):
>>>>>>> dlr_phot
    '''get aperture photometry for a single sn host'''
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    '''if zp_cut>0:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:'''
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Entered 'cap_phot.py' to do common aperture photometry for the host of %s"%sn_name)
    logger.info("Will search a radius of %s arcseconds around the SN location"%dist_thresh)
    # first let's get to the right directory and set up a stack class object for each band_dir
    bands = ['g','r','i','z']

    ra,dec,f,y,chip = get_sn_dat(sn_name)
    if y == 5:
        y = 'none'
    logger.info("Found transient in the database SNCAND")
    logger.info("It's in %s, in Season %s, on chip %s, at coordinates RA = %s, Dec = %s"%(f,y,chip,ra,dec))
    # Make a Stack instance for each band
    logger.info("Setting up Stack instances for each band")
<<<<<<< HEAD
    cuts = [get_cuts(f,b) for b in bands]
=======
    if autocuts:
        cuts = [get_cuts(f,b) for b in bands]
    else:
        cuts = [{'teff': 0.15, 'psf':None},{'teff': 0.15,'psf':None},{'teff': 0.25,'psf':None},{'teff': 0.25,'psf':None}]
>>>>>>> dlr_phot
    sg,sr,si,sz = [stack.Stack(f, b, y, chip ,wd,cuts[counter]) for counter,b in enumerate(bands)]

    # if there is no white image, make ones
    det_name = os.path.join(sg.out_dir,'CAP',sn_name,'%s_white_stamp.fits'%(sn_name))
    if not os.path.isfile(det_name):
        logger.info("Couldn't find a detection image, so going to make 300x300 pix stamps of each band plus white")
        det_name = make_cap_stamps(sg,sr,si,sz,chip,sn_name,ra,dec,300,300)
    # check to see if sexcats exist already
    existing_sexcats = glob.glob(os.path.join(sg.out_dir,'CAP',sn_name,'*.sexcat'))

    sexcats = {}
    '''for b in bands:
        sexcat_fn = '%s_%s_cap_sci.sexcat'%(sn_name,b)
        sexcat_path = os.path.join(sg.out_dir,'CAP',sn_name)
        full_sexcat_fn = os.path.join(sexcat_path,sexcat_fn)
        if full_sexcat_fn in existing_sexcats:
            sexcats[b]=os.path.join(sg.out_dir,'CAP',sn_name,full_sexcat_fn)
    if len(existing_sexcats)!=4:'''

        # do common aperture photometry
    logger.info("Going to cap_sex to do CAP on each band")

    sexcats =cap_sex_sn(sg,sr,si,sz,chip,sn_name)
    # set up an empty results dataframe
<<<<<<< HEAD
    rescols = ['SN_NAME','X_WORLD', 'Y_WORLD',
=======
    rescols = ['SN_NAME','X_WORLD', 'Y_WORLD','X_IMAGE','Y_IMAGE',
               'A_IMAGE','B_IMAGE','THETA_IMAGE','CXX_IMAGE','CYY_IMAGE','CXY_IMAGE',
>>>>>>> dlr_phot
                           'MAG_AUTO_g', 'MAGERR_AUTO_g','MAG_APER_g', 'MAGERR_APER_g',
                           'MAG_AUTO_r', 'MAGERR_AUTO_r','MAG_APER_r', 'MAGERR_APER_r',
                           'MAG_AUTO_i', 'MAGERR_AUTO_i','MAG_APER_i', 'MAGERR_APER_i',
                           'MAG_AUTO_z', 'MAGERR_AUTO_z','MAG_APER_z', 'MAGERR_APER_z',
                           'FWHM_WORLD_g','FWHM_WORLD_r','FWHM_WORLD_i','FWHM_WORLD_z',
<<<<<<< HEAD
                           'ELONGATION_g','ELONGATION_r','ELONGATION_i','ELONGATION_z',
                           'KRON_RADIUS_g','KRON_RADIUS_r','KRON_RADIUS_i','KRON_RADIUS_z',
                           'CLASS_STAR_g','CLASS_STAR_r','CLASS_STAR_i','CLASS_STAR_z',
                           'LIMMAG_g','LIMMAG_r','LIMMAG_i','LIMMAG_z',
                           'FLUX_RADIUS_g','FLUX_RADIUS_r','FLUX_RADIUS_i','FLUX_RADIUS_z',
                           'DLR_g','DLR_r','DLR_i','DLR_z']
    res_df = pd.DataFrame(columns=rescols)

    for s in [sg,sr,si,sz]:
        # load in the photometry from sextractor
        qualfiles = glob.glob(os.path.join(s.band_dir,str(chip),'ana','*_ana.qual'))
        quals= np.loadtxt(qualfiles[0])
        zp = float(quals[0])
        av_fwhm = float(quals[2])
        capcat = Table.read(os.path.join(sg.out_dir,'CAP',sn_name,'%s_%s_cap_sci.sexcat'%(sn_name,s.band))).to_pandas()
=======
                           'ELONGATION',
                           'KRON_RADIUS',
                           'CLASS_STAR_g','CLASS_STAR_r','CLASS_STAR_i','CLASS_STAR_z',
                           'LIMMAG_g','LIMMAG_r','LIMMAG_i','LIMMAG_z',
                           'FLUX_RADIUS_g','FLUX_RADIUS_r','FLUX_RADIUS_i','FLUX_RADIUS_z',
                           'DLR',
                           'DLR_RANK']
    res_df = pd.DataFrame(columns=rescols)
    for s in [sg,sr,si,sz]:
        # load in the photometry from sextractor

        quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        zp = float(quals[0])
        av_fwhm = float(quals[2])
        capcat_fn = os.path.join(sg.out_dir,'CAP',sn_name,'%s_%s_cap_sci.sexcat'%(sn_name,s.band))
        logger.info('Reading in the catalog from: %s'%capcat_fn)
        capcat = Table.read(capcat_fn).to_pandas()
>>>>>>> dlr_phot
        capcat['MAG_APER']=capcat['MAG_APER']+zp
        capcat['MAG_AUTO']=capcat['MAG_AUTO']+zp

        sncoord = SkyCoord(ra = ra*u.deg,dec = dec*u.deg)
        catalog = SkyCoord(ra = capcat.X_WORLD*u.deg,dec = capcat.Y_WORLD*u.deg)
        d2d= sncoord.separation(catalog)
        close_inds = d2d <dist_thresh*u.arcsec
        dists = d2d[close_inds]
        match = capcat.iloc[close_inds]
<<<<<<< HEAD
        with open(os.path.join(s.band_dir,str(chip),'ana','%s_%s_%s_%s_init.result'%(y,f,s.band,chip)),'r') as resheader:
            header = [next(resheader) for x in range(8)]
        limmag = header[-1].split(' ')[-1].strip('\n')
        if len(match)==0:

            logger.info("Didn't detect a galaxy within %s arcsec of the SN; reporting limit of %s in %s band"%(dist_thresh,limmag,s.band))

            init_lim_array = np.array([sn_name,ra,dec,limmag,-1,limmag,-1,-1,-1,-1,-1,limmag,-1,-1])
            
            init_lim_cols = ['SN_NAME','X_WORLD', 'Y_WORLD',
                   'MAG_AUTO_%s'%s.band, 'MAGERR_AUTO_%s'%s.band,'MAG_APER_%s'%s.band, 'MAGERR_APER_%s'%s.band,
                   'FWHM_WORLD_%s'%s.band,
                   'ELONGATION_%s'%s.band,
                   'KRON_RADIUS_%s'%s.band,
                   'CLASS_STAR_%s'%s.band,
                   'LIMMAG_%s'%s.band,
                   'FLUX_RADIUS_%s'%s.band,
                   'DLR_%s'%s.band]
            if s.band =='g': 
                logger.info('Filling in g band with limits') 
                res_df = pd.DataFrame([init_lim_array],columns=init_lim_cols)
                
            else:
                lim_cols = ['MAG_AUTO_%s'%s.band, 'MAGERR_AUTO_%s'%s.band,'MAG_APER_%s'%s.band, 'MAGERR_APER_%s'%s.band,
                   'FWHM_WORLD_%s'%s.band,
                   'ELONGATION_%s'%s.band,
                   'KRON_RADIUS_%s'%s.band,
                   'CLASS_STAR_%s'%s.band,
                   'LIMMAG_%s'%s.band,
                   'FLUX_RADIUS_%s'%s.band,
                   'DLR_%s'%s.band]
                lim_array = np.array([limmag,-1,limmag,-1,-1,-1,-1,-1,limmag,-1,-1])
                for counter,c in enumerate(lim_cols):
                   res_df[c] = ''
                   res_df[c].iloc[0] = lim_array[counter]
        else:
            match.index = ['%s_%s'%(sn_name,i) for i in range(len(match.index))]
            cols = match.columns
            band_cols = {}
            for col in cols:
                band_cols[col]=col+'_%s'%s.band
            match.rename(index=str,columns=band_cols,inplace=True)
            dlr = dists*3600/(match['FLUX_RADIUS_%s'%s.band]*0.27*u.arcsec)
            match['DLR_%s'%s.band] = np.array(dlr)

            if s.band =='g':

                res_df = res_df.append(match)

            for c in match.columns:
                res_df[c] = ''
                res_df[c].loc[match.index] = match[c]

            rank = res_df['DLR_%s'%s.band].rank().astype(int)
            res_df['DLR_RANK_%s'%s.band]=rank

            logger.info(os.path.join(s.band_dir,str(chip),'ana','%s_%s_%s_%s_init.result'%(y,f,s.band,chip)))

            res_df['LIMMAG_%s'%s.band]= limmag
            logger.info('Limiting magnitude in %s band: %s'%(s.band,limmag))


            # make region files for ds9
            reg = open(os.path.join(s.out_dir,'CAP',sn_name,'%s_%s.reg'%(sn_name,s.band)),'w')

            for i in range(len(capcat['X_WORLD'].values)):
                print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f}'%(capcat['X_WORLD'].iloc[i],capcat['Y_WORLD'].iloc[i],capcat['MAG_AUTO'].iloc[i],capcat['MAGERR_AUTO'].iloc[i]),file=reg)
            print ('fk5; point %s %s # point=cross text={%s} color=red'%(ra,dec,sn_name),file=reg)
            reg.close()

=======
        angsep = np.array([float(d2d[close_inds][j].to_string(unit=u.arcsec,decimal=True)) for j in range(len(d2d[close_inds]))])
        with open(os.path.join(s.band_dir,str(chip),'ana',
            '%s_%s_%s_%s_init.result'%(y,f,s.band,chip)),'r') as resheader:
            header = [next(resheader) for x in range(8)]
        limmag = header[-1].split(' ')[-1].strip('\n')
        logger.info("Found %s galaxies within %s arcseconds in %s band"%(len(match),dist_thresh,s.band))
        if len(match)==0:

            logger.info("Didn't detect a galaxy within 2 arcsec of the SN; reporting limit of %s in %s band"%(limmag,s.band))

            init_lim_array = np.array([sn_name,ra,dec,limmag,-1,limmag,-1,-1,-1,-1,-1,limmag,-1,-1])
            init_lim_cols = ['SN_NAME','X_WORLD', 'Y_WORLD',
                   'MAG_AUTO_%s'%s.band, 'MAGERR_AUTO_%s'%s.band,'MAG_APER_%s'%s.band, 'MAGERR_APER_%s'%s.band,
                   'FWHM_WORLD_%s'%s.band,
                   'ELONGATION_%s'%s.band,
                   'KRON_RADIUS_%s'%s.band,
                   'CLASS_STAR_%s'%s.band,
                   'LIMMAG_%s'%s.band,
                   'FLUX_RADIUS_%s'%s.band,
                   'DLR_%s'%s.band]
            if s.band =='g':
                res_df=pd.DataFrame([init_lim_array],
                columns=init_lim_cols)
            else:

                lim_cols = ['MAG_AUTO_%s'%s.band, 'MAGERR_AUTO_%s'%s.band,'MAG_APER_%s'%s.band, 'MAGERR_APER_%s'%s.band,
                   'FWHM_WORLD_%s'%s.band,
                   'ELONGATION_%s'%s.band,
                   'KRON_RADIUS_%s'%s.band,
                   'CLASS_STAR_%s'%s.band,
                   'LIMMAG_%s'%s.band,
                   'FLUX_RADIUS_%s'%s.band,
                   'DLR_%s'%s.band]
                lim_array = np.array([limmag,-1,limmag,-1,-1,-1,-1,-1,limmag,-1,-1])
                for counter,c in enumerate(lim_cols):
                    res_df[c] = ''

                    res_df[c].iloc[0] = lim_array[counter]
        else:
            #match.index = ['%s_%s'%(sn_name,i) for i in range(len(match.index))]

            band_col_keys = ['MAG_AUTO', 'MAGERR_AUTO', 'MAG_APER', 'MAGERR_APER','CLASS_STAR']
            band_cols = {}
            new_band_cols = []
            for col in band_col_keys:
                band_cols[col]=col+'_%s'%s.band
                new_band_cols.append(col+'_%s'%s.band)
            if s.band =='g':
                match =match.rename(index=str,columns=band_cols)

                res_df = res_df.append(match)
                res_df['SN_NAME']=sn_name
                dlr = get_DLR_ABT(ra,dec, match.X_WORLD, match.Y_WORLD, match['A_IMAGE'], match['B_IMAGE'],  match['THETA_IMAGE'], angsep)[0]

                res_df['DLR'] = np.array(dlr)
                rank = res_df['DLR'].rank().astype(int)
                for counter, r in enumerate(res_df['DLR'].values):
                    if r >4:
                        rank[counter]*=-1
                res_df['DLR_RANK']=rank
            else:
                match =match.rename(index=str,columns=band_cols)

                match = match[new_band_cols]
                for c in match.columns:
                    res_df[c]= match[c]
            res_df['LIMMAG_%s'%s.band]= limmag
            res_df = res_df[res_df['DLR']<10]
            # make region files for ds9
            reg = open(os.path.join(s.out_dir,'CAP',sn_name,'%s_%s.reg'%(sn_name,s.band)),'w')

            for i in range(len(capcat['X_WORLD'].values)):
                print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f}'%(capcat['X_WORLD'].iloc[i],capcat['Y_WORLD'].iloc[i],capcat['MAG_AUTO'].iloc[i],capcat['MAGERR_AUTO'].iloc[i]),file=reg)
            print ('fk5; point %s %s # point=cross text={%s} color=red'%(ra,dec,sn_name),file=reg)
            reg.close()

>>>>>>> dlr_phot
    all_sn_fn = os.path.join(sg.res_dir,savename)
    if os.path.isfile(all_sn_fn):
        all_sn = pd.read_csv(all_sn_fn,index_col=0)
    else:
        all_sn = pd.DataFrame(columns = rescols)
<<<<<<< HEAD
    all_sn = all_sn.append(res_df)
=======
    all_sn = all_sn.append(res_df.reset_index(drop=True)).reset_index(drop=True)
>>>>>>> dlr_phot
    print ('Saving result to %s'%all_sn_fn)
    all_sn.to_csv(all_sn_fn)

    logger.info("Done doing CAP for %s"%sn_name)
    return res_df

def cap_phot_all(y,f,chip,wd='coadding'):
    '''get aperture photometry for an entire chip'''
    logger = logging.getLogger(__name__)
    logger.handlers =[]
    ch = logging.StreamHandler()
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Entered 'cap_phot_all' to do common aperture photometry for MY%s, %s, chip %s"%(y,f,chip))
    # first let's get to the right directory and set up a stack class object for each band_dir
    bands = ['g','r','i','z']
    survey_flags = {
    'DES_AAOmega':['3','4','6'],
    'ZFIRE_UDS':['3'],
    'NOAO_0522':['4','6'],
    'NOAO_0334':['4','6'],
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
    'VIPERS':['2','3','4','9','22','23','24','29','12','13','14','19','212','213','214','219'],
    'DEEP2_DR4':['-1','3','4'],
    'VUDS_COSMOS':['3','4','13','14','23','24','43','44'],
    'VUDS_ECDFS':['3','4','13','14','23','24','43','44'],
    }
    sg,sr,si,sz = [stack.Stack(f, b, y, chip ,wd) for b in bands]

    # if there is no white image, make ones
    det_name = os.path.join(sg.out_dir,'MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_white.fits'%(y,f,chip))
    if not os.path.isfile(det_name):
        logger.info("Couldn't find a detection image, so going to resample each band plus white to the same pixels")
        det_name = resample_chip_for_cap(sg,sr,si,sz,chip)
    # do common aperture photometry
    logger.info("Going to cap_sex to do CAP on each band")
    sexcats =cap_sex_chip(sg,sr,si,sz,chip)
    # set up an empty results dataframe
    res_df = pd.DataFrame(columns=['X_WORLD', 'Y_WORLD', 'BAND','MAG_AUTO', 'MAGERR_AUTO',
     'MAG_APER', 'MAGERR_APER', 'FWHM_WORLD', 'ELONGATION', 'CLASS_STAR'])
    this_chip_lims = get_chip_vals(sg.field,chip,vals = 'lims')

    # find the galaxies that OzDES has redshifts for
    grc = Table.read(os.path.join(sg.cat_dir,'ozdes_grc.fits')).to_pandas()
    good_redshifts = pd.DataFrame()
    for survey,flags in survey_flags.items():
        if flags !=['Any']:
            for flag in flags:
                good_redshifts = good_redshifts.append(grc[(grc['source']==survey)&(grc['flag']==flag)])
        else:
            good_redshifts = good_redshifts.append(grc[grc['source']==survey])
    good_in_chip_inds = (good_redshifts['RA']< this_chip_lims[0][0])&(good_redshifts['RA']> this_chip_lims[2][0]) & (good_redshifts['DEC']> this_chip_lims[0][1]) & (good_redshifts['DEC']< this_chip_lims[1][1])
    gals_with_z = good_redshifts[good_in_chip_inds]
    z_gals = SkyCoord(ra=gals_with_z['RA']*u.degree,dec = gals_with_z['DEC']*u.degree)
    phot_plus_spec = pd.DataFrame()
    for s in [sg,sr,si,sz]:
        # load in the photometry from sextractor
        capcat = Table.read(sexcats[s.band]).to_pandas()
        quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        zp = float(quals[0])
        av_fwhm = float(quals[2])
        capcat = capcat.sort_values(by='X_WORLD')
        logger.info("Calibrating in %s band using zeropoint from result file: %.3f"%(s.band,zp))
        capcat['MAG_APER']=capcat['MAG_APER']+zp
        capcat['MAG_AUTO']=capcat['MAG_AUTO']+zp
        # get rid of clearly wrong values
        truth =capcat['MAG_AUTO']<35
        capcat = capcat.iloc[truth.values]
        capcat.to_csv(os.path.join(s.out_dir,'MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_%s_phot_galcat.result'%(y,f,chip,s.band)))
        catobjs = SkyCoord(ra = capcat['X_WORLD']*u.degree,dec = capcat['Y_WORLD']*u.degree)
        # match the cap catalog with the ozdes one
        idx,d2d,d3d = catobjs.match_to_catalog_sky(z_gals)
        init_good_zgals = gals_with_z.iloc[idx]
        good_match_inds = np.nonzero(d2d.arcsec < 2.0)[0]
        logger.info("In %s band, I found %s galaxies which match to within 3 arcseconds"%(s.band,len(good_match_inds)))
        good_phot_gals = capcat.iloc[good_match_inds]
        init_good_spec_gals = init_good_zgals.iloc[good_match_inds]
        good_spec_gals = init_good_spec_gals[~init_good_spec_gals.index.duplicated(keep='first')]
        good_phot_gals = good_phot_gals[~init_good_spec_gals.index.duplicated(keep='first')]
        # make region files for ds9
        photreg = open(os.path.join(s.out_dir,'MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_%s_phot.reg'%(y,f,chip,s.band)),'w')
        specreg = open(os.path.join(s.out_dir,'MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_%s_spec.reg'%(y,f,chip,s.band)),'w')
        for i in range(len(capcat['X_WORLD'].values)):
            print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f}'%(capcat['X_WORLD'].iloc[i],capcat['Y_WORLD'].iloc[i],capcat['MAG_AUTO'].iloc[i],capcat['MAGERR_AUTO'].iloc[i]),file=photreg)

        photreg.close()
        print ('global color=red',file=specreg)
        for i in range(len(gals_with_z['RA'].values)):
            print ('fk5; circle(%s,%s,1") # text={z = %.3f}'%(gals_with_z['RA'].iloc[i],gals_with_z['DEC'].iloc[i],gals_with_z['z'].iloc[i]),file=specreg)
        specreg.close()
        #write the phot and spec properties to a file
        if s.band == 'g' and chip == 1:
            phot_plus_spec =good_spec_gals

        '''elif len(phot_plus_spec) < len(good_phot_gals):

            logger.debug('New good_phot_gals is longer than old phot_plus_spec')
            duplicated = good_phot_gals[good_phot_gals.index.duplicated()]
            logger.debug(duplicated)
            good_phot_gals = good_phot_gals[~good_phot_gals.index.duplicated(keep='first')]
            logger.debug('It was duplicated, now it is the right length?')
            logger.debug(len(phot_plus_spec))
            logger.debug(len(good_phot_gals))'''
        if len(phot_plus_spec) < len(good_phot_gals):

            previous_phot_plus_spec = phot_plus_spec
            phot_plus_spec = good_spec_gals

            #find the already done columns
            done_cols =previous_phot_plus_spec.columns
            done_bands = []
            for col in done_cols:
                if 'MAG' in col:
                    done_bands.append(col)

            with open(os.path.join(s.band_dir,str(chip),'ana','%s_%s_%s_%s_init.result'%(y,f,s.band,chip)),'r') as res:
                header = [next(res) for x in range(8)]
            limmag = header[-1].split(' ')[-1].strip('\n')
            for col in done_bands:


                phot_plus_spec[col]=limmag
                try:
                    phot_plus_spec.loc[previous_phot_plus_spec.index,col]=previous_phot_plus_spec[col]
                except KeyError:
                    overlap_inds = []
                    for ind in previous_phot_plus_spec.index:
                        if ind in phot_plus_spec.index:
                            overlap_inds.append(ind)
                    phot_plus_spec.loc[overlap_inds,col]=previous_phot_plus_spec.loc[overlap_inds,col]
        if len(phot_plus_spec) > len(good_phot_gals):
            try:
                phot_plus_spec.loc[good_spec_gals.index,'MAG_AUTO_%s'%s.band],phot_plus_spec.loc[good_spec_gals.index,'MAGERR_AUTO_%s'%s.band] = good_phot_gals['MAG_AUTO'].values,good_phot_gals['MAGERR_AUTO'].values
            except:

                for ind in good_spec_gals.index:

                    try:
                       phot_plus_spec.loc[ind,'MAG_AUTO_%s'%s.band],phot_plus_spec.loc[ind,'MAGERR_AUTO_%s'%s.band] = good_phot_gals['MAG_AUTO'].values,good_phot_gals['MAGERR_AUTO'].values
                    except:
                       logger.warning("Failed in the %s band at index %s"%(s.band,ind))

        else:

            phot_plus_spec['MAG_AUTO_%s'%s.band],phot_plus_spec['MAGERR_AUTO_%s'%s.band] = good_phot_gals['MAG_AUTO'].values,good_phot_gals['MAGERR_AUTO'].values


    phot_plus_spec.to_csv(os.path.join(sg.out_dir,'MY%s'%y,f,'CAP',str(chip),'spec_phot_galcat_%s_%s_%s.result'%(sg.my,sg.field,chip)))
    return phot_plus_spec



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
