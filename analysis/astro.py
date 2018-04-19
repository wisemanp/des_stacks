
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
from des_stacks.utils.stack_tooks import make_white
from des_stacks.utils.sex_tools import cap_sex

sns.set_palette('Dark2')
sns.set_color_codes(palette='colorblind')

def calib(s,chip,sexcat,phot_type='AUTO'):
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
    star_inds = ((y3a1['MAG_AUTO_%s'%Band]<18) & (y3a1['CLASS_STAR_%s'%Band]>0.3) | ((y3a1['SPREAD_MODEL_%s'%Band]\
     + 3*y3a1['SPREADERR_MODEL_%s'%Band])<0.003)) & (y3a1['MAG_AUTO_%s'%Band]<22)

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
    resfile.write(reshead)
    resfile.write(psfstring)
    savestring = os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(s.my,s.field,s.band,chip))
    s.logger.info("Saved result file to: %s"%savestring)
    return (kr_lim,kr_lim2,skylim,np.mean([kr_lim,kr_lim2,skylim]))

#####################################################################################################
def cap_phot(y,f,chip,wd = 'coadding'):
    '''get aperture photometry'''
    s.logger.info("Entered 'cap_phot.py' to do common aperture photometry on MY%s, %s, %s"%(y,f,chip))
    # first let's get to the right directory and set up a stack class object for each band_dir
    bands = ['g','r','i','z']
    sg,sr,si,sz = [stack.Stack(f, b, y, chip ,wd) for b in bands]

    # if there is no white image, make ones
    white_name = os.path.join(sg.out_dir,'MY%s'%y,f,'cap',str(chip),'MY%s_%s_ccd_%s_white_sci.fits'%(y,f,chip))
    if not os.path.isfile(white_name):
        white_name = make_white(sg,sr,si,sz,y,f,chip)

    sexcats =cap_sex(sg,sr,si,sz,y,f,chip,white_name)
    for s in [sg,sr,si,sz]:

        capcat = pd.DataFrame(fits.getdata(sexcats[s.band],ext=1))
        quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        zp = float(quals[0])
        av_fwhm = float(quals[2])


        capcat = cat.sort_values(by='X_WORLD')
        capcat['MAG_APER']=cat['MAG_APER']+zp
        # get rid of clearly wrong values
        truth =capcat['MAG_AUTO']<35
        capcat = capcat.iloc[truth.values]

    pixscale=0.27

    qual = os.path.join(ana_dir,'%s_ana.qual'%s.cutstring)
    thresh = 5
    skyflux = skynoise*np.sqrt(np.pi*(av_fwhm/pixscale)**2)
    skymag = 2.5*np.log10(thresh*skyflux)
    zmag = zp_psf
    skylim = zmag -skymag
    #####################################
    'THIS IS WHERE YOU GOT TO'

    resfile = open(os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(s.my,s.field,s.band,chip)),'w')
    psf['FWHM_WORLD'] = psf['FWHM_WORLD']*3600
    for i in range(len(psf['FWHM_WORLD'].values)):
        psf['FWHM_WORLD'].values[i] = float(psf['FWHM_WORLD'].values[i])
    radec=psf[['X_WORLD','Y_WORLD']].applymap("{0:7.5f}".format)
    rest = psf[['MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','FWHM_WORLD','ELONGATION']].applymap("{0:4.3f}".format)
    rest[['X_WORLD','Y_WORLD']]=radec[['X_WORLD','Y_WORLD']]
    rest['CLASS_STAR']=psf['CLASS_STAR']
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
    reshead +='# Zeropoint based on PSF photometry: %s \n'%zp_psf
    reshead +='# Limiting Kron magnitude based on matched objects: %.3f\n'% kr_lim
    reshead +='# Limiting magnitude based on PSF photometry: %.3f\n'% psf_lim
    reshead +='# %s sigma limiting magnitude based on matched objects: %.3f\n'%(limsig,psf_lim2)
    reshead +='# %s sigma limiting magnitude using zeropoint %.3f: %.3f\n' %(thresh,zmag,skylim)
    reshead +='# Columns:\n'
    reshead +='# Dec (J2000)\n'
    reshead +='# Kron Magnitude\n'
    reshead +='# Kron Magnitude error\n'
    reshead +='# PSF Magnitude\n'
    reshead +='# PSF Magnitude error\n'
    reshead +='# FWHM of the source (arcsec)\n'
    reshead +='# Elongation of source\n'
    resfile.write(reshead)
    resfile.write(psfstring)
    savestring = os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(s.my,s.field,s.band,chip))
    s.logger.info("Saved result file to: %s"%savestring)
    return (kr_lim,psf_lim,psf_lim2,skylim,np.mean([kr_lim,psf_lim,psf_lim2,skylim]))
