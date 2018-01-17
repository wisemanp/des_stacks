
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import UnivariateSpline as spln

def astrometry(stack,chip,sexcat,phot_type='AUTO'):
    '''Load in the existing DES and the newly SExtracted catalogs'''
    stack.logger.info("Reading in catalog in order to do photometry")
    cmap = {'PSF':'red','AUTO':'green','cat':'blue'}
    old_cat = os.path.join(stack.cat_dir,'%s_All_filters_3.csv'%(stack.field[3]))
    old = pd.DataFrame.from_csv(old_cat)
    sexdat = fits.getdata(sexcat,ext=1)
    stack.logger.info("Successfully read in catalog: %s" %old_cat)
    stack.logger.info("Matching objects...")
    new =pd.DataFrame(sexdat)
    new_obj = SkyCoord(ra=new['X_WORLD']*u.degree,dec =new['Y_WORLD']*u.degree)
    old_obj = SkyCoord(ra=old['RA_%s'%stack.band]*u.degree,dec =old['DEC_%s'%stack.band]*u.degree)
    # match the catalogs
    idx, d2d, d3d = new_obj.match_to_catalog_sky(old_obj)
    match_ids = idx
    match_dists = d2d.arcsec
    stack.logger.info("Successfully matched %s objects!" %len(match_ids))
    # get old cat mags of the matched objects
    init_match_cat_mag =old['CLIPPED_MEAN_%s'%stack.band].iloc[match_ids]
    init_match_cat_magerr =old['CLIPPED_SIGMA_%s'%stack.band].iloc[match_ids]
    # get indices of the objects that are within a specified distance of their matches
    dist_cut =2.0
    good_inds = np.nonzero(match_dists < dist_cut)[0]

    good_new_ra = new['X_WORLD'].iloc[good_inds]
    good_new_dec = new['Y_WORLD'].iloc[good_inds]
    stack.logger.info("Using catalog magnitudes to calibrate photometry and get zeropoint")
    # find the new mags that correspond to that
    good_new_mag = new['MAG_%s'%phot_type].iloc[good_inds]
    good_new_magerr = new['MAGERR_%s'%phot_type].iloc[good_inds]
    # and the old ones
    good_cat_mag = init_match_cat_mag.iloc[good_inds]
    good_cat_magerr = init_match_cat_magerr.iloc[good_inds]
    # subtract to get the frame ZP
    zp = np.median(good_cat_mag.values - good_new_mag.values)
    psf = np.median(new['FWHM_WORLD']*3600)
    stack.logger.info("Successfully calbirated this DES stack of: %s, MY %s, %s band, CCD %s" %(stack.field,stack.my,stack.band,chip))
    return zp,psf

def init_phot(stack,chip,cat):
    stack.logger.info("Entered 'init_phot.py' to get Kron and PSF photometry and provide limiting magnitudes")

    try:
         ana_dir = stack.ana_dir
    except:
        ana_dir = os.path.join(stack.band_dir,chip,'ana')
    try:
        final = stack.final
    except AttributeError:
        final = True
    # first, get the raw magnitudes and add zero-points to make them proper magnitudes

    zp_cut,psf_cut = stack.zp_cut,stack.psf_cut
    if final ==True:
        imgname = os.path.join(stack.band_dir,'ccd_%s_%s_%.3f_%s.fits'%(chip,stack.band,zp_cut,psf_cut))
    else:
        imgname = stack.band_dir+'/ccd_%s_%s_%.3f_%s_temp.fits'%(chip,stack.band,zp_cut,psf_cut)
    cuts = imgname.split('_')
    q= open(os.path.join(ana_dir,'%s_%s_ana.qual'%(zp_cut,psf_cut)),'r')
    q = q.read()
    quals = q.split('\n')[-2]
    quals = quals.split(' ')
    av_fwhm = float(quals[2])
    zp_kr = float(quals[0])
    zp_psf = float(quals[1])
    cat = cat.sort_values(by='X_WORLD')
    cat['MAG_AUTO']=cat['MAG_AUTO']+zp_kr
    # get rid of clearly wrong values
    truth =cat['MAG_AUTO']<35
    cat = cat.iloc[truth.values]
    psftruth = cat['MAG_PSF']<98
    psf = copy.deepcopy(cat.iloc[psftruth.values])
    psf['MAG_PSF']=psf['MAG_PSF']+zp_psf
    # make region files for ds9
    krreg = open(os.path.join(ana_dir,'%s_%s_%s_%s_auto.reg'%(stack.my,stack.field,stack.band,chip)),'w')
    psfreg = open(os.path.join(ana_dir,'%s_%s_%s_%s_psf.reg'%(stack.my,stack.field,stack.band,chip)),'w')
    print ('global color=red',file=psfreg)
    for i in range(len(cat['X_WORLD'].values)):
        print ('fk5; circle(%s,%s,10p) # text={%.2f +/- %.2f}'%(cat['X_WORLD'].iloc[i],cat['Y_WORLD'].iloc[i],cat['MAG_AUTO'].iloc[i],cat['MAGERR_AUTO'].iloc[i]),file=krreg)
    for i in range(len(psf['X_WORLD'].values)):
        print ('fk5; circle(%s,%s,5p) # text={%.2f +/- %.2f}'%(psf['X_WORLD'].iloc[i],psf['Y_WORLD'].iloc[i],psf['MAG_PSF'].iloc[i],psf['MAGERR_PSF'].iloc[i]),file=psfreg)
    krreg.close()
    psfreg.close()
    stack.logger.info("Saved ds9 region files in /ana directory")
    sns.set_palette('Dark2')
    sns.set_color_codes(palette='colorblind')
    f,ax=plt.subplots()
    alp= 0.75
    cat.hist(column='MAG_AUTO',bins=150,normed=True,ax=ax,alpha=alp+0.25,label='Kron Magnitudes',color='r')
    psf.hist(column='MAG_PSF',bins=150,normed=True,ax=ax,alpha=alp,label='PSF Magnitudes',color='g')
    ax.set_xlabel('Mag')
    ax.set_ylabel('Frequency Density')
    ax.set_title('Magnitude Distribution in MY %s, %s, CCD %s, %s' %(stack.my,stack.field,chip,stack.band))
    #ax.set_yscale('log')
    hst,bin_edges = np.histogram(cat['MAG_AUTO'],bins=150,density=True)
    hstpsf,binspsf = np.histogram(psf['MAG_PSF'],bins=150,density=True)
    splkron = spln(bin_edges[1:],hst,s=0.02)
    splpsf = spln(binspsf[1:],hstpsf,s=0.02)
    x2 = np.linspace(bin_edges[0],bin_edges[-1],200)
    y2= splkron(x2)
    x3 = np.linspace(binspsf[0],binspsf[-1],200)
    y3 = splpsf(x3)
    ax.plot(x2,y2,c='r')
    ax.plot(x3,y3,c='g')
    ax.set_xlim(17,30)
    kr_lim = x2[np.argmax(y2)]
    psf_lim = x3[np.argmax(y3)]
    ax.vlines(kr_lim,0,1.1*np.max(y2),linestyle='--',label='Limiting Kron magnitude',color='r')
    ax.vlines(psf_lim,0,1.1*np.max(y3),linestyle='-.',label='Limiting PSF magnitude',color='g')
    ax.legend()
    f.savefig(os.path.join(ana_dir,'%s_%s_%s_%s_hist.jpg'%(stack.my,stack.field,stack.band,chip)))

    f2,ax2 = plt.subplots()
    cat.plot.scatter('MAG_AUTO','MAGERR_AUTO',s=5,ax=ax2,label='Kron Magnitudes',color='r')
    psf.plot.scatter('MAG_PSF','MAGERR_PSF',s=5,ax=ax2,label='PSF Magnitudes',color='g')
    ax2.set_xlabel('Magnitude')
    ax2.set_ylabel('Magnitude Error')
    limsig = 10
    errthresh = 2.5*np.log10(1+(1/limsig))
    ax2.hlines(errthresh,15,30,linestyle='--',color='#7570b3')
    ax2.set_xlim(17,30)
    ax2.set_ylim(-0.03,0.35)
    ax2.legend()
    f2.savefig(os.path.join(ana_dir,'%s_%s_%s_%s_mag_vs_err.jpg'%(stack.my,stack.field,stack.band,chip)))
    plt.close('all')
    b_hi = errthresh +(errthresh/500)
    b_lo = errthresh -(errthresh/500)

    c2 = cat[cat['MAGERR_AUTO']<b_hi]
    c2 = c2[c2['MAGERR_AUTO']>b_lo]
    kr_lim2 = c2['MAG_AUTO'].median()


    psf2 = psf[psf['MAGERR_PSF']<b_hi]
    psf2 = psf2[psf2['MAGERR_PSF']>b_lo]
    psf_lim2 = psf2['MAG_AUTO'].median()

    nclip=50

    out = iraf.imstat(imgname,fields='midpt,stddev',format=0,Stdout=1,nclip=nclip,usigma=2.8,lsigma=2.8)
    mean,skynoise = map(float, out[0].split())
    h = fits.getheader(imgname)
    exptime= h['EXPTIME']
    pixscale=0.27

    qual = os.path.join(ana_dir,'%s_%s_ana.qual'%(zp_cut,psf_cut))
    thresh = 5
    skyflux = skynoise*np.sqrt(np.pi*(av_fwhm/pixscale)**2)
    skymag = 2.5*np.log10(thresh*skyflux)
    zmag = zp_psf
    skylim = zmag -skymag
    stack.logger.info("Limiting Kron magnitude based on matched objects: %.3f\n"% kr_lim)
    stack.logger.info("Limiting magnitude based on PSF photometry: %.3f\n"% psf_lim)
    stack.logger.info("%s sigma limiting magnitude based on matched objects: %.3f\n"%(limsig,psf_lim2))
    stack.logger.info("%s sigma limiting magnitude using zeropoint %.3f: %.3f\n "%(thresh,zmag,skylim))

    resfile = open(os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(stack.my,stack.field,stack.band,chip)),'w')
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
    re.to_csv(os.path.join(stack.temp_dir,'temp_cat.csv'),index=False,sep=' ')
    stringthing = open(os.path.join(stack.temp_dir,'temp_cat.csv'),'r')
    psfstring = stringthing.read()
    stringthing.close()
    reshead = '# Result file for a stack of Dark Energy Survey data taken by DECam\n'
    reshead +='# Field: %s\n'% stack.field
    reshead +='# Minus year: %s\n'% stack.my
    reshead +='# Band: %s\n' % stack.band
    reshead +='# CCD Number: %s\n' % chip
    reshead +='# Total exposure time: %s s\n' %exptime
    reshead +='# Limiting Kron magnitude based on matched objects: %.3f\n'% kr_lim
    reshead +='# Limiting magnitude based on PSF photometry: %.3f\n'% psf_lim
    reshead +='# %s sigma limiting magnitude based on matched objects: %.3f\n'%(limsig,psf_lim2)
    reshead +='# %s sigma limiting magnitude using zeropoint %.3f: %.3f\n' %(thresh,zmag,skylim)
    reshead +='# Columns:\n'
    reshead +='# RA (J2000)\n'
    reshead +='# Dec (J2000)\n'
    reshead +='# Kron Magnitude\n'
    reshead +='# Kron Magnitude error\n'
    reshead +='# PSF Magnitude\n'
    reshead +='# PSF Magnitude error\n'
    reshead +='# FWHM of the source (arcsec)\n'
    reshead +='# Elongation of source\n'
    resfile.write(reshead)
    resfile.write(psfstring)
    savestring = os.path.join(ana_dir,'%s_%s_%s_%s_init.result'%(stack.my,stack.field,stack.band,chip))
    stack.logger.info("Saved result file to: %s"%savestring)
    return (kr_lim,psf_lim,psf_lim2,skylim)
