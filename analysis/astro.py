
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
from scipy.interpolate import UnivariateSpline as spln

def astrometry(stack,chip,sexcat,phot_type='AUTO'):
    '''Load in the existing DES and the newly SExtracted catalogs'''

    cmap = {'PSF':'red','AUTO':'green','cat':'blue'}
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

def init_phot(stack,chip,sexcat):
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
    q= open(os.path.join(stack.ana_dir,'%s_%s_ana.qual'%(zp_cut,psf_cut)),'r')
    q = q.read()
    quals = q.split('\n')[-2]
    quals = quals.split(' ')
    av_fwhm = quals[2]
    zp_kr = quals[0]
    zp_psf = quals[1]
    sexcat = Table.read(os.path.join(stack.ana_dir,'%s_%s_%s_%s.sexcat'%(stack.my,stack.field,stack.band,stack.chip)))
    sexcat = sexcat.to_pandas()
    cat['MAG_AUTO']=sexcat['MAG_AUTO']+zp_kr
    # get rid of clearly wrong values
    truth =cat['MAG_AUTO']<35
    cat = cat.iloc[truth.values]
    psftruth = sexcat['MAG_PSF']<98
    psf = sexcat.iloc[psftruth.values]
    psf['MAG_PSF']=psf['MAG_PSF']+zp_psf
    # make region files for ds9
    krreg = open(os.path.join(stack.ana_dir,'%s_%s_%s_%s_auto.reg'%(stack.my,stack.field,stack.band,stack.chip)),'w')
    psfreg = open(os.path.join(stack.ana_dir,'%s_%s_%s_%s_psf.reg'%(stack.my,stack.field,stack.band,stack.chip)),'w')
    print ('global color=red',file=psfreg)
    for i in range(len(cat['X_WORLD'].values)):
        print ('fk5; circle(%s,%s,10p) # text={%.2f +/- %.2f}'%(cat['X_WORLD'].iloc[i],cat['Y_WORLD'].iloc[i],cat['MAG_AUTO'].iloc[i],cat['MAGERR_AUTO'].iloc[i]),file=krreg)
    for i in range(len(psf['X_WORLD'].values)):
        print ('fk5; circle(%s,%s,5p) # text={%.2f +/- %.2f}'%(psf['X_WORLD'].iloc[i],psf['Y_WORLD'].iloc[i],psf['MAG_PSF'].iloc[i],psf['MAGERR_PSF'].iloc[i]),file=psfreg)
    krreg.close()
    psfreg.close()
    sns.set_palette('Dark2')
    sns.set_color_codes(palette='colorblind')
    f,ax=plt.subplots()
    alp= 0.75
    cat.hist(column='MAG_AUTO',bins=150,normed=True,ax=ax,alpha=alp+0.25,label='Kron Magnitudes',c='r')
    psf.hist(column='MAG_PSF',bins=150,normed=True,ax=ax,alpha=alp,label='PSF Magnitudes',c='g')
    x = np.linspace(0,30,10000)
    mu = 26.3
    sig = 2.1
    a = -3
    pdf = skewnorm.pdf(x,a,loc=mu,scale=sig)
    ax.set_xlabel('Mag')
    ax.set_ylabel('Frequency Density')
    ax.set_title('Magnitude Distribution in MY %s, %s, CCD %s' %(s.my,s.field,s.chips[0]))
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
    kr_lim = x2[np.argmax(y2)]
    psf_lim = x3[np.argmax(y3)]
    ax.vlines(kr_lim,0,1.1*np.max(y2),linestyle='--',label='Limiting Kron magnitude',c='r')
    ax.vlines(psf_lim,0,1.1*np.max(y3),linestyle='-.',label='Limiting Kron magnitude',c='g')
    ax.legend()
    f.savefig(os.path.join(stack.ana_dir,'%s_%s_%s_%s_hist.png'%(stack.my,stack.field,stack.band,stack.chip)))

    f2,ax2 = plt.subplots()
    cat.plot.scatter('MAG_AUTO','MAGERR_AUTO',s=5,ax=axscat,label='Kron Magnitudes',c='r')
    psf.plot.scatter('MAG_PSF','MAGERR_PSF',s=5,ax=axscat,label='PSF Magnitudes',c='g')
    ax2.set_xlabel('Magnitude')
    ax2.set_ylabel('Magnitude Error')
    errthresh = 2.5*np.log10(1+(1/limsig))
    ax2.hlines(errthresh,15,30,linestyle='--',c='#7570b3')
    ax2.set_xlim(16,30)
    ax2.set_ylim(-0.03,0.35)
    ax2.legend()
    f2.savefig(os.path.join(stack.ana_dir,'%s_%s_%s_%s_mag_vs_err.png'%(stack.my,stack.field,stack.band,stack.chip)))

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
    print (skynoise)
    h = fits.getheader(imgname)
    exptime= h['EXPTIME']
    pixscale=0.27
    qual = os.path.join(stack.ana_dir,'%s_%s_ana.qual'%(zp_cut,psf_cut))
    thresh = 5
    skyflux = skynoise*np.sqrt(np.pi*(av_fwhm/pixscale)**2)
    skymag = 2.5*np.log10(thresh*skyflux)
    zmag = zp_psf
    mlim = zmag -skymag


    resfile = open(os.path.join(stack.ana_dir,'%s_%s_%s_%s_init.result'%(stack.my,stack.field,stack.band,stack.chip)),'w')

    radec=psf[['X_WORLD','Y_WORLD']].applymap("{0:.5f}".format)
    rest = psf[['MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','FWHM_WORLD','ELONGATION']].applymap("{0:.3f}".format)

    for i in range(len(rest['FWHM_WORLD'].values)):
        rest['FWHM_WORLD'].values[i] = float(rest['FWHM_WORLD'].values[i])
    rest[['X_WORLD','Y_WORLD']]=radec[['X_WORLD','Y_WORLD']]
    cols = rest.columns.tolist()
    rearranged = cols[-2:]+cols[:-2]
    re = rest[rearranged]

    psfstring = re.to_csv(sep = ' ', index=False)
    reshead = '# Result file for a stack of Dark Energy Survey data taken by DECam\n'
    reshead +='# Field: %s\n'% stack.fields
    reshead +='# Minus year: %s\n'% stack.my
    reshead +='# Band: %s\n' % stack.band
    reshead +='# CCD Number: %s\n' % chip
    reshead +='# Total exposure time: %s\n' %exptime
    reshead +='# Limiting Kron magnitude: %.3f\n'% kr_lim
    reshead +='# Limiting PSF magnitude: %.3f\n'% psf_lim
    reshead +='%s sigma limiting magnitude based on matched objects: %.3f'%(limsig,psf_lim2)
    reshead +='%s sigma limiting magnitude using zeropoint %s: %.3f' %(thresh,zmag,mlim)
    resfile.write(reshead)
    resfile.write(re)
