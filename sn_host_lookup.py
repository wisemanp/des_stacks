# -*- coding: utf-8 -*-
'''Script to find the photometry for a certain supernova'''
# exectuable to run through the entire stack process
import numpy as np
import pandas as pd
import logging
import argparse
import easyaccess as ea
from time import gmtime, strftime
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.io.fits as fits
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import _pickle as cpickle

#Note: import this first else it crashes importing sub-modules
plot_locs={
'g':[0.2,0.53,0.39,0.42],
'r':[0.6,0.53,0.39,0.42],
'i':[0.2,0.1,0.39,0.42],
'z':[0.6,0.1,0.39,0.42]
}

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--sn_name',help='Full DES supernova ID, e.g. DES16C2nm (case sensitive)',default=None)
    parser.add_argument('-l','--namelist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-wd','--workdir',help='Path to directory to work in',required = False,default = '/media/data3/wiseman/des/coadding')
    parser.add_argument('-d','--distance',help='Maximum allowable match radius (arcsec) for galaxy association',required=False,default=5)
    parser.add_argument('-s','--stamp',help='Produce a stamp',action ='store_true')
    parser.add_argument('-ss','--size',help='Stamp size (arcsec)',default=30,type=float)
    parser.add_argument('-sh','--show',help='Show image now', action ='store_true')
    parser.add_argument('-p','-path',help='Full path to output image if you are making a stamp',default = 'sn_dir')
    return parser.parse_args()

def get_sn_dat(sn):
    f=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
    chiplims = cpickle.load(f)
    conn = ea.connect(section='desoper')
    q = 'select ra,dec,field,season,z_spec,z_spec_err from SNCAND \
    where transient_name = \'%s\''%sn
    dat = conn.query_to_pandas(q)
    ra,dec =dat[['RA','DEC']].iloc[0].values
    y = dat['SEASON'].values[0]
    obj_field = sn[5:7]
    the_field = chiplims[obj_field]
    for ccd in the_field.keys():
        if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
            if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                return (ra,dec,'SN-%s'%obj_field,y,ccd)

def main(args,logger):
    sn_name = args.sn_name
    match_dist = args.distance
    if sn_name:
        l = [sn_name]
    else:
        l = args.namelist
    for sn in l:
        logger.info("Locating result file for %s ..."%sn)
        logger.info("First querying the DES database to get the RA,Dec of %s."%sn)
        sn_ra,sn_dec,f,y,chip = get_sn_dat(sn)
        logger.info("%s has the following info:"%sn)
        logger.info("RA:     %s"%sn_ra)
        logger.info("Dec:    %s"%sn_dec)
        logger.info("Season: %s"%y)
        logger.info("Field:  %s"%f)
        logger.info("CCD:    %s"%chip)
        cap_chip_dir = os.path.join(args.workdir,'stacks','MY%s'%y,f,'CAP',str(chip))
        
        chip_res_fn = os.path.join(cap_chip_dir,'spec_phot_galcat_%s_%s_%s.result'%(y,f,chip))
        logger.info("Filename should look like %s."%chip_res_fn)
        bands = ['g','r','i','z']
        if os.path.isfile(chip_res_fn):
            logger.info("Found the result file for MY%s, %s, chip %s! Now trying to find the SN host"%(y,f,chip))

        else:
            logger.info("Result file doesn't exist yet; going to quickly do common aperture photometry on MY%s, %s, chip %s"%(y,f,chip))
            try:
                from des_stacks.analysis.astro import cap_phot_all
                sg,sr,si,sz = [stack.Stack(f, b, y, chip ,'coadding')for b in bands]
                cap_phot_all(sg,sr,si,sz,chip)
            except:
                if y=='none':
                    raise Exception("Unable to import the des_stacks module to do the photometry - Talk to Phil!")
                else:
                    logger.info("I'll try to do it in MYnone just to give you a look")
                    chip_res_fn = os.path.join(args.workdir,'stacks','MYnone',f,'CAP',str(chip),'spec_phot_galcat_%s_%s_%s.result'%('none',f,chip))
                    if os.path.isfile(chip_res_fn):
                        logger.info("Found a result file for MYnone, so I'll use that!")

                    else:
                        logger.info("%s, %s still not done even in MYnone, you should bug Phil to get it done."%(f,chip))
        # set up an empty results dataframe

        chip_res = pd.read_csv(chip_res_fn)
        res_objs = SkyCoord(ra=chip_res['RA']*u.degree,dec=chip_res['DEC']*u.degree)
        snloc = SkyCoord(ra=sn_ra*u.degree,dec = sn_dec*u.degree)
        idx,d2d,d3d = snloc.match_to_catalog_sky(res_objs)

        match = chip_res.iloc[int(idx)]
        logger.info(match)
        sn_dir = os.path.join(args.workdir,'stacks','CAP',sn)
        if not os.path.isdir(sn_dir):
            os.path.mkdir(sn_dir)

        if d2d.arcsec < (2*float(match_dist)):
            restype='spec'
            if d2d.arcsec < (float(match_dist)):
                logger.info("The SN lies within %s arcsec of a galaxy with a redshift!"%match_dist)
                logger.info("The galaxy details, including magnitudes from the deep DES SN stack, are: ")
                logger.info(match)
            else:
                logger.info("The SN lies within %s arcsec of a galaxy with a redshift, here are details, including magnitudes from the deep DES SN stack"%d2d.arcsec)
                
            for b in bands:
                reg = open(os.path.join(sn_dir,'%s_%s.reg'%(sn,b)),'w')
                for i in range(len(chip_res)):
                    print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f} color=green'%(chip_res['RA'].loc[ind],chip_res['DEC'].iloc[ind],chip_res['MAG_AUTO_%s'%b].iloc[ind],chip_res['MAGERR_AUTO_%s'%b].iloc[ind]),file=reg)
                for ind in match.index:
                    print ('fk5; circle(%s,%s,1.5") # text={%.2f +/- %.2f} color=blue width=3'%(match['RA'].loc[ind],match['DEC'].iloc[ind],match['MAG_AUTO_%s'%b].iloc[ind],match['MAGERR_AUTO_%s'%b].iloc[ind]),file=reg)
                print ('fk5; point %s %s # point=cross text={%s} color=red width=2'%(sn_ra,sn_dec,sn),file=reg)
                reg.close()
                logger.info("Saved region file to %s "%os.path.join(sn_dir,'%s_%s.reg'%(sn,b)))
        else:
            restype='phot'
            logger.info("No spectroscopically redshifted galaxies nearby; going to check the photometric catalog")
            for b in bands:
                phot_fn = os.path.join(args.workdir,'stacks','MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_%s_phot_galcat.result'%(y,f,chip,b))
                phot_res = pd.read_csv(phot_fn)
                res_objs = SkyCoord(ra=phot_res['X_WORLD']*u.degree,dec=phot_res['Y_WORLD']*u.degree)
                idx,d2d,d3d = snloc.match_to_catalog_sky(res_objs)
                logger.info("In %s band, I found the following source(s):"%b)
                match = phot_res.iloc[idx]
                logger.info(match)

                reg = open(os.path.join(sn_dir,'%s_%s.reg'%(sn,b)),'w')
                for i in range(len(phot_res)):
                    print ('fk5; circle(%s,%s,1") # text={%.2f +/- %.2f} color=green'%(phot_res['X_WORLD'].loc[ind],phot_res['Y_WORLD'].iloc[ind],phot_res['MAG_AUTO'].iloc[ind],phot_res['MAGERR_AUTO'].iloc[ind]),file=reg)
                for ind in match.index:
                    print ('fk5; circle(%s,%s,1.5") # text={%.2f +/- %.2f} color=blue width=3'%(match['X_WORLD'].loc[ind],match['X_WORLD'].iloc[ind],match['MAG_AUTO'].iloc[ind],match['MAGERR_AUTO'].iloc[ind]),file=reg)
                print ('fk5; point %s %s # point=cross text={%s} color=red width=2'%(sn_ra,sn_dec,sn),file=reg)
                reg.close()
                logger.info("Saved region file to %s "%os.path.join(sn_dir,'%s_%s.reg'%(sn,b)))
        if args.stamp:
            logger.info("You want a stamp of %s too. So I'm making one."%sn)
            logger.warning("You might need AplPy installed...")
            import aplpy
            fig,ax = plt.subplots()
            w = args.size/3600
            ax.set_xticks([])
            ax.set_yticks([])
            for loc in ['top','right','left','bottom']:
                ax.spines[loc].set_visible(False)
            ax.set_ylabel('Declination (J2000)',fontsize=12,labelpad = 30)
            hor_line = np.array([[sn_ra-0.00027,sn_ra+0.00027],[sn_dec,sn_dec]])
            ver_line = np.array([[sn_ra,sn_ra],[sn_dec-0.00027,sn_dec+0.00027]])
            for counter,band in enumerate(bands):
                img_fn = os.path.join(cap_chip_dir,'ccd_%s_%s_0.15_sci.resamp.fits'%(chip,band))

                if os.path.isfile(img) != 'Failed to load image':
                    if restype=='spec':
                        res = pd.read_csv(chip_res_fn)
                        res = res[res['RA']<sn_ra+(w/2)]
                        res = res[res['RA']>sn_ra-(w/2)]
                        res = res[res['DEC']<sn_dec+(w/2)]
                        res = res[res['DEC']>sn_dec-(w/2)]
                        ras,decs,mags,errs = res.RA.values,res.DEC.values,res['MAG_AUTO_%s'%band].values,res['MAGERR_AUTO_%s'%band].values

                    else:
                        print ("Loading the %s band region!" %band)
                        phot_fn = os.path.join(args.workdir,'stacks','MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_%s_phot_galcat.result'%(y,f,chip,b))
                        res = pd.read_csv(phot_fn)

                        res = res[res['X_WORLD']<sn_ra+(w/2)]
                        res = res[res['X_WORLD']>sn_ra-(w/2)]
                        res = res[res['Y_WORLD']<sn_dec+(w/2)]
                        res = res[res['Y_WORLD']>sn_dec-(w/2)]
                        ras,decs,mags,errs = res.X_WORLD.values,res.Y_WORLD.values,res.MAG_AUTO.values,res.MAGERR_AUTO.values
                    img = fits.open(img_fn)
                    fg = aplpy.FITSFigure(img,figure=fig,subplot=plot_locs[band])
                    fg.recenter(sn_ra,sn_dec,width=w,height=w)
                    fg.show_lines([ver_line,hor_line],color='r',linewidth=.5)
                    fg.show_grayscale(vmin=-1.,vmax=15.)
                    fg.axis_labels.hide()
                    fg.tick_labels.hide()
                    fg.set_theme('publication')
                    fg.ticks.set_length(-3)
                    fg.add_label(0.1,0.8,band,relative=True,color='r',fontsize=14)
                    if restype=='spec':
                         fg.show_circles(res.RA.values,res.DEC.values,radius=0.00027,edgecolor='g',facecolor='none',linewidth=.5,alpha=.8)
                    else:
                         fg.show_circles(res.X_WORLD.values,res.Y_WORLD.values,radius=0.00027,edgecolor='g',facecolor='none',linewidth=.5,alpha=.8)
                    for i in range(len(ras)):
                        fg.add_label(ras[i],decs[i]+0.0003,'%.2f +/- %.2f'%(mags[i],errs[i]),size=6,color='g')
                    print ("Finished loading the %s band region!" %band)

                    for index,row in match.iterrows():

                        if restype=='spec':
                            fg.show_circles(row['RA'],row['DEC'],radius = 0.00015,edgecolor='r',facecolor='none',linewidth=.8)
                            fg.add_label(row['RA']+0.00015,row['DEC']+0.00015,'%.2f +/- %.2f'%(row['MAG_AUTO_%s'%b],row['MAGERR_AUTO_%s'%b]),size=6,color='r')
                        else:
                            fg.show_circles(row['X_WORLD'],row['Y_WORLD'],radius = 0.00015,edgecolor='b',facecolor='none',linewidth=.8)
                            fg.add_label(row['X_WORLD']+0.00015,row['Y_WORLD']+0.00015,'%.2f +/- %.2f'%(row['MAG_AUTO'],row['MAGERR_AUTO']),size=6,color='b')

                    if counter in [0,2]:
                        fg.tick_labels.show_y()
                    if counter in [2,3]:
                        fg.tick_labels.show_x()
                else:
                    fg = aplpy.FITSFigure('/media/data3/wiseman/des/coadding/config/blank2.fits',figure=fig,subplot=plot_locs[band])
                    fg.axis_labels.hide()
                    fg.tick_labels.hide()
                    fg.add_label(0.5,0.5,'[Failed to load %s band image]'%band,relative=True,fontsize=12,color='black')
                if counter ==0:
                    fg.add_label(0.99,1.05,sn,relative=True,fontsize=14,color='black')

            plt.suptitle('Right Ascension (J2000)',x=0.57,y=0.04)
            if args.path =='sn_dir':
                savepath =os.path.join(sn_dir,'%s_stamp.pdf'%sn)
            else:
                savepath =os.path.join(args.path,'%s_stamp.pdf'%sn)

            plt.savefig(savepath)
            fig.close()
if __name__ == "__main__":
    logger = logging.getLogger('sn_host_lookup.py')
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(name)s - %(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("***********************************")
    logger.info("Initialising *** sn_host_lookup.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    logger.info("***********************************")
    args = parser()
    main(args,logger)
