#!/home/wiseman/anaconda3/bin/python
# -*- coding: utf-8 -*-
'''Script to find the photometry for a certain supernova'''
# exectuable to run through the entire stack process
import numpy as np
import pandas as pd
import logging
import argparse
from time import gmtime, strftime
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import astropy.io.fits as fits
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import _pickle as cpickle
import math
import glob
import seaborn as sns
import itertools


sns.set_color_codes(palette='colorblind')

plot_locs_labeled={
'g':[0.18,0.53,0.40,0.43],
'r':[0.59,0.53,0.40,0.43],
'i':[0.18,0.09,0.40,0.43],
'z':[0.59,0.09,0.40,0.43]
}

plot_locs_paper = {
'g':[0.02,0.51,0.47,0.47],
'r':[0.51,0.51,0.47,0.47],
'i':[0.02,0.02,0.47,0.47],
'z':[0.51,0.02,0.47,0.47]
}

bands = ['g','r','i','z']
pix_arcsec = 0.264
def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ra','--ra',help='RA of target',default=None)
    parser.add_argument('-de','--dec',help='DEC of target',default=None)
    parser.add_argument('-l','--coordlist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-wd','--workdir',help='Path to directory to work in',required = False,default = '/media/data3/wiseman/des/coadding')
    parser.add_argument('-p','--path',help='Full path to output image if you are making a stamp',default = 'sn_dir')
    parser.add_argument('-vm','--vmin',help='vmin for greyscale',default=-0.8)
    parser.add_argument('-vx','--vmax',help='vmax for greyscale',default=15.0)
    parser.add_argument('-fc','--finder',help = 'Is this a finder chart?',action='store_true')
    parser.add_argument('-pa','--paper',help = 'Is this for a paper?',action='store_true')
    parser.add_argument('-re','--resfile',help = 'File to find host phot results for this SN',default = None)
    parser.add_argument('-b','--band',help='Only use one band? If so, enter here',default='All')
    parser.add_argument('-f','--ftype',help='File type for save. Default = pdf',default='pdf')
    parser.add_argument('-s','--size',help='Size in arcsec',default=60,dtype='float')

    return parser.parse_args()

def find_chip(ra,dec):
    f=open('/media/data3/wiseman/des/coadding/config/chiplims.pkl','rb')
    chiplims = cpickle.load(f)

    #################
    for field in ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']:
        the_field = chiplims[field]
        for ccd in the_field.keys():
            if the_field[ccd][0][0] > ra > the_field[ccd][2][0]:
                if the_field[ccd][0][1] < dec < the_field[ccd][1][1]:
                    return ('SN-%s'%field,ccd)

def main(args,logger):
    if args.ra:
        l = [[args.ra,args.dec]]
    else:
        l = np.loadtxt(args.coordlist)
    for coords in l:
        ra,dec = coords[0],coords[1]
        f,ccd = find_chip(ra,dec)
        import aplpy
        fig,ax = plt.subplots() #figsize=(16,9)
        w = args.size/3600
        ax.set_xticks([])
        ax.set_yticks([])
        for loc in ['top','right','left','bottom']:
            ax.spines[loc].set_visible(False)
        #ax.set_ylabel('Declination (J2000)',fontsize=12,labelpad = 30)

        if args.band !='All':
            bands = [args.band]
        else:
            bands = ['g','r','i','z']
        palette = itertools.cycle(sns.color_palette(palette='colorblind',n_colors=5))
        fdir = '/media/data3/wiseman/des/coadding/5yr_stacks/MY1/%s/'%f
        for counter,b in enumerate(bands):
            bdir = os.path.join(fdir,b)
            capdir = os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks/CAP/%s_%s'%(f,ccd))
            if not os.path.isdir(capdir):
                os.mkdir(capdir)
            color = next(palette)
            if counter==1:
                color=next(palette)
            try:
                img_fn = glob.glob(os.path.join(capdir,'ccd_%s*%s*_sci.resamp.fits'%(ccd,b)))[0]
            except:
                from des_stacks import des_stack as stack
                from des_stacks.utils.stack_tools import make_cap_stamps,get_cuts
                cuts = [get_cuts(f,b) for b in bands]
                sg,sr,si,sz = [stack.Stack(f, b, y, [str(chip)] ,'coadding',cuts[counter]) for counter,b in enumerate(bands)]
                make_cap_stamps(sg,sr,si,sz,chip,'%s_%s'%(f,ccd),ra,dec,args.size/0.264,args.size/0.264)
                img_fn = glob.glob(os.path.join(capdir,'ccd_%s*%s*_sci.resamp.fits'%(ccd,b)))[0]

            if os.path.isfile(img_fn):


                img = fits.open(img_fn)
                if not args.paper:
                    plot_locs = plot_locs_labeled
                else:
                    plot_locs = plot_locs_paper

                if args.band != 'All':
                    fg = aplpy.FITSFigure(img,figure=fig)
                else:
                    fg = aplpy.FITSFigure(img,figure=fig,subplot=plot_locs[b])
                try:
                    fg.recenter(ra,dec,w)
                except:
                    logger.info('Could not recenter to outside the frame')


                fg.show_grayscale(vmin=float(args.vmin),vmax=float(args.vmax))
                fg.axis_labels.hide()
                fg.tick_labels.hide()
                fg.set_theme('publication')
                fg.ticks.set_length(0)
                if args.paper:
                    fg.ticks.hide()
                fg.add_label(0.1,0.8,b,relative=True,color=color,fontsize=24,weight='bold')
                fg.add_scalebar(1/360,color=color,linewidth=3,fontsize=20,weight='bold')
                fg.scalebar.set_label("1'")

            else:
                fg = aplpy.FITSFigure('/media/data3/wiseman/des/coadding/config/blank2.fits',figure=fig,subplot=plot_locs[b])
                fg.axis_labels.hide()
                fg.tick_labels.hide()
                fg.add_label(0.5,0.5,'[Failed to load %s band image]'%b,relative=True,fontsize=12,color='black')
            if counter ==0 and not args.paper:
                fg.add_label(0.99,1.05,sn,relative=True,fontsize=14,color='black')
        if args.paper:
            plt.subplots_adjust(left=0.02,right=0.98)
        #plt.suptitle('Right Ascension (J2000)',x=0.57,y=0.04)
        if args.path =='sn_dir':
            savepath =os.path.join(capdir,'%s_%s_stamp.%s'%(f,ccd,args.ftype))
        else:
            savepath =os.path.join(args.path,'%s_%s_stamp.%s'%(f,ccd,args.ftype))

        plt.savefig(savepath)
        plt.close(fig)
        logger.info("Figure saved to %s"%savepath)
        logger.info("************* Finished making stamp for %s, %s! *************"%(f,ccd))
if __name__ == "__main__":
    logger = logging.getLogger('des_deep_stamp.py')
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(name)s - %(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("***********************************")
    logger.info("Initialising *** des_deep_stamp.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    logger.info("***********************************")
    args = parser()
    main(args,logger)
    logger.info("************* Finished all coords you gave me! *************")
