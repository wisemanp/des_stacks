# -*- coding: utf-8 -*-
'''Stack Everything'''
# exectuable to run through the entire stack process
import numpy as np
import pandas as pd
import seaborn as sns
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import datetime
import configparser
import os
import logging
import argparse
import glob
from time import gmtime, strftime
from des_stacks import des_stack as stack
from des_stacks.utils.loop_stack import iterate_sex_loop, init_sex_loop


all_years = ['none','1','2','3','4']
all_fields = ['SN-X1','SN-X2','SN-X3','SN-C1','SN-C2','SN-C3','SN-E1','SN-E2','SN-S1','SN-S2']
all_chips = np.arange(1,62)

def parser():
    parser = argparse.ArgumentParser(description='Stack some DES SN images')
    parser.add_argument('-f','--field', help = 'Field(s) to stack. Separate with space or comma (e.g. X2 X3)',nargs='?',required=False,default='X2')
    parser.add_argument('-b', '--band', help = 'Bands(s) to stack. Separate with space or comma (e.g. g r)',nargs='?',required=False,default='r')
    parser.add_argument('-my','--minusyears', help = 'Which minus years to stack (e.g. 1,2,3,4,none)',nargs='?',required=False,default='1')
    parser.add_argument('-ch','--chips', help = 'Which chips to stack (e.g. [1,5] = 1,3,4)',nargs=1,required=False,default='All')
    parser.add_argument('-wd','--workdir', help = 'Working directory', default = './')
    parser.add_argument('-l','--looptype', help ='Type of loop (can be "psf", "zp", "b" or "both")',required = False, default = 'both')
    parser.add_argument('-ps','--psfstep', help ='The size of the cut step if using psf',required=False, default =0.25)
    parser.add_argument('-zs','--zpstep', help ='Size of the cut step for zeropoint cuts',required=False,default = 0.025)
    parser.add_argument('-zc','--zcut', help ='Define zp cut to do a stack with',required=False,nargs=1,default= None)
    parser.add_argument('-pc','--pcut', help ='Define psf cut to do a stack with',required=False,nargs=1,default= None)
    parser.add_argument('-tc','--tcut', help ='Define teff cut to do a stack with',required=False,nargs=1,default= None,type=float)

    parser.add_argument('-ic','--initcuts',help = 'Define starting cuts for a loop stack',required=False,default= [-0.150,2.5])
    parser.add_argument('-t','--tidy',help = 'Tidy up temporary files after? 1 = Yes, 0 = no. Default = 1, turn off when testing',default = True)
    args=parser.parse_args()
    parsed = {}
    try:
        fields = args.field.split(',')

    except:
        try:
            fields = args.field[0].split(' ')
        except:
            fields =args.field

    for i in range(len(fields)):
        try:
            field = fields[i]
            field = 'SN-'+field
            fields[i]=field
        except:
            fields = 'SN-'+fields[0]
    parsed['fields']=fields

    try:
        bands = args.band.split(',')
    except:
        try:
            bands = args.band[0].split(' ')
        except:
            bands = args.band
    parsed['bands']=bands
    try:
        mys = args.minusyears.split(',')
    except:
        try:
            mys = args.minusyears[0].split(' ')
        except:
            mys = args.minusyears
    parsed['mys']=mys

    if args.chips != 'All':
        try:
            chips = args.chips[0].split(',')
        except:
            if args.chips[0][0]== '[':
                chip_bounds = args.chips[0][1:-1].split(',')
                chips = np.arange(int(chip_bounds[0]), int(chip_bounds[-1]))
            else:
                chips = args.chips[0].split(' ')

    else:
        chips = args.chips
    parsed['chips']=chips
    print ('Parsed chips as: %s' % chips)

    if not args.workdir:
        workdir = 'current'
    else:
        workdir = args.workdir

    parsed['workdir']=workdir

    try:
        loop_type = args.looptype
        parsed['looptype']=loop_type
    except:
        pass
    try:
        parsed['psf_step']=args.psfstep
    except:
        parsed['psf_step']=0.25

    try:
        parsed['zp_step']=args.zpstep
    except:
        parsed['zp_step'] = 0.05
    try:
        parsed['zcut'] = float(args.zcut[0])
    except:
        parsed['zcut'] = None
    try:
        parsed['pcut'] = float(args.pcut[0])
    except:
        parsed['pcut'] = None
    try:
        parsed['tcut'] = float(args.tcut[0])
    except:
        parsed['tcut'] = None
    try:
        parsed['init_cuts'] = args.initcuts
    except:
        parsed['init_cuts'] = [-0.150,2.5]
    try:
        parsed['tidy'] = args.tidy
    except:
        parsed['tidy'] = 1
    parsed['tidy'] = bool(float(parsed['tidy']))
    return parsed

    parsed['tidy'] = bool(float(parsed['tidy']))
def simple_stack(logger,parsed):
    '''code to run the stacks'''
    #read in parameters from the command line
    fields = parsed['fields']
    bands = parsed['bands']
    mys = parsed['mys']
    chips = parsed['chips']
    workdir = parsed['workdir']
    cuts={'zp':parsed['zcut'],'psf':parsed['pcut'],'teff':parsed['tcut']}
    logger.info("Parsed command line and will work on:\n Fields: %s \n Bands: %s \n MYs: %s \n Chips: %s"%(fields,bands,mys,chips))
    for f in fields:
        for b in bands:
            for my in mys:
                s = stack.Stack(f,b,my,chips,workdir)
                s.do_my_stack(cuts=cuts,final=True)
                s.run_stack_sex(cuts=cuts,final=True)
                s.init_phot()
                #if parsed['tidy']in [1,True]:
                    #for temp_fn in glob.glob(os.path.join(s.temp_dir,'*resamp*')):
                        #if not os.path.isdir(temp_fn):
                            #os.remove(temp_fn)

                    #for temp_fn in glob.glob(s.band_dir):
                        #if temp_fn not in glob.glob(os.path.join(s.band_dir,'*sci*'))+glob.glob(os.path.join(s.band_dir,'*wgt*')):
                           # if not os.path.isdir(temp_fn):
                                #os.remove(temp_fn)
def looped_stack(logger,parsed):
    fields = parsed['fields']
    bands = parsed['bands']
    mys = parsed['mys']
    chips = parsed['chips']
    workdir = parsed['workdir']
    loop_type = parsed['looptype']
    psfstep = parsed['psf_step']
    init_cut,init_step = parsed['init_cuts'],[parsed['zp_step'],parsed['psf_step']]
    logger.info("Parsed command line and will work on:\n Fields: %s \n Bands: %s \n MYs: %s \n Chips: %s"%(fields,bands,mys,chips))
    for f in fields:
        for b in bands:
            for my in mys:
                logger.info("Doing the initial stack for %s, %s band, MY %s" %(f,b,my))
                norm,normlim,stringent,stringentlim,generous,generouslim= init_sex_loop(logger,f,b,my,chips,loop_type,init_cut,init_step,workdir)

                for chip in chips:
                    logger.info("Looping chip %s" %chip)
                    av_normlim= np.median(normlim[chip])
                    av_strinlim=np.median(stringentlim[chip])
                    av_genlim =np.median(generouslim[chip])
                    logger.info("Average limiting magnitude for generous stack: %.3f" %av_genlim)
                    logger.info("Average limiting magnitude for middle stack: %.3f" %av_normlim)
                    logger.info("Average limiting magnitude for stringent stack: %.3f" %av_strinlim)
                    logger.info("******************************************************")
                    logger.info("Now starting the loop to get the best stack for chip %s" %chip)

                    # find out which of the stacks is deeper
                    if av_genlim>av_normlim:
                        logger.info("Generous cut gives deeper stack than normal stack")
                        if av_genlim>av_strinlim:
                            logger.info("Generous cut gives deeper stack than stringent stack")
                            logger.info("Heading in the generous direction!")
                            # The generous stack is deepest, so go a bit further
                            cuts,qual = iterate_sex_loop(logger, f, b, my, chip, loop_type,init_cut,init_step,workdir,'ge',av_genlim)
                            logger.info("Saved stack and exiting loop for chip %s" %chip)
                        elif av_genlim<=av_strinlim:
                            logger.info("Generous cut not as deep as stringent stack.")
                            logger.info("Heading in the stringent direction!")
                            cuts,qual = iterate_sex_loop(logger, f, b, my, chip, loop_type,init_cut,init_step,workdir,'st',av_strinlim)
                            logger.info("Not been told to do anything, exiting")
                            # explore some options. Unlikely to get to this point.
                    elif av_genlim<av_normlim:
                        logger.info("Stack with a generous cut not as deep as the normal stack")
                        if av_normlim<av_strinlim:
                            # The stringent cut is deepest, so go more stringent
                            logger.info("Stack with a normal cut also not as deep as the stringent cut")
                            logger.info("Heading in the stringent direction!")
                            cuts,qual = iterate_sex_loop(logger, f, b, my, chip, loop_type,init_cut,init_step,workdir,'st',av_strinlim)
                            logger.info("Saved stack and exiting loop for chip %s" %chip)
                        elif av_normlim>av_strinlim:
                            logger.info("Stack with normal cut is deeper than the stringent cut. Normal is good.\n Doing it again as a final stack.")
                            s = stack.Stack(f,b,my,[chip],workdir)
                            cuts = {'zp':init_cut[0],'psf':init_cut[1]}
                            s.do_my_stack(cuts,final=True)
                            s.run_stack_sex(cuts)
                            s.init_phot()
                    if parsed['tidy']==True:

                        logger.info("********************* Tidying up *********************")
                        try:
                            for filename in glob.glob(os.path.join(s.band_dir,'*temp.fits'))+glob.glob(os.path.join(s.band_dir,chip,'*temp.fits')):
                                os.remove(filename)
                        except:
                            pass
    logger.info("Done! stack_all.py finished. Enjoy your stacked data!")




def check_done(proc,wd):
    sd = os.path.join(wd,'stacks')
    ld = os.path.join(wd,'log')
    for y in all_years:
        done_df = pd.DataFrame(index = all_bands,columns = fields)
        for f in all_fields:
            for b in all_bands:
                bd = os.path.join(sd,'MY%s'%y,f,b)
                if b in ['g','r']:
                    tc = 0.15
                elif b in ['i','z']:
                    tc = 0.25
                cn = 0
                for c in all_chips:
                    if proc =='stack':
                        sci_fn = os.path.join(bd,'ccd_%s_%s_%s_sci.fits'%(str(c),b,tc))
                    elif proc =='sex':
                        sci_fn = os.path.join(bd,str(c),'ana','MY%s_%s_%s_%s_%s_sci.sexcat'%(y,f,b,str(c),tc))
                    elif proc =='phot':
                        sci_fn = os.path.join(bd,str(c),'ana','%s_%s_%s_%s_init.result'%(y,f,b,str(c)))
                    if os.path.isfile(sci_fn):
                        cn +=1
                done_df.loc[b,[f]]=cn
        now = datetime.datetime.now()
        today = '%s_%s_%s'%(now.year,now.month,now.day)
        done_df.to_csv(os.path.join(ld,'%s_%s_progress_%s.csv'%(y,proc,today)))
        print ('For MY ',y,', the following %s has been done:' %proc)
        print (done_df)
        print ('*********************')

if __name__=="__main__":
    logger = logging.getLogger('stack_all.py')
    logger.setLevel(logging.DEBUG)
    formatter =logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("***********************************")
    logger.info("Initialising *** stack_all.py *** at %s UT" % strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    logger.info("***********************************")
    parsed = parser()
    if 'looptype' in parsed.keys():
        if parsed['looptype']in ['b','both','zp','z','psf']:
            looped_stack(logger,parsed)

        else:
            simple_stack(logger,parsed)
    else:
        simple_stack(logger,parsed)
    check_done('stack',parsed['workdir'])
    check_done('sex',parsed['workdir'])
    check_done('phot',parsed['workdir'])
