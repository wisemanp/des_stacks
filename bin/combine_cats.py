import pandas as pd
import numpy as np
import os
import argparse
from des_stacks.utils.gen_tools import get_good_des_chips
good_des_chips = get_good_des_chips()

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--field',default = 'all')
    parser.add_argument('-my','--year',default='none')
    parser.add_argument('-ch','--chip',default='all')
    parser.add_argument('-df','--df',default = 'none')
    parser.add_argument('-sf','--savename',default = 'combined.cat')
    return parser.parse_args()

def main(args):
    fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
    mys = ['none']

    if args.year !='none':
        mys = args.year
    if args.field != 'all':
        try:
            fields = args.field.split(',')

        except:
            try:
                fields = args.field[0].split(' ')
            except:
                fields =args.field

    if args.chip !='all':
        chips = args.chip.split(',')
    main_df = pd.DataFrame()
    if args.df !='none':
        main_df = pd.read_csv(args.df,index_col=0)
    for my in mys:

        for f in fields:
            f = 'SN-'+f
            for ch in good_des_chips:
                ch = int(ch)
                cap_chip_dir = '/media/data3/wiseman/des/coadding/5yr_stacks/MY%s/%s/CAP/%s'%(my,f,ch)
                cat = os.path.join(cap_chip_dir,'%s_%s_%s_obj_deep.cat'%(my,f,ch))
                cat_df = pd.read_csv(cat,index_col=0)
                print ('Adding cat: %s'%cat, ' of length ',len(cat_df))
                main_df = main_df.append(cat_df)
    main_df.to_csv(os.path.join('/media/data3/wiseman/des/coadding/results',args.savename))
    print ('Saved new file to %s'%os.path.join('/media/data3/wiseman/des/coadding/results',args.savename))
if __name__=="__main__":
    args=parser()
    main(args)
