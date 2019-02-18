import pandas as pd
import numpy as np
import os
import argparse
from des_stacks.utils.gen_tools import get_good_des_chips
good_des_chips = get_good_des_chips()

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--list',default = '/home/wiseman/code/des_stacks/source_lists/all_tranisents.txt')
    parser.add_argument('-sf','--savename',default = 'sngals_deep_v4.result')
    parser.add_argument('-df','--df',default = 'none')
    return parser.parse_args()

def main(args):

    #main_df = pd.DataFrame()
    #if args.df !='none':
        #main_df = pd.read_csv(args.df,index_col=0)

    resdir = '/media/data3/wiseman/des/coadding/5yr_stacks/CAP/'
    snlist = np.genfromtxt(args.list,dtype=str,delimiter='\n')
    for sn in snlist:
        main_f = open(os.path.join('/media/data3/wiseman/des/coadding/results',args.savename),'a')
        cat = os.path.join(resdir,sn,'%s.result'%sn)
        try:
            c = open(cat,'r')
            print ('Adding cat: %s'%cat, ' of length ',len(c.readlines()))
            for l in c.readlines():
                main_f.write(l)
            main_f.close()

        except:
            print ('Couldnt read cat: %s'%cat)
            main_f.close()
    print ('Saved new file to %s'%os.path.join('/media/data3/wiseman/des/coadding/results',args.savename))
if __name__=="__main__":
    args=parser()
    main(args)
