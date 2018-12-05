#!/home/wiseman/anaconda3/bin/python
'Tiny wrapper to do common aperture photometry on everything'
#Note: import this first else it crashes importing sub-modules
import argparse
import os
import multiprocessing
from multiprocessing import Process
import pathos.pools as pp

from des_stacks import des_stack as stack
from des_stacks.analysis.astro import cap_phot_all
from des_stacks.utils.gen_tools import get_good_des_chips


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--avoid',default=None)
    parser.add_argument('-s','--skipdone',action='store_true')
    parser.add_argument('-f','--field',default = 'all')
    parser.add_argument('-my','--year',default='all')
    parser.add_argument('-ch','--chip',default='all')
    return parser.parse_args()



def phot_worker(arg_pair):
    args, chip = arg_pair[0],arg_pair[1]
    print(args)
    my,f,parsed_args = [args[i]for i in range(len(args))]
    
    ch = int(chip)
    print ('In phot_worker to do cap_phot_all on chip %s'%ch)
    n_bad = 0
    for category in [my,f,ch]:
        if parsed_args.avoid:
            if category in [i for i in parsed_args.avoid.split(',')]:
               n_bad+=1
    if n_bad<2:
        capdir = '/media/data3/wiseman/des/coadding/5yr_stacks/MY%s/%s/CAP/%s'%(my,f,ch)
        done_phot = os.path.isfile(os.path.join(capdir,'%s_%s_%s_obj_deep.cat'%(my,f,ch)))
        if not parsed_args.skipdone:
            cap_phot_all(my,f,ch,autocuts=True)
            return 
        elif not done_phot:
            cap_phot_all(my,f,ch,autocuts=True)
            return 

def multi_phot(my,f,chips,parsed_args):

    args = [my,f,parsed_args]
    pool_size = multiprocessing.cpu_count()*2
    act = multiprocessing.active_children()
    pool = pp.ProcessPool(processes=pool_size,
                                maxtasksperchild=2,
                                )
    pool._clear()
    pool._serve()

    chips = list(chips)
    
    all_args = []
    for c in chips:
        all_args.append([args,c])
        #p = Process(target=worker,args=(args,c))
        #p.start()
        #p.join()
    
    results = pool.map(phot_worker,all_args)
    
    pool.close()
    pool.join()
    return results

def main(parsed_args):
    fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
    mys = ['1','2','3','4','5','none']
    chips = get_good_des_chips()
    if parsed_args.year !='all':
        mys = parsed_args.year.split(',')
    if parsed_args.field != 'all':
        try:
            fields = parsed_args.field.split(',')

        except:
            try:
                fields = parsed_args.field[0].split(' ')
            except:
                fields = parsed_args.field

    if parsed_args.chip !='all':
        chips = parsed_args.chip.split(',')
    for my in mys:

        for f in fields:
            f = 'SN-'+f
            res = multi_phot(my,f,chips,parsed_args)

if __name__=="__main__":
    parsed_args=parser()
    main(parsed_args)
