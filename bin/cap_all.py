#!/home/wiseman/anaconda3/bin/python
'Tiny wrapper to do common aperture photometry on everything'
#Note: import this first else it crashes importing sub-modules
import argparse
import os
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

def main(args):
    fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
    mys = ['1','2','3','4','5','none']
    chips = get_good_des_chips()
    if args.year !='all':
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
    for my in mys:

        for f in fields:
            f = 'SN-'+f
            for ch in chips:
                n_bad = 0
                for category in [my,f,ch]:
                    if args.avoid:
                        if category in [i for i in args.avoid.split(',')]:
                            n_bad+=1
                if n_bad<2:
                    capdir = '/media/data3/wiseman/des/coadding/5yr_stacks/MY%s/%s/CAP/%s'%(my,f,ch)
                    done_phot = os.path.isfile(os.path.join(capdir,'%s_%s_%s_obj_deep.cat'%(my,f,ch)))
                    if not args.skipdone:
                        cap_phot_all(my,f,ch,autocuts=True)
                    elif not done_phot:
                        cap_phot_all(my,f,ch,autocuts=True)

if __name__=="__main__":
    args=parser()
    main(args)
