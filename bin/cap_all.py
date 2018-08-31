#!/home/wiseman/anaconda3/bin/python
'Tiny wrapper to do common aperture photometry on everything'
#Note: import this first else it crashes importing sub-modules
import argparse
import os
from des_stacks import des_stack as stack
from des_stacks.analysis.astro import cap_phot_all

fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
mys = ['1','2','3','4','5','none']
good_des_chips = []
for c in range(1,63):
    if c not in [2,31,61]:
        good_des_chips.append(c)

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--avoid',default=None)
    parser.add_argument('-s','--skipdone',action='store_true')
    return parser.parse_args()
def main(args):
    for my in mys:
        for f in fields:
            f = 'SN-'+f
            for ch in good_des_chips:
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
