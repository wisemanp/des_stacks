import numpy as np
import pandas as pd
import os
from des_stacks import des_stack as stack
from des_stacks.utils.stack_tools import get_cuts
years = ['1','2','3','4','5','none']
fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
bands = ['g','r','i','z']

def add_lim(y,f,chip):
    res_df = pd.read_csv(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks','MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_obj_deep.cat'%(y,f,chip)),index_col=0)
    for b in bands:
        cuts = get_cuts(f,b)
        s = stack.Stack(f, b, y, [str(chip)] ,'coadding',cuts=cuts)
        quals= np.loadtxt(os.path.join(s.band_dir,str(chip),'ana','%s_ana.qual'%s.cutstring))
        zp,zp_sig,av_fwhm = (float(quals[i]) for i in [0,1,2])
        with open(os.path.join(s.band_dir,str(chip),'ana','%s_%s_%s_%s_init_wgtd.result'%(y,f,s.band,chip)),'r') as res:
            header = [next(res) for x in range(9)]
        limmag = header[-1].split(' ')[-1].strip('\n')
        limflux = 10**((float(limmag)-zp)/-2.5)
        res_df['LIMMAG_%s'%b] = limmag
        res_df['LIMFLUX_%s'%b] = limflux
    res_df.to_csv(os.path.join('/media/data3/wiseman/des/coadding/5yr_stacks','MY%s'%y,f,'CAP',str(chip),'%s_%s_%s_obj_deep.cat'%(y,f,chip)),index_col=0)
    print ('Done',y,f,chip)
def main():
    good_des_chips = []
    for c in range(1,63):
        if c not in [2,31,61]:
            good_des_chips.append(c)


    for y in years:
        for f in fields:
            f = 'SN-'+f
            for chip in good_des_chips:
                add_lim(y,f,chip)

if __name__ == '__main__':
    main()
