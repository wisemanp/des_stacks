import os
import glob

from des_stacks.utils.gen_tools import get_des_bands,get_good_des_chips

fields = ['X1','X2','X3','C1','C2','C3','E1','E2','S1','S2']
mys = ['1','2','3','4','5','none']
chips = get_good_des_chips()
bands = get_des_bands()
wd = '/media/data3/wiseman/des/coadding/5yr_stacks'
for my in mys:
    for f in fields:
        for b in bands:
            for ch in chips:
                clipweighted_sci_fns = glob.glob(os.path.join(wd,'MY%s'%my,'SN-%s'%f,b,'*%s*clipweighted_sci*'%ch))
                for fn in clipweighted_sci_fns:
                    fn_start = os.path.split(fn)[0]
                    fn_end = os.path.split(fn)[-1]
                    os.rename(fn,os.path.join(fn_start,'MY%s_SN-%s'%(my,f)+fn_end))
