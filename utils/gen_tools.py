"""gen_tools.py: A collection of generic functions to assist with processing and analysing stack data"""

import numpy as np
###################
### DES stuff ###
def get_des_bands():
    return ['g','r','i','z']

def get_good_des_chips():
    good_des_chips = []
    for c in range(1,63):
        if c not in [2,31,61]:
            good_des_chips.append(c)
    return good_des_chips
###################
# stats

def mc_robust_median(dist,
                     return_sigma = False,
                     return_indices = False,
                     nmad_cut = 4.0,
                     niter_max = 10):
    '''Credit: Mike Childress'''
    c = 0.6745
    data = np.copy(dist)
    ind0 = np.arange(len(data))
    ind = np.arange(len(data))
    prev_npts = -1
    niter = 0
    while(len(ind) != prev_npts and niter<niter_max):
        prev_npts = len(ind)
        m = np.median(data)
        MAD = np.median(np.abs(data - m)/c)
        ind = np.where(abs(data-m)<nmad_cut*MAD)[0]
        data = data[ind]
        ind0 = ind0[ind]
        niter += 1
    final_median = np.median(data)
    final_sigma = np.median(np.abs(data - final_median)/c)
    if return_sigma == True and return_indices == True:
        return final_median, final_sigma, ind0
    elif return_sigma == True:
        return final_median, final_sigma
    elif return_indices == True:
        return final_median, ind0
    else:
        return final_median
