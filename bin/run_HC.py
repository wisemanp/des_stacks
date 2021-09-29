import numpy as np
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import seaborn as sns
from astropy.coordinates import SkyCoord
import logging
from astropy.table import Table
import astropy.io.fits as fits
import os
import pickle
import scipy.stats as stats
from shutil import copyfile
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM 
from astropy import wcs
from astropy.nddata import NDData
import aplpy
from astroimtools import nddata_stats
from des_stacks import des_stack as stack
from des_stacks.bin import stack_all
from des_stacks.utils import stack_tools,source_tools,gen_tools
from des_stacks.analysis import astro
from des_stacks.utils.gen_tools import mc_robust_median as r_median
bands = gen_tools.get_des_bands()
sns.set_color_codes(palette='colorblind')
import time
import _pickle as cpickle
import itertools
import progressbar

plt.rcParams['errorbar.capsize']=4
from IPython.display import display, HTML
display(HTML(data="""
<style>div#notebook-container    { width: 95%; }
 div#menubar-container     { width: 65%; }  
div#maintoolbar-container { width: 99%; }
</style>
"""))
def compute_HC(DLR,SEP):
    dlr = np.sort(DLR) # first sort DLR array
    e = 1e-5 # define epsilon, some small number
    if len(dlr)==1: # only 1 object in radius
        HC = -99.0
    else:
        delta12 = dlr[1] - dlr[0] + e
        D1D2 = dlr[0]**2/dlr[1] + e
        Dsum = 0
        for i in range(0, len(dlr)):
            for j in range(i+1, len(dlr)):
                didj = dlr[i]/dlr[j] + e
                delta_d = dlr[j] - dlr[i] + e
                Dsum += didj/((i+1)**2*delta_d)
        HC = np.log10(D1D2*Dsum/delta12)
    return HC

v17 = pd.read_hdf('/media/data3/wiseman/des/coadding/results/deep/sngals_deep_v17.h5',index_col=0)
feature_cols = ['HC']
v17[feature_cols] = pd.DataFrame(columns=feature_cols,index=v17.index)
snidgroups = v17.groupby('SNID')
from tqdm import tqdm
for (n,g) in tqdm(snidgroups):
        
    if len(g[(g['DLR']>0)&(g['Z_RANK']<2)])>2:
        host_candidates = g[(g['DLR']>0)&(g['Z_RANK']<2)]

        v17.loc[g.index,feature_cols] = compute_HC(host_candidates['DLR'].values,host_candidates['SEPARATION'].values)
    elif len(g[g['DLR']==0])==1:
        host_candidates = g[(g['DLR']>0)&(g['Z_RANK']<2)]
        if len(g[(g['DLR']>0)&(g['Z_RANK']<2)])>2:
            v17.loc[g.index,feature_cols] = compute_HC(host_candidates['DLR'].values,host_candidates['SEPARATION'].values)
        else:
            v17.loc[g.index,feature_cols]=99         

    else:
        v17.loc[g.index,feature_cols]=99    
    
v17.to_hdf('/media/data3/wiseman/des/mismatch/v17_features.h5',index=True,key="main")      
