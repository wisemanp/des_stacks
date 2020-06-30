# des_stacks
==========
# v0.1.1
==========

des_stacks is a simple Python module based on SWarp, SExtractor, and PSFEx
available for free at https://www.astromatic.net/software
des_stacks is intended to be used on data from the Dark Energy Survey Supernova Programme (DES-SN). A description of the process is available in Wiseman et al. 2020 https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.4040W.

For details about using the software, please contact Phil Wiseman (pacelweb@gmail.com)

For a simple stack, use stack_all.py with command line arguments:

  -f: field (e.g. -f X2,X3)
  -b: band (e.g. -b g,r,i,z, default = all)
  -my: minusyear, i.e. the year to leave out of the stack (e.g. -my 1 or default = none)
  -ch:  chips -  (e,g [1,3,4,5], default = all
  -wd: working directory - This should be the highest level of directory for the stacks, in which sub-directories will be created automatically
  -pc: psf cut - if you want to stack excluding objects above a certain seeing value, define tit here. If nothing is given, default is 2.5. 
  -tc: teff cut - same as above but for t_effective. Default is 0.15
    OR
  -o: optimized - uses the optimized values for each field/band as defined in Wiseman et al. 2020. 
  -t: tidy - tidies up temporary files (1/0 for y/n)

For an optimized stack, some or all of the following parameters are needed in addition:

  -l: looptype - If you want to do an optimized stack, this determines the parameter you want to optimize. Current options: zeropoint (zp), seeing (psf), or both.
  -ps: psf_step - If optimizing the psf, this is the step size of the iterations (in arcsec).    default = 0.25
  -zs: zp_step - If optimizing the zeropoint, this is the step size of the iterations (in mag).    default = 0.025
  -ic: init_cuts - Initial cuts in (residual) zeropoint and seeing. Default is [-0.150,2.5].

# Requirements:

## python packages:
 * numpy
 * matplotlib
 * scipy
 * astropy
 * pandas
 * configparser
 * subprocess
 * multiprocessing
 * pathos

## AstroMatic software:
 * SWarp
 * Source Extractor
 * PSFex
