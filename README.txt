des_stacks
==========

des_stacks is a simple Python module based on SWarp and SExtractor,
available for free at https://www.astromatic.net/software

For a simple stack, use stack_all.py with command line arguments:

  -f: field (e.g. -f X2,X3)
  -b: band (e.g. -b g,r,i,z, default = all)
  -my: minusyear, i.e. the year to leave out of the stack (e.g. -my 1 or default = none)
  -ch:  chips -  (e,g [1,3,4,5], default = all
  -wd: working directory - This should be the highest level of directory for the stacks, in which sub-directories will be created automatically
  -c: cuts - if you want to stack based on zeropoint and seeing cuts, define them here in a list like [zpcut,psfcut]. If nothing is given, defaults are [-0.15,2.5] and seem to give pretty good stacks. The zp cut is given as a residual magnitude of the zeropoint of any given frame compared to the median for that chip,field,band,season.
  -t: tidy - tidies up temporary files (1/0 for y/n)

For an optimized stack, some or all of the following parameters are needed in addition:

  -l: looptype - If you want to do an optimized stack, this determines the parameter you want to optimize. Current options: zeropoint (zp), seeing (psf), or both.
  -ps: psf_step - If optimizing the psf, this is the step size of the iterations (in arcsec).    default = 0.25
  -zs: zp_step - If optimizing the zeropoint, this is the step size of the iterations (in mag).    default = 0.025
  -ic: init_cuts - Initial cuts in (residual) zeropoint and seeing. Default is [-0.150,2.5].
