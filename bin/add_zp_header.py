#!/home/wiseman/anaconda3/bin/python
outdir = '/media/data3/wiseman/des/coadding/stacks/'
years = ['1','2','3','4','none']
fields = ['E1','E2','S1','S2','C1','C2','C3','X1','X2','X3']
bands = ['g','r','i','z']
chips = range(1,63)
for y in years:
    my = 'MY'+y
    for f in fields:
        f='SN-'+f
        for b in bands:

            for chip in chips:
                if chip not in [2,31,61]:
                    sci_fns = glob.glob(os.path.join(outdir,my,f,b,'ccd_%s_%s_????_sci.fits'%(chip,b)))
                    print (sci_fns)
                    fn = sci_fns[0]
                    qual = np.loadtxt(os.path.join(outdir,my,f,b,str(chip),'ana','%s_ana.qual'%fn[-13:-9]))
                    zp = qual[0]
                    hdulist=fits.open(sci_fns[0])
                    h = hdulist[0].header
                    h['MAG_ZP'] = zp
                    hdulist[0].header=h
                    hdulist.writeto(sci_fns[0],overwrite=True)
