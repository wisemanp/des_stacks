from des_stacks import des_stack as stack
def init_sex_loop(logger,f,b,my,chips,loop_type,init_cut,init_step,workdir):
    #do initial stack
    s = stack.Stack(f,b,my,chips,workdir)
    zpcut = init_cut[0]
    psfcut =init_cut[1]
    cuts = {'zp':zpcut,'psf':psfcut}
    logger.info("Starting by going in each direction in %s" %loop_type)
    if loop_type == 'zp':
        init_step[1] =0
    elif loop_type == 'psf':
        init_step[0] = 0
    zp_stringe,zp_gen = zpcut+init_step[0],zpcut-init_step[0]
    psf_stringe,psf_gen= psfcut-init_step[1],psfcut+init_step[1]
    s.do_my_stack(cuts,final=False)
    qual = s.run_stack_sex(cuts)
    # do stack with more generous cut
    cuts['zp'],cuts['psf']=zp_gen,psf_gen
    s.do_my_stack(cuts, final=False)
    qual_gen = s.run_stack_sex(cuts)
    # do stack with more stringent cut
    cuts['zp'],cuts['psf']=zp_stringe,psf_stringe
    s.do_my_stack(cuts,final=False)
    qual_stringe = s.run_stack_sex(cuts)
    #test which is better
    #try going further in that direction
    #until we get no further.
    logger.info("Quality of normal stack %s" %qual)
    logger.info("Quality of generous stack %s" %qual_gen)
    logger.info("Quality of stringent stack %s" %qual_stringe)
    return (qual, qual_stringe, qual_gen)

def iterate_sex_loop(logger,f,b,my,chip,loop_type,init_cut,init_step,workdir,direction,quals):
    s = stack.Stack(f,b,my,[chip],workdir)
    if direction == 'ge':
        f = -1
    elif direction == 'st':
        f = 1
    if loop_type == 'zp':
        init_step[1] ==0
    elif loop_type == 'psf':
        init_step[0] == 0
    step=[f*init_step[0],-1*f*init_step[1]]
    carryon = True
    n = 2
    nstep = 1
    while carryon == True:
        zp_start = init_cut[0]+(n*step[0])
        psf_start = init_cut[1]+(n*step[1])
        if zp_start >0:
            logger.info("This cut is going to cost you the majority of your frames, so the code will probably crash. Aborting now to save face")
            n-=nstep
            cuts = {'zp':zp_start-(n*step[0]),'psf':psf_start-(n*step[1])}
            carryon = False
        else: 
           cuts = {'zp':zp_start,'psf':psf_start}
        s.do_my_stack(cuts,final=False)
        new_qual = s.run_stack_sex(cuts)
        impr = new_qual.loc[chip,'zp']-quals.loc[chip,'zp']
       
        if impr>0.01:
            logger.info("Iteration %s improved the zp by %s mag; carrying on" %(n-1,impr))
            n+=nstep
        elif impr >0.005:
            logger.info("Iteration %s improved the zp by %s mag; carrying on" %(n-1,impr))
            n+=nstep/2
        elif impr < 0:
            logger.info("Iteration %s didn't improve the zp at all; stopping."% (n-1))
            n=2
            carryon = False
            
        else:
            logger.info("Converged to a best zp of %s mag after %s iterations" %(new_qual.loc[chip,'zp'],n-1))
            s.do_my_stack(cuts,final = True)
            carryon = False
    if n==2:
        zp_start = init_cut[0]+(n*step[0])
        psf_start = init_cut[1]+(n*step[1])
        cuts = {'zp':zp_start,'psf':psf_start}
        s.do_my_stack(cuts,final =True)
        new_qual = s.run_stack_sex(cuts)
    return (new_qual,cuts)
