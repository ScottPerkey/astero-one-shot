'''
fits using lsq
'''

class Params(object):
    '''
    various parameters for fitting dnu
    '''
    # use the maximum likelihood if True. if False, use the mean likelihood solution
    use_max = True
    fixed = False
    do_mcmc = False
    # DO Maxmium Likelihood Estimation. IF TRUE, SET AN OPTION BELOW (do_sa (simulated annealing), do_de (differential evolution), do_pso (particle swarm), do_nelder_mead (Nelder-Mead/amoeba/simplex)
    do_mle = True
    do_bfgs = False
    do_nelder_mead = False
    do_newton_cg = False
    do_sa = False
    do_de = False
    do_pso = False
    # JCZ 210916 
    # misnomer.
    # really, whether or not to use the least-squares fitting algorithm, curve_fit
    # whether or not Gaussian (least-squares) stats are used for \ln\mathcal{L} depends solely on whether or not collapsed == True. If True, then lsq stats are used. but fitting algorithm can still be any of these options.
    do_lsq = True
    # how many MLE runs should be tried to try to get the global maximum not
    # just a local one... if only do_mle is run, the std among the <do_mle_nloop>
    # runs will be taken to be the error. you can increase the number of these runs and
    # maintain a reasonable runtime by increasing do_mle_nloop (default 4) and increasing minfunc, (default 0.0001), which dictates the minimum difference bewteen succssive swarm fitnesses below which the swarm run is terminated.
    do_mle_nloop = 1
    minfunc = 0.0001
    # if True, will fix fwhm and use its parameter slot in the parameter
    # vector as an error for the chi2 part of logl to marginalize over the unknown
    # error of the collapsed echelle diagram
    # defer to the bayes_harvey_pypar.py file for whether or not to be
    # in debugging mode.
    debug = False
    verbose = False
    debug_plot = False
    # counter for debugging plots, called 'numax_region_w_model'+_params.plot_count. only plotting things for
    # debugging if above <debug_plot> is True.
    plot_count = 0
    fix_fwhm = False
    # don't allow the different orders to have independent positions of the radial position
    fixed_epsilon = False
    plot = True
    # if True, run the fiducial test object only
    fiducial = False
    ln = True
    entropy = False
    # if False, will fit the entire region around numax.
    # if True, will fit the collapsed echelle diagram
    collapsed = True
    # compute the error on dnu using the second derivatve of likelihood function wrt dnu. also makes a plot of the second derivative in the region around dd.
    # JCZ 210916
    # should really not be used any more...
    dd = False
    # do lsq fitting across a grid of dnu, and compute dd (if dd == True, above) using this grid of maximum lnL values as a function of dnu. also plots this.
    grid = False
    # JCZ 280916
    # whether or not to apply numax prior
    numax_prior = True
    # whether or not to apply small freq sep prior
    # JCZ 100717
    # CHANGED THIS TO FALSE TO TRY TO REMOVE REMAINING DIFFERENCE BETWEEN DANIEL AND MY 02 SEPARATIONS
    deltanu_prior = False
    
    # whether or not to apply visibility prior
    amp_prior = True
    # search around not the expected nu from numax but rather the dnu that has
    # highest correlation with the expected three-spike pattern. otherwise, use
    # the expected dnu from numax as the starting dnu for mcmc/MLE
    # JCZ 141118
    # this used to be False at some point...
    use_corr_dnu = True
params = Params()
