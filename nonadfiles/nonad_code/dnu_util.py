import pandas as pd
import numpy as np
from .plot_funcs import Plot
from .util import muhz2pday, pday2muhz
# from memory_profiler import profile
#condorized @profile
def lorentzian(x, loc, amp, width):
    return amp/(1.0 + (x - loc)**2/width)
#condorized @profile
def lorentzians(f, xs, amps, dnu, fwhm, fixed=False):
    '''
    returns the sum of three lorentzians up to <ell>, with epsilon describing the shift of the radial peak from mod nu/dnu = 0.0., and can go up to 1.0.
    the radial position is epsilon.
    the free parameters are dnu, epsilon, xs, and amps
    Inputs
    f : frequency array, in muhz
this will be transformed into mod nu/dnu within this program
    xs : ndarray
    locations of the peaks
    amps : ndarray
    amplitudes of the peaks, relative to the highest
    power : power array normalized to the guess as to the radial peak
    fwhm : float
    fwhm of the smoothed gaussian excess. this determines the width of the lorentzian profiles of each mode
    x0 : int
    location in the f array of the HIGHEST peak (l = 0)
    [ ell : int ]
    how many modes to fit. only 3 for K2
    Outputs
    l : ndarray
     the total lorentzian profile (sum of <p.ell> components)
    [ fixed : bool ]
    if true, the absolute difference between peaks is fixed
    Notes
       
    '''

    ell = len(xs)
    # first we want to normalize the peaks to the zeroth peak, which we will label as the second highest peak
    # width of the lorentzians as a fraction of the FWHM
    # JCZ 160218
    # !!! DEBUG

    if type(fwhm) == type(np.float64(1.0)) or type(fwhm) == type(int):
        w = np.array([0.05, 0.1, 0.05, 0.0])*fwhm
    else:
        w = np.array([0.0]*4)
        w[0] = fwhm[0]
        w[1] = fwhm[1]
        # JCZ 170218
        w[2] = fwhm[2]

       #v = np.array([1.0, 0.48, 1.54, 0.043])
    # this uses dennis' values
    # v = np.array([1.0, 0.5, 0.8, 0.0])
    lor = f*0.
    epsilon = xs[0]
    # deltanu_const = np.array([0.0, 0.0, 0.047, 0.16])[0:ell]
    if fixed :
        diffs = xs*0.
        if ell >= 1:
            diffs[1] = (dnu/2. - dnu) + 0.025*dnu
        if ell >= 2:
            diffs[2] = -(0.121)*dnu - 0.047

    if fixed:
        # JCZ 081117
        # this should not be the case. fixed should only mean that the relative locations of the modes are fixed, not the radial amplitude.
        # epsilon = amps[0]
        # _amps = amps[:]
        # _amps[0] = 1.0
        _amps = amps
        pass
    else:
        # prepend the radial amplitude, which is fixed to 1 in the case where one allows the small freq
        # separations to vary (i.e., fixed = False)
        # JCZ 170517
        # think this was wrong. amps was actually [a1, a2, dnu] instead of [a0, a1, a2]
        # _amps = np.hstack(([1.],amps))
        _amps = amps

    for l in range(ell):
        if fixed:
            pos = (diffs[l]+epsilon) % dnu
            lor += lorentzian(f, pos, _amps[l], w[l])
            # wrap around to the other side
            #print diffs[l], epsilon[0], dnu
            if pos > dnu/2:
                lor += lorentzian(f, pos-dnu, _amps[l], w[l])
            else:
                lor += lorentzian(f, pos+dnu, _amps[l], w[l])
        else:
            pos = xs[l] % dnu
            lor += lorentzian(f, pos, _amps[l], w[l])
            # wrap around to the other side
            if pos > dnu/2:
                lor += lorentzian(f, pos-dnu, _amps[l], w[l])
            else:
                lor += lorentzian(f, pos+dnu, _amps[l], w[l])




    # lor += const
    return lor


#condorized @profile
def lorentzians_full(f, xs, amps, dnu, fwhm, fixed=False, epsilon=None, fixed_epsilon=False, wrap=False):
    '''
    returns the sum of three lorentzians up to <ell> for multiple orders, with epsilon describing the shift of the radial peak from mod nu/dnu = 0.0., and can go up to 1.0. 
    all the peaks' locations and amplitudes are expressed wrt the radial order.
    the free parameters are dnu, epsilon, xs, and amps
    Inputs
    f : frequency array, in muhz
this will be transformed into mod nu/dnu within this program
    xs : ndarray [n, ell]
    locations of the peaks
    amps : ndarray [n, ell]
    amplitudes of the peaks, relative to the highest, of size [n, ell]
    power : power array normalized to the guess as to the radial peak
    fwhm : ndarray [n, ell]
    fwhm of the smoothed gaussian excess. this determines the width of the lorentzian profiles of each mode
    epsilon : float | None
     guess for the offset of the radial mode. If None, will use the radial mode locations of each order (or the first order + dnu*n, if fixed_epsilon is True). Default None.
    x0 : int
    location in the f array of the HIGHEST peak (l = 0)
    [ ell : int ]
    how many modes to fit. only 3 for K2
    [ wrap : bool ]
    if True, assumes that <f> has a maximum value of input dnu and will wrap the model for a collapsed echelle diagram model. if False, no wrapping is done.
    Outputs
    l : ndarray
     the total lorentzian profile (sum of <p.ell> components)
    [ fixed : bool ]
    if true, the absolute difference between peaks is fixed
    [ fixed_epsilon : bool ]
     if True, only consider the first order radial position as the epsilon, rather than allowing each order
    to have its own epsilon. only valid if <fixed> is True.
    Notes
       
    '''

    # how many ms are being fit
    ell = xs.shape[1]
    # how many orders, n, are being fit
    n = xs.shape[0]

    # first we want to normalize the peaks to the zeroth peak, which we will label as the second highest peak
    # width of the lorentzians as a fraction of the FWHM
    # JCZ 160218
    # DEBUG !!!
    if type(fwhm) == type(np.float64(1.0)):
        w = np.array([0.05, 0.1, 0.05, 0.0])*fwhm
    else:
        w = np.array([0.0]*4)
        w[0] = fwhm[0]
        w[1] = fwhm[1]
        # JCZ 170218
        w[2] = fwhm[2]

    lor = f*0.
    # JCZ 150817
    epsilon_0 = epsilon
    # 
    epsilons = xs*0.0

    # sum over orders, n
    for _n in range(n):

        # sum over each mode in a given order, _n
        if fixed :
            diffs = xs[_n,:]*0.
            if ell >= 1:
                diffs[1] = (dnu/2. - dnu) + 0.025*dnu
            if ell >= 2:
                diffs[2] = -(0.121)*dnu - 0.047

            if fixed_epsilon:
                # JCZ 150817
                # epsilon = amps[0,0] + _n*dnu
                # JCZ 201017
                # added if and else clause. to undo, remove else clause entirely and keep the commented contents of the if clause but without the if statement.
                if epsilon_0 is not None:
                    epsilon = epsilon_0 + _n*dnu
                    # epsilon = epsilon_0 + _n*dnu
                else:
                    epsilon = xs[0,0] + _n*dnu
            else:
                # JCZ 201017
                # this is the radial amplitude, not used for epsilon anymore...
                # epsilon = amps[_n,0]
                # JCZ 201017
                # changing to xs[_n,0]
                # changing to only the first mode's position because this is fixed so there is only one free position parameter!
                epsilon = xs[0,0]
                # JCZ 201017
                # making it a vector
                # epsilons[_n] = xs[_n,0]
                
            _amps = amps[_n,list(range(ell))]
            # JCZ 081117
            # removing setting radial amplitude to 1.0.
            # _amps[0] = 1.0
            # JCZ 081117
            # added + epsilon here instead of below, where i had
            #             if fixed:
                # pos = (_xs[l] + epsilon)

            _xs = diffs*0. + epsilon
            _xs[1] = dnu/2. - diffs[1]
            _xs[2] = -diffs[2]

        else:
            # prepend the radial amplitude, which is fixed to 1 in the case where one allows the small freq
            # separations to vary (i.e., fixed = False)
            # JCZ 170517
            # think this was wrong. amps was actually [a1, a2, dnu] instead of [a0, a1, a2]
            # _amps = np.hstack(([1.], amps[_n,:]))
            _amps = amps[_n,:]
            # added epsilon info for the non-fixed case.
            if epsilon_0 is not None:
                epsilon = epsilon_0 + _n*dnu
                # epsilon = epsilon_0 + _n*dnu
            else:
                epsilon = xs[0,0] + _n*dnu
            # print epsilon
            # JCZ 141117
            # don't add epsilon...
            # xs += epsilon
        # JCZ 191217
        # commented out the % dnu parts 
        for l in range(ell):
            if fixed:
                # pos = (_xs[l])
                # lor += lorentzian(f, pos, _amps[l], w[l])
                pos = (diffs[l]+epsilon) #% dnu
                lor += lorentzian(f, pos, _amps[l], w[l])
                # wrap around to the other side
                #print diffs[l], epsilon[0], dnu
                # if pos > dnu/2:
                    # lor += lorentzian(f, pos-dnu, _amps[l], w[l])
                # else:
                    # lor += lorentzian(f, pos+dnu, _amps[l], w[l])

            else:
                pos = xs[_n,l] #% dnu
                lor += lorentzian(f, pos, _amps[l], w[l])
            # JCZ 211017
            # !!! DEBUG changed from wrap to True
            # JCZ 191217
            # changing back to wrap
            if wrap:
                # wrap around to the other side
                if pos > dnu/2:
                    lor += lorentzian(f, pos-dnu, _amps[l], w[l])
                else:
                    lor += lorentzian(f, pos+dnu, _amps[l], w[l])

                
    return lor

# for every position, make a thing
# and then scale their amplitudes to be in the model ratio
# these expressions are from Chaplin+ 1997
def xt(t, A, dxdt_0, omega_0, eta, x_0):
    '''
    Inputs
    dxdt_0 : float
     initial velocity (scalar)
    omega_0 : float
     central freqeuncy of the mode (scalar).
    eta : float
     damping rate (scalar).
    x_0 : float
     initial displacement (scalar).
    A : float
     amplitude of the driving term (scalar).
    t : ndarray | float
     time.
    Outputs
    x(t) from Chaplin+ 1997
    See also
    lorentzians_fake()
    '''
    exp = np.exp(-eta*t)
    arg = np.sqrt(omega_0**2 - eta**2)*t
    # this is a constant (t dependence is divided out)
    C = (A + dxdt_0 + eta*x_0)/(arg/t)    
    return C*exp*np.sin(arg) + x_0*exp*np.cos(arg)

def dxdt(t, A, dxdt_0, omega_0, eta, x_0):
    '''
    yields the velocitiy of the modes as a function of time
    Inputs
    dxdt_0 : float
     initial velocity (scalar)
    omega_0 : float
     central freqeuncy of the mode (scalar).
    eta : float
     damping rate (scalar).
    x_0 : float
     initial displacement (scalar).
    A : float
     amplitude of the driving term (scalar).
    t : ndarray | float
     time.
    Outputs
    dxdt(t) from Chaplin+ 1997
    See also
    lorentzians_fake()
    '''
    exp = np.exp(-eta*t)
    arg = np.sqrt(omega_0**2 - eta**2)*t
    # this is a constant (t dependence is divided out)
    C = (A + dxdt_0 + eta*x_0)/(arg/t)
    
    return (C*exp*arg - x_0*eta*exp)*np.cos(arg) -(x_0*exp*arg + C*eta*exp)*np.sin(arg)
def make_x(t, A, dxdt_0, omega_0, eta, x_0):
    '''
    will actually randomly excite the mode that is requested to be generated 
    '''
    # what is the median sampling time
    ts_sample = np.median(np.diff(t))
    # choose excitation interval to be smaller than the damping timescale
    ts_damp = 1./eta
    # excitation sampling timescale is <fac>x smaller than the damping timescale
    fac = 300.
    interval = int(ts_damp/fac/ts_sample)
    print(ts_damp)
    print(ts_sample)
    print(interval)
    print(omega_0)
    print(eta)
    print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    # make sure that there are enough time points to do excite that often
    if interval < 10:
        raise Exception('sampling time is {:4.2f} days but the damping timescale is {:4.1f} days -- leaving each excitation event to last only {} points. Make the time series finer!'.format(ts_sample, ts_damp, interval))
    if interval > len(t)-1:
        raise Exception('excitation time of {:4.2f} days actually exceeds the baseline of {:6.1f} days Make the time series longer!'.format(ts_damp/fac, t[-1]-t[0]))
    x = t*0.0
    _dxdt = t*0.0
    i = 0
    # start off with the requested intial conditions
    x[0] = x_0
    _dxdt[0] = dxdt_0
    np.random.seed(1)
    print((i+1)*interval, len(t)-1)
    while (i+1)*interval < len(t)-1:
        # choose random damping, with mean equal to the passed amplitude, A
        _A = np.random.normal(A, 1.0)
        _dxdt[(i)*interval+1:(i+1)*interval+1] = dxdt(t[(i)*interval+1:(i+1)*interval+1], _A, _dxdt[i*interval], omega_0, eta, x[i*interval])
        x[(i)*interval+1:(i+1)*interval+1] = xt(t[(i)*interval+1:(i+1)*interval+1], _A, _dxdt[i*interval], omega_0, eta, x[i*interval])
        i += 1
    # fill in the last remaining bit of the time series
    if (i)*interval < len(t) - 1:
        
        # choose random damping, with mean equal to the passed amplitude, A
        _A = np.random.normal(A, 1.0)
        _dxdt[(i)*interval+1:] = dxdt(t[(i)*interval+1:], _A, _dxdt[(i)*interval], omega_0, eta, x[(i)*interval])
        x[(i)*interval+1:] = xt(t[(i)*interval+1:], _A, _dxdt[(i)*interval], omega_0, eta, x[(i)*interval])
    # only excite at the beginning
    # x = xt(t, _A, _dxdt[0], omega_0, eta, x[0])
    print(x)
    return x

def lorentzians_fake(p_mcmc, _params, dnu, data, freq_space=False):
    '''
    Will create fake data in the same way that lorentzians_full will make a model of the data
    Inputs
    freq_space : bool 
     If True, will generate fake data by adding the expected noise to the limit spectrum from lorentzians_full() (Anderson, Duvall, and Jefferies 1990). Otherwise, will create time series according to Chaplin+ 1994... Default False.
    '''

    params = (make_params(p_mcmc, _params, dnu, data)).reshape(p_mcmc.n+1,-1)

    ell = p_mcmc.ell
    if freq_space:
        x = data[0].__array__()
        # below taken from dnu_full_fit.get_model():
        if _params.fix_fwhm:
            m = lorentzians_full(x, params[:-1,0:ell], params[:-1,ell:2*ell], params[-1,-2], params[-1,-2]*0.1, fixed=_params.fixed, fixed_epsilon=_params.fixed_epsilon, wrap=_params.collapsed)
        else:
            m = lorentzians_full(x, params[:-1,0:ell], params[:-1,ell:2*ell], params[-1,-2], params[-1,-1], fixed=_params.fixed, fixed_epsilon=_params.fixed_epsilon, wrap=_params.collapsed)
        return x, -m*np.log(np.random.uniform(size=len(m)))
    else:
        # first make the time series -- assuming K2-like stats
        # 80-day length
        # 30-min cadence

        t = np.arange(0.0, 80., 1./(48.))
        x = t*0.0
        # now choose a random amplitude for all the modes to generate -- will scale them according to the
        

        
        print(params)
        for j in range(p_mcmc.n):
            for pos, amp, fwhm in zip(params[j,0:ell], params[j,ell:2*ell], np.array([1.0, 1.0, 1.0])*params[-1,-1]):

                print(pos)
                print(fwhm)
                print(amp)
                print('****************************************************************')

                eta = muhz2pday(fwhm)
                A = amp
                x_0 = 1.0
                # !!! DEBUG fac of 8 --- can't get correct frequency for some reason...
                omega_0 = muhz2pday(pos)*8.
                dxdt_0 = 0.0
                # add this mode's contribution to the light curve. make_x uses x(t) to randomly excite the mode at regular intervals
                x += make_x(t, A, dxdt_0, omega_0, eta, x_0)


        # now compute a power spectrum from this and return it
        # this is ligted directly from bayes_harvey_test_pymc.py
        # want it so that Nf*2 is power of 2
        Nf = len(t)
        # !!! this will not result in a criticall-sampled spectrum.........
        num = 2.**(np.ceil(np.log2(Nf)))
        # !!! DEBUG setting number of freqs to be those for critically-sampled spec
        num = Nf
        f_nyq = (48.)/2.0
        # make the flux be in parts
        x = x/np.median(x)
        import pylab as plt
        # print t
        # print x
        # plt.clf()
        # plt.scatter(t, x)
        # plt.savefig('ts.png', format='png')
        # plt.clf()
        from .bayes_harvey_test_pymc import psd_dft
        __f, __p = psd_dft(t, x-1.0, f_nyq, int(Nf), {})
        __f /= muhz2pday(1.0)
        # put back into parts^2/muhz instead of parts^2/pday
        # put power in parts per million
        __p *= 1.e12*muhz2pday(1.0)

        return __f, __p
def make_params(p_mcmc, _params, dnu, data, params=None, debug=False):
    '''
    returns a parameter array of shape p_mcmc.n, p_mcmc.ell that can be used by dnu_full_fit.py
    Inputs
    p_mcmc : MCMC class
     from dnu_full_fit.py
    _params : Pypar class
     from params_class.py, e.g.
    dnu : float
     dnu
    data : 2, N array
     contains time ([0,:]) and power ([1,:]) for the numax +/- 3\times dnu region being fitted.
    params : ndarray 
     result of best_dnu_entropy(). will only be used if _params.use_corr_dnu. Default None.
    '''
    print(_params)
    # fiducial values are given behind the comment, with rounding to the nearest muhz in order of 0, 1, 2
    ell = p_mcmc.ell
    deltanu = np.array([0.0, -0.025, 0.121, 0.282])[0:ell]
    deltanu_const = np.array([0.0, 0.0, 0.047, 0.16])[0:ell]
    diffs =    deltanu*dnu + deltanu_const #np.array([5.0, 2.1, 4.1])#*0.2
    # find the lowest frequency that corresponds to f mod dnu = zero, which is the reference point for epsilon.
    # expected frequency for that to occur at is...
    # don't want to actually start modeling the outermost orders, but find the lowest order symmetric around numax

    print(_params)
    if not _params.collapsed:
        # JCZ 201017
        # if want to be competely consistent with the collapsed method, the radial order position should correspond to the lowerst frequency for which f % dnu is 0... NOT the below line, which is something that I don't quite understand
        f0 = data[0].iloc[len(data[0])//2]
        f_zero = f0 - ((p_mcmc.n//2 == 0) + p_mcmc.n//2)*(dnu - (f0 % dnu ))
        f_zero = data[0].iloc[np.argmin(data[0].__array__() % dnu)]
        i = 0
        while f_zero - i*dnu > np.min(data[0]):
            f_zero -= i*dnu
            i += 1

        # DEBUG !!!
        # setting to be a constant so that the test_grad grid doesn't have ariticial jumps because dnu has changed slightly and now the first (f % dnu) == 0 point has moved by dnu...
        # + 0.5 makes it align with the test in ~/stello/pipeline_paper/compare_logls.py
        # f_zero = data[0].iloc[0] + 0.5

    else:
        # make f_zero be epsilon for the collapsed case.
        # if _params.collpased:
        # 
            # if not _params.use_corr_dnu:
                # f_zero = 0.634 + 0.63*np.log10(dnu)
            # else:
                # f_zero = params[0]
        f_zero = 0.
    # try to fit up to three orders
    n = p_mcmc.n
    xs = np.zeros(shape=(n, ell))
    amps = np.zeros(shape=(n,ell))
    # print params[0]

    for _n in range(n):
        amps[_n,:] =  np.array([1.0, 0.5, 0.8, 0.0])[0:ell] #np.array([1.0, 0.5, 0.8])#*0.5
        # the zeroth amp is by definition 1 so it is used for epsilon, this is the expected from cosaro+2012

        if not _params.use_corr_dnu:
            amps[_n,0] = 0.634 + 0.63*np.log10(dnu)
            xs[_n,0] =  0.634 + 0.63*np.log10(dnu)
        else:
            # use the epsilon from the correlation solution, if using it for dnu. otherwise the initial guess could be very off if usin ghte dnu from correlation fit and the epsilon from cosaro+2012 scaling relations, above
            # JCZ 210916
            # added the small factor because an error will be thrown if it so happens that the grid chooses epsilon of 0.
            amps[_n,0] = params[0] + 1e-2
            xs[_n,0] = params[0] + 1e-2
            # print params[0]
            # ggg
        
        diffs = xs[_n,:]*0.
        if ell >= 1:
            diffs[1] = -0.025*dnu
        if ell >= 2:
            diffs[2] = (0.121)*dnu + 0.047
        xs[_n,1] = dnu/2. + xs[_n,0] - diffs[1]
        xs[_n,2] = xs[_n,0] - diffs[2]

        
        # xs are the absolute positions of the peaks, so need to go from diffs to positions
        # xs[_n,1] = diffs[1] + xs[_n,0] - dnu/2.

        # if xs[_n,2] < 0.:
        #     xs[_n,1] = xs[_n,0] + dnu/2. - diffs[1]
        # xs[_n,2] = xs[_n,0] - diffs[2]
        # if xs[_n,2] < 0:
        #     xs[_n,2] = xs[_n,0] + dnu - diffs[2]

        xs[_n,:] = xs[_n,:] + dnu*_n
    if _params.collapsed:
        xs = xs % dnu


    # JCZ 201017
    # the problem here is that in the case you want to do fixed, the radial amplitude slot is used for the radial position, and that leads to the amplitude blowing up... because i've started using the radial amplitude, allowing it to float...
    # so removed the if condition and made it the case for everything... !!! still need to figure out how to get epsilon to be a free thing
    xs += f_zero

    # print params[0]
    
    # if not _params.fixed:
    amps[:,0] = 1.0
    # else:
        # pass
    # !!! DEBUG
        # amps[:,0] = [dnu*a for a in range(len(amps[:,0]))] + f_zero + amps[:,0]
    # print xs
    # ggg

    ################################################################

    # first we need to set up 
    
    # we flatten out the parameters for emcee's sake
    start = np.zeros(shape=(p_mcmc.nvariables,))

    start[0:p_mcmc.nvariables - p_mcmc.ell*(2)] =     np.hstack((xs, amps)).reshape(-1)
    start[-p_mcmc.ell*(2):-2] =start[0:p_mcmc.ell*(2)-2]


    start[-2] = dnu
    start[-1] = dnu**2*0.0002
    print('starting with following params:')

    # JCZ 160218
    # the first of the extra position arguments is now used as a fwhm
    # JCZ 170218
    # and the second one as the quad. fwhm
    
    start[-p_mcmc.ell*2] = start[-1]*2.0
    start[-p_mcmc.ell*2+1] = start[-1]
    # JCZ 180218
    # constant bg to model the noise.
    start[-p_mcmc.ell*2+2] = 0.05
    # JCZ 191217
    # the following is for making sure that the f_zero thing makes it so that both the unfolded and folded cases actually ine up with that model.
    if debug:
        import pylab as plt
        from .dnu_full_fit import get_model
        _f, _p = lorentzians_fake(p_mcmc, _params, start.reshape(p_mcmc.n+1,-1)[-1,-2], data, freq_space=True)
        
        plt.clf()
        x = data[0]
        plt.plot(_f, _p)
        plt.plot(x, get_model(x, (start.reshape(p_mcmc.n+1, -1)), p_mcmc, fixed=_params.fixed, fixed_epsilon=_params.fixed_epsilon, wrap=False), label='fit')
        plt.savefig('testing_no_wrap.png', format='png')
        plt.clf()
        x = data[0] % dnu
        _params.collapsed = True
        start =  make_params(p_mcmc, _params, dnu, data, params=params, debug=False)
        _params.collapsed = False
        plt.plot(x, _p)
        plt.plot(x, get_model(x, (start.reshape(p_mcmc.n+1, -1)), p_mcmc, fixed=_params.fixed, fixed_epsilon=_params.fixed_epsilon, wrap=True), label='fit')
        plt.savefig('testing_wrap.png', format='png')
    return start
#condorized @profile
def lorentzians_full_ell(f, xs, amps, dnu, fwhm, l, fixed=False, epsilon=0., fixed_epsilon=False, wrap=False):
    '''
    returns the model for just one mode, across all orders. otherwise the same as lorentzians_full(), with epsilon describing the shift of the radial peak from mod nu/dnu = 0.0., and can go up to 1.0. 
    all the peaks' locations and amplitudes are expressed wrt the radial order.
    the free parameters are dnu, epsilon, xs, and amps
    Inputs
    f : frequency array, in muhz
this will be transformed into mod nu/dnu within this program
    xs : ndarray [n, ell]
    locations of the peaks
    amps : ndarray [n, ell]
    amplitudes of the peaks, relative to the highest, of size [n, ell]
    power : power array normalized to the guess as to the radial peak
    fwhm : ndarray [n, ell]
    fwhm of the smoothed gaussian excess. this determines the width of the lorentzian profiles of each mode
    epsilon : float
    guess for the offset of the radial mode
    x0 : int
    location in the f array of the HIGHEST peak (l = 0)
    l : int
     which mode to return -- 0 for radial, 1 for dipole, 2 for quadrupole
    [ ell : int ]
    how many modes to fit. only 3 for K2
    [ wrap : bool ]
    if True, the model will be that of a collapsed echelle diagram and <f> is assumed to take a maximum value of dnu.
    Outputs
    l : ndarray
     the total lorentzian profile (sum of <p.ell> components)
    [ fixed : bool ]
    if true, the absolute difference between peaks is fixed
    [ fixed_epsilon : bool ]
     if True, only consider the first order radial amplitude as the epsilon, rather than allowing each order
    to have its own epsilon. only valid if <fixed> is True.
    Notes
       
    '''

    # how many ms are being fit
    ell = xs.shape[1]
    # how many orders, n, are being fit
    n = xs.shape[0]

    # first we want to normalize the peaks to the zeroth peak, which we will label as the second highest peak
    # width of the lorentzians as a fraction of the FWHM
    # JCZ 160218
    # DEBUG !!!
    if type(fwhm) == type(np.float64(1.0)):
        w = np.array([0.05, 0.1, 0.05, 0.0])*fwhm
    else:
        w = np.array([0.0]*4)
        w[0] = fwhm[0]
        w[1] = fwhm[1]
        # JCZ 170218
        w[2] = fwhm[2]

    lor = f*0.
    # sum over orders, n
    for _n in range(n):
        # sum over each mode in a given order, _n
        if fixed :
            diffs = xs[_n,:]*0.
            if ell >= 1:
                diffs[1] = 0.025*dnu
            if ell >= 2:
                diffs[2] = -(0.121)*dnu - 0.047
            if fixed_epsilon:
                epsilon = amps[0,0] + _n*dnu
            else:
                epsilon = amps[_n,0]
            _amps = amps[_n,list(range(ell))]
            _amps[0] = 1.0
            _xs = diffs*0.
            _xs[1] = -dnu/2. + diffs[1]
            _xs[2] =  - diffs[2]

        else:
            # prepend the radial amplitude, which is fixed to 1 in the case where one allows the small freq
            # separations to vary (i.e., fixed = False)
            # JCZ 170517
            # think this was wrong. amps was actually [a1, a2, dnu] instead of [a0, a1, a2]
            # _amps = np.hstack(([1.], amps[_n,:]))
            _amps = amps[_n,:]

        # print l
        if fixed:        
            pos = (_xs[l]+epsilon)
            lor += lorentzian(f, pos, _amps[l], w[l])
        else:
            pos = xs[_n,l]
            # print xs[_n,l]
            # print l
            # print '^^^'
            lor += lorentzian(f, pos, _amps[l], w[l])

        if wrap:
            # wrap around to the other side
            if pos > dnu/2:
                lor += lorentzian(f, pos-dnu, _amps[l], w[l])
            else:
                lor += lorentzian(f, pos+dnu, _amps[l], w[l])

    return lor



#condorized @profile
def dnu_e(numax):
    '''
    given a numax, will return the expected value of dnu
    '''

    # dnu = A\times_max^{dnu_numax_slope}
    # the fitting is done in log-log (base 10) so this becomes:
    # \log dnu = dnu_numax_slope\times\log \nu_{max} + \log A
    # from stello+ 2009 : slope = 0.772 and intercept = \log A = \log 0.263
    dnu_numax_slope = 0.772
    dnu_numax_intercept = np.log10(0.263) 
    dnu_e = 10.**(dnu_numax_slope*np.log10(numax) + dnu_numax_intercept)
    return dnu_e

#condorized @profile
def f2echelle_df(data, dnu, SNR=None, median=False, mean=False, sum=True, error=False):
    '''
    takes frequency and power arrays and turns it into a collapsed echelle diagram
    USES PANDAS TO DO THIS
    Inputs
    f : ndarray
    frequency array
    p : ndarray
    [ median : bool ]
     if true, use median to collapase the echelle diagram
    [ mean : bool ]
     if true, use mean to collapse the echelle diagram
    [ sum : bool ]
     if true, use sum of amplitudes in each bin to collapse the echelle diagram
    [ error : bool ]
     if true, returns a third element of the returned tuple being the std of the amplitude bins
    power array

    
    '''
    # first calculate the mod thing
    
    fres = (data[0].iloc[1]-data[0].iloc[0])
    data['mod'] = data[0] % dnu
    n = int(dnu/fres)
    cut_data = pd.cut(data['mod'], int((np.max(data['mod']) - np.min(data['mod']))/fres))
    group = data.groupby(by = cut_data)
    if mean:
        func = np.mean
    if median:
        func = np.median
    if sum:
        func = np.sum
    middle = group.aggregate(np.median)
    result = group.aggregate(func)
    maxamp = np.max(result[1])
    if error:
        error = group.aggregate(np.std)
        return middle['mod'], result[1]/np.max(result[1]), error[1]/maxamp
    # normalize the collapsed echelle diagram
    return middle['mod'], result[1]/maxamp


#condorized @profile
def f2echelle_slice(data, dnu, SNR=None, median=False, mean=False, sum=True, error=False):
    '''
    takes frequency and power arrays and turns it into a collapsed echelle diagram
    this is an alernate to f2echelle_df to make the calculation faster by using
    ndarray slices, ASSUMING FREQUENCIES ARE EQUALLY-SPACED
    USES PANDAS TO DO THIS
    Inputs
    f : ndarray
    frequency array
    p : ndarray
    [ median : bool ]
     if true, use median to collapase the echelle diagram
    [ mean : bool ]
     if true, use mean to collapse the echelle diagram
    [ sum : bool ]
     if true, use sum of amplitudes in each bin to collapse the echelle diagram
    [ error : bool ]
     if true, returns a third element of the returned tuple being the std of the amplitude bins
    power array

    
    '''
    amps = data[1].__array__()
    freqs = data[0].__array__()

    fres = (freqs[1] - freqs[0])
    # bin size in f mod dnu space
    bin_mod = fres

    # how many indices to go in frequency space to go dnu
    onewrap = int(dnu / fres)
    # how many frequency elements to lump into each bin_mod bin
    stride = int(bin_mod/fres)
    # the first element in frequency will probably not evenly divide
    # into dnu. so this indicates how to order the for loop below such
    # that the first element in the collapsed echelle diagram histogram
    # corresponds to frequency mod dnu = 0
    offset = int((freqs[0] % dnu)/float(bin_mod))
    
    

    # how many bins in the collapsed echelle diagram there will be
    n_bins = int(dnu/bin_mod)

    # middle values of the collapsed echelle diagram bins
    middle = np.linspace(0., dnu, num=n_bins)# + bin_mod*0.5

    # the metric computed in each bin to be the collapsed echelle diagram
    result = np.zeros(shape=(n_bins,))
    if mean:
        pass
            
    if sum:
        for j in range(n_bins):
            result[j] += np.sum(amps[stride*j::onewrap])
    if median:
        pass
    
    maxresult = np.max(result)
    result = np.roll(result, offset)
    
    ind_arr = list(range(len(amps)))
    if error:
        std = np.zeros(shape=(n_bins,))
        for j in range(n_bins):
            std[j] = np.std(amps[stride*j::onewrap])
        std = np.roll(std, offset)
        return middle, result/maxresult, std/maxresult
    return middle, result/maxresult
#condorized @profile
def f2echelle2d(f, mod, amp, fdirect='./', fname='echelle' ):
    '''
    returns a 2D array with (f mod dnu, f) that is smoothed to be pretty.
    Inputs
    f : ndarray
     frequencies
    mod : ndarray
     frequency mod dnu
    amp : ndarray
     power at each frequency
    '''

    # for now just do a dumb iteration because i don't know how best to do the 2d hist
    # for 

    from matplotlib import image
    pl = Plot(fdirect=fdirect, fname=fname, ext='eps', save=True)
    ax = pl.ax
    #im = image.NonUniformImage(ax)#, interpolation='bilinear')




    f_res = np.median(f[1:] - f[0:-1])
    # make histogram resolution be <fac_f>th the frequency resolution of the array, f
    # JCZ 260717
    # makes the resolution equal to dnu, which is the smallest it can be...
    fac_f = float(int(np.max(mod)/f_res))
    # resolution in f mod dnu 
    # how much less should mod resolution be than frequency resolution
    fac_mod = 1.0
    N_mod = np.max(mod)/(f_res)/fac_mod
    yedges = np.linspace(np.min(mod), np.max(mod), num=int(N_mod))

    
    N_f = (f[-1] - f[0])/f_res/fac_f
    _f = np.linspace(f[0], f[-1], num=int(N_f))
    xedges = _f[:-1] - (_f[1:] - _f[0:-1])*0.5

    if len(xedges) < 2 or len(yedges) < 2:
        print('dnu was too small. not making 2D echelle diagram.')
        return -1
    # ignore all things within factor, <fac_noise> of the noise
    amp_noise = np.median(amp)
    fac_noise = 0.#1.5
    print(amp)
    print(amp_noise*fac_noise)
    _amp = amp.copy()
    _amp[_amp < amp_noise*fac_noise] = 0.
    print(_amp)
    H, xedges, yedges = np.histogram2d(f, mod, bins=(xedges, yedges), weights=_amp)
    print(np.max(H))
    xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])

    extent = (ycenters[0], ycenters[-1], xcenters[0], xcenters[-1])
    import pylab as plt    
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=np.min(H), vmax=np.max(H))
    
    cmap = mpl.cm.get_cmap('viridis')
    # JCZ 260717
    # oversample by <ofac>
    ofac = 5.
    _f = np.linspace(f[0], f[-1], num=int(N_f)*ofac)
    xedges_hi =  _f[:-1] - (_f[1:] - _f[0:-1])*0.5
    yedges_hi = np.linspace(np.min(mod), np.max(mod), num=int(N_mod)*ofac)
    xcenters_hires = xedges_hi[:-1] + 0.5 * (xedges_hi[1:] - xedges_hi[:-1])
    ycenters_hires = yedges_hi[:-1] + 0.5 * (yedges_hi[1:] - yedges_hi[:-1])
    print(xcenters_hires.shape, ycenters_hires.shape)
    y, x = np.meshgrid(ycenters_hires, xcenters_hires, indexing='ij')

    from scipy.interpolate import RegularGridInterpolator

    print(np.max(y), np.max(x), np.max(ycenters), np.max(xcenters))
    print(np.min(y), np.min(x), np.min(ycenters), np.min(xcenters))
    a = np.vstack((x.flatten(), y.flatten())).T.reshape(-1,2)
    print(a)
    print(x.flatten())
    print(y.flatten())
    
    print(a.shape)
    print(H.shape)
    print(ycenters.shape)
    print(xcenters.shape)
    H_hires = RegularGridInterpolator((xcenters, ycenters), H, bounds_error=False, fill_value=0.0)(a)
    print(H_hires.shape)
    extent = (ycenters_hires[0], ycenters_hires[-1], xcenters_hires[0], xcenters_hires[-1])
    H_hires = H_hires.reshape((len(ycenters_hires), len(xcenters_hires)),order='F')
    ax.imshow(H_hires, extent=extent,aspect=0.5, cmap=cmap, norm=norm)
    # this is slower
    ax.set_xlabel(r'$\nu$ mod $\Delta \nu$ [$\mu$Hz]'); ax.set_ylabel(r'Frequency [$\mu$Hz]')
    #plt.pcolormesh(H, cmap=plt.cm.Blues)
    pl.wrapup()

    return
