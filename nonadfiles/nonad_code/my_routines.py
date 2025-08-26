
__version__ = "0.5"
'''
0.1
-- added imports() to check what and what has not been imported
-- then use this in new function sigma_to_percent which converts z-scores into percentages.
0.2
-- find_nearest_val and find_nearest_ind were switched
0.3
-- added pass statements so that when condorized, will behave
0.4
-- added smooth function
0.5
-- added fullprint from http://stackoverflow.com/questions/1987694/print-the-full-numpy-array
'''
import pandas as pd
import numpy as np
import os
from astropy.io import fits
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized from ipdb import set_trace as pause

class empty(object):
    def __init__(self):
        pass


import types
def fullprint(*args, **kwargs):
  from pprint import pprint
  import numpy
  opt = numpy.get_printoptions()
  numpy.set_printoptions(threshold='nan')
  pprint(*args, **kwargs)
  numpy.set_printoptions(**opt)
def fits2pd(fitsfile, hdu=1):
    '''
    given an inputs fits file, will return pandas dataframe of that data.
    Inputs
    hdu : int
     which HDU to access. Default 1.
    fitsfile : str
     fits file with at least a column called <match_on> (default 'epic'). Assumes only has on HUD
    Outputs
    fits_df : pd DataFrame
     pd DataFrame
    '''
    fits = fits.open(fitsfile)[hdu].data
    dict = {}
    for key in fits.names:
        # JCZ 070817
        # should not affect normal operations, but will successfully add data with more than one dimension.
        try:
            entry = fits.field(key).byteswap().newbyteorder()
        except:
            print('couldnt do {}'.format(key))
            continue
        entry_dim = len(entry.shape)
        n_dim = entry.shape[np.argsort(entry.shape)[0]]
        dim_suffixes = np.arange(n_dim).astype(str)
        # JCZ 251017
        # for some reason now that i have added lengths to string columns, there are some that have no dimension? so need to skip over those...
        try:
            dim_suffixes[0] = ''
        except:            
            continue

        # add the dimensions one-by-one, along the short dimension -- for the case of two dimensions. for more dimensions, the array is simply flattened
        if entry_dim > 1 and entry_dim < 3:

            for dim, suffix in zip(range(n_dim), dim_suffixes):
                print('could not add {} to the DataFrame because it had multiple dimensions. breaking it up in to {} entries, instead...'.format(key, np.min(entry.shape)))
                print('adding entry {}'.format(key+suffix))
                if np.argsort(entry.shape)[0] == 0:
                    dict[key+suffix] = entry[dim, :]
                else:
                    dict[key+suffix] = entry[:, dim]
        elif entry_dim == 1:
            dict[key] = entry
    fits_df = pd.DataFrame(dict)
    return fits_df

    
def incremental_std(arr):
    # incrementally calculated std by way of keeping track of 
    # meansq : \sum_i \langle x_i \rangle ^2
    # sqmean : \sum_i \langle x_i^2 \rangle
    # and then use
    # \sigma = \sqrt{sqmean - meansq}

    
    std = []
    meansq = 0.
    sqmean = 0.
    c = 0.
    
    for j in range(n_bins):
        pass
        # c += 1.
        # meansq[j+offset] = meansq[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap])/c
        # sqmean[j+offset] = sqmean[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap]**2)/c
            
    std = np.sqrt(np.array(sqmean - meansq))
    return 9999

def imports():
    '''
    returns imported modules
    '''
    for name, val in globals().items():
        if isinstance(val, types.ModuleType):
            yield val.__name__

def str_ind(array, st, get_bool=False):
    '''
    returns the indices in <array> that contains <st>                                                                     
    Inputs                                                                                                                
    array : ndarray                                                                                                       
    st : str                         
    [ get_bool : bool ]
     if True, returns a boolean array of the same length as <array>, instead of an index array with a length not necessarily equal to the length of <array>. Default False.
    Outputs                                                                                                               
    idx : ndarray                                                                                                         
    containing indices corresponding to <array> that have <st> in them                                                    
    Notes
     will interpret elements in the passed array that are not strings as strings.
    '''
    if get_bool:
        indices = [True if st in str(x) else False for i,x in enumerate(array)]
    else:
        indices = [i for i,x in enumerate(array) if st in str(x)]
    # JCZ 160118
    # added the conversion to array
    return np.array(indices)

# def strip_ext(file):
#     '''
#     Inputs
#     file : str
#      path that you want to strip of the extension (e.g., /this/is/cool.notcool)
#     Outputs
#     stripped filed : str
#      e.g., /this/is/cool
#     Notes
#     uses os.path to do this
#     does not work on multiple files!
#     '''
#     return os.path.join(os.
def ang2p(ang):
    return 2.*np.pi/ang

def p2ang(p):
    return 2.*np.pi/p

def f2ang(f):
    return f*np.pi*2.

def ang2f(ang):
    return ang/2./np.pi

def f2p(f):
    return 1./f


def p2f(p):
    return 1./p

def isinfnan(x):

    return np.logical_or((np.isinf(x)), (np.isnan(x)))
def delete_inf_nan(x, y=None, ind=False, fill=None, full_ind=False):
    '''
    Returns an array without the values that are NaN or Inf. Optionally can replace these entries with a fill factor, <fill>.
    full ind : bool
     if true, returns a boolean array of the same shape as input array, to be contrasted with the ind option, which returns an int ndarray of
     length equal to or less than input array.
    '''
    if fill is not None:
        if y is not None:
            a = np.where((np.isinf(y)) | (np.isnan(y)))[0]
        else:
            a = np.where((np.isinf(x)) | (np.isnan(x)))[0]
        if ind:
            return a
        x[a] = fill
        
        if y is not None:
            y[a] = fill
            return (x, y)
        else:
            return x
        

    if y is not None:
        a = np.where((~np.isinf(y)) & (~np.isnan(y)))[0]
    else:
        a = np.where((~np.isinf(x)) & (~np.isnan(x)))[0]
    if ind:
        return a


    if full_ind:
        return np.isinf(x) | np.isnan(x)

    if y is not None:
        return (x[a], y[a])
    else:
        return x[a]

def n_elements(array):
    '''
    Will count the number of elements along the zeroth dimension of a list or numpy array.
    '''
    if type(array) == type([]):
        return len(array)
    if type(array) == type(np.array([])):
        return array.shape[0]
    return 0

def plot(x, y=None):
    '''
    Automatically shows a plot, with interactive window on
    '''
    #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized import pylab as plt
    if y is not None:
        #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized plt.plot(x,y); plt.ion(); plt.show()
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized pause()
        pass
    else:
        #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized plt.plot(x); plt.ion(); plt.show()
        pass
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized pause()

def find_nearest_ind(array, value):
    '''
    taken from unutbu http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    '''

    ind = (np.abs(array-value)).argmin()
    return ind

def find_nearest_val(array, value):
    '''
    taken from unutbu http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    '''

    return array[find_nearest_ind(array, value)]

def smooth(raw, w, type='boxcar'):
    """
    Inputs
     w : float or int
      triangle:
       w is the FWHM of the triangle
      boxcar:
       w is the width of the boxcar - 2
     
    """

    N = raw.shape[0]

    sm = raw.copy()

    if type == 'boxcar':
        # print w % 2 > 0
        if not w % 2 > 0:
            w = w - 1
            print('changed w from {} to {}.'.format(w+1, w))
        inds = np.linspace(w, N-w, num=N-2*w+1)
        w = int(w)-2

        w_i = np.ones(w+2)

    if type == 'triangle':
        # base of the triangle in units of indices
        inds = np.linspace(w, N-w-1, num=N-2*w-1)
        b = w*2.
        w = int(b) - 1
        center = int((w+1)/2.)
        w_i = np.zeros(shape=(w+2))
        for i in range(w+1+1):

            w_i[i] = (b - 2*np.abs(center - i))/b

        
    for ind in inds:

        sm[ind] = np.sum(w_i*raw[ind - (w+1)/2 : ind + (w+1)/2+1])/np.sum(w_i)

        

    return sm


def sigma_to_percent(sigma):
    '''
    uses z-score to convert to fraction
    '''
    if 'norm' not in imports():
        import scipy.stats
        
    return (scipy.stats.norm.cdf(sigma)*2.) - 1.0
