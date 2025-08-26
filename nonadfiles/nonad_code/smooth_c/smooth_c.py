import pandas as pd
from pudb import set_trace as pause
#from memory_profiler import profile
__all__ = ['smooth_c']
#@profile
def smooth_c(raw, w, type):
    '''
    Inputs 
    raw : ndarray
     array to be smoothed
    w : int
    width of smoothing
     = std for type == 'gaussian'
     = width of boxcar for type == 'boxcar'
     = FWHM of triangle for type == 'triangle'
    type : str
     method of smoothing : gaussian | boxcar | triangle
    Outputs
    ndarray
     smoothed array, of samem length as input, with the edges not doing weird things.     
    '''
    # pandas doesn't recognize 'triangle', only 'triang'
    if type == 'triangle': type = 'triang'
    series = pd.Series(raw)
    # special case because it requires std

    # JCZ 170321
    # this was deprecated in python
    # if type == 'gaussian':
    #     return pd.rolling_window(series, win_type=type, window=10*len(raw), center=True, min_periods=0, std=w).__array__()
    # else:
    #     return pd.rolling_window(series, win_type=type, window=w, center=True, min_periods=0).__array__()

    if type == 'gaussian':

        return series.rolling(win_type=type, window=10*len(raw), center=True, min_periods=0).mean(std=w)
    else:
        return series.rolling(win_type=type, window=w, center=True, min_periods=0).mean()

#@profile
def smooth_std(raw, w, type):
    '''
    Inputs 
    raw : ndarray
     array to be smoothed
    w : int
    width of smoothing
     = std for type == 'gaussian'
     = width of boxcar for type == 'boxcar'
     = FWHM of triangle for type == 'triangle'
    type : str
     method of smoothing : gaussian | boxcar | triangle
    Outputs
    ndarray
     smoothed array, of samem length as input, with the edges not doing weird things.     
    '''
    # pandas doesn't recognize 'triangle', only 'triang'
    if type == 'triangle': type = 'triang'
    series = pd.Series(raw)
    if type == 'gaussian':
        # FIX THIS DOESN"T ACCEPT WIN_TYPE ARG LIKE ROLLING_WINDOW DOES SO THIS ISN"T A GAUSSIAN WINDOW
        return pd.rolling_std(series,  window=w, center=True, min_periods=0).__array__()
    else:
        return pd.rolling_std(series, window=w, center=True, min_periods=0).__array__()

