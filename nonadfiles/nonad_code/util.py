from __future__ import print_function
from pudb import set_trace as pause
import numpy as np
import re
import pandas as pd
from shapely import geometry
from .setup_graphics import *
from astroML.plotting import hist
from .my_routines import sigma_to_percent
import pylab as plt

from astropy.io import ascii,fits
# import matplotlib as mpl

# changing to latex
# JCZ 290121
# plt.style.use('jcz_paper_latex')

# mpl.rc('font', family='serif', serif='cm10')
# mpl.rc('text', usetex=True)
# mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

from astropy.io import ascii,fits

def remove_dredge_up(ls, teffs):
    # require a pad of 5 evo. points
    pad = 1
    ind_old = 0
    inds = [0]
    assert len(ls) == len(teffs), 'ls and teffs must be the same length for remove_dredge_up to work.'
    for i in np.arange(len(ls)-pad)+1:

        if (ls[i] <= ls[i:i+pad]).all() and (teffs[i] >= teffs[i:i+pad]).all(): # and (ls[i] > ls[inds[-1]]):
            # require that the temperature always be less than all the previous temperatures, not just the immediately preceding teff.
            # this is to prevent a double-valued bump part.
            if (teffs[0:i] >= teffs[i:i+pad]).all():
                ind_old_ = i + pad - 1
                inds.append(ind_old_)
        # this is to allow the dip in the hook
        elif (teffs[i] >= teffs[i:i+pad]).all():# and (ls[i] > ls[inds[-1]]):
            # require that the temperature always be less than all the previous temperatures, not just the immediately preceding teff.
            # this is to prevent a double-valued bump part.
            if (teffs[0:i] >= teffs[i:i+pad]).all():

                ind_old_ = i + pad - 1
                inds.append(ind_old_)


    return np.array(inds)


def weighted_avg(x, xerr, get_err=True, min_n=0):
    '''
    Inputs
    [ min_n ] : int
     will only return averages for things with at least min_n entries. Default 0.
    x : ndarray
    xerr : ndarray
    weighted avg and optinoally the standard unc. on that. will automatically deal with inf and nan entries in x and xerr and only use valid pairs of (x, xerr) to compute xavg and xavg_err.
    get_err : bool 
     will return (xavg, xavg_err) if True (otherwise will return just xavg). Default True. 
    '''
    print('in wfun')
    # cannot be lists
    assert ~isinstance(x, list)
    assert ~isinstance(xerr, list)
    # input must be arrays
    assert (~np.isscalar(x))*(~np.isscalar(xerr))
    assert (len(x) > 0)
    assert (len(xerr) > 0)
    # can't be ints
    assert (type(x[0]) != np.int)
    assert (type(xerr[0]) != np.int)
    assert (type(x[0]) != np.int64)
    assert (type(xerr[0]) != np.int64)
    
    _x = x[np.where( (~np.isnan(x))*(~np.isnan(xerr)))[0]]
    _xerr = xerr[np.where( (~np.isnan(x))*(~np.isnan(xerr)))[0]]
    # print(_x)
    # print(_xerr)
    # print(x)
    # print(xerr)
    xavg = np.sum(_x/_xerr**2)/np.sum(1./_xerr**2)
    # print(xavg)
    # print('out of wfun')






    if len(_x) < min_n:
        if get_err:
            return(np.nan, np.nan)
        xavg = np.nan
    if get_err:
        xavg_err = 1./np.sqrt(np.sum(1./_xerr**2)  )
        return(xavg, xavg_err)
    else:
        return xavg


def weighted_std(x, xerr, min_n=0):
    '''
    Inputs
    [ min_n ] : int
     will only return averages for things with at least min_n entries. Default 0.
    x : ndarray
    xerr : ndarray
    weighted std that. will automatically deal with inf and nan entries in x and xerr and only use valid pairs of (x, xerr) to compute xstd
    '''
    print('in wfun')
    # cannot be lists
    assert ~isinstance(x, list)
    assert ~isinstance(xerr, list)
    # input must be arrays
    assert (~np.isscalar(x))*(~np.isscalar(xerr))
    assert (len(x) > 0)
    assert (len(xerr) > 0)
    # can't be ints
    assert (type(x[0]) != np.int)
    assert (type(xerr[0]) != np.int)
    assert (type(x[0]) != np.int64)
    assert (type(xerr[0]) != np.int64)
    
    _x = x[np.where( (~np.isnan(x))*(~np.isnan(xerr)))[0]]
    _xerr = xerr[np.where( (~np.isnan(x))*(~np.isnan(xerr)))[0]]
    # print(_x)
    # print(_xerr)
    # print(x)
    # print(xerr)
    xavg = np.sum(_x/_xerr**2)/np.sum(1./_xerr**2)
    xstd = np.sqrt(np.average((_x - xavg)**2, weights=1./_xerr**2))
    
    if len(_x) < min_n:
        xstd = np.nan

    return xstd


def add_contours(x, y, color='grey', label='Galaxia', ax=None, side=10, range=None, contour_labels=True, reverse=True, zorder=None, sigmas=[1.0,2.0,3.0], plot_legend=False, reverse_y=False, smooth=False, lw=2.0, alpha=1.0, filled=False, get_outside=False, get_inside=False):
    '''
    [ get_outside : bool ]
     if True, will return ax, (x, y), where x and y are the points lying outside the contours. Default False. Cannot have both get_inside and get_outside.
    [ get_inside : bool ]
     if True, will return ax, (x, y), where x and y are the points lying inside the contours. Default False. Cannot have both get_inside and get_outside.
    # [ contour_labels : bool ]

     if True, will put on the contours the labels corresponding to the percentage of points within each contour. Default True.
    NB: even if range is None, if <ax> is passed, range will be taken from ax so that the percentiles plotted in the contours will represent the percent of the population IN THE RANGE PLOTTED IN AX AS add_countours() SEES IT. If ax is modified after add_contours, then the percentiles will no longer represent the percentiles for the new plotted range.
    [ alpha : float ]
     alpha for the cotours. Default 1.0.
    [ lw : float ]
     max linewidth for the smallest sigma plotted.
    [ reverse_y : bool ]
     will only reverse y-axis, not both (reverse does both).
    [ reverse : bool ]
     if True, the axis limits are considered to be reversed; like in an HR diagram or Kiel diagram.
    [ zorder : int | float ]
     very negative will make the contours appear under everything else in the plot; very positive will make the contours appear over everything else in the plot.
    [ sigmas : list of floats ]
     sigmas corresponding to which contours will be plotted, from the inside out (e.g., [1.0, 2.0, 3.0] will plot 68, 95, and 99.7% contours). Default [1.0, 2.0, 3.0].
    [ filled : bool ]
     If True, fill the contours. should also provide a colors thing. Default False. If True, color can be a list of length sigmas to speify the filling colors.
    [ color : str or list of str ]
     If list, then the lines/fillings (depending on if filled is False or True) will be different colors. Otherwise, they will all be the same color.
    
    
    '''
    # JCZ 150920
    # removing this
    # mplt.style.use('jcz_paper_latex')    
    if get_inside and get_outside:
        raise Exception("cannot set both get_inside and get_outside. try doing two different calls.")
    if ax is None:
        ax = plt.figure().gca()
        resize = False
    else:
        xmin, xmax, ymin, ymax = ax.get_xlim()[0], ax.get_xlim()[1], ax.get_ylim()[0], ax.get_ylim()[1]
        # JCZ 140118
        # need to re-apply the original limits on the plot because the contours will slightly change the limits.
        resize = True
        # NOTE THAT EDGES ARE MIXED UP BECAUSE THIS IS KIEL DIAGRAM AND THINGS ARE BACKWARD
        if reverse:
            range = [[ymax, ymin], [xmax, xmin]]
        elif reverse_y:
            range = [[ymax, ymin], [xmin, xmax]]
        else:
            range = [[ymin, ymax], [xmin, xmax]]

    H, xedges, yedges = np.histogram2d(y, x, range=range, bins=(side, side))
    # NOTE THAT EDGES ARE MIXED UP BECAUSE THIS IS KIEL DIAGRAM AND THINGS ARE BACKWARD
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    # JCZ 220117
    # MUST BE FLOATS
    # the following levels work well for 10x10 grid for galaxia
    # levels = (50., 200., 250., 300., 500.)

    # JCZ 080217
    # automatically choose levels such that they represent fixed fractions
    levels = []
    # JCZ 191217
    # changing the levels to correspond to sigma
    fractions_desired = np.array(sigmas)
    fractions_desired = sigma_to_percent(fractions_desired)
    fractions_desired = np.sort(np.array(fractions_desired))[::-1]

    
    n_add = 0.05*int(np.sum(H)) + 1

    if smooth:
        from scipy import interpolate
        from scipy import ndimage
        from scipy.ndimage.filters import gaussian_filter
        # x, y = np.meshgrid(xedges, yedges)
        
        H_fine = ndimage.zoom(H, 2)
        # then smooth w gaussian
        sigma = np.min(np.abs(np.array([xedges[1] - xedges[0], yedges[1] - yedges[0]])))*4.0
        # print('smoothing w sigma of {}'.format(sigma))
        H = gaussian_filter(H, sigma)

    
    for f in fractions_desired:
        f0 = 1.
        # changed from 0 to -2.0 to try to get a sigma = 0 level...
        level = -2.0
        while f < f0:
            level += 0.1#n_add
            f0 = np.sum((H[np.where(H > level)]))/np.sum(H)
        levels.append(level)

    # print('levels:')
    # print( levels)
    # print('H:')
    # print (H)
    fracs = [np.sum((H[np.where(H > level)]))/np.sum(H) for level in levels]

    fmt = {}
    # express contour labels as percentage of the total range of points in this part of the HR diagram that are wihtin the contour
    strs = ['{:3.1f}%'.format(frac*100.) for level, frac in zip(levels, fracs)]
    scales = np.linspace(lw/2.0, lw, num=len(levels))
    if filled:
        cont = ax.contourf
        linewidths=None
    else:
        cont = ax.contour
        linewidths=4*scales

    a = cont(H, levels, origin='lower',extent=extent, colors=color, label=label, linewidths=linewidths, zorder=zorder, alpha=alpha, extend='max')

    # this code is from
    # https://gis.stackexchange.com/questions/99917/converting-matplotlib-contour-objects-to-shapely-objects
    # answer of kpenner

    if get_outside or get_inside:
    #     for i, collection in enumerate(a.collections):
    #         for path in collection.get_paths():
    #             if path.to_polygons():
    #                 for npoly, polypoints in enumerate(path.to_polygons()):
    #                     poly_lons = polypoints[:, 0]
    #                     poly_lats = polypoints[:, 1]
    #                     # be careful with the following---check to make sure
    #                     # your coordinate system expects lat first, lon second
    #                     poly_init = geometry.Polygon([coords for coords in zip(poly_lats, poly_lons)])
    #                     if poly_init.is_valid:
    #                         poly_clean = poly_init
    #                         if npoly == 0:
    #                             poly = poly_clean
    #                     else:
    #                         poly_clean = poly_init.buffer(0.)
    #                         if npoly == 0:
    #                             poly = poly_clean
    #                         else:
    #                             poly = poly.difference(poly_clean)
        # the above one doesn't work for some reason... but the below does. from a diff. answer in the same url above, but i added the poly.union thing. 
        for i in np.arange(len(a.collections)):
            p = a.collections[i].get_paths()[0]
            v = p.vertices
            _x = v[:,0]
            _y = v[:,1]
            if i == 0:
                
                poly = geometry.Polygon([(i[0], i[1]) for i in zip(_x,_y)])
            else:
                poly = poly.union(geometry.Polygon([(i[0], i[1]) for i in zip(_x,_y)]))
        
        points = geometry.MultiPoint(tuple(zip(x,y)))
        out = []
        count = 0
        for p in points:
            if (not (poly.contains(p))) and get_outside:
                out.append(count)
            elif ((poly.contains(p))) and get_inside:
                out.append(count)
            count += 1
            
        out = (x.iloc[out], y.iloc[out])
    
    fmt = {}
    # express contour labels as percentage of the total range of points in this part of the HR diagram that are wihtin the contour

    for l, s in zip(levels, strs):
        fmt[l] = s
    ax.plot([],[], linestyle='solid', color=color, linewidth=scales[-1]*4, label=label)
    # plt.clabel(a, inline=1, fontsize=12, fmt=fmt)
    if resize:
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        
    if plot_legend:
        ax.legend(loc='upper left')
    if contour_labels:
        plt.clabel(a, inline=1, fontsize=12, fmt=fmt)

    for c in a.collections:
        c.set_linestyle('solid')
    if get_outside or get_inside:
        return ax, out
    return ax

def rename_df(df, bam_dnu=True, bam_numax=True):
    '''
    Inputs
    df : pd DataFrame
     WILL NOT change <df> in place. Only returned pd DataFrame has updated names.
    [ bam_dnu : bool ]
     Default True.
    [ bam_numax : bool ]
     Default True.
    will change the names of the columns of the passed dataframe like:
    a1_sig --> a1_err
    whitenoise_sig --> whitenoise_err
    maxam --> maxamp (if maxam is defined and not maxamp)
    maxamp_sig/maxam_sig --> maxamp_err
    maxam(gaus) --> maxamp(gaus)
    maxam(gaus)sig/maxamp(gaus)sig --> maxamp(gaus)_err
    FWHM(gaus) --> FWHM
    FWHM(gaus)sig --> FWHM_err
    and 
    numax(gaus) --> numax (if bam_numax is True; otherwise, numax --> numax)
    numax(gaus)sig --> numax_err (if bam_numax is True; otherwise, numax_sig --> numax_gaus)
    dnu(gaus) --> dnu "
    dnu(gaus)sig --> dnu_err
    etc.
    Outputs
    df : pd DataFrame
     with new naming scheme.
    '''
    names_old = ['a1_sig', 'b1_sig', 'a2_sig', 'b2_sig', 'maxam', 'maxam(gaus)', 'maxamp(gaus)sig', 'maxam_sig', 'maxam(gaus)sig', 'maxamp_sig', 'whitenoise_sig', 'FWHM(gaus)sig', 'FWHM(gaus)']
    names_new = ['a1_err', 'b1_err', 'a2_err', 'b2_err', 'maxamp', 'maxamp(gaus)', 'maxamp(gaus)_err', 'maxamp_err', 'maxamp(gaus)_err', 'maxamp_err', 'whitenoise_err', 'FWHM_err', 'FWHM']
    if bam_numax:
        names_old += ['numax(gaus)', 'numax(gaus)sig']
    else:
        names_old += ['numax', 'numax_sig']
    names_new += ['numax', 'numax_err']

    if bam_dnu:
        # JCZ 210917
        # sometimes I have named the column dnu_sig(gaus)
        names_old += ['dnu(gaus)', 'dnu(gaus)sig', 'dnu_sig(gaus)']
        names_new += ['dnu', 'dnu_err', 'dnu_err']
    else:
        names_old += ['dnu', 'dnu_sig']
        names_new += ['dnu', 'dnu_err']
    # JCZ 210917
    # need to replace names one-by-one, in case a renaming creates a thing that is then trying to be replaced by another re-naming. this shouldn't be the case in theory, but just to make sure it doesn't happen, i do the for loop and check if the name already exists. if it does, the existing name, <name>, is renamed to <name>_renamed.
    for old,new in zip(names_old, names_new):
        # only rename the existing key if the old one attempting to be replaced exists... (hence the second and clause)
        if new in df.keys() and old in df.keys():
            df = df.rename(columns={new:new+'_renamed'})
            df[new] = df[old].copy()
        else:
            df = df.rename(columns={old:new})
    return df

def python2fortran_format(data, verbose=False):
    '''
    will return the code for the data type for the array, <data>, needed for astropy.io.fits
    if an error is encountered in getting the data type, format is returned as None.
    Examples
    >>> code = python2fortran_format([1., 2.], verbose=False)
    >>> print code
    >>> 'D'
    
    >>> code = python2fortran_format(np.array([1., 2.]), verbose=False)
    >>> print code
    >>> 'D'
    
    '''
    try:
        s = str(data.dtype)
    except:
        try:
            # this clause is if providing a list, in which case need to take an element out to get the type, otherwise the type is 'list'.
            s = str(type(data[0]))
        except:

            try:
                s = str(type(data))
            except:
                if verbose:
                    print('could not get data type by the dtype method or the type() built-in function. returning None.')
                s = ''
                f = None
    
    # from astropy doc:
    # L                        logical (Boolean)               1
    # X                        bit                             *
    # B                        Unsigned byte                   1
    # I                        16-bit integer                  2
    # J                        32-bit integer                  4
    # K                        64-bit integer                  4
    # A                        character                       1
    # E                        single precision floating point 4
    # D                        double precision floating point 8
    # C                        single precision complex        8
    # M                        double precision complex        16
    # P                        array descriptor                8
    # Q                        array descriptor                16

    if 'bool' in s:
        f = 'L'
    elif 'uint8' in s:
        f = 'B'
    elif 'int' in s and '16' in s:
        f = 'I'
    elif 'int' in s and '32' in s:
        f = 'J'
    elif ('int' in s and '64' in s) or 'long' in s:
        f = 'K'
    # JCZ 231017
    # added object as a condition 
    elif 'str' in s or 'S' in s or 'object' in s:
        # JCZ 241017
        # adding the length of the string!!! otherwise, cut off after one string...
        try:
            s = str(len(data[0]))
            print(data)
        except:
            'could not determine length of string. setting to 1.'
            s = 1
        f = s + 'A'
    elif 'float' in s and '32' in s:
        f = 'E'
    elif 'float' in s and '64' in s:
        f = 'D'
    elif 'complex' in s and '32' in s:
        f = 'C'
    elif 'complex' in s and '64' in s:
        f = 'M'
    elif 'complex' in s:
        if verbose:
            print('interpreting type complex as complex128')
        f = 'M'
    elif 'int' in s:
        if verbose:
            print('interpreting type int as int64')
        f = 'K'
    elif 'float' in s:
        if verbose:
            print('interpreting type complex as float64')
        f = 'D'
    else:
        if verbose:
            print('the {} format could not be translated to a fortran one... returning None.'.format(s))
        f = None
    return f
def pd2fits(df, get_cols=False):
    '''
    returns an astropy fits object from a panda dataframe that can then be written to file, e.g.:
    df = pd.DataFrame({'a':[1,2]})
    tbhdu = pd2fits(df)
    tbhdu.writeto('test.fits')
    Inputs
    [ get_cols : bool ]
     if True, return the columns instead of the actual Table object made from those columns. Useful if want to add columns to an existing list of columns, which you will then make into a Table and save as a fits. Default False.
    Notes
    uses python2fortran_format() (defined above) to get the correct Fortran data type to write to file.
    '''
    cols = []

    for key in df.keys():
        cols.append(fits.Column(name=key,array=df[key], format=python2fortran_format(df[key])))
    if get_cols:
        return cols
    cols=fits.ColDefs(cols)
    tbhdu=fits.BinTableHDU.from_columns(cols)
    return tbhdu
def write_fits(filename, tbhdu):
    '''
    will write info in tbhdu to a file, <filename>.
    
    '''
    tbhdu.writeto(filename)

# JCZ 311018
# took out fits2pd and put it instead in my_routines
# !!! debug doesn't work right now because SYDOSU has its own version of my_routines which differs from the one on mac (~/packages/my_routines.py) where I put fits2pd. so need to bring them into alignment.
from .my_routines import fits2pd as _fits2pd
def fits2pd(*args, **kwargs):
    return _fits2pd(*args, **kwargs)

def join_fits_df(fitsfile, df, match_on='epic'):
    '''
    given an inputs fits file and a dataframe, will join the two on <match_on> (default 'epic'). if match_on='epic', 
    then will convert epic to an integer before matching. this could cause problems if epic cannot be cast to an int.
    Inputs
    fitsfile : str
     fits file with at least a column called <match_on> (default 'epic'). Assumes only has on HUD
    df : pd DataFrame
     pd DataFrame with at least a column called <match_on> (default 'epic')
    Outputs
    merged : pd DataFrame
     pd DataFrame via a merge (so columns in common will be named col_x and col_y, where x is the fits one and y is the df one.
    '''
    fits = fits.open(fitsfile)[1].data
    dict = {}
    for key in fits.names:
        dict[key] = fits.field(key).byteswap().newbyteorder()
    fits_df = pd.DataFrame(dict)
    if match_on == 'epic' or match_on == 'KEPLER_ID':
        fits_df[match_on] = fits_df[match_on].astype(int)
    merged = fits_df.merge(df, how='inner', on=match_on)
    return merged

#from memory_profiler import profile
# from http://stackoverflow.com/questions/26783719/efficiently-get-indices-of-histogram-bins-in-python
def mydate():
    '''
    returns a date of the form
    010217
    (1st of February 2017)
    '''
    import time
    return time.strftime('%d%m%y')


def read_ts(file, ignore_err=False):
    '''
    reads a timeseries, no matter if has format time flux error or time flux
    Inputs
    file : str
     time series file to be read. ASCII ONLY.
    [ ignore_err : bool ]
     if True, will only return (time, flux) EVEN IF THERE IS AN ERROR ARRAY TO BE READ IN. Ignored if no error array can be found
    Outputs
    (time, flux, error) | (time, flux)
     tuple of ndarrays
    '''
    try:
        t, f, e = np.loadtxt(file, unpack=True)
        if ignore_err:
            return (t, f)
        else:
            return (t, f, e)
    except:
        t, f = np.loadtxt(file, unpack=True)
        return (t, f)


def write_ts(file, t, f, err=None):
    '''
    write a timeseries, no matter if to have format time flux error or time flux
    Inputs
    file : str
     time series file to be read.
    [ t : ndarray ]
     time series
    [ f : ndarray ]
     flux series
    [ err : ndarray ]
     error on flux
    
    
    '''

    if err is not None:
        np.savetxt(file, np.vstack((t, f, err)).T)
    else:
        np.savetxt(file, np.vstack((t, f)).T)

def muhz2pday(muhz):
    '''
    Inputs
    muhz : float
     frequency in muhz
    Outputs
    pday : float
     frequency in per day
    Examples
    >>> print(muhz2pday(47.))
    >>> 4.0608
    '''
    return 1.e-6*3600.*24.*muhz


def pday2muhz(pday):
    '''
    Inputs
    pday : float
     frequency in per day
    Outputs
    muhz : float
     frequency in muhz
    Examples
    >>> print(pday2muhz(1./20.))
    >>>
    '''
    return 1./(1.e-6*3600.*24.)*pday

def pandas_bin_weighted(df, index_field, value_field, value_err_field, min, max, nbins, min_n):
    '''
    will return a series of weighted averages and their errors for things with at least min_n things in the bin, otherwise that bin will be nan.
    min : float
     the leftmost bin edge to consider
    max : float
     the rightmost bin edge to consider
    nbins : int
     how many bins to make within the bin range defined by min and max. this is done using bins = np.linspace(min, max, nbins)
    value_field : str
     the name of the column that you want the weighted avg of
    value_err_field : str
     uncertainty column corresponding to value_field
    index_field : str
     what to cut on using the bins
    min_n : int
     only bins with min_n or more entires will return non-nan entries for the weighted avg
    Returns
    (avg, rr_on_avg) : tuple of ndarray(nbins)
    
    '''
    ws = []
    es = []
    bins = np.linspace(min, max, nbins+1)
    data_cut = pd.cut(df[index_field].__array__(), bins)
    group = df.groupby(by=data_cut)               
    for s in group:
        # print (s[1][value_field].values)
        # print  (s[1][value_err_field].values)
        # pause()
        # print(weighted_avg(s[1][value_field].values, s[1][value_err_field].values, min_n=min_n, get_err=False))
        try:
            ws.append(weighted_avg(s[1][value_field].values, s[1][value_err_field].values, min_n=min_n, get_err=False))
            es.append(weighted_avg(s[1][value_field].values, s[1][value_err_field].values, get_err=True, min_n=min_n)[1])
        except:
            # pass
            ws.append(np.nan)
            es.append(np.nan)
    return np.array(ws), np.array(es)


def pandas_bin_weighted_std(df, index_field, value_field, value_err_field, min, max, nbins, min_n):
    '''
    will return a series of weighted std for things with at least min_n things in the bin, otherwise that bin will be nan.
    min : float
     the leftmost bin edge to consider
    max : float
     the rightmost bin edge to consider
    nbins : int
     how many bins to make within the bin range defined by min and max. this is done using bins = np.linspace(min, max, nbins)
    value_field : str
     the name of the column that you want the weighted std of
    value_err_field : str
     uncertainty column corresponding to value_field
    index_field : str
     what to cut on using the bins
    min_n : int
     only bins with min_n or more entires will return non-nan entries for the weighted std
    Returns
    std : ndarray(nbins)
    
    '''
    stds = []
    bins = np.linspace(min, max, nbins+1)
    data_cut = pd.cut(df[index_field].__array__(), bins)
    group = df.groupby(by=data_cut)               
    for s in group:
        # print (s[1][value_field].values)
        # print  (s[1][value_err_field].values)
        # pause()
        # print(weighted_avg(s[1][value_field].values, s[1][value_err_field].values, min_n=min_n, get_err=False))
        

        try:
            std = weighted_std(s[1][value_field].values, s[1][value_err_field].values, min_n=min_n)

            stds.append(std)

        except:
            # pass
            stds.append(np.nan)

    return np.array(stds)


def pandas_bin_median(df, index_field, value_field, value_err_field, min, max, nbins, min_n):
    '''
    will return a series of medians and their errors for things with at least min_n things in the bin, otherwise that bin will be nan.
    min : float
     the leftmost bin edge to consider
    max : float
     the rightmost bin edge to consider
    nbins : int
     how many bins to make within the bin range defined by min and max. this is done using bins = np.linspace(min, max, nbins)
    value_field : str
     the name of the column that you want the median of
    value_err_field : str
     uncertainty column corresponding to value_field
    index_field : str
     what to cut on using the bins
    min_n : int
     only bins with min_n or more entires will return non-nan entries
    Returns
    (avg, rr_on_avg) : tuple of ndarray(nbins)
    
    '''
    ws = []
    es = []
    bins = np.linspace(min, max, nbins+1)
    data_cut = pd.cut(df[index_field].__array__(), bins)
    
    group = df.groupby(by=data_cut)

    for s in group:
        # print (s[1][value_field].values)
        # print  (s[1][value_err_field].values)
        # pause()
        # print(weighted_avg(s[1][value_field].values, s[1][value_err_field].values, min_n=min_n, get_err=False))
        try:

            l = len(np.where(~np.isnan(s[1][value_err_field].values))[0])
            assert l >= min_n
            ws.append(np.nanmedian(s[1][value_field].values))
            es.append(np.sqrt(np.nansum(s[1][value_err_field].values**2)/(float(l))**2))

                    
        except:
            # pass
            ws.append(np.nan)
            es.append(np.nan)
    return np.array(ws), np.array(es)

def pandas_bin(df, index_field, value_field, min, max, nbins, func):                                          
    '''             
    applies <func> to the binned values of <value_field>, binned on <index_field>                             
    Examples        
    pandas_hist()   
    df : pandas DataFrame                         
     must have fields <index_field>, <value_field>
    index_field : str                             
     field on which binning is to be done         
    value_field :  str                            
     field that you want binned                   
    min : float     
     min index value to consider when binning     
    max : float     
     max index value to consider when binning     
    nbins : int     
     how many bins to break in to                 
    func : function 
     function to apply to binned values           
    '''             
    bins = np.linspace(min, max, nbins+1)
    data_cut = pd.cut(df[index_field].__array__(), bins)
    group = df[[value_field]].groupby(by=data_cut)               
    result = group.aggregate(func)             
    return result[value_field]

# def pandas_bin(df, index_field, value_field, min, max, nbins, func):
#     '''
#     applies <func> to the binned values of <value_field>, binned on <index_field>
#     Examples
#     pandas_hist()
#     df : pandas DataFrame
#      must have fields <index_field>, <value_field>
#     index_field : str
#      field on which binning is to be done 
#     value_field :  str
#      field that you want binned
#     min : float
#      min index value to consider when binning
#     max : float
#      max index value to consider when binning
#     nbins : int
#      how many bins to break in to
#     func : function
#      function to apply to binned values
#     '''
#     bins = np.linspace(min, max, num=nbins)
#     data_cut = pd.cut(df[index_field], bins)    
#     group = df.groupby(by=data_cut)
#     result = group.aggregate(func)
#     return result

def process_one(file, verbose=False, kepler=False, drop_duplicates=True, match_on='epic', sep=None, drop_no_epic=True, make_index=False):
    '''
    read in the table and rename the file column to make it 'EPIC' so that this can be later matched on.
    Inputs
    [ make_index : bool ]
     If True, will make epic to be the index of the dataframe.Default False. NOT COMPATIBLE YE TWTH COMBINE_TABLES, SO IT'S CURRENTLY NOT EVEN POSSIBLE FOR THIS TO BE DONE FROM A CALL TO COMBINE_TABLES. WOULD HAVE TO MANUALLY CHANGE THE COMBINE_TABLES CODE TO MAKE THAT POSOSIBLE. !!!
    [ drop_no_epic : bool ]
     if True, will not return rows for which no EPIC could be found.
    
    '''
    # put the whitespace separator last in case there are issues with some columns not being filled in
    # JCZ 221017
    # added nested try/except structure
    if sep is None:
        try:
            one = pd.read_table(file, comment='#', sep='\s+', header=0)
            if len(one.keys()) > 1:
                if verbose:
                    print('{} read successfully'.format(file))
            elif len(one.keys()) == 1:  
                one = pd.read_table(file, comment='#', sep=',', header=0)
            elif len(one.keys()) == 1:
                one = pd.read_table(file, comment='#', sep='\|', header=0)
            else:
                print('{} could not be read using , | or whitespace delimiters'.format(file))

        except:
            try:
                one = pd.read_table(file, comment='#', sep=',', header=0)
                if len(one.keys()) > 1:
                    if verbose:
                        print('{} read successfully'.format(file))
                elif len(one.keys()) == 1:
                    one = pd.read_table(file, comment='#', sep='\|', header=0)
                else:
                    print('{} could not be read using , | or whitespace delimiters'.format(file))
            except:
                one = pd.read_table(file, comment='#', sep='\|', header=0)
                if len(one.keys()) > 1:
                    if verbose:
                        print('{} read successfully'.format(file))
                else:
                    print('{} could not be read using , | or whitespace delimiters'.format(file))

            
                
    else:
        try:
            one = pd.read_table(file, comment='#', sep=sep, header=0)
        except:
            print('{} could not be read using the \'{}\' delimiter'.format(file,sep))
    # JCZ 070317
    # added this as an option to be kepler-friendly or not. jie's KICs don't prepend 00, so there are issues...
    
    if kepler:
        func = lambda x : (re.search('([0-9]{6,10})', x).group(1))
    else:
        func = lambda x : (re.search('([0-9]{6,10})', x).group(1))

    if verbose:
        print('found the following keys:')
        print(one.keys())
    if match_on == 'epic' or match_on == 'epic_orig':
        found = False
    else:
        found = True
    # JCZ 120317
    # aded kepid to deal with MAST files without having to rename the KIC col
    # JCZ 071118
    # 'epic' HAS TO BE FIRST !!!

    for k in ['epic', 'file', 'filename', 'EPIC', 'KICID', 'ID', 'KIC', 'KID', 'files', 'kepid', 'epic_orig','Kepler_ID', 'K2_id', 'id']:
        if k in one.keys() and found is False:
            # now cut on the pipeline_rating:
            one['__epic_orig'] = one[k].copy()
            try:
                one[k] = one[k].astype(int)
            except:
                pass
            if ((one[k]).dtype != 'int64' and (one[k]).dtype != 'float64') and found is False:
                
                one[k] = one[k].apply(func).astype(int)
            # merge one and pipeline output table (flags_df) so that if they are different sizes, one = one[one['pipeline_rating'] > 0] will work
            # JCZ 071118
            # added k != 'epic', because otherwise there will be two epic columns
            if found is False and k != 'epic':
                one = one.rename(columns={k:'epic'})
            # JCZ 040317
            # added because will have two epic_orig cols if epic and epic_orig are already present in table
            if 'epic_orig' in one.keys():
                one = one.rename(columns={'epic_orig':'epic_orig_orig'})
            one = one.rename(columns={'__epic_orig':'epic_orig'})
            print(one.keys())

            found = True
    if not found:
        if verbose:
            print('couldnt find a column named file epic or filename so assuming the first unnamed column is the correct one. assuming no header in the input file')
        one = pd.read_table(file, comment='#', sep='\s+|,', engine='python', header=None)
        # JCZ 270121
        # updating because of deprecated behvior in recent pandas
        # one = one.convert_objects(convert_numeric=True)
        one = one.apply(pd.to_numeric, errors='ignore')

        
        print(one)

        one['epic_orig'] = one[0].copy()
        if not((one[0]).dtype == 'int64' or one[0].dtype == 'float64'):
            one[0] = one[0].astype(str).apply(func).astype(int)
        if one[0].dtype == 'float64':
            one[0] = one[0].astype(int)
        one = one.rename(columns={0:'epic'})
        one['epic'] = one['epic'].__array__().astype(int)
    # except:
        
    #     if 'epic' not in one.keys():
        
    #         print 'could not find file or EPIC column for {}...'.format(file)
    #     # if epic does exist as a field then it must be converted to an int
    #     else:
    #         print 'found epic field -- converting to int'
    #         one['epic'] = one['epic'].astype(int)
    # JCZ 290916
    # drop any duplicates.

    if drop_duplicates:
        if verbose:
            print('dropping duplicates')
        one = one.drop_duplicates(subset=[match_on])
    # JCZ 151117
    # drop any objects that don't have an EPIC
    if drop_no_epic:
        one = one.dropna(subset=[match_on])

        # print one[match_on]
    if one[match_on].dtype == 'float64':
        one[match_on] = one[match_on].astype(int)
    if make_index:
        one.set_index(match_on, inplace=True, drop=False)
        # needed to add this for recent version of pandas
        one.rename_axis(None,inplace=True)
    return one

def combine_tables(*files, **kwargs):
    '''
    do separate = True and drop_duplicates=False to process many files with the same EPIC, but different values for their parameters
    [ match_on : str ]
     what field to match on. Should usually be 'epic' (Default). Using 'epic_orig' will match on whatever column is interpreated as the filename column, but before the EPIC digits are extracted from it. e.g., matching on kplr012504054.txt.clean.hipass.rebin.psd.sp.gif instead of the default ('epic') of matching on the EPIC embedded within (012504054, in the example case).
    [ deprecated : bool ]
     Default True. This will keep the old way of having a duplicate version of the first merged file. Otherwise, if False, will make it so there is no duplication. Required if how='append' is going to return the unique entries in the merged files --- otherwise the entries of the first file willb e duplicated once... Default is False if how=='append'.
    '''
    
    deprecated = kwargs.pop('deprecated', True)

    separate = kwargs.pop('separate',True)
    drop_duplicates = kwargs.pop('drop_duplicates',True)
    verbose = kwargs.pop('verbose',False)
    kepler = kwargs.pop('kepler',False)
    match_on = kwargs.pop('match_on', 'epic')
    sep = kwargs.pop('sep', None)
    how = kwargs.pop('how', 'merge')
    if how == 'append':
        deprecated = False
    drop_no_epic = kwargs.pop('drop_no_epic', True)
    # JCZ 070317
    # added kepler option
    # old = process_one(files[0], verbose=verbose, kepler=kepler, drop_duplicates=drop_duplicates, match_on=match_on, sep=sep)
    count = 0
    if len(files) < 2:
        print('only one file was provided, so returning info from that one file using process_one().')
        return process_one(files[0], verbose=verbose, kepler=kepler, drop_duplicates=drop_duplicates, match_on=match_on, sep=sep, drop_no_epic=drop_no_epic, make_index=True)
    if not deprecated:
        for file in files:
            files = files[1:]
            try:
                old = process_one(file, verbose=verbose, kepler=kepler, drop_duplicates=drop_duplicates, match_on=match_on, sep=sep,  drop_no_epic=drop_no_epic, make_index=True)
                old[file] = 1
                break
            except:
                print('could not read {}. moving on to next file.'.format(file))
                continue
        

        
    for file in files:
        if verbose:
            print('################################################################')
        # JCZ 070317
        # kepler option
        if count == 0 and deprecated:
            try:
                old = process_one(file, verbose=verbose, kepler=kepler, drop_duplicates=drop_duplicates, match_on=match_on, sep=sep,  drop_no_epic=drop_no_epic, make_index=True)
            except:
                print('could not read {}. moving on to next file.'.format(file))
                continue
        try:
            new = process_one(file, verbose=verbose, kepler=kepler, drop_duplicates=drop_duplicates, match_on=match_on, sep=sep, drop_no_epic=drop_no_epic,make_index=True)
            # to figure out which file the info came from, make a new column that is named the file name and that is 1 if there are columns from that file in a given row.
            new[file] = 1
        except:
            print('could not read {}. moving on to next file.'.format(file))
            continue
        if verbose:
            print('################################################################################################')
        if separate:

            new = old.merge(new, on=match_on, how='inner', suffixes=('', '_'+str(count)))
        else:
            if verbose:
                print('not separating')
            # update the old frame with values from the new. match on epic still
            if how == 'merge':
                new = old.merge(new, on=match_on, how='outer', suffixes=('', '_'+str(count)))
            # new = old.merge(new, how='outer', on='epic')
            
            # new = pd.concat([old, new])
            # new = old.join(new, how='outer')
            # this operation works in place!!
            # old.update(new)
            if how == 'combine_first':
                new = old.combine_first(new)
            if how == 'append':
                new = old.append(new)
            if drop_duplicates:
                print('dropping duplicates')
                new = new.drop_duplicates(subset=[match_on], keep='last')
        l = len(new)
        new.dropna(how='all', inplace=True, subset=[col for col in new.keys() if str(col) not in match_on])
        l_ = len(new)
        if l != l_:
            if verbose:
                print('dropped {} row due to NaNs being found in it'.format(l - l_))
            else:
                pass
        if verbose:
            print('length of table is now {}'.format(l_))
        if l == 0:
            print('could not combine the table {} with the others'.format(file))
        old = new
        count += 1
    global N
    N = count
    if verbose:
        print('combined a total of {} files'.format(N))
        print('and combined table has the following keys:')
        print(new.keys())

    return new


def file2epic_int(file):
    '''
    takes a file, extracts the KIC, and returns it as an int
    Inputs
    file : str
    e.g., hlspktwo123456789_blah
    Ouputs
    epic : int
    '''
    func = lambda x : (re.search('([0-9]{6,10})', file).group(1))


    try:
        epic = func(file)
    except:
        print('could not do {}'.format(file))
        epic = 0
    # i'm not sure why this was added, so i'm getting rid of it...
    # if type(epic) == type(''):
        # epic = np.array([epic])
    # epic = epic.astype(int)
    return int(epic)
#@profile
def binned_statistic(x, values, func, nbins, range):
    '''The usage is approximately the same as the scipy one''' 
    from scipy.sparse import csr_matrix
    N = len(values)
    r0, r1 = range

    digitized = (float(nbins) / (r1-r0) * (x-r0)).astype(int)
    S = csr_matrix((values, [digitized, np.arange(N)]), shape=(nbins, N))

    return np.array([func(group) for group in np.split(S.data, S.indptr[1:-1])])

def add_gauss(c,ax=None, color='black', add_zero=True, **kwargs):
    '''
    [ add_zero : bool ]
     If True, put on a black dashed line showing the median of the Gaussian at 0.0. Will extend tot he top of the bell curve
    '''
    c = c.__array__()
    if ax is None:
        ax = plt.figure().gca()
        # xmin = np.nanpercentile(c, 5)
        # xmax= np.nanpercentile(c, 95)
        xmin = np.nanmin(c)
        xmax = np.nanmax(c)
    else:
        xmin, xmax = ax.get_xlim()[0], ax.get_xlim()[1]

    N = 500

    x = np.linspace(xmin, xmax, num=N)
    mean = 0.0#np.mean(c)
    ns_normed =np.exp(-(x - mean)**2/2.)/np.sqrt(2.*np.pi)
    max_y = ns_normed[np.argmin(np.abs(x -  0.0))]
    ax.vlines(0., 0., max_y, linestyle='dashed', color='black')
    ax.plot(x, ns_normed, color=color, linestyle='dashed')

    return ax

def hist_w_gauss(*args, **kwargs):
    get_bins = kwargs.pop('get_bins', False)
    ax = plot_hist(*args, **kwargs)
    if get_bins:
        kwargs['get_bins'] = get_bins
        out = plot_hist(*args, **kwargs)
        return out
    kwargs['ax'] = ax
    ax = add_gauss(args[0], **kwargs)
    if kwargs.get('save', True):
        ext = kwargs.get('ext', 'jpg')
        ax.figure.tight_layout()
        ax.figure.savefig(args[1] + '.' + ext, format=ext)

    else:
        return ax
    
def plot_hist(c, plotname, ax=None, xlabel='', ylabel='', color='black', label='', plot_legend=False, ext='jpg', nbins=None, xlims=[-5, 5], ylims=[None, None], normed=True, loc='upper left', yscale='linear', xscale='linear', plot_median=True, plot_zero=True, save=True, plot_err=True, caption=None,linestyle='solid', alpha=0.5, get_bins=False,**kwargs):
    '''
    plots histogram of the quantity, with a Gaussian overlaid.
    Inputs
    [ get_bins : bool ]
     If True, return (bins, vals, errors) and NOT ax. Otherwise, return ax.
    [ alpha : float ]
     transparency.
    [ plot_err : bool ]
     If True, will plot error bars on each histogram bin, according to Poisson statistics. Default True. Will make them alpha of 1 no matter what, and will also make them relatively thick...
    [ ylims : list of 2 elements ]
     Could be [None, None] or [0.0, 1.0], e.g..
    [ plot_median : bool ]
     if True, put in vertical line corresponding to the median of the sample
    [ plot_zero : bool ]
     if True, put in vertical line at x = 0. Only valid if plot_median is True.
    [ yscale : str ['log', 'linear'] ]
     whether or not to put y-scale in log. Default 'linear'.
    [ xscale : str ['log', 'linear'] ]
     whether or not to put x-scale in log. Default 'linear'.
    [ loc : str ]
     legend location. Default 'upper left'.
    [ xlims : list (2,) | tuple (2,) ]
     only <c> values between xlimits will be plotted. Can be None.
    [ nbins : int | None ]
     if None, uses sqrt(N(c))
    c : ndarray
     array that will be plotted
    plotname : str
     <plotname>.eps and *.png will be saved
    [ ax : axis instance ]
     use this to overplot things
    [ xlabel : str ]
    [ ylabel : str ]
    [ color : str ]
    [ label : str ]
     if legend is turned on, this will be the label for the line. useful if overplotting by passing <ax>
    [ caption : str | None ]
     if not None, will be something like 'a)' and will plot that on the plot 
    [ ext : str ]
     can be 'eps', 'png', 'jpg', etc. a 'png' will be made no matter what.
    Examples
    c = (tab[x] - tab[y])/np.sqrt(tab[x_err]**2 + tab[y_err]**2)
    hist_w_gauss(c)
    
    '''


    c = c.__array__()
    # JCZ 300819
    # comemnting this out
    # plt.style.use('jcz_paper_latex')
    if ax is None:
        ax = plt.figure().gca()

    if nbins is None:
        nbins = int((len(c))**(1./2.))
        
    med = np.nanmedian(c)
    print('min/max : {} {}'.format(np.nanmin(c), np.nanmax(c)))
    __xlims = [xlims[0], xlims[1]]
    if __xlims[0] is None:
        __xlims[0] = -np.inf
    if __xlims[1] is None:
        __xlims[1] = np.inf
    l_c0 = len(c)
    print(__xlims)
    if 'lw' not in kwargs.keys():
        lw = mplt.rcParams['lines.linewidth']
    else:
        lw = kwargs['lw']
    n_outliers = len(np.where(np.logical_and(~np.isnan(c), np.logical_or(c < __xlims[0], c > __xlims[1])))[0])
    # n_outliers = len(np.where(n(np.logical_and(c < __xlims[0], c > __xlims[1])))[0])
    # print np.where(np.logical_and(~np.isnan(c), np.logical_and(c < __xlims[0], c > __xlims[1])))[0]
    # print ~np.isnan(c)
    c = c[np.where(np.logical_and(c > __xlims[0], c < __xlims[1]))[0]]
    l_c = len(c)
    print('{} things being plotted'.format(l_c))
    print('{} outliers not being plotted.'.format(n_outliers))
    print('{} nans not being plotted.'.format(l_c0 - l_c - n_outliers))
    # JCZ 200320
    # changing this. this is a major change, meaning that when you overplot histograms on the same plot, they will actually be normalized the same.
    # _xlims = [np.nanmin(c), np.nanmax(c)]
    _xlims = [__xlims[0], __xlims[1]]
    if xscale == 'log' or xscale == 'symlog':
        __xlims = np.log(_xlims)
    
        bins = np.exp(np.linspace(__xlims[0], __xlims[1], nbins))
    else:
        bins = (np.linspace(_xlims[0], _xlims[1], nbins))
    ns_normed, bins, _ = hist(c, histtype='step', color='white', edgecolor=color, ax=ax, label=label, density=normed,alpha=alpha, bins=bins, linewidth=mplt.rcParams['lines.linewidth']*0.5, linestyle=linestyle)
    if ylims[0] is None:
        ylims[0] = list(ax.get_ylim())[0]
        
    if ylims[1] is None:
        ylims[1] = list(ax.get_ylim())[1]
        
    if xscale == 'log' or xscale == 'symlog':
        bins = np.exp(np.linspace(__xlims[0], __xlims[1], nbins))
    else:
        bins = (np.linspace(_xlims[0], _xlims[1], nbins))
    ################################################################
    # add Poisson errors
    ################################################################
    # re-compute a non-normalized histogram to get normalization for each bin and the real number for each bin
    ns, bins, _ = hist(c, histtype='stepfilled', color='white', edgecolor=color, ax=ax, density=False,alpha=0, bins=bins, linewidth=mplt.rcParams['lines.linewidth']*0.5, linestyle=linestyle)
    if xscale == 'log' or xscale == 'symlog':
        xlims = np.exp(__xlims)
    
    bins = np.array(bins)
    # these are mid-point bin values
    _bins = np.diff(bins)
    bins = np.cumsum(_bins) + bins[0] - _bins/2.0
    print(bins, ns_normed)
    # \sigma = \sqrt{N_{\mathrm{bin}}}
    # and then additionally scaled down to be on the same scale as the normalized histogram values
    if normed:
        errors = np.sqrt(ns)*(ns_normed/ns)
        # JCZ 241217
        # made None because i figure if ylabel is set to '' it was done manually and the intention was to have no label. Otherwise, a smart label should be chosen if ylabel was not set to anything (None)
        if ylabel is None:
            ax.set_ylabel('Normalized count per bin')
    else:
        errors = np.sqrt(ns)
        if ylabel is None:
            ax.set_ylabel('Count per bin')

    if ylabel != '' and ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.set_yscale(yscale)
    if plot_err:
        # JCZ 090620
        # plot only every fourth one

            
        ax.errorbar(bins[::4], ns_normed[::4], xerr=None, yerr=errors[::4], linestyle='none', color=color, elinewidth=mplt.rcParams['lines.linewidth']*1.5, capsize=0)
    ################################################################
    ax.set_xlim(xlims)

    ax.set_xscale(xscale)
    
    # JCZ 220118
    # commenting this out because labelsize should already be OK.
    # ax.tick_params('x', which='major', labelsize=15.)

    # ax.set_xlabel(r'{}'.format(xlabel))
    ax.set_xlabel(xlabel)
    # JCZ 191016
    # remove legend for nicer plot
    if plot_legend:
        ax.legend(loc=loc, fontsize=25*0.6)
    # increase upper edge of y-limit to take into account the errorbars that have been added on since the definition of the y-limits.
    if True:
        # ylims = list(ax.get_ylim())
        ylims[1] = np.max([ylims[1], (np.nanmax(ns_normed) + np.nanmax(errors))*1.13])
    # ylims = ax.get_ylim()
    # plot median of the sample

    if plot_median:
        max_y = ns_normed[np.argmin(np.abs(bins - med))]
        ax.vlines(med, 0.,max_y, linestyle='solid', color=color, lw=lw*0.5)
        # overplot median of zero
        if plot_zero:
            max_y = ns_normed[np.argmin(np.abs(bins -  0.0))]
            ax.vlines(0., 0., max_y, linestyle='dashed', color='black', lw=lw*0.5)

    ax.set_ylim(ylims)
    # !!! DEBUG removed this temporarily.
    # plt.tight_layout()
    # JCZ 220118
    # added the a) b) c) label option
    if caption is not None:
        ax.text(0.9, 0.9, caption, fontsize=25, transform=ax.transAxes)

    if save:
        ax.figure.tight_layout()
        ax.figure.savefig(plotname + '.' + ext, format=ext)

        # import os
        # filename = plotname + '.' + ext
        # basename = filename.split('.png')[0]
        # os.system('convert {} {}.eps'.format(filename, basename))
    if get_bins:
        return (bins,ns_normed, errors)
    else:
        return ax

#@profile
def freq_cut(freqs, values, range):
    '''
    will excise freqs and <values> corresponding to frequencies outside of range
    Inputs
    freqs : ndarray
    Inputs 
    values : ndarray
    range : tuple
     (min_f, max_f)
    Outputs
    freqs, values
     cut outside of [<range>]

    '''

    inds = np.where(np.logical_and(freqs < range[1], freqs > range[0]))[0]

    return freqs[inds], values[inds]
    
