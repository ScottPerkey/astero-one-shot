from __future__ import print_function
import os
import numpy as np
# if os.environ.get("SSH_TTY"):
print('setting up graphics for ssh connection')
import matplotlib as mplt
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec
mplt.use("Agg")
import pylab as plt
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# else:
    # import matplotlib as mplt
    # import matplotlib.gridspec as gridspec
    # import pylab as plt
    # from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    # JCZ 280721
    # commented these out
# fig = plt.figure()
# ax = fig.gca()
# # empty handle
# empty_handle = mplt.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none', visible=False) 


def plot_sequence(x, y, z, cmap='Blues', norm=None, ax=None, lw=1.0, alpha=1.0):
    '''
    plots the x and y points as a spectrum with cmap with colors corresponding to z
    code taken from https://matplotlib.org/examples/pylab_examples/multicolored_line.html
    norm : LineCollection thing
     passed to LineCollection
    cmap : str
     passed to plt.get_cmap() and then LineCollection
    ax : pylab axis object
     is None, one will be made
    lw : float
     line with for the segments of the spectrum
    alpha : float
     the transparency of the segments
    '''
    
    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be numlines x points per line x 2 (x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Create the line collection object, setting the colormapping parameters.
    # Have to set the actual values used for colormapping separately.
    # if norm is None:
        # vmin = np.nanmin(x)
        # vmax = np.nanmax(z)
        # norm = plt.Normalize(vmin=vmin,vmax=vmax)
    lc = LineCollection(segments, cmap=plt.get_cmap(cmap), norm=norm)
    lc.set_array(z)
    lc.set_linewidth(lw)
    # lc.set_alpha(alpha)
    print (lc)
    ax.add_collection(lc)
    return ax
