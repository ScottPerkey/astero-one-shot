#!/home/zinn.44/.conda/envs/stan2/bin/python
# helper functions and classes for plotting in python
# Joel C. Zinn 2015
from __future__ import print_function
# from setup_graphics import *

#from pudb import set_trace as pause
import pylab as plt
import os
import numpy as np


class Plot(object):
    '''
    plot class to make it easier to see things and manipulate them.
    '''
    def __init__(self, **kwargs):
        self.show = kwargs.pop('show', False)
        self.save = kwargs.pop('save', False)
        self.mail = kwargs.pop('mail', False)
        self.ext = kwargs.pop('ext', 'eps')
        self.tight_layout = kwargs.pop('tight_layout', True)
        self.figsize = kwargs.pop('figsize', (9,9))
        try:
            from python_params_pipeline import pypar
            default_plot_dir = os.path.join(pypar.astero_plot_data_dir, 'plots')
        except:

            default_plot_dir = './'
        print('setting default plot output directory to {}'.format(default_plot_dir))
        self._fdirect = kwargs.pop('fdirect', default_plot_dir)
        self._fname = kwargs.pop('fname', 'tmp')
        self.fig = kwargs.pop('fig', plt.figure(figsize=self.figsize))
        self.ax = kwargs.pop('ax', self.fig.gca())
    @property
    def fname(self):
        # JCZ 310317
        # force the directory to be demarked by a '/', just in case the user doesn't make fdirect end in '/'
        return self._fdirect + '/' + self._fname + '.' + self.ext
    @fname.setter
    def fname(self, value):
        self._fname = value

    def wrapup(self, hspace=1.08):
        # show the plot
        if self.show:
            plt.ion()
            plt.show()
            # pause()
            pass
        if self.save:
            # save the plot to current directory
            print('saving plot to {}'.format(self.fname))
            # plt.tight_layout(h_pad=hspace)
            if self.tight_layout:
                plt.tight_layout()
            self.fig.savefig(self.fname, format=self.ext)
            pass

        if self.mail:
            self.fig.savefig(self.fname, format=self.ext)
            print('saved {}'.format(self.fname))
            #os.system('/bin/touch tmp_mail')
            #os.system('/bin/mail -a ' + self.fname + ' -s "' + self.fname + '" joelczinn@gmail.com < tmp_mail')
            #os.system('/bin/rm tmp_mail')


def makeTempHistogramPlot(xdata,ydata,filename='test',xlims=-99, ylims =-99 , \
                          nxbins = 10,nybins=10, bw=True, nbins=20,contours=1,sigma=1,line=1, xlabel='', ylabel='', title='', color=None, fig=None, scatter=True, hist=True, bottomx=None, bottomy=None, get_bottoms=False, **scatter_args):
    '''
    Combination scatter plot and histogram
    Inputs
    xdata : np array
    ydata : np array
    color : str (None)
     if None (default), then the colors of the scattered points and the lines of the histogram will be colored
    fig : object
     figure object, e.g., as returned by a previous call to this function
    filename : str
     saved plot will be <filename>.eps
    bw : bool
     black and white or not
    bottomx, bottomy : ndarray or list or None
     if not None these will be the bottoms for the x and y histograms (useful for stacking multiple data)
    Outputs
     fig
     saves plot <filename>.eps
    
    Notes
    See also
    will set the limits to autmatically be 
    
    '''
    from matplotlib import ticker
    #bw = 0 for color, = 1 for black and white
    #line = 0 for no line, =1 for line
    #sigma = 1 for display % below line, =0 for not
    #contours = 1 for display 1,2,3 sigma contours, = 0 for not.

    # Define the x and y data
    small = 1e-5
    x = xdata
    y = ydata
    if fig is None:
        fig_ = None
    else:
        fig_ = fig
    bad_y = (y < small)
    bad_x = (x < small)
    y = y[~np.logical_or((bad_y), (bad_x))]
    x = x[~np.logical_or((bad_y), (bad_x))]

    # Set up the size of the figure
    if fig is None:
        fig = plt.figure(1, figsize=(9,9))
        fig.clf()
    # if there are multiple plots being appended and limits are not explicitly stated, just use the limits from teh first plot
    elif xlims == -99 and ylims == -99:
        # the below lines don't work because the previous axes' limits are not changed...
        # ylims = [np.min([fig.axes[0].get_ylim()[0], min(y)]), np.max([fig.axes[0].get_ylim()[1], max(y)])]
        # ylims = [np.min([fig.axes[0].get_xlim()[0], min(x)]), np.max([fig.axes[0].get_xlim()[1], max(x)])]
        ylims = fig.axes[0].get_ylim()
        xlims = fig.axes[0].get_xlim()

    # Set up default x and y limits
    if (xlims == -99): xlims = [max(x),min(x)]
    if (ylims == -99): ylims = [min(y),max(y)]

    # Set up your x and y labels
    xlabel = xlabel
    ylabel = ylabel


    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02

    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height] # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram


    # Make the three plots
    if fig_ is None:
        axTemperature = plt.axes(rect_temperature) # temperature plot
        axHistx = plt.axes(rect_histx) # x histogram
        axHisty = plt.axes(rect_histy) # y histogram
    else:
        axTemperature = fig.axes[0]
        axHistx = fig.axes[1]
        axHisty = fig.axes[2]

   # Remove the inner axes numbers of the histograms
    nullfmt = ticker.NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)

    # Make the 'main' temperature plot
    xbins = np.linspace(start = ymax, stop = xmax, num = nxbins)
    ybins = np.linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    aspectratio = 1.0*(xmax - 0)/(1.0*ymax - 0)
#    H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins))
#    X = xcenter
#    Y = ycenter
#    Z = H

    # Plot the temperature data
    if scatter:
        print(scatter_args)
        axTemperature.scatter(x,y, color=color, **scatter_args)
        axTemperature.text(0.5, 0.8,title, fontsize=15, transform=axTemperature.transAxes, horizontalalignment='center', verticalalignment='center')
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=20)
    axTemperature.set_ylabel(ylabel,fontsize=20)

    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(15)
        label.set_family('serif')

    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(15)
        label.set_family('serif')

    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)

    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)

    #Plot the histograms

    bottomx_, dum1, dum2 = axHistx.hist(x, bins=xbins, color = color, histtype='bar', bottom=bottomx)
    bottomy_ , dum1, dum2= axHisty.hist(y, bins=ybins, orientation='horizontal', color = color, histtype='bar', bottom=bottomy)
    print('set ', bottomx_, ' on top of ', bottomx)
    print('set ', bottomy_, ' on top of ', bottomy)
    bottomy = bottomy_
    bottomx = bottomx_

    #Set up the histogram limits
    axHistx.set_xlim(xlims )
    axHisty.set_ylim(ylims)

    #Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family('serif')

    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(plt.MaxNLocator(4))
    axHistx.yaxis.set_major_locator(plt.MaxNLocator(4))

    if filename:
        plt.savefig('./' + filename + '.eps',format = 'eps')
        # savefig(filename + '.pdf',format = 'pdf', transparent=True)
        # savefig(filename + '.png',format = 'png', transparent=True)
    if get_bottoms:
        return fig, bottomx, bottomy
    return fig


def log_neg(_x, _y, fig, ylog=True, xlog=False, debug=False, **scatter_kwds):
    '''
    takes a figure and makes a plot composed of two axes with a shared axis that can support negative values
    Inputs
    _x : ndarray
     x-values to be plotted
    _y : ndarray
     y-values to be plotted. can contain negative values even though this will be a log plot
    [ xlog : bool ]
     whether or not to make the xscale log, as well.
    [ debug : bool]
     option to be turned on when editing log_neg
    OUtputs
    fig : object
     the updated figure
    Notes
    THIS ROUTINE WILL CLEAR THE FIGURE THAT IS INPUT before updating it and returning it
    See also
    
    '''
    x = _x.copy()
    y = _y.copy()
    mask = np.where(np.logical_and(~np.isnan(x), ~np.isnan(y)))[0]
    x = np.array(x[mask])
    y = np.array(y[mask])
    if 'c' in scatter_kwds.keys():
        c = scatter_kwds['c'][mask]
        del scatter_kwds['c']
    else:
        c = np.zeros(shape=x.shape)
    if debug:
        y = np.random.normal(loc=0,scale=0.1,size=100)
        x = np.random.normal(loc=0,scale=0.1, size=100)
        # x = x**2
    xneg = False
    
    yneg=False
    if len(np.where(x < 0.)[0]) > 0:
        xneg = True
    if len(np.where(y < 0.)[0]) > 0:
        yneg=True

    fig.clf()
    fig.subplots_adjust(hspace=0., wspace=0.)

    # NOTE ON ORDERING
    # azx1 is Q2 (top left)
    # ax2 is Q3 (bottom left)
    # ax3 is Q1 (top right)
    # ax4 is Q4 (bottom right)
    lax1 = np.logical_and(y > 0. , x > 0.).any()
    lax2 = np.logical_and(y < 0. , x >0.).any()
    lax3 = np.logical_and(y > 0. , x < 0.).any()

    lax4 = np.logical_and(y<0. , x < 0.).any()

    if lax1 and lax2:
        gs1 = gridspec.GridSpec(2,1)
        gs1.update(wspace=0.0, hspace=0.0)
        ax1 = fig.add_subplot(gs1[0]) #221
        ax2 = fig.add_subplot(gs1[1]) #223

    if lax3 and lax4:
        gs1 = gridspec.GridSpec(2,1)
        gs1.update(wspace=0.0, hspace=0.0)
        ax3 = fig.add_subplot(gs1[0]) #221
        ax4 = fig.add_subplot(gs1[1]) #223

    if (lax1 and lax3) or (lax2 and lax4):
        gs1 = gridspec.GridSpec(2, 2)
        gs1.update(wspace=0.0, hspace=0.0)
        ax1 = fig.add_subplot(gs1[0]) #221
        ax2 = fig.add_subplot(gs1[2]) #223
        ax3 = fig.add_subplot(gs1[1]) #222
        ax4 = fig.add_subplot(gs1[3]) # 224
    if lax4 and not lax1 and not lax2 and not lax3:
        ax4 = fig.add_subplot(111)
    if lax3 and not lax1 and not lax2 and not lax4:
        ax3 = fig.add_subplot(111)
    if lax2 and not lax1 and not lax3 and not lax4:
        ax2 = fig.add_subplot(111)
    if lax1 and not lax2 and not lax3 and not lax4:
        ax1 = fig.add_subplot(111)

    if lax1:
        ax1.scatter(x[np.logical_and(y > 0. , x > 0.)], y[np.logical_and(y > 0., x>0.)], c=c[np.logical_and(y > 0., x>0.)], **scatter_kwds)
    if lax2:
        ax2.scatter(x[np.logical_and(y < 0. ,x > 0.)], -y[np.logical_and(y < 0. , x >0.)],c=c[np.logical_and(y < 0. , x >0.)], **scatter_kwds)
    if lax3:
        ax3.scatter(-x[np.logical_and(y > 0. , x < 0.)], y[np.logical_and(y > 0., x<0.)], c=c[np.logical_and(y > 0., x<0.)],**scatter_kwds)
    if lax4:
        ax4.scatter(-x[np.logical_and(y<0. , x < 0.)], -y[np.logical_and(y< 0. ,x <0.)],c=c[np.logical_and(y< 0. ,x <0.)],**scatter_kwds)

    if ylog:
        if lax2:
            ax2.set_yscale('log')
        if lax1:
            ax1.set_yscale('log')
        if lax3:
            ax3.set_yscale('log')
        if lax4:
            ax4.set_yscale('log')
    if xlog:
        if lax2:
            ax2.set_xscale('log')
        if lax1:
            ax1.set_xscale('log')
        if lax3:
            ax3.set_xscale('log')
        if lax4:
            ax4.set_xscale('log')

    all_axes = []
    if lax1:
        all_axes.append(ax1)
    if lax2:
        all_axes.append(ax2)
    if lax3:
        all_axes.append(ax3)
    if lax4:
        all_axes.append(ax4)
    # force the seam to the lower of the two limits

    seam_low = min([ax.get_ylim()[0] for ax in all_axes])
    if lax1:
        ax1.set_ylim([seam_low,ax1.get_ylim()[1]])
    if lax2:
        ax2.set_ylim([seam_low, ax2.get_ylim()[1]])
        ax2.set_ylim(ax2.get_ylim()[::-1])
    if lax3 and lax1:
        ax3.set_ylim([seam_low, max([ax3.get_ylim()[1], ax1.get_ylim()[1]])])
    if lax4 and lax2:
        ax4.set_ylim([seam_low, max([ax4.get_ylim()[1], ax2.get_ylim()[0]])])
        ax4.set_ylim(ax4.get_ylim()[::-1])
    seam_hi = min([ax.get_xlim()[0] for ax in all_axes])
    if lax1 and lax2:
        ax1.set_xlim([max([ax2.get_xlim()[1], ax1.get_xlim()[1]]), seam_hi])
        ax2.set_xlim([max([ax2.get_xlim()[1], ax1.get_xlim()[0]]), seam_hi])

    if lax3 and lax4:
        ax3.set_xlim([seam_hi, max([ax3.get_xlim()[1], ax4.get_xlim()[1]])])
        ax4.set_xlim([seam_hi, max([ax3.get_xlim()[1], ax4.get_xlim()[1]])])

    # NB these labels will only make sense when both axes are log-scaled...
    if xlog and lax2:
        ax2.set_xticklabels([r'$10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax2.get_xticks()[:-1]])
    if lax2 and not xlog:
        ax2.set_xticklabels([r'${{{:3.2f}}}$'.format((_x)) for _x in ax2.get_xticks()[:-1]])
        ax2.xaxis.set_tick_params(labelsize=15)
    if lax2:
        ax2.set_yticklabels([r'-$10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax2.get_yticks()])
    if lax1:
        ax1.set_yticklabels([r'$10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax1.get_yticks()[:-1]])
        ax1.set_xticklabels([])
    if lax3:
        ax3.set_yticklabels([r'$10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax3.get_yticks()[:-1]])
        ax3.set_xticklabels([])
    if lax4:
        ax4.set_yticklabels([r'-$10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax4.get_yticks()])


    if xlog and lax4:
        ax4.set_xticklabels([r'$-10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax4.get_xticks()])
            # ax3.set_xticklabels([r'$-10^{{{:1.0f}}}$'.format(np.log10(_x)) for _x in ax3.get_yticks()])
    if lax4 and not xlog:
        ax4.set_xticklabels([r'$-{{{:3.2f}}}$'.format((_x)) for _x in ax4.get_xticks()[:-2]])
        ax4.xaxis.set_tick_params(labelsize=15)
    if lax3:
        ax3.yaxis.set_ticks_position('right')
    if lax4:
        ax4.yaxis.set_ticks_position('right')

    fig.tight_layout()

    if debug:
        fig.savefig('tmp.eps', format='eps')
    # fig.show()
    
    return fig

def new_color(ax=None):
    if ax is None:
        ax = fig.gca()
    return ax._get_lines.color_cycle.next()
