# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Module for routines common to many plotting scripts

import numpy as np
from matplotlib import colors, ticker
import matplotlib.pyplot as plt
from common import *

color_order = ['b', 'orange', 'r', 'm', 'c', 'k', 'y', 'g']*2
style_order = ['-', '--', '-.', ':']*2
marker_order = [".", "o", "v","s", "*", "x", "^", "<", ">"]*2
default_lw = 1.0
default_s = 0.2 # markersize
default_labelsize = 12
default_titlesize = 12
default_ticksize = 12
default_margin = 1/8
default_line_height = 7/32 # height of a line of text
default_margin_xlabel = 1/2
default_margin_ylabel = 3/4
# ylabels take up more space because floating
# point numbers are longer than they are tall
default_margin_title = 3/4

# default figure sizes 
sub_width_inches_default = 3.5
sub_height_inches_default = 2.5
sub_aspect_default = sub_height_inches_default/sub_width_inches_default

# lineplots 
lineplot_fig_dimensions = dotdict(dict({'sub_width_inches': sub_width_inches_default, 'sub_height_inches': sub_height_inches_default, 'sub_margin_top_inches': 1/4, 'sub_margin_bottom_inches': default_margin_xlabel, 'sub_margin_left_inches': default_margin_ylabel, 'margin_top_inches': 1.0}))

def axis_range(ax): # gets subplot coordinates on a figure in "normalized"
        # coordinates
    pos = plt.get(ax, 'position')
    bottom_left = pos.p0
    top_right = pos.p1
    xmin, xmax = bottom_left[0], top_right[0]
    ymin, ymax = bottom_left[1], top_right[1]
    return xmin, xmax, ymin, ymax

def xy_grid(x, y):
    """
    plt.pcolormesh() takes arguments X, Y, and C; treats the X, Y 
    arrays as VERTICES of quadrilaterals (not center points)
    For plotting, we usually have X, Y, and C of shape (m, n) and
    (X[i, j], Y[i, j]) represents the CENTER of a quadrilateral
    Thus, a call to plt.pcolormesh will ignore the last row/column
    of C, using only the smaller array C[:m-1, :n-1]
    "xy_grid" takes 1D arrays x, y and makes new arrays 
    # X_new, Y_new, which have dimension
    (m+1, n+1) and represent the vertices of quadrilaterals with centers
    at X[i, j] (if X, Y = meshgrid(x, y, indexing='ij'))
    """
    m, n = len(x), len(y)
    xmid = (0.5*(x[:m-1] + x[1:])).tolist()
    ymid = (0.5*(y[:n-1] + y[1:])).tolist()
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    x_new = [xmid[0] - dx] + xmid + [xmid[-1] + dx]
    y_new = [ymid[0] - dy] + ymid + [ymid[-1] + dy]

    return np.meshgrid(x_new, y_new, indexing='ij')

def testtex(label):
    plt.plot(range(10))
    plt.title(label)
    plt.show()

make_figure_kwargs_default =\
dict({'nplots': None, 'nrow': None, 'ncol': None,\

    'sub_width_inches': None,\
    'sub_height_inches': None,\
    'sub_aspect': None,\
    
    'sub_margin_left_inches': default_margin_ylabel, \
    'sub_margin_right_inches': default_margin,\
    'sub_margin_bottom_inches': default_margin_xlabel,\
    'sub_margin_top_inches': default_line_height,\

    'margin_left_inches': default_margin,\
    'margin_right_inches': default_margin,\
    'margin_bottom_inches': default_margin,\
    'margin_top_inches': default_margin_title,\

    'width_inches': None,\
    'height_inches': None,\
    'aspect': None})

def make_figure(**kwargs):
    kw = update_dict(make_figure_kwargs_default, kwargs)
    find_bad_keys(make_figure_kwargs_default, kwargs, 'make_figure')

    # unpack everything (annoying)
    nplots = kw.nplots; nrow = kw.nrow; ncol = kw.ncol

    sub_width_inches = kw.sub_width_inches
    sub_height_inches = kw.sub_height_inches
    sub_aspect = kw.sub_aspect
        
    sub_margin_left_inches = kw.sub_margin_left_inches
    sub_margin_right_inches = kw.sub_margin_right_inches
    sub_margin_bottom_inches = kw.sub_margin_bottom_inches
    sub_margin_top_inches = kw.sub_margin_top_inches

    margin_left_inches = kw.margin_left_inches
    margin_right_inches = kw.margin_right_inches
    margin_bottom_inches = kw.margin_bottom_inches
    margin_top_inches = kw.margin_top_inches

    width_inches = kw.width_inches
    height_inches = kw.height_inches
    aspect = kw.aspect

    # figure out nplots, nrow, ncol
    # default values are nplots = nrow = ncol = 1

    # assume two parameters are specified, otherwise go to default
    nspec = 0
    for val in [nplots, nrow, ncol]:
        nspec += not val is None
    if nspec == 2: 
        if nplots is None:
            nplots = nrow*ncol
        elif ncol is None:
            ncol = int(np.ceil(nplots/nrow))
        else:
            nrow = int(np.ceil(nplots/ncol))
    else: # all unspecified
        nplots = ncol = nrow = 1

    # OK, now the subplot structure is specified

    # next we need two parameters to get the subplot size
    nspec = (not (height_inches is None and sub_height_inches is None)) +\
            (not (width_inches is None and sub_width_inches is None)) +\
            (not (aspect is None and sub_aspect is None))

    if nspec == 2: # determine one parameter from the other two
        # first, take care of the fact that if either one of the "figure" lengths is specified,
        # the corresponding "sub_" length is specified
        if not width_inches is None:
            sub_width_inches = (width_inches - margin_left_inches - margin_right_inches -\
                    ncol*(sub_margin_left_inches + sub_margin_right_inches))/ncol
        elif not sub_width_inches is None:
            width_inches = margin_left_inches + margin_right_inches +\
                    ncol*(sub_margin_left_inches + sub_width_inches + sub_margin_right_inches)

        if not height_inches is None:
            sub_height_inches = (height_inches - margin_bottom_inches - margin_top_inches -\
                    nrow*(sub_margin_bottom_inches + sub_margin_top_inches))/nrow
        elif not sub_height_inches is None:
            height_inches = margin_bottom_inches + margin_top_inches +\
                    nrow*(sub_margin_bottom_inches + sub_height_inches + sub_margin_top_inches)

        # maybe get the aspects from the lengths
        if not width_inches is None and not height_inches is None: 
            aspect = height_inches/width_inches
            sub_aspect = sub_height_inches/sub_width_inches
        elif height_inches is None: # use the (sub_)aspect to get heights
            if not aspect is None:
                height_inches = aspect*width_inches
                sub_height_inches = (height_inches - margin_bottom_inches - margin_top_inches -\
                        nrow*(sub_margin_bottom_inches + sub_margin_top_inches))/nrow
                sub_aspect = sub_height_inches/sub_width_inches
            else:
                sub_height_inches = sub_aspect*sub_width_inches
                height_inches = margin_bottom_inches + margin_top_inches +\
                    nrow*(sub_margin_bottom_inches + sub_height_inches + sub_margin_top_inches)
                aspect = height_inches/width_inches
        else: # width_inches is None
            if not aspect is None:
                width_inches = height_inches/aspect
                sub_aspect = sub_height_inches/sub_width_inches
                sub_width_inches = (width_inches - margin_left_inches - margin_right_inches -\
                    ncol*(sub_margin_left_inches + sub_margin_right_inches))/ncol

            else:
                sub_width_inches = sub_height_inches/sub_aspect
                width_inches = margin_left_inches + margin_right_inches +\
                    ncol*(sub_margin_left_inches + sub_width_inches + sub_margin_right_inches)
                aspect = height_inches/width_inches
    else: # use the default
        sub_width_inches = sub_width_inches_default
        sub_height_inches = sub_height_inches_default
        sub_aspect = sub_height_inches/sub_width_inches
        width_inches = margin_left_inches + margin_right_inches +\
                ncol*(sub_margin_left_inches + sub_width_inches + sub_margin_right_inches)
        height_inches = margin_bottom_inches + margin_top_inches +\
                nrow*(sub_margin_bottom_inches + sub_height_inches + sub_margin_top_inches)
        sub_aspect = sub_height_inches/sub_width_inches

    #  now we're totally specified

    # collect all figure parameters in dictionary 
    fpar = dotdict()
    fpar['nplots'] = nplots
    fpar['nrow'] = nrow
    fpar['ncol'] = ncol
    fpar['width_inches'] = width_inches
    fpar['height_inches'] = height_inches
    fpar['aspect'] = aspect
    fpar['margin_left'] = margin_left_inches/width_inches
    fpar['margin_right'] = margin_right_inches/width_inches
    fpar['margin_bottom'] = margin_bottom_inches/height_inches
    fpar['margin_top'] = margin_top_inches/height_inches

    fpar['sub_width'] = sub_width_inches/width_inches
    fpar['sub_height'] = sub_height_inches/height_inches
    fpar['sub_aspect'] = sub_aspect
    fpar['sub_margin_left'] = sub_margin_left_inches/width_inches
    fpar['sub_margin_right'] = sub_margin_right_inches/width_inches
    fpar['sub_margin_bottom'] = sub_margin_bottom_inches/height_inches
    fpar['sub_margin_top'] = sub_margin_top_inches/height_inches

    # Generate the figure + axes
    fig = plt.figure(figsize=(width_inches, height_inches))
    #axs = np.zeros((nrow, ncol), dtype=object)
    axs = []
    for iplot in range(fpar['nplots']):
        icol = iplot%fpar['ncol']
        irow = iplot//fpar['ncol']
        ax_left = fpar['margin_left'] + fpar['sub_margin_left'] + icol*(fpar['sub_width'] + fpar['sub_margin_left'] + fpar['sub_margin_right'])
        ax_bottom = 1.0 - fpar['margin_top'] - fpar['sub_height'] - fpar['sub_margin_top'] - irow*(fpar['sub_height'] + fpar['sub_margin_top'] + fpar['sub_margin_bottom'])
        axs.append(fig.add_axes((ax_left, ax_bottom, fpar['sub_width'], fpar['sub_height'])))

        #axs[irow,icol] = fig.add_axes((ax_left, ax_bottom, fpar['sub_width'], fpar['sub_height']))
    axs  = np.reshape(axs, (nrow, ncol))

    return fig, axs, fpar

lineplot_minmax_kwargs_default = dict({'logscale': False, 'buff_ignore': None, 'plotleg': False, 'legfrac': None, 'symmetrize': False, 'domain_bounds': None, 'ixcut': 0})
def lineplot_minmax(xx, profiles, **kwargs):
    kw = update_dict(lineplot_minmax_kwargs_default, kwargs)
    find_bad_keys(lineplot_minmax_kwargs_default, kwargs, 'lineplot_minmax')

    # possibly ignore nastiness around domain bounds
    if not kw.buff_ignore is None:
        if kw.domain_bounds is None: # just look at the two ends
            kw.domain_bounds = np.min(xx), np.max(xx)
        delta_x = np.max(xx) - np.min(xx)
        profiles_old = profiles.copy()
        profiles = []
        for profile_old in profiles_old:
            tmp = []
            for ix in range(len(xx)):
                x_loc = xx[ix]
                # check if x_loc is in a "bad" location 
                # (near the domain_bounds)
                add_it = True
                for domain_bound in kw.domain_bounds:
                    if abs(x_loc - domain_bound) < kw.buff_ignore*delta_x:
                        add_it = False
                if add_it:
                    tmp.append(profile_old[ix])
            profiles.append(np.array(tmp))
                
    mmin = np.infty
    mmax = -np.infty
    for profile in profiles:
        mmin = min(np.min(profile[kw.ixcut:]), mmin)
        mmax = max(np.max(profile), mmax)
    if kw.plotleg:
        if kw.legfrac is None: # legfrac is how much of plot (y dimensions)
        # the legend should take up
            kw.legfrac = 1/3
        buff_frac_min = kw.legfrac/(1 - kw.legfrac) + buff_frac
    else:
        buff_frac_min = buff_frac

    if kw.logscale:
        yratio = mmax/mmin
        ymin = mmin/(yratio**buff_frac_min)
        ymax = mmax*(yratio**buff_frac)
    else:
        ydiff = mmax - mmin
        ymin = mmin - buff_frac_min*ydiff
        ymax = mmax + buff_frac*ydiff

    if kw.symmetrize:
        maxabs = max(abs(ymin), abs(ymax))
        ymin, ymax = -maxabs, maxabs

    # need to check for singular transormations (ymin = ymax)
    tol = 1.0e-100
    if np.abs(ymax - ymin) < tol: # it's a singular transformation!
        # first check for zero = zero singularity
        if np.abs(ymin) < tol and np.abs(ymax) < tol: # it's zero = zero
            ymin, ymax = -1.0, 1.0
        else: # ymin = ymax (but finite) singularity
            maxabs = np.abs(ymin)
            ymin, ymax = ymin - 0.5*maxabs, ymin + 0.5*maxabs

    return ymin, ymax

def get_symlog_params(field, field_max=None, sgnlog=False):
    if field_max is None:
        #maxabs = np.max(np.abs(field))
        #maxabs_exp = np.floor(np.log10(maxabs))
        #field_max = 10.**maxabs_exp
        field_max = np.max(np.abs(field))
    sig = np.std(field)
    linthresh = 0.15*sig
    dynamic_range = field_max/linthresh
    dynamic_range_decades = np.log10(dynamic_range)
    linscale = dynamic_range_decades
    if sgnlog:
        linscale *= 0.1
    return linthresh, linscale
    
def saturate_array(arr, my_min, my_max):
    arr[np.where(arr < my_min)] = my_min
    arr[np.where(arr > my_max)] = my_max

def sci_format(num, ndec=1, compact=False, nomant=False):
    exponent = get_exp(num)
    mantissa = num/10.**exponent
    if nomant:
        out = r'$10^{%i}$' %exponent
        if mantissa < 0:
            out = r'$-$' + out
        return out
    elif compact:
        return ( ('%1.' + ('%i' %ndec) +'fe%i')\
            %(mantissa, exponent))
    else:
        return ((r'$%1.' + (r'%i' %ndec) + r'f\times10^{%i}$')\
                %(mantissa, exponent))

lineplot_kwargs_default = dict({'xlabel': None, 'ylabel': None, 'title': None, 'xvals': np.array([]), 'yvals': np.array([]), 'labels': None, 'xlogscale': False, 'xminmax': None, 'minmax': None, 'xcut': None, 'minmax2': None, 'scatter': False, 'colors': color_order, 'linestyles': style_order[0], 'markers': marker_order[0], 'lw': default_lw, 's': default_s, 'ncolleg': 3, 'fontsize': default_labelsize, 'nosci': False, 'noscix': False, 'nosciy': False, 'legloc': 'lower left'})
lineplot_kwargs_default.update(lineplot_minmax_kwargs_default)

def lineplot(xx, profiles, ax, **kwargs):
    kw = update_dict(lineplot_kwargs_default, kwargs)
    kw_lineplot_minmax = update_dict(lineplot_minmax_kwargs_default, kwargs)
    find_bad_keys(lineplot_kwargs_default, kwargs, 'lineplot')

    # need to have a list of axes
    axs = [ax]

    # convert xvals, yvals to lists
    kw.xvals = make_array(kw.xvals, tolist=True)
    kw.yvals = make_array(kw.yvals, tolist=True)

    # deal with scalar 1-D arrays (everything must be "lists")
    if not isinstance(profiles, list):
        if isinstance(profiles, np.ndarray):
            if profiles.ndim == 1:
                profiles = [profiles]

    nprofiles = len(profiles)
    if np.isscalar(kw.markers):
        kw.markers = [kw.markers]*nprofiles
    if np.isscalar(kw.linestyles):
        kw.linestyles = [kw.linestyles]*nprofiles
    if np.isscalar(kw.lw):
        kw.lw = [kw.lw]*nprofiles
    if kw.labels is None:
        kw.labels = [None]*nprofiles

    if not kw.xcut is None:
        kw_lineplot_minmax.symmetrize = True
        ax2 = ax.twinx()
        axs.append(ax2)
        ixcut = np.argmin(np.abs(xx - kw.xcut))
        if xx[0] < xx[-1]: # x axis goes in "correct" direction
            ax_left = ax
            ax_right = ax2
        else: # x axis is reversed, so the second one is on the left
            ax_left = ax2
            ax_right = ax
        profiles2 = []
        for profile in profiles:
            profiles2.append(profile[ixcut:])
    else:
        kw_lineplot_minmax.symmetrize = False
        ixcut = len(xx)
        ax_left = ax

    # reset profiles based on ixcut
    profiles_old = profiles.copy()
    profiles = []
    for profile in profiles_old:
        profiles.append(profile[:ixcut])

    # get ylimits 
    if kw.minmax is None:
        kw.minmax = lineplot_minmax(xx[:ixcut], profiles, **kw_lineplot_minmax)
    ax.set_ylim(kw.minmax)

    if not kw.xcut is None:
        if kw.minmax2 is None:
            kw.minmax2 = lineplot_minmax(xx[ixcut:], profiles2, **kw_lineplot_minmax)
        ax2.set_ylim(kw.minmax2)
  
    # loop over profiles and make plots
    for iprof in range(nprofiles):
        profile = profiles[iprof]
        kw_scatter = dict({'label': kw.labels[iprof], 'marker': kw.markers[iprof], 'color': kw.colors[iprof], 's': kw.s})
        kw_plot = dict({'label': kw.labels[iprof], 'linestyle': kw.linestyles[iprof], 'color': kw.colors[iprof], 'linewidth': kw.lw[iprof]})
        if kw.xcut is None:
            if kw.scatter:
                ax.scatter(xx, profile, **kw_scatter)
            else:
                ax.plot(xx, profile, **kw_plot)
        else:
            profile2 = profiles2[iprof]
            if kw.scatter:
                ax.scatter(xx[:ixcut], profile2, **kw_scatter)
                ax2.scatter(xx[ixcut:], profile2, **kw_scatter)
            else:
                ax.plot(xx[:ixcut], profile, **kw_plot)
                ax2.plot(xx[ixcut:], profile2, **kw_plot)

    # ticks (mostly everywhere, deal with split axes)
    if kw.xcut is None:
        plt.sca(ax)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in', which='both')
        if not kw.ylabel is None:
            plt.ylabel(kw.ylabel, fontsize=kw.fontsize)
    else:
        plt.sca(ax_right)
        plt.minorticks_on()
        plt.tick_params(top=True, left=False, right=True, direction='in', which='both')
        ax_right.yaxis.tick_right()
        plt.sca(ax_left)
        plt.minorticks_on()
        plt.tick_params(top=True, left=True, right=False, direction='in', which='both')
        if not kw.ylabel is None:
            ax_left.set_ylabel(kw.ylabel, fontsize=kw.fontsize)
            ax_left.yaxis.set_label_position('left')
        ax_left.yaxis.tick_left()

    if kw.xlogscale:
        ax.set_xscale('log')
    if kw.logscale:
        for ax in axs:
            ax.set_yscale('log')

    if kw.xminmax is None:
        kw.xminmax = np.min(xx), np.max(xx)
    ax.set_xlim(kw.xminmax)

    if not kw.xcut is None:
        kw.xvals.append(xx[ixcut])
        kw.yvals.append(0)
 
    npoints = 100
    xpoints = np.linspace(kw.xminmax[0], kw.xminmax[1], npoints)

    # possibly mark x/y - values
    for xval in kw.xvals:
        if kw.xcut is None:
            ax_loc = ax
        else:
            if xval <= kw.xcut:
                ax_loc = ax_left
            else:
                ax_loc = ax_right
        y1, y2 = ax_loc.get_ylim()
        ypoints = np.linspace(y1, y2, npoints)
        ax_loc.plot(xval + np.zeros(npoints), ypoints, 'k--', linewidth=kw.lw[0])
    for yval in kw.yvals:
        for ax_loc in axs:
            # only plot if the line is within the range
            y1, y2 = ax_loc.get_ylim()
            if y1 < yval < y2:
                ax_loc.plot(xpoints, yval + np.zeros(npoints), 'k--', linewidth=kw.lw[0])

    if not kw.xlabel is None:
        ax.set_xlabel(kw.xlabel, fontsize=kw.fontsize)

    if not kw.title is None:
        ax.set_title(kw.title, fontsize=kw.fontsize)

    # set the tick label size
    for ax_loc in axs:
        plt.sca(ax_loc)
        plt.xticks(fontsize=kw.fontsize)
        plt.yticks(fontsize=kw.fontsize)

        # Get the non-log in scientific notation
        if not kw.nosci:
            if not kw.noscix:
                if not kw.xlogscale:
                    plt.ticklabel_format(useMathText=True, axis='x', scilimits=(0,0))
            if not kw.nosciy:
                if not kw.logscale:
                    plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

    # make the legend
    if kw.plotleg:
        ax.legend(loc=kw.legloc, ncol=kw.ncolleg, fontsize=0.8*default_labelsize)

add_cbar_kwargs_default = dict({'cbar_thick': 1/8, 'cbar_aspect': 1/20, 'cbar_prec': 2, 'cbar_no': 1, 'cbar_offset': None, 'cbar_pos': 'bottom', 'cbar_total_width': 1/2, 'units': '', 'nosci': False, 'cbar_fs': default_labelsize, 'tickvals': None, 'ticklabels': None, 'exp': 0, 'logscale': False, 'posdef': False, 'fullrange2': False, 'symlog': False, 'sgnlog': False, 'tol': 0.75, 'no0': False, 'cbar_label': None})
def add_cbar(fig, ax, im, **kwargs):
    # deal with kwargs
    kw = update_dict(add_cbar_kwargs_default, kwargs)
    find_bad_keys(add_cbar_kwargs_default, kwargs, 'add_cbar')

    # get fig dimensions
    fig_width_inches, fig_height_inches = fig.get_size_inches()
    fig_aspect = fig_height_inches/fig_width_inches
    # get ax dimensions
    ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
    ax_width = ax_right - ax_left
    ax_height = ax_top - ax_bottom

    if kw.sgnlog: # sgnlog is a special subcase of symlog
        kw.symlog = True

    if kw.cbar_pos == 'bottom':
        orientation = 'horizontal'
        cbar_height = kw.cbar_thick/fig_height_inches
        cbar_width = cbar_height/kw.cbar_aspect*fig_aspect
        cbar_width = min(cbar_width, kw.tol*ax_width) # don't let cbar be thicker than plot!
        # centrally position colorbar underneath the axes
        cbar_total_height = kw.cbar_total_width/fig_height_inches 
        # needs to contain
        # the colorbar ticklabels and little buffer space,
        # default 1/2
        if kw.cbar_offset is None:
            kw.cbar_offset = 1/16
        lilbit = kw.cbar_offset/fig_height_inches
        cbar_left = ax_left + 0.5*ax_width - 0.5*cbar_width
        cbar_bottom = ax_bottom - lilbit - cbar_height -\
                (kw.cbar_no - 1)*cbar_total_height
    elif kw.cbar_pos == 'right':
        orientation = 'vertical'
        cbar_width = kw.cbar_thick/fig_width_inches
        cbar_height = cbar_width/kw.cbar_aspect/fig_aspect
        cbar_height = min(cbar_height, kw.tol*ax_height) # don't let cbar be longer than plot!

        # centrally position colorbar to right of axes
        label_buff = 3/4/fig_width_inches # needs to contain
        # the colorbar ticklabels and little buffer space
        lilbit = 1/16/fig_width_inches
        cbar_bottom = ax_bottom + 0.5*ax_height - 0.5*cbar_height
        cbar_left = ax_right + lilbit + (kw.cbar_no - 1)*(label_buff + cbar_width)
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))

    cbar = plt.colorbar(im, cax=cax, orientation=orientation)

    # deal with labels
    if kw.cbar_label is None:
        if kw.nosci or kw.logscale or kw.symlog:
            kw.cbar_label = kw.units
        else:
            kw.cbar_label = (r'$\times10^{%i}\ $' %kw.exp) + kw.units

    # font size for the tick labels
    cax.tick_params(labelsize=kw.cbar_fs)
    #cbar.ax.tick_params(labelsize=fontsize)   

    # ticklabel format
    if kw.logscale: # set tickvals and labels through "smart locator"
        locator = ticker.LogLocator(subs='all')
        cbar.set_ticks(locator)
    else: # set tickvals and ticklabels separately, depending on norm
        arr = im.get_array()
        spacing = np.mean(np.diff(arr))
        # almost positive this behavior changed recently...
        # not sure what to do for symlog
        minn, maxx = get_minmax(arr)
        minn, maxx = minn - spacing/2., maxx + spacing/2.
        nneeded = len(arr) + 1
        levelsfield = np.linspace(minn, maxx, nneeded)
        nlevelsfield = len(levelsfield) - 1
        # first, tickvals
        if kw.tickvals is None:
            # just thin out actual field levels
            # hopefully field levels are a "nice" number:
            # 3 * (multiple of 4) for symlog
            # 2 * (multiple of 4) otherwise
            if kw.symlog:
                maxabs = levelsfield[-1]
                if kw.sgnlog: # ignore linthresh range ticks
                    # except invisible "zero" tick
                    linthresh = levelsfield[nlevelsfield//2]
                    kw.tickvals = np.hstack((\
                        np.linspace(-maxabs, -linthresh, 5), np.array([0.]),
                        np.linspace(linthresh, maxabs, 5) ))
                else:
                    linthresh = -levelsfield[nlevelsfield//3]
                    kw.tickvals = np.hstack((\
                        np.linspace(-maxabs, -linthresh, 4, endpoint=False),
                        np.linspace(-linthresh, linthresh, 4, endpoint=False),
                        np.linspace(linthresh, maxabs, 5) ))
            elif kw.fullrange2:
                kw.tickvals = np.array([levelsfield[0], 0., levelsfield[-1]])
            else:
                nskip = nlevelsfield//8
                kw.tickvals = levelsfield[::nskip]
        # then, ticklabels
        if kw.ticklabels is None:
            nticks = len(kw.tickvals)
            kw.ticklabels = ['']*nticks
            if kw.symlog: 
                if kw.sgnlog: # want "plus/minus linthresh" offset by 1
                    # from actual linthresh
                    ind = (nticks-1)//2 # location of zero tick
                    kw.ticklabels[ind] = r'$\pm$' + sci_format(kw.tickvals[ind+1], ndec=kw.cbar_prec, compact=True)
                    for ind in [0, nticks - 1]:
                        kw.ticklabels[ind] = sci_format(kw.tickvals[ind], ndec=kw.cbar_prec, compact=True)
                else: # -linthresh and maxabs
                    indvals = [(nticks-1)//3, nticks-1]
                    for ind in indvals:
                        kw.ticklabels[ind] = sci_format(kw.tickvals[ind], ndec=kw.cbar_prec, compact=True)
            else:
                if kw.posdef or kw.no0:
                    # just min/max
                    indvals = [0, nticks-1]
                else:
                    # min/max and zero
                    indvals = [0, (nticks-1)//2, nticks-1]

                # loop over indices and set values
                for ind in indvals:
                    fmt = '%.' + str(kw.cbar_prec) + 'f'
                    kw.ticklabels[ind] = fmt %kw.tickvals[ind]
                if kw.fullrange2: # remove zero tick
                    kw.ticklabels[1] = ''
        cbar.set_ticks(kw.tickvals)
        cbar.set_ticklabels(kw.ticklabels)

    if kw.cbar_pos == 'bottom':
        fig.text(cbar_left + cbar_width + 1/16/fig_width_inches,\
                cbar_bottom + 0.5*cbar_height, kw.cbar_label,\
                ha='left', va='center', fontsize=kw.cbar_fs) 
    elif kw.cbar_pos == 'right':
        #fig.text(cbar_left + cbar_width + lilbit/fig_aspect,\
        #        cbar_bottom + 0.5*cbar_height, kw.cbar_label,\
        #        ha='left', va='center', fontsize=kw.fontsize) 
        cax.set_title(kw.cbar_label, ha='left', fontsize=kw.cbar_fs)

contourf_minmax_kwargs_default = dict({'posdef': False, 'no0': False, 'logscale': False, 'symlog': False, 'sgnlog': False, 'fullrange': False, 'fullrange2': False, 'buff_ignore1': buff_frac, 'buff_ignore2': buff_frac}) 

def contourf_minmax(field, **kwargs):
    # Get good boundaries to saturate array [field], assuming either
    # posdef (True or False) and/or logscale (True or False)
    # first, possibly cut the array (ignore boundary vals)
    kw = update_dict(contourf_minmax_kwargs_default, kwargs)
    find_bad_keys(contourf_minmax_kwargs_default, kwargs, 'contourf_minmax')

    if not kw.buff_ignore1 is None:
        n1, dummy = np.shape(field)
        icut = int(n1*kw.buff_ignore1)
        field = np.copy(field[icut:n1-icut, :])
    if not kw.buff_ignore2 is None:
        dummy, n2 = np.shape(field)
        icut = int(n2*kw.buff_ignore2)
        field = np.copy(field[:, icut:n2-icut])
    # purge the field of nans
    field = field[np.where(1 - np.isnan(field))]
    if kw.fullrange or kw.symlog:
        maxabs = np.max(np.abs(field))
        minmax = -maxabs, maxabs       
    elif kw.fullrange2:
        mmin = np.min(field[np.where(field != 0.0)])
        mmax = np.max(field[np.where(field != 0.0)])
        minmax = mmin, mmax
    elif kw.logscale:
        logfield = np.log(field[np.where(field != 0)])
        meanlog = np.mean(logfield)
        stdlog = np.std(logfield)

        minexp = meanlog - 3.*stdlog
        maxexp = meanlog + 3.*stdlog
        minmax = np.exp(minexp), np.exp(maxexp)        
    elif kw.posdef:
        sig = rms(field)
        minmax = 0., 3.*sig        
    elif kw.no0:
        mean = np.mean(field)
        sig = np.std(field)
        minmax = mean - 3*sig, mean + 3*sig
    else:
        sig = np.std(field)
        minmax = -3.*sig, 3.*sig
    # Make sure minmax isn't 0, 0
    tinybit = 1.0e-100
    minmax = minmax[0] - tinybit, minmax[1] + tinybit
    return minmax

my_contourf_kwargs_default = dict({
        # basic flags:
         'plotfield': True,
         'plotcontours': True, 'ncontours': 8, 'contourlevels': None, 'contourstyles': '--', 'contourcolors': 'k', 'contourwidths': default_lw, 
        # colorbar stuff
        'plotcbar': True, 'cmap': None, 'norm': None, 'units': '', 'nlevels': 64,       
        # only need this for time-lat plots or such, since need ticks there
        'allticksoff': True,\
        # symlog stuff (again)
        'minmax': None, 'linthresh': None, 'linscale': None, 'scaleby': 1.0})

# color map stuff: symlog, logscale, posdef, are here:
my_contourf_kwargs_default.update(contourf_minmax_kwargs_default)
my_contourf_kwargs_default.update(add_cbar_kwargs_default)

def my_contourf(xx, yy, field, fig, ax, **kwargs):
    kw = update_dict(my_contourf_kwargs_default, kwargs)
    kw_contourf_minmax = update_dict(contourf_minmax_kwargs_default, kwargs)
    kw_add_cbar = update_dict(add_cbar_kwargs_default, kwargs)
    find_bad_keys(my_contourf_kwargs_default, kwargs, 'my_contourf')

    # make sure Python does not modify any of the arrays it was passed
    field = np.copy(field)/kw.scaleby

    if kw.sgnlog: # sgnlog is a special subcase of symlog
        kw.symlog = True

    # get default bounds if not specified
    if kw.minmax is None:
        kw.minmax = contourf_minmax(field, **kw_contourf_minmax)

    # plot the field, maybe
    # Factor out the exponent on the field and put it on the color bar
    # can turn this behavior off with "nosci=True"
    if not (kw.nosci or kw.logscale or kw.symlog):
        maxabs = max(np.abs(kw.minmax[0]), np.abs(kw.minmax[1]))
        kw_add_cbar.exp = get_exp(maxabs)
        divisor = 10.0**kw_add_cbar.exp
        field /= divisor
        kw.minmax = kw.minmax[0]/divisor, kw.minmax[1]/divisor

    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, kw.minmax[0], kw.minmax[1])

    # deal with norm
    if kw.norm is None and kw.logscale:
        kw.norm = colors.LogNorm(vmin=kw.minmax[0], vmax=kw.minmax[1])

    if kw.symlog:
        # levels determined by linscale and linthresh
        linthresh_default, linscale_default =\
            get_symlog_params(field, field_max=kw.minmax[1], sgnlog=kw.sgnlog)
        if kw.linthresh is None:
            kw.linthresh = linthresh_default
        if kw.linscale is None:
            kw.linscale = linscale_default
        log_thresh = np.log10(kw.linthresh)
        log_max = np.log10(kw.minmax[1])*0.999 # reduce by a bit
        # (otherwise creates white-filled space outside max contours)

        # special symlog norm
        kw.norm = colors.SymLogNorm(linthresh=kw.linthresh,\
            linscale=kw.linscale, vmin=kw.minmax[0], vmax=kw.minmax[1])

        if kw.sgnlog: # ignore linear area (make it one contour level)
            levels_neg = -np.logspace(log_max, log_thresh, kw.nlevels+1)
            levels_pos = np.logspace(log_thresh, log_max, kw.nlevels+1)
            levels = np.hstack((levels_neg, levels_pos))
        else: # include linear range (multiple contour levels)
            levels_neg = -np.logspace(log_max, log_thresh, kw.nlevels, endpoint=False)
            levels_mid = np.linspace(-kw.linthresh, kw.linthresh, kw.nlevels, endpoint=False)
            levels_pos = np.logspace(log_thresh, log_max, kw.nlevels+1)
            levels = np.hstack((levels_neg, levels_mid, levels_pos))
    elif kw.logscale:
        levels = np.logspace(np.log10(kw.minmax[0]), np.log10(kw.minmax[1]), 2*kw.nlevels+1)
    else:
        levels = np.linspace(kw.minmax[0], kw.minmax[1], 2*kw.nlevels+1)
        if kw.fullrange2: # need equally spaced levels on each side of zero
            kw.norm = colors.TwoSlopeNorm(vmin=kw.minmax[0], vcenter=0, vmax=kw.minmax[1])

    if kw.cmap is None:
        if kw.posdef:
            kw.cmap = 'plasma'
        elif kw.logscale:
            kw.cmap = 'Greys'
        else:
            kw.cmap = 'RdYlBu_r'

    # finally we make the contour plot!
    if kw.plotfield:
        im = ax.contourf(xx, yy, field, cmap=kw.cmap, levels=levels, norm=kw.norm)

    # now deal with color bar, if one is desired
    if kw.plotfield and kw.plotcbar:
        add_cbar(fig, ax, im, **kw_add_cbar)

    # Plot contours if desired
    if kw.plotcontours:
        # Determine the contour levels
        if kw.contourlevels is None:
            # just thin out the field levels
            nskip = kw.nlevels//kw.ncontours
            kw.contourlevels = levels[::nskip]
            if kw.fullrange2: # need equally spaced levels on 
                # each side of zero separately
                levels_neg = np.linspace(kw.minmax[0], 0, kw.ncontours//2 - 1, endpoint=False)
                levels_pos = np.linspace(0, kw.minmax[1], kw.ncontours//2, endpoint=False)
                kw.contourlevels = np.hstack((levels_neg, levels_pos))

        # plot the contours
        ax.contour(xx, yy, field, kw.contourlevels, norm=kw.norm,\
                colors=kw.contourcolors, linewidths=kw.contourwidths, linestyles=kw.contourstyles)


    if kw.allticksoff:
        # Set ax ranges to be just outside the boundary lines
        # avoid weird whitespace cutoffs
        # then turn off ticks
        lilbit = 0.01
        xmin, xmax = np.min(xx), np.max(xx)
        ymin, ymax = np.min(yy), np.max(yy)
        Dx = xmax - xmin
        Dy = ymax - ymin
        ax.set_xlim((xmin - lilbit*Dx, xmax + lilbit*Dx))
        ax.set_ylim((ymin - lilbit*Dy, ymax + lilbit*Dy))
        ax.axis('off') 
    
    return  fig, ax

# my_pcolormesh: pixellated my_contourf
my_pcolormesh_kwargs_default = dict({'x': None, 'y': None, 'xvals': [], 'yvals': [], 'linewidth': default_lw, 'minmax': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'yminmax': None, 'ymin': None, 'ymax': None, 'xymin': None, 'xymax': None, 'xyminmax': None,\
    # more cbar stuff
    'plotcbar': True, 'cmap': None, 'norm': None, 'units': '', 'scaleby': 1., 'fontsize': default_labelsize})

my_pcolormesh_kwargs_default.update(contourf_minmax_kwargs_default)
my_pcolormesh_kwargs_default.update(add_cbar_kwargs_default)
def my_pcolormesh(field, fig, ax, **kwargs):
    kw = update_dict(my_pcolormesh_kwargs_default, kwargs)
    kw_add_cbar = update_dict(add_cbar_kwargs_default, kwargs)
    kw_contourf_minmax = update_dict(contourf_minmax_kwargs_default, kwargs)
    find_bad_keys(my_pcolormesh_kwargs_default, kwargs, 'my_pcolormesh')

    # make sure Python does not modify any of the arrays it was passed
    field = np.copy(field)/kw.scaleby

    # use full integer grid if no specific axis was specified
    nx, ny = np.shape(field)
    if kw.x is None:
        kw.x = np.arange(nx)
    if kw.y is None:
        kw.y = np.arange(ny)

    # by default plot whole spectrum
    ix1, ix2 = 0, nx - 1
    iy1, iy2 = 0, ny - 1

    # might set both axis boundaries at once with "xy" min/max
    if not kw.xymin is None:
        kw.xmin = kw.xymin
        kw.ymin = kw.xymin
    if not kw.xymax is None:
        kw.xmax = kw.xymax
        kw.ymax = kw.xymax
    if not kw.xyminmax is None:
        kw.xminmax = kw.xyminmax
        kw.yminmax = kw.xyminmax

    # pick and choose part of spectrum to plot
    if not kw.xminmax is None:
        ix1 = np.argmin(np.abs(kw.x - kw.xminmax[0]))
        ix2 = np.argmin(np.abs(kw.x - kw.xminmax[1]))
    if not kw.xmin is None:
        ix1 = np.argmin(np.abs(kw.x - kw.xmin))
    if not kw.xmax is None:
        ix2 = np.argmin(np.abs(kw.x - kw.xmax))

    if not kw.yminmax is None:
        iy1 = np.argmin(np.abs(kw.y - kw.yminmax[0]))
        iy2 = np.argmin(np.abs(kw.y - kw.yminmax[1]))
    if not kw.ymin is None:
        iy1 = np.argmin(np.abs(kw.y - kw.ymin))
    if not kw.ymax is None:
        iy2 = np.argmin(np.abs(kw.y - kw.ymax))

    # now adjust everything by the (x, y) range we want
    kw.x = kw.x[ix1:ix2+1]
    kw.y = kw.y[iy1:iy2+1]
    field = field[ix1:ix2+1, iy1:iy2+1]

    # get the squares to correspond to modes
    xx, yy = xy_grid(kw.x, kw.y)

    # get default bounds if not specified
    if kw.minmax is None:
        kw.minmax = contourf_minmax(field, **kw_contourf_minmax)

    # Factor out the exponent on the field and put it on the color bar
    # can turn this behavior off with "nosci=True"
    if not (kw.nosci or kw.logscale):
        maxabs = max(np.abs(kw.minmax[0]), np.abs(kw.minmax[1]))
        kw_add_cbar.exp = get_exp(maxabs)
        divisor = 10**kw_add_cbar.exp
        field /= divisor
        kw.minmax = kw.minmax[0]/divisor, kw.minmax[1]/divisor
    else:
        divisor = 1.0

    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, kw.minmax[0], kw.minmax[1])
  
    # deal with norm and colormap
    if kw.norm is None and kw.logscale:
        kw.norm = colors.LogNorm(vmin=kw.minmax[0], vmax=kw.minmax[1])
        vmin, vmax = None, None
    else:
        vmin, vmax = kw.minmax

    if kw.fullrange2: # need equally spaced levels on each side of zero
        kw.norm = colors.TwoSlopeNorm(vmin=kw.minmax[0], vcenter=0, vmax=kw.minmax[1])

    if kw.cmap is None:
        #kw.cmap = 'jet'
        kw.cmap = 'RdYlBu_r'

    # make color plot
    im = ax.pcolormesh(xx, yy, field, cmap=kw.cmap, norm=kw.norm, vmin=vmin, vmax=vmax)  

    # now deal with color bar, if one is desired
    if kw.plotcbar:
        add_cbar(fig, ax, im, **kw_add_cbar)

    # set bounds
    space_x = kw.x[1] - kw.x[0]
    space_y = kw.y[1] - kw.y[0]
    ax.set_xlim(kw.x[0] - 0.5*space_x, kw.x[-1] + 0.5*space_x)
    ax.set_ylim(kw.y[0] - 0.5*space_y, kw.y[-1] + 0.5*space_y)

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # possibly plot x and y vals
    xvals = make_array(kw.xvals)
    yvals = make_array(kw.yvals)
    count = 1
    for axisvals in xvals, yvals:
        for axisval in axisvals:
            if count == 1:
                xline, yline = axisval + np.zeros_like(kw.y), kw.y
            if count == 2:
                xline, yline = kw.x, axisval + np.zeros_like(kw.x)
            ax.plot(xline, yline, 'k--')
        count += 1
    return kw.minmax[0]*divisor, kw.minmax[1]*divisor

def mark_axis_vals(ax, which='x', vals=[0.], style='k-', lw=0.5):
    n = 100
    zero = np.zeros(n)
    if which == 'x':
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax)
        yrange = np.linspace(ymin, ymax, n)
        for val in vals:
            ax.plot(val + zero, yrange, style, linewidth=lw)
    elif which == 'y':
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmin, xmax)
        xrange = np.linspace(xmin, xmax, n) # is this a built-in?
        for val in vals:
            ax.plot(xrange, val + zero, style, linewidth=lw)
