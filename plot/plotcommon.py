# Module for routines common to many plotting scripts
# Created: 02/08/2019
import numpy as np
from matplotlib import colors, ticker
import matplotlib.pyplot as plt
from common import *

color_order = ['b', 'g', 'r', 'c', 'm', 'y', 'k']*2
style_order = ['-', '--', '-.', ':']
marker_order = [".", "o", "v","s", "*", "x", "^", "<", ">"]
default_lw = 1.0
default_s = 0.2 # markersize
default_labelsize = 12
default_titlesize = 12
default_ticksize = 12
default_margin = 1/16
default_line_height = 1/4 # height of a line of text
default_margin_xlabel = 1/2
default_margin_ylabel = 3/4
# ylabels take up more space because floating
# point numbers are longer than they are tall
default_margin_title = 3/4

# default figure sizes 
sub_width_inches_default, sub_height_inches_default = 3.5, 2.5

def axis_range(ax): # gets subplot coordinates on a figure in "normalized"
        # coordinates
    pos = plt.get(ax, 'position')
    bottom_left = pos.p0
    top_right = pos.p1
    xmin, xmax = bottom_left[0], top_right[0]
    ymin, ymax = bottom_left[1], top_right[1]
    return xmin, xmax, ymin, ymax

def xy_grid(X, Y):
    """
    plt.pcolormesh() takes arguments X, Y, and C; treats the X, Y 
    arrays as VERTICES of quadrilaterals (not center points)
    For plotting, we usually have X, Y, and C of shape (m, n) and
    (X[i, j], Y[i, j]) represents the CENTER of a quadrilateral
    Thus, a call to plt.pcolormesh will ignore the last row/column
    of C, using only the smaller array C[:m-1, :n-1]
    "xy_grid" makes new arrays X_new, Y_new, which have dimension
    (m+1, n+1) and represent the vertices of quadrilaterals with centers
    at X[i, j]
    """
    m, n = np.shape(X)
    X_mid = 0.5*(X[:m-1, :n-1] + X[1:, :n-1])
    Y_mid = 0.5*(Y[:m-1, :n-1] + Y[:m-1, 1:])
    X_new, Y_new = np.zeros((m+1, n+1)), np.zeros((m+1, n+1))
    X_new[1:m, 1:n] = X_mid
    Y_new[1:m, 1:n] = Y_mid

    X_new[1:m, 0] = X_new[1:m, 1] - (X_new[1:m, 2] - X_new[1:m, 1])
    X_new[1:m, n] = X_new[1:m, n-1] + (X_new[1:m, n-1] - X_new[1:m, n-2])
    Y_new[1:m, 0] = Y_new[1:m, 1] - (Y_new[1:m, 2] - Y_new[1:m, 1])
    Y_new[1:m, n] = Y_new[1:m, n-1] + (Y_new[1:m, n-1] - Y_new[1:m, n-2])

    X_new[0, :] = X_new[1, :] - (X_new[2, :] - X_new[1, :])
    X_new[m, :] = X_new[m-1, :] + (X_new[m-1, :] - X_new[m-2, :])
    Y_new[0, :] = Y_new[1, :] - (Y_new[2, :] - Y_new[1, :])
    Y_new[m, :] = Y_new[m-1, :] + (Y_new[m-1, :] - Y_new[m-2, :])
    return (X_new, Y_new)

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

    # if nplots is specified but others aren't, ncol = min(nplots, 3)
    # and nrow follows
    # 
    # if nplots nrow or ncol specified but others aren't, nplots = (nrow or ncol) and ncol or nrow = 1
    # 
    # if two are specified, the other follows
    # 
    # if all three are specified, user overdetermined the problem
    # and the combination (nplots, ncol) determines nrow

    nspec = 0
    for val in [nplots, nrow, ncol]:
        nspec += not val is None
    if nspec == 0: # all unspecified
        nplots = ncol = nrow = 1
    if nspec == 1: # 1 is specified
        if not nplots is None:
            ncol = min(3, nplots)
            nrow = int(np.ceil(nplots/ncol))
        elif not nrow is None:
            nplots = nrow
            ncol = 1
        else:
            nplots = ncol
            nrow = 1
    if nspec >= 2: # two or more are specified
        if nplots is None:
            nplots = nrow*ncol
        elif ncol is None:
            ncol = int(np.ceil(nplots/nrow))
        else:
            nrow = int(np.ceil(nplots/ncol))

    sub_nspec = 0
    for val in [sub_width_inches, sub_height_inches, sub_aspect]:
        sub_nspec += not val is None
    nspec = 0
    for val in [width_inches, height_inches, aspect]:
        nspec += not val is None

    # logic similar to above for nplots, nrow, ncol
    if nspec == sub_nspec == 0:
        which_spec = 'sub'
        sub_width_inches = sub_width_inches_default
        sub_height_inches = sub_height_inches_default
        sub_aspect = sub_height_inches/sub_width_inches
    elif sub_nspec == 1:
        which_spec = 'sub'
        if not sub_width_inches is None:
            sub_height_inches = sub_height_inches_default
            sub_aspect = sub_height_inches/sub_width_inches
        elif not sub_height_inches is None:
            sub_width_inches = sub_width_inches_default
            sub_aspect = sub_height_inches/sub_width_inches
        elif not sub_aspect is None:
            sub_width_inches = sub_width_inches_default
            sub_height_inches = sub_aspect*sub_width_inches

    elif sub_nspec >= 2:
        which_spec = 'sub'
        if sub_width_inches is None:
            sub_width_inches = sub_height_inches/sub_aspect
        elif sub_height_inches is None:
            sub_height_inches = sub_width_inches*sub_aspect
        else: # width and height --> aspect in overdetermined problem
            sub_aspect = sub_height_inches/sub_width_inches

    elif nspec >= 2:
        which_spec = 'fig'
        if width_inches is None:
            width_inches = height_inches/aspect
        elif height_inches is None:
            height_inches = width_inches*aspect
        else:
            aspect = height_inches/width_inches

    # set the "other" parameters based on if fig or sub dimensions were set
    if which_spec == 'fig':
        sub_width_inches = (width_inches - margin_left_inches - margin_right_inches - ncol*(sub_margin_left_inches + sub_margin_right_inches))/ncol
        sub_height_inches = (height_inches - margin_bottom_inches - margin_top_inches - nrow*(sub_margin_bottom_inches + sub_margin_top_inches))/nrow
        sub_aspect = sub_height/sub_width

    if which_spec == 'sub':
        width_inches = ncol*(sub_width_inches + sub_margin_left_inches + sub_margin_right_inches) + margin_left_inches + margin_right_inches
        height_inches = nrow*(sub_height_inches + sub_margin_bottom_inches + sub_margin_top_inches) + margin_bottom_inches + margin_top_inches
        aspect = height_inches/width_inches

    # collect all figure parameters in dictionary 
    fpar = dict({})
    fpar['nplots'] = nplots
    fpar['nrow'] = nrow
    fpar['ncol'] = ncol
    fpar['width_inches'] = width_inches
    fpar['height_inches'] = height_inches
    fpar['aspect'] = height_inches/width_inches
    fpar['margin_left'] = margin_left_inches/width_inches
    fpar['margin_right'] = margin_right_inches/width_inches
    fpar['margin_bottom'] = margin_bottom_inches/height_inches
    fpar['margin_top'] = margin_top_inches/height_inches

    fpar['sub_width'] = sub_width_inches/width_inches
    fpar['sub_height'] = sub_height_inches/height_inches
    fpar['sub_aspect'] = sub_height_inches/sub_width_inches
    fpar['sub_margin_left'] = sub_margin_left_inches/width_inches
    fpar['sub_margin_right'] = sub_margin_right_inches/width_inches
    fpar['sub_margin_bottom'] = sub_margin_bottom_inches/height_inches
    fpar['sub_margin_top'] = sub_margin_top_inches/height_inches

    # Generate the figure + axes
    fig = plt.figure(figsize=(width_inches, height_inches))
    axs = np.zeros((nrow, ncol), dtype=object)
    for iplot in range(fpar['nplots']):
        icol = iplot%fpar['ncol']
        irow = iplot//fpar['ncol']
        ax_left = fpar['margin_left'] + fpar['sub_margin_left'] + icol*(fpar['sub_width'] + fpar['sub_margin_left'] + fpar['sub_margin_right'])
        ax_bottom = 1.0 - fpar['margin_top'] - fpar['sub_height'] - fpar['sub_margin_top'] - irow*(fpar['sub_height'] + fpar['sub_margin_top'] + fpar['sub_margin_bottom'])
        axs[irow,icol] = fig.add_axes((ax_left, ax_bottom, fpar['sub_width'], fpar['sub_height']))

    return fig, axs, fpar

contourf_minmax_kwargs_default = dict({'posdef': False, 'logscale': False, 'fullrange': False, 'buff_ignore1': buff_frac, 'buff_ignore2': buff_frac}) 
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
    if kw.logscale:
        logfield = np.log(field)
        medlog = np.median(logfield)
        shiftlog = logfield - medlog
        std_plus =\
            np.std(shiftlog[np.where(shiftlog > 0.)].flatten())
        std_minus =\
            np.std(shiftlog[np.where(shiftlog <= 0.)].flatten())
        av_std = (std_plus + std_minus)/2.

        minexp = medlog - 5.*av_std
        maxexp = medlog + 5.*av_std
        minmax = np.exp(minexp), np.exp(maxexp)        
    elif kw.posdef:
        sig = rms(field)
        minmax = 0., 3.*sig        
    elif kw.fullrange:
        maxabs = np.max(np.abs(field))
        minmax = -maxabs, maxabs       
    else:
        sig = np.std(field)
        minmax = -3.*sig, 3.*sig
    # Make sure minmax isn't 0, 0
    tinybit = 1.0e-100
    minmax = minmax[0] - tinybit, minmax[1] + tinybit
    return minmax

lineplot_minmax_kwargs_default = dict({'logscale': False, 'buff_ignore': buff_frac, 'legfrac': None, 'symmetrize': False, 'domain_bounds': None, 'xx': None})
def lineplot_minmax(profiles, **kwargs):
    kw = update_dict(lineplot_minmax_kwargs_default, kwargs)
    find_bad_keys(lineplot_minmax_kwargs_default, kwargs, 'lineplot_minmax')

    # possibly ignore nastiness around domain bounds
    # x axis (xx) must also be provided
    if not kw.domain_bounds is None:
        delta_x = np.max(kw.xx) - np.min(kw.xx)
        profiles_old = profiles.copy()
        profiles = []
        for profile_old in profiles_old:
            tmp = []
            for ix in range(len(kw.xx)):
                x_loc = kw.xx[ix]
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
        mmin = min(np.min(profile), mmin)
        mmax = max(np.max(profile), mmax)
    if not kw.legfrac is None: # legfrac is how much of plot (y dimensions)
        # the legend should take up
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
    return ymin, ymax

def get_symlog_params(field, field_max=None):
    if field_max is None:
        maxabs = np.max(np.abs(field))
        maxabs_exp = np.floor(np.log10(maxabs))
        field_max = 10.**maxabs_exp
    sig = np.std(field)
    linthresh = 0.15*sig
    dynamic_range = field_max/linthresh
    dynamic_range_decades = np.log10(dynamic_range)
    linscale = dynamic_range_decades
    return linthresh, linscale
    
def saturate_array(arr, my_min, my_max):
    arr[np.where(arr < my_min)] = my_min
    arr[np.where(arr > my_max)] = my_max

def get_exp(num):
    if num != 0.:
        return int(np.floor(np.log10(np.abs(num))))
    else:
        return 1

def sci_format(num, ndec=1):
    exponent = get_exp(num)
    mantissa = num/10.**exponent
    return ((r'$%1.' + (r'%i' %ndec) + r'f\times10^{%i}$')\
            %(mantissa, exponent))

lineplot_kwargs_default = dict({'xlabel': None, 'ylabel': None, 'title': None, 'xvals': np.array([]), 'yvals': np.array([]), 'labels': None, 'xlogscale': False, 'xminmax': None, 'minmax': None, 'xcut': None, 'minmax2': None, 'scatter': False, 'colors': color_order, 'linestyles': style_order[0], 'markers': marker_order[0], 'lw': default_lw, 's': default_s, 'plotleg': True})
lineplot_kwargs_default.update(lineplot_minmax_kwargs_default)

del lineplot_kwargs_default['xx']

def lineplot(xx, profiles, ax, **kwargs):
    kw = update_dict(lineplot_kwargs_default, kwargs)
    # make room for legend by default
    if kw['plotleg']:
        lineplot_minmax_kwargs_default['legfrac'] = 1/3
    kw_lineplot_minmax = update_dict(lineplot_minmax_kwargs_default, kwargs)
    find_bad_keys(lineplot_kwargs_default, kwargs, 'lineplot')

    # need to have a list of axes
    axs = [ax]

    # convert xvals, yvals to lists
    kw.xvals = make_array(kw.xvals, tolist=True)
    kw.yvals = make_array(kw.yvals, tolist=True)

    # deal with scalar "lists"
    nprofiles = len(profiles)
    if np.isscalar(kw.markers):
        kw.markers = [kw.markers]*nprofiles
    if np.isscalar(kw.linestyles):
        kw.linestyles = [kw.linestyles]*nprofiles
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
        kw_lineplot_minmax.xx = xx[:ixcut]
        kw.minmax = lineplot_minmax(profiles, **kw_lineplot_minmax)
    ax.set_ylim(kw.minmax)

    if not kw.xcut is None:
        if kw.minmax2 is None:
            kw_lineplot_minmax.xx = xx[ixcut:]
            kw.minmax2 = lineplot_minmax(profiles2, **kw_lineplot_minmax)
        ax2.set_ylim(kw.minmax2)
  
    # loop over profiles and make plots
    for iprof in range(nprofiles):
        profile = profiles[iprof]
        kw_scatter = dict({'label': kw.labels[iprof], 'marker': kw.markers[iprof], 'color': kw.colors[iprof], 's': kw.s})
        kw_plot = dict({'label': kw.labels[iprof], 'linestyle': kw.linestyles[iprof], 'color': kw.colors[iprof], 'linewidth': kw.lw})
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
        plt.tick_params(top=True, right=True, direction='in',\
                which='both')
        if not kw.ylabel is None:
            plt.ylabel(kw.ylabel, fontsize=kw.fontsize)
    else:
        plt.sca(ax_right)
        plt.minorticks_on()
        plt.tick_params(top=True, left=False, right=True, direction='in', which='both')
        ax_right.yaxis.tick_right()
        plt.sca(ax_left)
        plt.minorticks_on()
        plt.tick_params(top=True, left=True, right=False, direction='in',\
                which='both')
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
        ax_loc.plot(xval + np.zeros(npoints), ypoints, 'k--', linewidth=kw.lw)
    for yval in kw.yvals:
        for ax_loc in axs:
            # only plot if the line is within the range
            y1, y2 = ax_loc.get_ylim()
            if y1 < yval < y2:
                ax_loc.plot(xpoints, yval + np.zeros(npoints), 'k--', linewidth=kw.lw)

    if not kw.xlabel is None:
        ax.set_xlabel(kw.xlabel, fontsize=kw.fontsize)

    if not kw.title is None:
        ax.set_title(kw.title, fontsize=kw.fontsize)

    # set the tick label size
    for ax_loc in axs:
        plt.sca(ax_loc)
        plt.xticks(fontsize=kw.fontsize)
        plt.yticks(fontsize=kw.fontsize)

        # Get the y-axis in scientific notation
        plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

    # make the legend
    if kw.plotleg:
        ax.legend(loc='lower left', ncol=3, fontsize=0.8*default_labelsize)

my_contourf_kwargs_default = dict({
        # saturation of field values stuff
        'minmax': None,    
        # basic flags:
         'plotfield': True,\
        'plotcontours': True, 'ncontours': 8, 'contourlevels': None,\
        # colorbar stuff
        'plotcbar': True, 'cbar_thick': 1/16, 'cbar_aspect': 1/20, 'cbar_prec': 2, 'cbar_no': 1, 'cbar_pos': 'bottom', 'cmap': None, 'units': '', 'nosci': False, 'fontsize': default_labelsize,\
        # coordinate line stuff; do up to two "types"
        'vals1': np.array([]), 'func1': None, 'vals2': np.array([]), 'func2': None,\
                'plotboundary': True, 'lw': 1.,\
        # only need this for time-lat plots or such, since need ticks there
        'allticksoff': True})
my_contourf_kwargs_default.update(contourf_minmax_kwargs_default)

def my_contourf(xx, yy, field, fig, ax, **kwargs):
    kw = update_dict(my_contourf_kwargs_default, kwargs)
    kw_contourf_minmax = update_dict(contourf_minmax_kwargs_default, kwargs)
    find_bad_keys(my_contourf_kwargs_default, kwargs, 'my_contourf')

    # make sure Python does not modify any of the arrays it was passed
    field = np.copy(field)

    # get default bounds if not specified
    if kw.minmax is None:
        kw.minmax = contourf_minmax(field, **kw_contourf_minmax)
    elif kw.minmax == 'fullrange':
        kw_contourf_minmax.fullrange = True
        kw.minmax = contourf_minmax(field, **kw_contourf_minmax)

    # plot the field, maybe
    # Factor out the exponent on the field and put it on the color bar
    # can turn this behavior off with "nosci=True"
    if not kw.nosci:
        maxabs = max(np.abs(kw.minmax[0]), np.abs(kw.minmax[1]))
        exp = get_exp(maxabs)
        divisor = 10.0**exp
        field /= divisor
        kw.minmax = kw.minmax[0]/divisor, kw.minmax[1]/divisor

    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, kw.minmax[0], kw.minmax[1])

    nlevelsfield = 160
    levels = np.linspace(kw.minmax[0], kw.minmax[1], nlevelsfield + 1)
    if kw.cmap is None:
        if kw.posdef:
            kw.cmap = 'plasma'
        else:
            kw.cmap = 'RdYlBu_r'

    # finally we make the contour plot!
    if kw.plotfield:
        im = ax.contourf(xx, yy, field, cmap=kw.cmap, levels=levels)

    # now deal with color bar, if one is desired
    if kw.plotfield and kw.plotcbar:
        # get fig dimensions
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        fig_aspect = fig_height_inches/fig_width_inches
        # get ax dimensions
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
        ax_width = ax_right - ax_left
        ax_height = ax_top - ax_bottom
        if kw.cbar_pos == 'bottom':
            orientation = 'horizontal'
            cbar_height = kw.cbar_thick/fig_height_inches
            cbar_width = cbar_height/kw.cbar_aspect*fig_aspect
            if cbar_width > ax_width: # don't let cbar be thicker than plot!
                cbar_width = ax_width

            # centrally position colorbar underneath the axes
            label_buff = 3/8/fig_height_inches # needs to contain
            # the colorbar ticklabels and little buffer space
            lilbit = 1/16/fig_height_inches
            cbar_left = ax_left + 0.5*ax_width - 0.5*cbar_width
            cbar_bottom = ax_bottom - lilbit - cbar_height -\
                    (kw.cbar_no - 1)*(label_buff + cbar_height)
        elif kw.cbar_pos == 'right':
            orientation = 'vertical'
            cbar_width = kw.cbar_thick/fig_width_inches
            cbar_height = cbar_width/kw.cbar_aspect/fig_aspect
            if cbar_height > ax_height: 
                cbar_height = ax_height

            # centrally position colorbar to right of axes
            label_buff = 3/4/fig_width_inches # needs to contain
            # the colorbar ticklabels and little buffer space
            lilbit = 1/16/fig_width_inches
            cbar_bottom = ax_bottom + 0.5*ax_height - 0.5*cbar_height
            cbar_left = ax_right + lilbit + (kw.cbar_no - 1)*(label_buff + cbar_width)
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width,\
                cbar_height))        
        cbar = plt.colorbar(im, cax=cax, orientation=orientation)
            
        # font size for the tick labels
        cax.tick_params(labelsize=kw.fontsize)
        #cbar.ax.tick_params(labelsize=fontsize)   

        if kw.nosci:
            cbar_label = kw.units
        else:
            cbar_label = (r'$\times10^{%i}\ $' %exp) + kw.units
        # ticklabel format
        fmt = '%.' + str(kw.cbar_prec) + 'f'
        if kw.posdef:
            tickvals = [kw.minmax[0], kw.minmax[1]]
        else:
            tickvals = [kw.minmax[0], 0., kw.minmax[1]]
        ticklabels = []
        for tickval in tickvals:
            ticklabels.append(fmt %tickval)
        cbar.set_ticks(tickvals)
        cbar.set_ticklabels(ticklabels)
        if kw.cbar_pos == 'bottom':
            fig.text(cbar_left + cbar_width + lilbit*fig_aspect,\
                    cbar_bottom + 0.5*cbar_height, cbar_label,\
                    ha='left', va='center', fontsize=kw.fontsize) 
        elif kw.cbar_pos == 'right':
            #fig.text(cbar_left + cbar_width + lilbit/fig_aspect,\
            #        cbar_bottom + 0.5*cbar_height, cbar_label,\
            #        ha='left', va='center', fontsize=kw.fontsize) 
            cax.set_title(cbar_label, ha='left', fontsize=kw.fontsize)

    # Plot contours if desired
    linestyle = '--'
    if kw.plotcontours:
        # Determine the contour levels
        if kw.contourlevels is None:
            # just thin out the field levels
            nskip = nlevelsfield//kw.ncontours
            kw.contourlevels = levels[::nskip]
                    # contourf whitespace, I think. Not sure why)

        # Determine how to color the contours
        if kw.posdef:
            contourcolor = 'w'
        else:
            contourcolor = 'k'
        contourlw = kw.lw

        # plot the contours
        ax.contour(xx, yy, field, kw.contourlevels,\
                colors=contourcolor, linewidths=contourlw, linestyles=linestyle)

    # finally, plot some lines!

    # need to check if user provided the appropriate functions
    # if plotboundary == True
    if kw.plotboundary:
        if kw.func1 is None or kw.func2 is None:
            print ("my_contourf(): plotboundary = True, but either ")
            print ("func1 or func2 was not provided. Setting plotboundary=False")
            kw.plotboundary = False

    for ind in [1, 2]:
        if ind == 1:
            vals = list(kw.vals1)
            func = kw.func1
        if ind == 2:
            vals = list(kw.vals2)
            func = kw.func2

        linewidths = [kw.lw]*len(vals)
        linestyles = [linestyle]*len(vals)
        if kw.plotboundary:
            vals = [np.min(func)] + vals + [np.max(func)]
            vals = np.array(vals)
            linewidths = [kw.lw] + linewidths + [kw.lw] # make boundary
            # lines a a bit thicker... maybe
            linestyles = ['-'] + linestyles + ['-'] # make boundary line 
            # solid

        if len(vals) > 0:
            # check to make sure values are within (strictly) 
            # the func's range (otherwise they can't be plotted)
            vmin, vmax = np.min(func), np.max(func)
            maxabs = np.max(np.abs(func))
            for i in range(len(vals)):
                val = vals[i]
                if val <= vmin:
                    vals[i] = vmin + maxabs*1.0e-15
                if val >= vmax:
                    vals[i] = vmax - maxabs*1.0e-15 
                    # need to do this (same issue with
                    # contourf whitespace, I think. Not sure why)
            ax.contour(xx, yy, func, vals, colors='k', linestyles=linestyles, linewidths=linewidths)

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
