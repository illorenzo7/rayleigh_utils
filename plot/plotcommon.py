# Module for routines common to many plotting scripts
# Created: 02/08/2019
import numpy as np
from matplotlib import colors, ticker
import matplotlib.pyplot as plt
from common import *

color_order = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
style_order = ['-', '--', '-.', ':']
marker_order = [".", "o", "v","s", "*", "x", "^", "<", ">"]
default_lw = 1.0
default_s = 0.2 # markersize
default_labelsize = 12
default_titlesize = 12
default_ticksize = 12
default_margin = 1/16
default_margin_label = 3/8
default_margin_xlabel = 1/2
default_margin_ylabel = 3/4
# ylabels take up more space because floating
# point numbers are longer than they are tall
default_margin_title = 3.0/4.0

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

def make_figure(nplots=None, sub_width_inches=None, sub_aspect=None, sub_height_inches=None, nrow=None, ncol=None, margin_left_inches=None, margin_right_inches=None, margin_bottom_inches=None, margin_top_inches=None, sub_margin_left_inches=None, sub_margin_right_inches=None, sub_margin_bottom_inches=None, sub_margin_top_inches=None, width_inches=None, height_inches=None, aspect=None):

    # first, if any margin is unspecified, then it equals the default
    # "default_margin"
    # subplot margins are above and beyond figure margins
    if margin_left_inches is None:
        margin_left_inches = default_margin
    if margin_right_inches is None:
        margin_right_inches = default_margin
    if margin_bottom_inches is None:
        margin_bottom_inches = default_margin
    if margin_top_inches is None:
        margin_top_inches = default_margin_title
    if sub_margin_left_inches is None:
        sub_margin_left_inches = default_margin
    if sub_margin_right_inches is None:
        sub_margin_right_inches = default_margin
    if sub_margin_bottom_inches is None:
        sub_margin_bottom_inches = default_margin
    if sub_margin_top_inches is None:
        sub_margin_top_inches = default_margin_label

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

    # figure out aspect, width, height
    # need to specify TWO of these (other gets determined) for EITHER
    # the figure or the subplots; check subplots first (defaults set here)
    sub_width_inches_default, sub_height_inches_default = 3.5, 2.5
    width_inches_default, height_inches_default = 7.0, 4.0

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
    elif sub_nspec == 2:
        which_spec = 'sub'
        if sub_width_inches is None:
            sub_width_inches = sub_height_inches/sub_aspect
        if sub_height_inches is None:
            sub_height_inches = sub_width_inches*sub_aspect
        if sub_aspect is None:
            sub_aspect = sub_height_inches/sub_width_inches

    elif nspec == 1:
        which_spec = 'fig'
        if not width_inches is None:
            height_inches = height_inches_default
            aspect = height_inches/width_inches
        elif not height_inches is None:
            width_inches = width_inches_default
            aspect = height_inches/width_inches
        elif not aspect is None:
            width_inches = width_inches_default
            height_inches = aspect*width_inches
    elif nspec == 2:
        which_spec = 'fig'
        if width_inches is None:
            width_inches = height_inches/aspect
        if height_inches is None:
            height_inches = width_inches*aspect
        if aspect is None:
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

def lineplot(xx, yy, ax, xlabel=None, ylabel=None, title=None, xvals=np.array([]), yvals=np.array([]), label=None, xlog=False, logscale=False, xminmax=None, minmax=None, minmax2=None, scatter=False, xcut=None, color=color_order[0], linestyle=style_order[0], marker=marker_order[0], lw=default_lw, s=default_s, showplot=False, domain_bounds=None):

    axs = [ax]

    # convert xvals, yvals to lists
    xvals = make_array(xvals, tolist=True)
    yvals = make_array(yvals, tolist=True)

    if not xcut is None:
        symmetrize = True
        ax2 = ax.twinx()
        axs.append(ax2)
        ixcut = np.argmin(np.abs(xx - xcut))
        if xx[0] < xx[-1]: # x axis goes in "correct" direction
            ax_left = ax
            ax_right = ax2
        else: # x axis is reversed, so the second one is on the left
            ax_left = ax2
            ax_right = ax
    else:
        symmetrize = False
        ixcut = len(xx)
        ax_left = ax
   
    if xcut is None:
        if scatter:
            ax.scatter(xx, yy, label=label, marker=marker, s=s, color=color)
        else:
            ax.plot(xx, yy, label=label, linestyle=linestyle, linewidth=lw, color=color)
    else:
        if scatter:
            ax.scatter(xx[:ixcut], yy[:ixcut], marker=marker, s=s, color=color)
            ax2.scatter(xx[ixcut:], yy[ixcut:], marker=marker, s=s, color=color)
        else:
            ax.plot(xx[:ixcut], yy[:ixcut], linestyle=linestyle, linewidth=lw, color=color)
            ax2.plot(xx[ixcut:], yy[ixcut:], linestyle=linestyle, linewidth=lw, color=color)

    # ticks (mostly everywhere, deal with split axes)
    if xcut is None:
        plt.sca(ax)
        plt.minorticks_on()
        plt.tick_params(top=True, right=True, direction='in',\
                which='both')
        if not ylabel is None:
            plt.ylabel(ylabel, fontsize=default_labelsize)
    else:
        plt.sca(ax_right)
        plt.minorticks_on()
        plt.tick_params(top=True, left=False, right=True, direction='in', which='both')
        ax_right.yaxis.tick_right()
        plt.sca(ax_left)
        plt.minorticks_on()
        plt.tick_params(top=True, left=True, right=False, direction='in',\
                which='both')
        if not ylabel is None:
            ax_left.set_ylabel(ylabel, fontsize=default_labelsize)
            ax_left.yaxis.set_label_position('left')
        ax_left.yaxis.tick_left()

    if xlog:
        ax.set_xscale('log')
    if logscale:
        for ax in axs:
            ax.set_yscale('log')

    if xminmax is None:
        xminmax = np.min(xx), np.max(xx)
    ax.set_xlim(xminmax)

    if minmax is None:
        minmax = lineplot_minmax([yy[:ixcut]], logscale=logscale, symmetrize=symmetrize, domain_bounds=domain_bounds, xx=xx[:ixcut])
    ax.set_ylim(minmax)

    if not xcut is None:
        if minmax2 is None:
            minmax2 = lineplot_minmax([yy[ixcut:]], logscale=logscale, symmetrize=symmetrize, domain_bounds=domain_bounds, xx=xx[ixcut:])
        ax2.set_ylim(minmax2)

    if not xcut is None:
        xvals.append(xx[ixcut])
        yvals.append(0)
 
    npoints = 100
    xpoints = np.linspace(xminmax[0], xminmax[1], npoints)

    # possibly mark x/y - values
    for xval in xvals:
        if xcut is None:
            ax_loc = ax
        else:
            if xval <= xcut:
                ax_loc = ax_left
            else:
                ax_loc = ax_right
        y1, y2 = ax_loc.get_ylim()
        ypoints = np.linspace(y1, y2, npoints)
        ax_loc.plot(xval + np.zeros(npoints), ypoints, 'k--', linewidth=lw)
    for yval in yvals:
        for ax_loc in axs:
            # only plot if the line is within the range
            y1, y2 = ax_loc.get_ylim()
            if y1 < yval < y2:
                ax_loc.plot(xpoints, yval + np.zeros(npoints), 'k--', linewidth=lw)

    if not xlabel is None:
        ax.set_xlabel(xlabel, fontsize=default_labelsize)

    if not title is None:
        ax.set_title(title, fontsize=default_titlesize)

    # set the tick label size
    for ax_loc in axs:
        plt.sca(ax_loc)
        plt.xticks(fontsize=default_ticksize)
        plt.yticks(fontsize=default_ticksize)

        # Get the y-axis in scientific notation
        plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

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
            fig.text(cbar_left + cbar_width + lilbit/fig_aspect,\
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
