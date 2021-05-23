# Module for routines common to many plotting scripts
# Created: 02/08/2019
import numpy as np
from matplotlib import colors, ticker
import matplotlib.pyplot as plt
from common import *
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
#from matplotlib import ticker
#sys.path.append(os.environ['raco'])

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
default_margin_xlabel = 3/8
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

def my_set_title(fig, ax, title, fontsize=default_titlesize, ha='center', va='bottom', pad=1./8.):
    ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
    ax_width = ax_right - ax_left
    ax_height = ax_top - ax_bottom
    fig_width_inches, fig_height_inches = fig.get_size_inches()
    margin_y = pad/fig_height_inches
    fig.text(ax_left + 0.5*ax_width, ax_top + margin_y, title, fontsize=fontsize, ha=ha, va=va, **csfont)

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

def make_cbar(fig, ax, im, aspect=1./20., length=0.75, orientation='horizontal', units='', fs=10., posdef=False, logscale=False, symlog=False, linthresh=None, linscale=None):

    # get figure and axis dimensions
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    width_inches, height_inches = fig.get_size_inches()
    aspect = height_inches/width_inches

    # make colorbar axes and colorbar itself
    default_margin = 1./8.
    margin_x = default_margin/width_inches
    margin_y = default_margin/height_inches
    if orientation == 'horizontal':
        cbar_width = length*ax_delta_x
        cbar_height = cbar_width*aspect/aspect
        cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
        cbar_bottom = ax_ymin - margin_y - cbar_height
    if orientation == 'vertical':
        cbar_height = length*ax_delta_y
        cbar_width = cbar_height*aspect*aspect
        cbar_left = ax_xmax + margin_x
        cbar_bottom = ax_ymin + 0.5*ax_delta_y - 0.5*cbar_height
    cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))   
    cbar = plt.colorbar(im, cax=cax, orientation=orientation)
        
    # font size for the ticks
    cax.tick_params(labelsize=fs)
    cbar.ax.tick_params(labelsize=fs)   

    # location of the ticks
    if logscale:
        locator = ticker.LogLocator(subs='all')
        cbar.set_ticks(locator)
        cbar_label = units
    elif posdef:
        cbar_label = (r'$\times10^{%i}\ $' %exp) + units
        cbar.set_ticks([minmax[0], minmax[1]])
        cbar.set_ticklabels(['%1.1f' %minmax[0],\
                '%1.1f' %minmax[1]])
    elif symlog:
        cbar_label = units
        nlin = 5
        nlog = 6
        lin_ticks = np.linspace(-linthresh, linthresh, nlin)
        log_ticks1 = np.linspace(minmax[0], -linthresh, nlog,\
                endpoint=False)
        log_ticks2 = -log_ticks1[::-1]
        ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
        nticks = nlin + 2*nlog
        cbar.set_ticks(ticks)
        ticklabels = []
        for i in range(nticks):
            ticklabels.append(r'')
#                ticklabels[0] = sci_format(minmax[0])
        ticklabels[nlog] = sci_format(-linthresh)
#                ticklabels[nticks//2] = r'$0$'
#                ticklabels[nlog + nlin - 1] = sci_format(linthresh)
        ticklabels[nticks - 1] = sci_format(minmax[1])
        cbar.set_ticklabels(ticklabels)
    else:
        if nosci:
            cbar_label = units
        else:
            cbar_label = (r'$\times10^{%i}\ $' %exp) + units
        cbar.set_ticks([minmax[0], 0, minmax[1]])
        cbar.set_ticklabels(['%.1f' %minmax[0], '0', '%.1f'\
                %minmax[1]])

    # Title the colorbar based on the field's units
    line_height = 1./4./height_inches
    fig.text(cbar_left + 0.5*cbar_width, cbar_bottom - line_height,\
             cbar_label, ha='center', va='top', **csfont,\
             fontsize=cbar_fs) 
    #fig.text(cbax_left - 0.3*cbax_width, cbax_center_y,\
    #        cbar_label, ha='right', va='center', rotation=90,\
    #        fontsize=cbar_fs)

    if not rbcz is None: # Make a colorbar for the RZ
        cbar_bottom = ax_ymin - 2.5*cbar_height - 3*line_height
        cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width,\
                cbar_height))        
        cbar = plt.colorbar(imrz, cax=cax, orientation='horizontal')
            
        cax.tick_params(labelsize=cbar_fs)
        cbar.ax.tick_params(labelsize=cbar_fs)   
        # font size for the ticks

        if logscale:
            locator = ticker.LogLocator(subs='all')
            cbar.set_ticks(locator)
            cbar_label = units
        elif posdef:
            cbar_label = (r'$\times10^{%i}\ $' %exprz) + units
            cbar.set_ticks([minmaxrz[0], minmaxrz[1]])
            cbar.set_ticklabels(['%1.1f' %minmaxrz[0],\
                    '%1.1f' %minmaxrz[1]])
        elif symlog:
            cbar_label = units
            nlin = 5
            nlog = 6
            lin_ticks = np.linspace(-linthreshrz, linthreshrz, nlin)
            log_ticks1 = np.linspace(minmaxrz[0], -linthreshrz,\
                    nlog, endpoint=False)
            log_ticks2 = -log_ticks1[::-1]
            ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
            nticks = nlin + 2*nlog
            cbar.set_ticks(ticks)
            ticklabels = []
            for i in range(nticks):
                ticklabels.append(r'')
#                ticklabels[0] = sci_format(minmax[0])
            ticklabels[nlog] = sci_format(-linthreshrz)
#                ticklabels[nticks//2] = r'$0$'
#                ticklabels[nlog + nlin - 1] = sci_format(linthresh)
            ticklabels[nticks - 1] = sci_format(minmaxrz[1])
            cbar.set_ticklabels(ticklabels)
        else:
            cbar_label = (r'$\times10^{%i}\ $' %exprz) + units
            cbar.set_ticks([minmaxrz[0], 0, minmaxrz[1]])
            cbar.set_ticklabels(['%1.1f' %minmaxrz[0], '0', '%1.1f'\
                    %minmaxrz[1]])

        # Title the colorbar based on the field's units
        line_height = 1./4./height_inches
        fig.text(cbar_left + 0.5*cbar_width, cbar_bottom -\
                line_height, cbar_label, ha='center', va='top',\
                **csfont, fontsize=cbar_fs) 

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

def lineplot(xx, yy, ax=None, axtwin=None, xlabel=None, ylabel=None, title=None, xvals=None, yvals=None, label=None, xlog=False, ylog=False, xminmax=None, yminmax=None, yminmax2=None, scatter=False, xcut=None, color=color_order[0], linestyle=style_order[0], marker=marker_order[0], lw=default_lw, s=default_s, showplot=False):

    if ax is None: # probably called from command line
        fig, axs, fpar = make_figure()
        ax = axs[0,0]
        showplot = True

    axs = [ax]

    if not xcut is None:
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
            plt.ylabel(ylabel, **csfont, fontsize=default_labelsize)
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
            ax_left.set_ylabel(ylabel, **csfont, fontsize=default_labelsize)
            ax_left.yaxis.set_label_position('left')
        ax_left.yaxis.tick_left()

    if xlog:
        ax.set_xscale('log')
    if ylog:
        for ax in axs:
            ax.set_yscale('log')

    if xminmax is None:
        xminmax = ax.get_xlim()
    else:
        ax.set_xlim(xminmax)

    if xcut is None:
        if yminmax is None:
            yminmax = ax.get_ylim()
    else:
        if yminmax is None:
            yminmax = ax.get_ylim()
            maxabs = max(np.abs(yminmax[0]), np.abs(yminmax[1]))
            yminmax = -maxabs, maxabs
        if yminmax2 is None:
            yminmax2 = ax2.get_ylim()
            maxabs = max(np.abs(yminmax2[0]), np.abs(yminmax2[1]))
            yminmax2 = -maxabs, maxabs
    
    ax.set_ylim(yminmax)
    if not xcut is None:
        ax2.set_ylim(yminmax2)

    # possibly append to xvals and yvals (0 point and xcut)
    xvals_to_add = []
    if xminmax[0] < 0.0 < xminmax[1]:
        xvals_to_add.append(0.0)
    if not xcut is None:
        xvals_to_add.append(xx[ixcut])
    xvals_to_add = np.array(xvals_to_add)
    if xvals is None:
        xvals = xvals_to_add
    else:
        xvals = np.array(xvals)
        xvals = np.hstack((xvals, xvals_to_add))
    if yvals is None:
        yvals = np.array([0.0])
    else:
        yvals = np.hstack((yvals, np.array([0.0])))
 
    npoints = 100
    xpoints = np.linspace(xminmax[0], xminmax[1], npoints)

    # possibly mark x/y - values
    if not xvals is None:
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
    if not yvals is None:
        for yval in yvals:
            for ax_loc in axs:
                # only plot if the line is within the range
                y1, y2 = ax_loc.get_ylim()
                if y1 < yval < y2:
                    ax_loc.plot(xpoints, yval + np.zeros(npoints), 'k--', linewidth=lw)

    if not xlabel is None:
        ax.set_xlabel(xlabel, **csfont, fontsize=default_labelsize)


    if not title is None:
        ax.set_title(title, **csfont, fontsize=default_titlesize)

    # set the tick label size
    for ax_loc in axs:
        plt.sca(ax_loc)
        plt.xticks(fontsize=default_ticksize)
        plt.yticks(fontsize=default_ticksize)

        # Get the y-axis in scientific notation
        plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

    if showplot:
        plt.show()
        return fig, ax

def my_contourf(xx, yy, field, fig=None, ax=None,\
        # saturation of field values stuff
        minmax=None, fullrange=False, posdef=False, logscale=False, symlog=False, linthresh=None, linscale=None,\
        # basic flags:
        showplot=False, plotfield=True, nlevelsfield=160, levels=None,\
        plotcontours=True, ncontours=8, contourlevels=None, contourlw=1.0, contourcolor=None,\
        # colorbar stuff
        plotcbar=True, cbar_thick=1./8., cbar_aspect=1/15, cbar_prec=2, cbar_no=1, cmap=None, units='', nosci=False, fontsize=default_labelsize,\
        # coordinate line stuff; do up to two "types"
        vals1=np.array([]), func1=None, linestyles1=[], colors1=[], lw1=1.0,\
        vals2=np.array([]), func2=None, linestyles2=[], colors2=[], lw2=1.0,\
        plotboundary=True):

    # Create a default set of figure axes if they weren't specified
    if fig is None or ax is None:
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line
        
        # get default axes with correct axis ratio
        delta_x = np.max(xx) - np.min(xx)
        delta_y = np.max(yy) - np.min(yy)
        nplots = 1
        sub_width_inches = 4.0
        sub_aspect = delta_y/delta_x
        sub_margin_bottom_inches = 1.0
        fig, axs, fpar = make_figure(nplots, sub_width_inches, sub_aspect, sub_margin_bottom_inches=sub_margin_bottom_inches)
        ax = axs[0, 0]

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed
    field = np.copy(field)

    # If using logscale, you better have positive values!
    if logscale:
        posdef = True

    # Get default bounds if not specified
    if minmax is None:
        minmax = get_satvals(field, posdef=posdef, logscale=logscale, symlog=symlog, fullrange=fullrange)

    # plot the field, maybe
    # By default (for the linear-scaled default and posdef),
    # Factor out the exponent on the field and put it on the color bar
    # can turn this behavior off with "nosci=True"
    if not (logscale or symlog or nosci):
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = get_exp(maxabs)
        divisor = 10.0**exp
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor

    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, minmax[0], minmax[1])

    norm = None
    if logscale:
        if cmap is None:
            cmap = 'Greys'
        log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
        if levels is None:
            levels = np.logspace(log_min, log_max, nlevelsfield + 1)
        norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1])
    elif posdef:
        if cmap is None:
            cmap = 'plasma'
        if levels is None:
            levels = np.linspace(minmax[0], minmax[1], nlevelsfield + 1)
    elif symlog:
        if cmap is None:
            cmap = 'RdYlBu_r'
        linthresh_default, linscale_default =\
            get_symlog_params(field, field_max=minmax[1])
        if linthresh is None:
            linthresh = linthresh_default
        if linscale is None:
            linscale = linscale_default
        if levels is None:
            log_thresh = np.log10(linthresh)
            log_max = np.log10(minmax[1])
            nlevels_per_interval = nlevelsfield//4
            levels_neg = -np.logspace(log_max, log_thresh,\
                    nlevels_per_interval,\
                    endpoint=False)
            levels_mid = np.linspace(-linthresh, linthresh,\
                    2*nlevels_per_interval, endpoint=False)
            levels_pos = np.logspace(log_thresh, log_max,\
                    nlevels_per_interval + 1)
            levels = np.hstack((levels_neg, levels_mid, levels_pos))
        norm = colors.SymLogNorm(linthresh=linthresh,\
            linscale=linscale, vmin=minmax[0], vmax=minmax[1], base=10)
    else: # this is the default...just do linear scaled pos/neg values
        if cmap is None:
            cmap = 'RdYlBu_r'
        if levels is None:
            levels = np.linspace(minmax[0], minmax[1], nlevelsfield + 1)

    # finally we make the contour plot!
    if plotfield:
        im = ax.contourf(xx, yy, field, cmap=cmap, levels=levels, norm=norm)

    # now deal with color bar
    if plotfield and plotcbar:
        # get cbar dimensions
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        fig_aspect = fig_height_inches/fig_width_inches
        cbar_height = cbar_thick/fig_height_inches
        cbar_width = cbar_height/cbar_aspect*fig_aspect

        # centrally position colorbar underneath the axes
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
        ax_width = ax_right - ax_left
        ax_height = ax_top - ax_bottom
        line_height = 3/8/fig_height_inches # needs to contain
        # the colorbar ticklabels and little buffer space
        lilbit = 1.0/16.0/fig_height_inches
        cbar_left = ax_left + 0.5*ax_width - 0.5*cbar_width
        cbar_bottom = ax_bottom - lilbit - cbar_height -\
                (cbar_no - 1)*(line_height + cbar_height)
        cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width,\
                cbar_height))        
        cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            
        # font size for the tick labels
        cax.tick_params(labelsize=fontsize)
        #cbar.ax.tick_params(labelsize=fontsize)   

        if logscale:
            locator = ticker.LogLocator(subs='all')
            cbar.set_ticks(locator)
            cbar_label = units
        elif posdef:
            cbar_label = (r'$\times10^{%i}\ $' %exp) + units
            cbar.set_ticks([minmax[0], minmax[1]])
            fmt = '%.' + str(cbar_prec) + 'f'
            cbar.set_ticklabels([fmt %minmax[0],\
                    fmt %minmax[1]])
        elif symlog:
            cbar_label = units
            nlin = 5
            nlog = 6
            lin_ticks = np.linspace(-linthresh, linthresh, nlin)
            log_ticks1 = np.linspace(minmax[0], -linthresh, nlog,\
                    endpoint=False)
            log_ticks2 = -log_ticks1[::-1]
            ticks = np.hstack((log_ticks1, lin_ticks, log_ticks2))
            nticks = nlin + 2*nlog
            cbar.set_ticks(ticks)
            ticklabels = []
            for i in range(nticks):
                ticklabels.append(r'')
#                ticklabels[0] = sci_format(minmax[0])
            ticklabels[nlog] = sci_format(-linthresh, ndec=cbar_prec)
#                ticklabels[nticks//2] = r'$0$'
#                ticklabels[nlog + nlin - 1] = sci_format(linthresh)
            ticklabels[nticks - 1] = sci_format(minmax[1], ndec=cbar_prec)
            cbar.set_ticklabels(ticklabels)
        else:
            if nosci:
                cbar_label = units
            else:
                cbar_label = (r'$\times10^{%i}\ $' %exp) + units
            cbar.set_ticks([minmax[0], 0, minmax[1]])
            fmt = '%.' + str(cbar_prec) + 'f'
            cbar.set_ticklabels([fmt %minmax[0], '0', fmt\
                    %minmax[1]])
            #cbar.set_ticklabels(['%.1f' %minmax[0], '0', '%.1f'\
            #        %minmax[1]])

        # Title the colorbar based on the field's units
        fig.text(cbar_left + cbar_width + lilbit/fig_aspect,\
                cbar_bottom + 0.5*cbar_height, cbar_label,\
                ha='left', va='center', **csfont, fontsize=fontsize) 

    # Plot contours if desired
    if plotcontours:
        # Determine the contour levels
        if contourlevels is None:
            # just thin out the field levels
            nskip = nlevelsfield//ncontours
            contourlevels = levels[::nskip]
                    # contourf whitespace, I think. Not sure why)

        # Determine how to color the contours
        if contourcolor is None:
            if logscale:
                contourcolor = 'r'
            elif posdef:
                contourcolor = 'w'
            else:
                contourcolor = 'k'

        # plot the contours
        ax.contour(xx, yy, field, contourlevels,\
                colors=contourcolor, linewidths=contourlw)

    # finally, plot some lines!

    # need to check if user provided the appropriate functions
    # if plotboundary == True
    if plotboundary:
        if func1 is None or func2 is None:
            print ("my_contourf(): plotboundary = True, but either ")
            print ("func1 or func2 was not provided. Setting plotboundary=False")
            plotboundary = False

    for ind in [1, 2]:
        if ind == 1:
            vals = vals1
            func = func1
            linestyles = linestyles1
            colors_loc = colors1 # need to be careful with this one!
            lw = lw1
        if ind == 2:
            vals = vals2
            func = func2
            linestyles = linestyles2
            colors_loc = colors2
            lw = lw2

        vals = list(vals)
        linestyles = list(linestyles)
        colors_loc = list(colors_loc)
        if linestyles == []:
            linestyles = ['--']*len(vals)
        if colors_loc == []:
            colors_loc = ['k']*len(vals)
        linewidths = [lw]*len(vals)

        if plotboundary:
            vals = [np.min(func)] + list(vals) + [np.max(func)]
            vals = np.array(vals)
            linestyles = ['-'] + linestyles + ['-']
            colors_loc = ['k'] + colors_loc + ['k']
            linewidths = [2*lw] + linewidths + [2*lw]

        if len(vals) > 0:
            # check to make sure values are within (strictly) the func's range
            # (otherwise they can't be plotted)
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
            ax.contour(xx, yy, func, vals,\
                    linestyles=linestyles, colors=colors_loc, linewidths=linewidths)

    # Set ax ranges to be just outside the boundary lines
    # avoid weird whitespace cutoffs
    lilbit = 0.01
    xmin, xmax = np.min(xx), np.max(xx)
    ymin, ymax = np.min(yy), np.max(yy)
    Dx = xmax - xmin
    Dy = ymax - ymin
    ax.set_xlim((xmin - lilbit*Dx, xmax + lilbit*Dx))
    ax.set_ylim((ymin - lilbit*Dy, ymax + lilbit*Dy))
    ax.axis('off') 

    return fig, ax
