# Module for routines common to many plotting scripts
# Created: 02/08/2019
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import ticker

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

def integerticks(lmax):
    ''' This is for plotting l time-traces
    Gives O(5) ticklabels between -lmax and lmax, picked sensibly
    '''
    delta_l = lmax/3
    # round this to nearest even number
    delta_l = 2*round(delta_l/2)
    positive_ticks = np.arange(0, lmax, delta_l)
    negative_ticks = np.arange(-1, -lmax - 1, -delta_l)
    if delta_l > 3:
        negative_ticks = negative_ticks[1:]
    tickvals = np.hstack((positive_ticks, negative_ticks))
    st = []
    for tickval in tickvals:
        tickval = int(tickval)
        st.append(str(np.abs(tickval)))
    ticklabels = np.array(st)
    return tickvals, ticklabels

def default_axes_1by1(width=6.):
    # Good for orthographic projections and AZ_Avgs in HALF meridional plane 
    # and equatorial slices
    sub_width_inches = width
    sub_height_inches = width
    margin_inches = 1./8.
    margin_bottom_inches = 3./4.

    width_inches = sub_width_inches + 2.*margin_inches
    height_inches = sub_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/width_inches
    margin_bottom = margin_bottom_inches/height_inches
    sub_width = sub_width_inches/width_inches
    sub_height = sub_height_inches/height_inches

    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, sub_width,\
            sub_height))
    return fig, ax

def default_axes_2by1(width=8.):
    # Good for Mollweide projections
    sub_width_inches = width
    sub_height_inches = 0.5*width
    margin_inches = 1./8.
    margin_bottom_inches = 3./4. # leave room for colorbar

    width_inches = sub_width_inches + 2.*margin_inches
    height_inches = sub_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/width_inches
    margin_bottom = margin_bottom_inches/height_inches
    sub_width = sub_width_inches/width_inches
    sub_height = sub_height_inches/height_inches

    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, sub_width,\
            sub_height))
    return fig, ax

def default_axes_1by2(width=3.75):
    # Good for AZ_Avgs in full meridional plane
    sub_width_inches = width
    sub_height_inches = 2*width
    margin_inches = 1./8.

    width_inches = sub_width_inches + 2.*margin_inches
    height_inches = sub_height_inches + 2.*margin_inches

    margin_x = margin_inches/width_inches
    margin_y = margin_inches/height_inches
    sub_width = sub_width_inches/width_inches
    sub_height = sub_height_inches/height_inches

    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.add_axes((margin_x, margin_y, sub_width, sub_height))
    return fig, ax

def default_axes_tl(width=7.):
    # Good for time-latitude plots
    sub_width_inches = width
    sub_height_inches = 1.5
    margin_inches = 1./8.
    margin_left_inches = 1./2.
    margin_right_inches = 1.25
    margin_bottom_inches = 1./2.

    width_inches = sub_width_inches + margin_left_inches +\
            margin_right_inches
    height_inches = sub_height_inches + margin_inches +\
            margin_bottom_inches

    margin_bottom = margin_bottom_inches/height_inches
    margin_left = margin_left_inches/width_inches
    sub_width = sub_width_inches/width_inches
    sub_height = sub_height_inches/height_inches

    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.add_axes((margin_left, margin_bottom, sub_width,\
            sub_height))
    return fig, ax

def testtex(label):
    plt.plot(range(10))
    plt.title(label)
    plt.show()

def generate_cbar(fig, ax, im, aspect=1./20., length=0.75, orientation='horizontal', units='', fs=10., posdef=False, logscale=False, symlog=False, linthresh=None, linscale=None):

    # get figure and axis dimensions
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    width_inches, height_inches = fig.get_size_inches()
    aspect = height_inches/width_inches

    # make colorbar axes and colorbar itself
    margin_inches = 1./8.
    margin_x = margin_inches/width_inches
    margin_y = margin_inches/height_inches
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

def make_figure(nplots=None, sub_width_inches=None, sub_aspect=None, sub_height_inches=None, margin_inches=0.125, nrow=None, ncol=None, margin_left_inches=None, margin_right_inches=None, margin_bottom_inches=None, margin_top_inches=None, sub_margin_left_inches=None, sub_margin_right_inches=None, sub_margin_bottom_inches=None, sub_margin_top_inches=None, width_inches=None, height_inches=None, aspect=None):

    # first, if any margin is unspecified, then it equals the default
    # "margin_inches"
    # subplot margins are above and beyond figure margins
    if margin_left_inches is None:
        margin_left_inches = margin_inches
    if margin_right_inches is None:
        margin_right_inches = margin_inches
    if margin_bottom_inches is None:
        margin_bottom_inches = margin_inches
    if margin_top_inches is None:
        margin_top_inches = margin_inches
    if sub_margin_left_inches is None:
        sub_margin_left_inches = margin_inches
    if sub_margin_right_inches is None:
        sub_margin_right_inches = margin_inches
    if sub_margin_bottom_inches is None:
        sub_margin_bottom_inches = margin_inches
    if sub_margin_top_inches is None:
        sub_margin_top_inches = margin_inches

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
