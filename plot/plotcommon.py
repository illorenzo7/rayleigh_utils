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
    subplot_width_inches = width
    subplot_height_inches = width
    margin_inches = 1./8.
    margin_bottom_inches = 3./4.

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/fig_width_inches
    margin_bottom = margin_bottom_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    return fig, ax

def default_axes_2by1(width=8.):
    # Good for Mollweide projections
    subplot_width_inches = width
    subplot_height_inches = 0.5*width
    margin_inches = 1./8.
    margin_bottom_inches = 3./4. # leave room for colorbar

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + margin_inches +\
            margin_bottom_inches

    margin_x = margin_inches/fig_width_inches
    margin_bottom = margin_bottom_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    return fig, ax

def default_axes_1by2(width=3.75):
    # Good for AZ_Avgs in full meridional plane
    subplot_width_inches = width
    subplot_height_inches = 2*width
    margin_inches = 1./8.

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + 2.*margin_inches

    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))
    return fig, ax
