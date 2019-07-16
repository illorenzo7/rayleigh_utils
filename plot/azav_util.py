##########################################################################################
# This file contains utilities useful for plotting the data from azimuthal averages,
# files in the directory AZ_Avgs.
# Written by Nick Featherstone
# First modified by Loren Matilsky, 03/09/2018

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
plt.rcParams['contour.negative_linestyle'] = 'solid'
from common import rms, get_satvals, get_exp, rsun, trim_field
from binormalized_cbar import MidpointNormalize

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def default_axes_2by1():
    # Create plot
    subplot_width_inches = 3.75
    subplot_height_inches = 7.5
    margin_inches = 1./8.

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + 2.*margin_inches
    fig_aspect = fig_height_inches/fig_width_inches

    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))
    return fig, ax

def default_axes_1by1():
    # Create plot
    subplot_width_inches = 3.75
    subplot_height_inches = 3.75
    margin_inches = 1./8.

    fig_width_inches = subplot_width_inches + 2.*margin_inches
    fig_height_inches = subplot_height_inches + 2.*margin_inches
    fig_aspect = fig_height_inches/fig_width_inches

    margin_x = margin_inches/fig_width_inches
    margin_y = margin_inches/fig_height_inches
    subplot_width = subplot_width_inches/fig_width_inches
    subplot_height = subplot_height_inches/fig_height_inches

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))
    return fig, ax

def plot_azav(field, rr, cost, sint, fig=None, ax=None, cmap='RdYlBu_r',\
    units='', minmax=None, posdef=False, logscale=False,\
    plotcontours=True, plotfield=True, nlevs=10, levels=None,\
	plotlatlines=False, rvals=None, rvals_norm=None, fsize=8,\
    showplot=False, plot_cbar=True):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x2
    and adds a plot of [field] in the meridional plane to the axes,
    with colorbar in the "cavity" of the meridional plane.'''

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed (shouldn't fucking have to do this)
    field = np.copy(field)
    rr = np.copy(rr)
    cost = np.copy(cost)
    sint = np.copy(sint)

    # Derivative grid info
    ri, ro = np.min(rr), np.max(rr)
    shell_depth = ro - ri
    rr_depth = (np.max(rr) - rr)/shell_depth
    nr = len(rr)
    nt = len(cost)

    # If using logscale, you better have positive values!
    if logscale:
        posdef = True
	
    # Get default bounds if not specified
    if minmax is None:
        # Cut away the data near the domain boundaries
        trimmed_field = trim_field(field, rr, cost)
        minmax = get_satvals(trimmed_field, posdef, logscale)

    # Need these if logscale is True; made need them for other stuff later
    minexp, maxexp = get_exp(minmax[0]), get_exp(minmax[1])

    # Get the exponent to use for scientific notation
    if not logscale:
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = float(np.floor(np.log10(maxabs)))
        divisor = 10.**exp
        
        # Normalize field by divisor
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor

    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_2by1()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot
   
    # Get the position of the axes on the figure
    pos = ax.get_position().get_points()
    ax_left, ax_bottom = pos[0]
    ax_right, ax_top = pos[1]
    ax_width = ax_right - ax_left
    ax_height = ax_top - ax_bottom
    ax_aspect = ax_height/ax_width
  
    if plot_cbar:
        # Set the colorbar ax to be in the "cavity" of the meridional plane
        # The colorbar height is set by making sure it "fits" in the cavity
        chi = np.min(rr)/np.max(rr)
        cavity_height = ax_height*chi
        cbax_center_x = ax_left + 0.3*ax_width
        cbax_center_y = ax_bottom + ax_height/2.
        cbax_aspect = 10.
        cbax_height = 0.5*cavity_height
        cbax_width = cbax_height/cbax_aspect/ax_aspect
        
        cbax_left = cbax_center_x - cbax_width/2.
        cbax_bottom = cbax_center_y - cbax_height/2.
    
    # Calculate the grid on which to plot
    rr2d = rr.reshape((1, nr))
    cost2d = cost.reshape((nt, 1))
    sint2d = sint.reshape((nt, 1))
    xx = rr2d*sint2d/ro
    zz = rr2d*cost2d/ro

    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    lw = 1.
    contour_lw = 0.2
    
    if (plotfield):
        plt.sca(ax)
        levs = np.linspace(minmax[0], minmax[1], 100)
#        im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
#        levels=levs)
#        plt.sca(ax)
        if logscale:
            plt.pcolormesh(xx, zz, field, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]))
        else:
            if posdef:
                cmap = 'plasma'
            else:
                cmap = 'RdYlBu_r'
            plt.pcolormesh(xx, zz, field, vmin=minmax[0], vmax=minmax[1],\
                    cmap=cmap)
        if plot_cbar:
            cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                           cbax_width, cbax_height])
            cbar = plt.colorbar(cax=cbaxes)

            cbaxes.tick_params(labelsize=fsize)
            cbar.ax.tick_params(labelsize=fsize)   #font size for the ticks

            if not logscale:
                if posdef:
                    mid = (minmax[0] + minmax[1])/2.
                    ticks = np.array([minmax[0], mid, minmax[1]])
                else:
                    ticks = np.array([minmax[0], 0., minmax[1]])
                ticklabels = []
                for i in range(len(ticks)):
                    ticklabels.append(str(round(ticks[i], 1)))
                ticks = np.array(ticks)
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticklabels)
                cbar_label = (r'$\times10^{%i}\ $' %exp) + units
            else:
                cbar_label = units
            # Put the units and exponent to left of colorbar
            fig.text(cbax_left - 0.3*cbax_width, cbax_center_y, cbar_label,\
                ha='right', va='center', rotation=90, fontsize=fsize)

    # Plot contours in the meridional plane, if desired
    if plotcontours:
        # Determine the contour levels
        if levels is None:
            if logscale:
                levels = np.logspace(minexp, maxexp, nlevs)
            else:
                levels = np.linspace(minmax[0], minmax[1], nlevs)
        else: # the caller specified specific contour levels to plot!
            levels = np.array(levels)
        # Determine how to color the contours
        if logscale:
            contour_color = 'r'
        elif posdef:
            contour_color = 'w'
        else:
            contour_color = 'k'
        # plot the contours
        plt.sca(ax)
        plt.contour(xx, zz, field, colors=contour_color, levels=levels,\
                linewidths=contour_lw)

    # Plot the boundary of the meridional plane
    plt.sca(ax)
    plt.plot(rr[0]/ro*sint, rr[0]/ro*cost, 'k', linewidth=lw)
    plt.plot(rr[-1]/ro*sint, rr[-1]/ro*cost, 'k', linewidth=lw)
    plt.plot([0.,0.], [rr[-1]/ro, rr[0]/ro], 'k', linewidth=lw)
    plt.plot([0.,0.], [-rr[-1]/ro, -rr[0]/ro], 'k', linewidth=lw)

    # Plot latitude lines, if desired
    if plotlatlines: 
        lats = (np.pi/2. - np.arccos(cost))*180./np.pi
        lats_to_plot = np.arange(-75., 90., 15.) 
        for lat in lats_to_plot:
            it = np.argmin(np.abs(lats - lat))
            xx_in, zz_in = rr[-1]/ro*sint[it], rr[-1]/ro*cost[it]
            xx_out, zz_out = rr[0]/ro*sint[it], rr[0]/ro*cost[it]
            plt.sca(ax)
            plt.plot([xx_in, xx_out], [zz_in, zz_out], 'k',\
                    linewidth=contour_lw)

    # Plot various radii, if desired
    # rvals must be given normalized to outer rr
    if not rvals is None:
        if rvals_norm is None:
            rvals_norm = rsun
        for rval in rvals: 
            plt.sca(ax)
            rval_dim = rval*rvals_norm/ro
            # "dimensional" rval (in dimensions of ro!)
            plt.plot(rval_dim*sint, rval_dim*cost, 'k--', linewidth=.7)

    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

    if showplot:
        plt.show()

def plot_azav_half(field, rr, cost, sint, sym='even',\
        fig=None, ax=None, cmap='RdYlBu_r', units='', minmax=None,\
        posdef=False, logscale=False, plotcontours=True, plotfield=True,\
        nlevs=10, levels=None, plotlatlines=False, rvals=None, norm=None,\
        fsize=8, showplot=False):
	
    '''Takes a figure with a subplot (axis) of aspect ratio 1x1 (or
    generates default axes if they are not provided) and adds
    a plot of the upper meridional plane to the axis, averaging the full
    plane with either even or odd symmetry
    '''

    # First, "fold" the field in half, using even symmetry for now
    # Please don't try this with odd N_theta!
    nt = len(cost)
    it_half = int(nt/2)
    
    # Field is ordered from theta=pi (south pole) to theta=0 (north pole)
    # Average the Northern and Southern hemispheres together (must flip 
    # the Southern hemisphere)
    if sym=='even':
        field = 0.5*(field[it_half:, :] +\
                np.flip(field[:it_half, :], axis=0))
    elif sym=='odd':
        field = 0.5*(field[it_half:, :] -\
                np.flip(field[:it_half, :], axis=0))
    cost = cost[it_half:]
    lats = (np.pi/2. - np.arccos(cost))*180./np.pi
    sint = sint[it_half:]
    nt = len(sint)

    # Derivative grid info
    nr = len(rr)
    ri, ro = np.min(rr), np.max(rr)
    shell_depth = ro - ri
    rr_depth = (np.max(rr) - rr)/shell_depth

    # Get default bounds if not specified
    if minmax is None:
        # Compute the indices beyond +/- 75 degrees 
        # latitude, which usually shouldn't be included in any sort 
        # bounds estimates.
        lat_cutoff = 75.
        it_cut = np.argmin(np.abs(lats - lat_cutoff))

        # Also stay away from within 5 percent of top and bottom!
        ir_cuttop = np.argmin(np.abs(rr_depth - 0.05))
        ir_cutbot = np.argmin(np.abs(rr_depth - 0.95))
        field_cut = field[:it_cut, ir_cuttop:ir_cutbot + 1]
        minmax = get_satvals(field_cut, posdef, logscale)
    
    # Need these if logscale is True; made need them for other stuff later
    minexp, maxexp = get_exp(minmax[0]), get_exp(minmax[1])

    # Get the exponent to use for scientific notation
    # and normalize the field by 10**exp
    if not logscale:
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = float(np.floor(np.log10(maxabs)))
        divisor = 10.**exp
        
        # Normalize field by divisor
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor
    
    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by1()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot

    # Get the position of the axes on the figure
    pos = ax.get_position().get_points()
    ax_left, ax_bottom = pos[0]
    ax_right, ax_top = pos[1]
    ax_width = ax_right - ax_left
    ax_height = ax_top - ax_bottom
    ax_aspect = ax_height/ax_width
   
    # Set the colorbar axis to be in the "cavity" of the meridional plane
    # The colorbar height is set by making sure it "fits" in the cavity
    chi = ri/ro # aspect ratio of the shell
    cavity_height = ax_height*chi
    cbax_aspect = 10
    cbax_height = 0.7*cavity_height
    cbax_width = cbax_height/cbax_aspect / ax_aspect
    
    cbax_left = ax_left + 0.1*ax_width
    cbax_bottom = ax_bottom + 0.1*ax_height
   
    # Calculate the grid on which to plot
    rr2d = rr.reshape((1, nr))
    cost2d = cost.reshape((nt, 1))
    sint2d = sint.reshape((nt, 1))
    xx = rr2d*sint2d/ro
    zz = rr2d*cost2d/ro
    
    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    lw = 1
    contour_lw = .2
    
    if (plotfield):
        plt.sca(ax)
#        levs = np.linspace(minmax[0], minmax[1], 100)
#        im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
#        levels=levs)
#        plt.sca(ax)
        if logscale:
            plt.pcolormesh(xx, zz, field, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]))
        else:
            if posdef:
                cmap = 'plasma'
            else:
                cmap = 'RdYlBu_r'
            plt.pcolormesh(xx, zz, field, vmin=minmax[0], vmax=minmax[1],\
                    cmap=cmap)
        cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                       cbax_width, cbax_height])
        cbar = plt.colorbar(cax=cbaxes)

        cbaxes.tick_params(labelsize=fsize)
        cbar.ax.tick_params(labelsize=fsize)   #font size for the ticks

        if not logscale:
            if posdef:
                mid = (minmax[0] + minmax[1])/2.
                ticks = np.array([minmax[0], mid, minmax[1]])
            else:
                ticks = np.array([minmax[0], 0., minmax[1]])
            ticklabels = []
            for i in range(len(ticks)):
                ticklabels.append(str(round(ticks[i], 1)))
            ticks = np.array(ticks)
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticklabels)
            cbar_label = (r'$\times10^{%i}\ $' %exp) + units
        else:
            cbar_label = units
        # Put the units and exponent to left of colorbar
        fig.text(cbax_left - 0.3*cbax_width, cbax_bottom + cbax_height/2.,\
                cbar_label, ha='right', va='center', rotation=90,\
                fontsize=fsize)

    # Plot the boundary of the meridional plane
    plt.sca(ax)
    # outer boundary
    plt.plot(rr[0]/ro*sint, rr[0]/ro*cost, 'k', linewidth=lw)
    # inner boundary
    plt.plot(rr[-1]/ro*sint, rr[-1]/ro*cost, 'k', linewidth=lw) 
    # polar "edge"
    plt.plot([0,0], [rr[-1]/ro, rr[0]/ro], 'k', linewidth=lw) 
    # equatorial "edge"
    plt.plot([rr[-1]/ro, rr[0]/ro], [0,0], 'k', linewidth=lw) 

    # Plot latitude lines, if desired
    if plotlatlines: 
        lats = (np.pi/2. - np.arccos(cost))*180./np.pi
        lats_to_plot = np.arange(0., 90., 15.) 
        for lat in lats_to_plot:
            it = np.argmin(np.abs(lats - lat))
            xx_in, zz_in = rr[-1]/ro*sint[it], rr[-1]/ro*cost[it]
            xx_out, zz_out = rr[0]/ro*sint[it], rr[0]/ro*cost[it]
            plt.sca(ax)
            plt.plot([xx_in, xx_out], [zz_in, zz_out], 'k',\
                    linewidth=contour_lw)

    # Set axis ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-lilbit, 1 + lilbit))
    ax.axis('off') 

    # Plot contours in the meridional plane, if desired
    if plotcontours:
        # Determine the contour levels
        if levels is None:
            if logscale:
                levels = np.logspace(minexp, maxexp, nlevs)
            else:
                levels = np.linspace(minmax[0], minmax[1], nlevs)
        else: # the caller specified specific contour levels to plot!
            levels = np.array(levels)
        # Determine how to color the contours
        if logscale:
            contour_color = 'r'
        elif posdef:
            contour_color = 'w'
        else:
            contour_color = 'k'
        # plot the contours
        plt.sca(ax)
        plt.contour(xx, zz, field, colors=contour_color, levels=levels,\
                linewidths=contour_lw)
    if showplot:
        plt.show()

def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and cost and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the rr and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 cost from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), cost(N_Theta)
    order      = If greater than zero, integration begins at the
                 outer shell and the north pole and proceeds
                 inward and southward.  If less than zero,
                 integration begins at the inner shell and 
                 south pole and proceeds upward and northward.
                 If equal to zero, both are done and an average
                 is taken.
   ------------------------------------------------------------
    OUTPUTS:
   
    psi = the streamfunction
   ------------------------------------------------------------
    """

    n_t,n_r=np.shape(vr)
    dtheta = np.zeros(n_t)
    dr     = np.zeros(n_r)

    psi = np.zeros((n_t,n_r))

    dpsi_dr = np.zeros((n_t,n_r))
    dpsi_dt = np.zeros((n_t,n_r))

    theta = np.arccos(cost)
    sint  = np.sqrt(1.0-cost**2)

    for i in range(n_t):
        dpsi_dr[i,:] = -r*sint[i]*vt[i,:]
        dpsi_dt[i,:] = r*r*sint[i]*vr[i,:]

    if (order >= 0):
        # double precision accumulation
        dtheta[1:n_t] = theta[1:n_t]-theta[0:n_t-1]
        dr[1:n_r] = r[1:n_r]-r[0:n_r-1]

        dtheta[0]=0 
        dr[0]=0

        for i in range(1,n_r):
            psi[1:n_t,i] = psi[1:n_t,i-1] + dpsi_dr[1:n_t,i]*dr[i]
        for i in range(1,n_t):
            psi[i,1:n_r] = psi[i-1,1:n_r] + dpsi_dt[i,1:n_r]*dtheta[i]

    if (order <= 0):
        psi2=np.zeros((n_t,n_r))
        
        dtheta[0:n_t-1] = theta[0:n_t-1]-theta[1:n_t]
        dr[0:n_r-1] = r[0:n_r-1]-r[1:n_r]
        
        dtheta[n_t-1]=0 
        dr[n_r-1]=0
        
        for i in range (n_r-2, -1, -1): 
            psi[0:n_t-1,i] = psi[0:n_t-1,i+1] + dpsi_dr[0:n_t-1,i]*dr[i]
        for i in range(n_t-2, -1, -1):
            psi[i,0:n_r-1] = psi[i+1,0:n_r-1] + dpsi_dt[i,0:n_r-1]*dtheta[i]
        
        if (order < 0):
            return psi2
        else:
            psi=0.5*(psi+psi2)
            
    return psi
