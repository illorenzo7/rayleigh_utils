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
from common import get_satvals

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def default_axes_2by1():
    # Create plot
    subplot_width_inches = 2.5
    subplot_height_inches = 5.
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

def plot_azav(field, radius, costheta, sintheta, fig=None, ax=None,\
	cmap='RdYlBu_r', units='', minmax=None, posdef=False,\
	plotcontours=True, plotfield=True, nlevs=10, levels=None,\
	plotlatlines=False, norm=None, fsize=8, showplot=False):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x2
    and adds a plot of [field] in the meridional plane to the axes,\
    with colorbar in the "cavity" of the meridional plane fuck around'''
	
    # Compute the indices beyond +/- 75 degrees 
    # latitude, which usually shouldn't be included in any sort 
    # bounds estimates.
    lats = (np.pi/2. - np.arccos(costheta))*180./np.pi
    lat_cutoff = 75.
    it_cutm = np.argmin(np.abs(lats + lat_cutoff))
    it_cutp = np.argmin(np.abs(lats - lat_cutoff))

    # Get default bounds if not specified
    if minmax is None:
        field_cut = field[it_cutm:it_cutp+1, :]
        if posdef:
            mini, maxi = get_satvals(field_cut, posdef=True)
        else:
            sig = np.std(field_cut, :])
            mini, maxi = -3*sig, 3*sig
    else:
        mini, maxi = minmax

    if not posdef:
        # Get the exponent to use for scientific notation
        maxabs = max(np.abs(mini), np.abs(maxi))
        exp = int(np.floor(np.log10(maxabs)))
        divisor = 10**exp
        
        # Normalize field by divisor
        field /= divisor
        mini /= divisor
        maxi /= divisor

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
   
    # Set the colorbar ax to be in the "cavity" of the meridional plane
    # The colorbar height is set by making sure it "fits" in the cavity
    chi = np.min(radius)/np.max(radius)
    cavity_height = ax_height*chi
    cbax_center_x = ax_left + 0.3*ax_width
    cbax_center_y = ax_bottom + ax_height/2.
    cbax_aspect = 10.
    cbax_height = 0.5*cavity_height
    cbax_width = cbax_height/cbax_aspect / ax_aspect
    
    cbax_left = cbax_center_x - cbax_width/2.
    cbax_bottom = cbax_center_y - cbax_height/2.
    
    # Calculate the grid on which to plot
    rr = radius/np.max(radius)
    nr = len(rr)
    nt = len(costheta)
    rr2 = rr.reshape((1, nr))
    cost2 = costheta.reshape((nt, 1))
    sint2 = sintheta.reshape((nt, 1))
    xx = rr*sint2
    zz = rr2*cost2

    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    lw = 1
    contour_lw = .2
    
    if (plotfield):
        plt.sca(ax)
# norm=colors.LogNorm(vmin=my_min, vmax=my_max),\
        if posdef:
            im = ax.contourf(xx, zz, field, cmap='Greys',\
            norm=colors.LogNorm(vmin=mini, vmax=maxi),\
            levels=np.logspace(minexp, maxexp, 50, base=np.exp(1.)))
        else:
            im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
            levels=np.linspace(my_min, my_max, 50),\
            norm=MidpointNormalize(0))
#        plt.pcolormesh(xx, zz, field, vmin=mini, vmax=maxi, cmap=cmap,\
#                norm=norm)
        cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                       cbax_width, cbax_height])
        cbar = plt.colorbar(cax=cbaxes)

        #fsize = 8 # fontsize for colorbar ticks and labels
        cbaxes.tick_params(labelsize=fsize)
        cbar.ax.tick_params(labelsize=fsize)   #font size for the ticks
        ticks = np.array([mini, 0, maxi])
        ticklabels = []
        for i in range(len(ticks)):
            ticklabels.append(str(round(ticks[i], 1)))
        ticks = np.array(ticks)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticklabels)

        # Put the units and exponent to left of colorbar
        cbar_label = (r'$\times10^{%i}\ $' %exp) + units
        fig.text(cbax_left - 0.3*cbax_width, cbax_center_y, cbar_label,\
            ha='right', va='center', rotation=90, fontsize=fsize)

    # Plot the boundary of the meridional plane
    plt.sca(ax)
    plt.plot(rr[0]*sintheta, rr[0]*costheta, 'k', linewidth=lw)
    plt.plot(rr[-1]*sintheta, rr[-1]*costheta, 'k', linewidth=lw)
    plt.plot([0.,0.], [rr[-1], rr[0]], 'k', linewidth=lw)
    plt.plot([0.,0.], [-rr[-1], -rr[0]], 'k', linewidth=lw)

    # Plot latitude lines, if desired
    if plotlatlines:
        lats_to_plot = np.arange(-75, 90, 15)
        for lat in lats_to_plot:
            theta_val = (90 - lat)*np.pi/180
            x_in, z_in = r[-1]*np.sin(theta_val), r[-1]*np.cos(theta_val)
            x_out, z_out = r[0]*np.sin(theta_val), r[0]*np.cos(theta_val)
            plt.sca(ax)
            plt.plot([x_in, x_out], [z_in, z_out], 'k',\
                    linewidth=contour_lw)


    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

    if (plotcontours):
        if levels is None:
            levs = np.linspace(mini, maxi, nlevs)
        else: # the caller specified specific contour levels to plot!
            levs=np.array(levels)
        plt.contour(xx, zz, field, colors='k', levels=levs,\
                linewidths=contour_lw)

    if showplot:
        plt.show()

def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and costheta and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the radius and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 costheta from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), costheta(N_Theta)
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
