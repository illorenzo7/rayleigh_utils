############################################################################ 
# This file contains utilities useful for plotting the data from 
# Equatorial_Slices
# Written by Loren Matilsky
# Created: 10/02/2019

import numpy as np
import matplotlib as mpl
from matplotlib import ticker
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
plt.rcParams['contour.negative_linestyle'] = 'solid'
from common import *
from plotcommon import axis_range, default_axes_1by2, default_axes_1by1

def plot_eqslice(field, rr, phi, fig=None, ax=None, cmap='RdYlBu_r',\
    units='', minmax=None, posdef=False, logscale=False, symlog=False,\
    plotcontours=True, plotfield=True, nlevs=None, levels=None,\
	plotlonlines=False, rvals=None, rvals_norm=None, cbar_fs=10,\
    showplot=False, plot_cbar=True, lw=1., linthresh=None, linscale=None,\
    plotboundary=True):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x1
    and adds a plot of [field] in the equatorial plane to the axes,
    with colorbar beneath the equatorial plane.'''

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed (shouldn't fucking have to do this)
    field = np.copy(field)
    rr = np.copy(rr)
    phi = np.copy(phi)

    # Derivative grid info
    ri, ro = np.min(rr), np.max(rr) # inner/ outer radii
    nr = len(rr)
    nphi = len(phi) + 1 # I predict I will extend phi...

    # now extend phi ...
    phi = np.hstack((phi, phi[0]))
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    # Extend the field
    field = np.vstack((field, field[0].reshape((1, nr))))

    rr_2d = rr.reshape((1, nr))
    phi_2d = phi.reshape((nphi, 1))

    xx = rr_2d*np.cos(phi_2d)/ro
    yy = rr_2d*np.sin(phi_2d)/ro

    # Deal with saturation values for the field
    
    # If using logscale, you better have positive values!
    if logscale:
        posdef = True
	
    # Get default bounds if not specified
    if minmax is None:
        minmax = get_satvals(field, posdef=posdef,\
                logscale=logscale, symlog=symlog)

    # Factor out the exponent on the field and put it on the color bar
    # for the linear-scaled color bars (default and posdef)
    if not (logscale or symlog) and plotfield:
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = float(np.floor(np.log10(maxabs)))
        divisor = 10.**exp
        
        # Normalize field by divisor
        field /= divisor
        minmax = minmax[0]/divisor, minmax[1]/divisor

    # Saturate the array (otherwise contourf will show white areas)
    saturate_array(field, minmax[0], minmax[1])

    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by1()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot

    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    contour_lw = 0.2
    
    if (plotfield):
        if logscale:
            log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
            levs = np.logspace(log_min, log_max, 150)
            im = ax.contourf(xx, yy, field, cmap='Greys',\
                norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
                levels=levs)  
        elif posdef:
            levs = np.linspace(minmax[0], minmax[1], 150)
            im = ax.contourf(xx, yy, field, cmap='plasma', levels=levs)
        elif symlog:
            linthresh_default, linscale_default =\
                get_symlog_params(field, field_max=minmax[1])
            if linthresh is None:
                linthresh = linthresh_default
            if linscale is None:
                linscale = linscale_default
            log_thresh = np.log10(linthresh)
            log_max = np.log10(minmax[1])
            nlevs_per_interval = 100
            levels_neg = -np.logspace(log_max, log_thresh,\
                    nlevs_per_interval,\
                    endpoint=False)
            levels_mid = np.linspace(-linthresh, linthresh,\
                    nlevs_per_interval, endpoint=False)
            levels_pos = np.logspace(log_thresh, log_max,\
                    nlevs_per_interval)
            levs = np.hstack((levels_neg, levels_mid, levels_pos))
            im = ax.contourf(xx, yy, field, cmap='RdYlBu_r',\
                norm=colors.SymLogNorm(linthresh=linthresh,\
                linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
                levels=levs)
        else:
            im = ax.contourf(xx, yy, field, cmap='RdYlBu_r',\
                    levels=np.linspace(minmax[0], minmax[1], 150))                
        if plot_cbar:
            # Get the position of the axes on the figure
            ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
            ax_width = ax_right - ax_left
            ax_height = ax_top - ax_bottom
            fig_width_inches, fig_height_inches = fig.get_size_inches()
            fig_aspect = fig_height_inches/fig_width_inches
          
            # Set the colorbar ax to be just below equatorial plane
            cbax_aspect = 1./20.
            cbax_width = 0.5*ax_width
            if symlog:
                cbax_width *= 1.5
                cbax_aspect /= 1.5

            cbax_height = cbax_width*cbax_aspect/fig_aspect
            
            cbax_left = ax_left + 0.5*(ax_width - cbax_width)
            cbax_bottom = ax_bottom - 1/8/fig_height_inches - cbax_height
            
            cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                           cbax_width, cbax_height])
            cbar = plt.colorbar(im, cax=cbaxes, orientation='horizontal')

    
            cbaxes.tick_params(labelsize=cbar_fs)
            cbar.ax.tick_params(labelsize=cbar_fs)   #font size for the ticks

            if logscale:
                locator = ticker.LogLocator(subs='all')
                cbar.set_ticks(locator)
                cbar_label = units
            elif posdef:
                cbar_label = (r'$\times10^{%i}\ $' %exp) + units
                cbar.set_ticks([minmax[0], minmax[1]])
                cbar.set_ticklabels(['%1.1f' %minmax[0], '%1.1f' %minmax[1]])
            elif symlog:
                cbar_label = units
                cbar.set_ticks([-minmax[1], -linthresh, 0, linthresh,\
                        minmax[1]])
                cbar.set_ticklabels([sci_format(-minmax[1]),\
                        sci_format(-linthresh), '0', sci_format(linthresh),\
                        sci_format(minmax[1])])
        #            cax.minorticks_on()
            else:
                cbar_label = (r'$\times10^{%i}\ $' %exp) + units
                cbar.set_ticks([minmax[0], 0, minmax[1]])
                cbar.set_ticklabels(['%1.1f' %minmax[0], '0', '%1.1f'\
                        %minmax[1]])
    
            # Put the units (and possibly the exponent) to left of colorbar
            fig.text(cbax_left + cbax_width + 1/16/fig_width_inches,\
                    cbax_bottom + 0.5*cbax_height, cbar_label, ha='left',\
                    va='center', fontsize=cbar_fs)

    # Plot contours in the equatorial plane, if desired
    if plotcontours:
        # Determine the contour levels
        if levels is None:
            if nlevs is None:
                nlevs = 11 # Default nlevs to 11
            if logscale:
                min_log = np.log10(minmax[0])
                max_log = np.log10(minmax[1])
                levels = np.logspace(min_log, max_log, nlevs)
            elif symlog:
                log_thresh = np.log10(linthresh)
                log_max = np.log10(minmax[1])
                nlevs_per_interval = nlevs//3
                levels_neg = -np.logspace(log_max, log_thresh,\
                        nlevs_per_interval, endpoint=False)
                levels_mid = np.linspace(-linthresh, linthresh,\
                        nlevs_per_interval, endpoint=False)
                levels_pos = np.logspace(log_thresh, log_max,\
                        nlevs_per_interval)
                levels = np.hstack((levels_neg, levels_mid, levels_pos))
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
        ax.contour(xx, yy, field, colors=contour_color, levels=levels,\
                linewidths=contour_lw)

    # Plot the boundary of the meridional plane, if desired
    if plotboundary:
        plt.sca(ax)
        plt.plot(rr[0]/ro*cosphi, rr[0]/ro*sinphi, 'k', linewidth=lw)
        plt.plot(rr[-1]/ro*cosphi, rr[-1]/ro*sinphi, 'k', linewidth=lw)

    # Plot longititude lines, if desired
    if plotlonlines: 
        lons = phi*180./np.pi
        lons_to_plot = np.arange(0., 360., 15.) 
        for lon in lons_to_plot:
            iphi = np.argmin(np.abs(lons - lon))
            xx_in, yy_in = rr[-1]/ro*cosphi[iphi], rr[-1]/ro*sinphi[iphi]
            xx_out, yy_out = rr[0]/ro*cosphi[iphi], rr[0]/ro*sinphi[iphi]
            plt.sca(ax)
            plt.plot([xx_in, xx_out], [yy_in, yy_out], 'k',\
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
            plt.plot(rval_dim*cosphi, rval_dim*sinphi, 'k--', linewidth=.7)

    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-1 - lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

    if showplot:
        plt.show()
