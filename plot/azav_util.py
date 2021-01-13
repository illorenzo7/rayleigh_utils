# This file contains utilities useful for plotting the data from azimuthal
# averages,
# files in the directory AZ_Avgs.
# Written by Nick Featherstone
# First modified by Loren Matilsky, 03/09/2018

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

def plot_azav(field, rr, cost, fig=None, ax=None, cmap='RdYlBu_r',\
    units='', minmax=None, posdef=False, logscale=False, symlog=False,\
    plotcontours=True, plotfield=True, nlevs=10, levels=None,\
	plotlatlines=False, rvals=[], cbar_fs=10,\
    showplot=False, plot_cbar=True, lw_scaling=1.,\
    linthresh=None, linscale=None, plotboundary=True, rbcz=None,\
    minmaxrz=None, linthreshrz=None, linscalerz=None):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x2
    and adds a plot of [field] in the meridional plane to the axes,
    with colorbar in the "cavity" of the meridional plane.'''

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed (shouldn't fucking have to do this)
    field = np.copy(field)
    rr = np.copy(rr)
    cost = np.copy(cost)
    sint = np.sqrt(1. - cost**2.)

    # Derivative grid info
    ri, ro = np.min(rr), np.max(rr) # inner/ outer radii
    nr = len(rr)
    nt = len(cost)

    rr_2d = rr.reshape((1, nr))
    cost_2d = cost.reshape((nt, 1))
    sint_2d = sint.reshape((nt, 1))

    if rbcz is None: # just plotting 1 domain
        xx = rr_2d*sint_2d/ro
        zz = rr_2d*cost_2d/ro
    else:
        irbcz = np.argmin(np.abs(rr - rbcz))
        fieldcz = field[:, :irbcz+1]
        fieldrz = field[:, irbcz+1:]
        rrcz = rr[:irbcz+1]
        rrrz = rr[irbcz+1:]
        xxcz = (rr_2d*sint_2d)[:, :irbcz+1]/ro
        zzcz = (rr_2d*cost_2d)[:, :irbcz+1]/ro
        xxrz = (rr_2d*sint_2d)[:, irbcz+1:]/ro
        zzrz = (rr_2d*cost_2d)[:, irbcz+1:]/ro

    # Deal with saturation values for the field
    
    # If using logscale, you better have positive values!
    if logscale:
        posdef = True
	
    # Get default bounds (for CZ) if not specified
    if minmax is None:
        # Cut away the data near the domain boundaries
        if rbcz is None:
            trimmed_field = trim_field(field, rr, cost)
            minmax = get_satvals(trimmed_field, posdef=posdef,\
                    logscale=logscale, symlog=symlog)
        else:
            trimmed_fieldcz = trim_field(fieldcz, rrcz, cost)
            minmax = get_satvals(trimmed_fieldcz, posdef=posdef,\
                    logscale=logscale, symlog=symlog)

    # Get default bounds (for RZ) if not specified
    if not rbcz is None:
        if minmaxrz is None:
            # Cut away the data near the domain boundaries
            trimmed_fieldrz = trim_field(fieldrz, rrrz, cost)
            minmaxrz = get_satvals(trimmed_fieldrz, posdef=posdef,\
                    logscale=logscale, symlog=symlog)

    # Factor out the exponent on the field and put it on the color bar
    # for the linear-scaled color bars (default and posdef)
    if not (logscale or symlog) and plotfield:
        if rbcz is None:
            maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
            exp = get_exp(maxabs)
            divisor = 10.**exp
            
            # Normalize field by divisor
            field /= divisor
            minmax = minmax[0]/divisor, minmax[1]/divisor
        else:
            maxabscz = max(np.abs(minmax[0]), np.abs(minmax[1]))
            maxabsrz = max(np.abs(minmaxrz[0]), np.abs(minmaxrz[1]))
            exp = expcz = get_exp(maxabscz)
            exprz = get_exp(maxabsrz)
            divisorcz = 10.**expcz
            divisorrz = 10.**exprz
            
            # Normalize fields by divisor
            fieldcz /= divisorcz
            fieldrz /= divisorrz
            minmax = minmax[0]/divisorcz, minmax[1]/divisorcz
            minmaxrz = minmaxrz[0]/divisorrz, minmaxrz[1]/divisorrz

    if rbcz is None:
        # Saturate the array (otherwise contourf will show white areas)
        saturate_array(field, minmax[0], minmax[1])
    else:
        saturate_array(fieldcz, minmax[0], minmax[1])
        saturate_array(fieldrz, minmaxrz[0], minmaxrz[1])

    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by2()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot

    # Specify linewidths to be used in the meridional plane, one for the 
    # boundary (lw) and one for the contours (contour_lw)
    contour_lw = 0.3*lw_scaling
    
    if plotfield:
        if rbcz is None:
            if logscale:
                log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
                levs = np.logspace(log_min, log_max, 150)
                im = ax.contourf(xx, zz, field, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)  
            elif posdef:
                levs = np.linspace(minmax[0], minmax[1], 150)
                im = ax.contourf(xx, zz, field, cmap='plasma', levels=levs)
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
                im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
                    norm=colors.SymLogNorm(linthresh=linthresh,\
                    linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)
            else:
                im = ax.contourf(xx, zz, field, cmap='RdYlBu_r',\
                        levels=np.linspace(minmax[0], minmax[1], 150))                
        else:
            if logscale:
                # First plot field in CZ
                log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
                levs = np.logspace(log_min, log_max, 150)
                im = ax.contourf(xxcz, zzcz, fieldcz, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)  
                # Now plot field in RZ
                log_min, log_max = np.log10(minmaxrz[0]),\
                        np.log10(minmaxrz[1])
                levs = np.logspace(log_min, log_max, 150)
                imrz = ax.contourf(xxrz, zzrz, fieldrz, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmaxrz[0],\
                    vmax=minmaxrz[1]), levels=levs)  
            elif posdef:
                # First plot field in CZ
                levs = np.linspace(minmax[0], minmax[1], 150)
                im = ax.contourf(xxcz, zzcz, fieldcz, cmap='plasma',\
                        levels=levs)
                # Then plot field in RZ
                levs = np.linspace(minmaxrz[0], minmaxrz[1], 150)
                imrz = ax.contourf(xxrz, zzrz, fieldrz, cmap='plasma',\
                        levels=levs)
            elif symlog:
                # First plot field in CZ
                linthresh_default, linscale_default =\
                    get_symlog_params(fieldcz, field_max=minmax[1])
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
                im = ax.contourf(xxcz, zzcz, fieldcz, cmap='RdYlBu_r',\
                    norm=colors.SymLogNorm(linthresh=linthresh,\
                    linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)

                # Then plot field in RZ
                linthresh_default, linscale_default =\
                    get_symlog_params(fieldrz, field_max=minmaxrz[1])
                if linthreshrz is None:
                    linthreshrz = linthresh_default
                if linscalerz is None:
                    linscalerz = linscale_default

                log_thresh = np.log10(linthreshrz)
                log_max = np.log10(minmaxrz[1])
                nlevs_per_interval = 100

                levels_neg = -np.logspace(log_max, log_thresh,\
                        nlevs_per_interval,\
                        endpoint=False)
                levels_mid = np.linspace(-linthreshrz, linthreshrz,\
                        nlevs_per_interval, endpoint=False)
                levels_pos = np.logspace(log_thresh, log_max,\
                        nlevs_per_interval)
                levs = np.hstack((levels_neg, levels_mid, levels_pos))

                imrz = ax.contourf(xxrz, zzrz, fieldrz, cmap='RdYlBu_r',\
                    norm=colors.SymLogNorm(linthresh=linthreshrz,\
                    linscale=linscalerz, vmin=minmaxrz[0],\
                    vmax=minmaxrz[1]), levels=levs)
            else:
                # First plot field in CZ
                im = ax.contourf(xxcz, zzcz, fieldcz, cmap='RdYlBu_r',\
                        levels=np.linspace(minmax[0], minmax[1], 150)) 
                # Then plot field in RZ
                imrz = ax.contourf(xxrz, zzrz, fieldrz, cmap='RdYlBu_r',\
                        levels=np.linspace(minmaxrz[0], minmaxrz[1], 150)) 

        if plot_cbar:
            ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
            ax_delta_x = ax_xmax - ax_xmin
            cbar_aspect = 1./20.
            cbar_width = 0.75*ax_delta_x # make cbar a fraction as long as\
                    # plot is wide            
            fig_width_inches, fig_height_inches = fig.get_size_inches()
            fig_aspect = fig_height_inches/fig_width_inches
            cbar_height = cbar_width*cbar_aspect/fig_aspect
            cbar_bottom = ax_ymin - 2.5*cbar_height
            cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
            cax = fig.add_axes((cbar_left, cbar_bottom, cbar_width,\
                    cbar_height))        
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
                
            cax.tick_params(labelsize=cbar_fs)
            cbar.ax.tick_params(labelsize=cbar_fs)   
            # font size for the ticks

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
                cbar_label = (r'$\times10^{%i}\ $' %exp) + units
                cbar.set_ticks([minmax[0], 0, minmax[1]])
                cbar.set_ticklabels(['%1.1f' %minmax[0], '0', '%1.1f'\
                        %minmax[1]])
    
            # Title the colorbar based on the field's units
            line_height = 1./4./fig_height_inches
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
                line_height = 1./4./fig_height_inches
                fig.text(cbar_left + 0.5*cbar_width, cbar_bottom -\
                        line_height, cbar_label, ha='center', va='top',\
                        **csfont, fontsize=cbar_fs) 

    # Plot contours in the meridional plane, if desired
    if plotcontours:
        if rbcz is None:
            # Determine the contour levels
            if levels is None:
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
            ax.contour(xx, zz, field, colors=contour_color, levels=levels,\
                    linewidths=contour_lw)
        else:
            # Determine the contour levels
            if levels is None:
                if logscale:
                    min_log = np.log10(minmax[0])
                    max_log = np.log10(minmax[1])
                    levelscz = np.logspace(min_log, max_log, nlevs)
                    min_log = np.log10(minmaxrz[0])
                    max_log = np.log10(minmaxrz[1])
                    levelsrz = np.logspace(min_log, max_log, nlevs)
                elif symlog:
                    log_thresh = np.log10(linthresh[0])
                    log_max = np.log10(minmax[1])
                    nlevs_per_interval = nlevs//3
                    levels_neg = -np.logspace(log_max, log_thresh,\
                            nlevs_per_interval, endpoint=False)
                    levels_mid = np.linspace(-linthresh, linthresh,\
                            nlevs_per_interval, endpoint=False)
                    levels_pos = np.logspace(log_thresh, log_max,\
                            nlevs_per_interval)
                    levelscz = np.hstack((levels_neg, levels_mid,\
                            levels_pos))
                    log_thresh = np.log10(linthreshrz)
                    log_max = np.log10(minmaxrz[1])
                    nlevs_per_interval = nlevs//3
                    levels_neg = -np.logspace(log_max, log_thresh,\
                            nlevs_per_interval, endpoint=False)
                    levels_mid = np.linspace(-linthresh, linthresh,\
                            nlevs_per_interval, endpoint=False)
                    levels_pos = np.logspace(log_thresh, log_max,\
                            nlevs_per_interval)
                    levelsrz = np.hstack((levels_neg, levels_mid,\
                            levels_pos))
                else:
                    levelscz = np.linspace(minmax[0], minmax[1], nlevs)
                    levelsrz = np.linspace(minmaxrz[0], minmax[1], nlevs)
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
            ax.contour(xxcz, zzcz, fieldcz, colors=contour_color,\
                    levels=levelscz, linewidths=contour_lw)
            ax.contour(xxrz, zzrz, fieldrz, colors=contour_color,\
                    levels=levelsrz, linewidths=contour_lw)

    # Plot the boundary of the meridional plane, if desired
    if plotboundary:
        plt.sca(ax)
        plt.plot(rr[0]/ro*sint, rr[0]/ro*cost, 'k', linewidth=lw_scaling)
        plt.plot(rr[-1]/ro*sint, rr[-1]/ro*cost, 'k', linewidth=lw_scaling)
        plt.plot([0.,0.], [rr[-1]/ro, rr[0]/ro], 'k', linewidth=lw_scaling)
        plt.plot([0.,0.], [-rr[-1]/ro, -rr[0]/ro], 'k', linewidth=lw_scaling)
        # Look for consecutive repeated values in radius array: 
        # these will correspond to domain divisions, which should be plotted
        for ir in range(1, nr):
            if rr[ir] == rr[ir - 1]:
                plt.plot(rr[ir]/ro*sint, rr[ir]/ro*cost, 'k',\
                        linewidth=0.5*lw_scaling)

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
    # If rbcz has been provided, it should be one of the radii plotted
    # crude way to ensure "rvals" is passed by value, not reference
    tmp = list(rvals)
    del rvals
    rvals = tmp
    if not rbcz is None:
        rvals.append(rbcz)
    
    for rval in rvals: 
        plt.sca(ax)
        rval /= ro
        # "dimensional" rval (in dimensions of ro!)
        plt.plot(rval*sint, rval*cost, 'k--',\
                linewidth=.7*lw_scaling)

    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

    if showplot:
        plt.show()
    if plotfield:
        return im

def plot_azav_half(field, rr, cost, sym='even', fig=None, ax=None,\
        cmap='RdYlBu_r', units='', minmax=None, posdef=False, logscale=False,\
        symlog=False, linthresh=None, linscale=None, plotcontours=True,\
        plotfield=True, nlevs=10, levels=None, plotlatlines=False, rvals=[],\
        norm=None, cbar_fs=10, showplot=False, plot_cbar=True):
	
    '''Takes a figure with a subplot (axis) of aspect ratio 1x1 (or
    generates default axes if they are not provided) and adds
    a plot of the upper meridional plane to the axis, averaging the full
    plane with either even or odd symmetry
    '''

    # Get default bounds if not specified
    # Do this before "folding" the field
    if minmax is None:
        field_cut = trim_field(np.copy(field), rr, cost)
        minmax = get_satvals(field_cut, posdef=posdef, logscale=logscale,\
                symlog=symlog)
        
    # Grid info -- make "folded" grid
    nr = len(rr)
    nt = len(cost)
    ri, ro = np.min(rr), np.max(rr)
    
    # "Fold" the field in half
    # Please don't try this with odd N_theta!
    it_half = int(nt/2)
    
    # "Half grid" quantities
    cost = cost[it_half:]
    lats = (np.pi/2. - np.arccos(cost))*180./np.pi
    sint = np.sqrt(1. - cost**2.)
    
    # Calculate the grid on which to plot
    rr_2d = rr.reshape((1, nr))
    cost_2d = cost.reshape((int(nt/2), 1))
    sint_2d = sint.reshape((int(nt/2), 1))
    xx = rr_2d*sint_2d/ro
    zz = rr_2d*cost_2d/ro    
    # Field is ordered from theta=pi (south pole) to theta=0 (north pole)
    # Average the Northern and Southern hemispheres together (must flip 
    # the Southern hemisphere with respect to latitude, then add or subtract it
    if sym=='even':
        field = 0.5*(field[it_half:, :] +\
                np.flip(field[:it_half, :], axis=0))
    elif sym=='odd':
        field = 0.5*(field[it_half:, :] -\
                np.flip(field[:it_half, :], axis=0))

    # Get the exponent to use for scientific notation
    # and normalize the field by 10**exp
    if not logscale:
        maxabs = max(np.abs(minmax[0]), np.abs(minmax[1]))
        exp = get_exp(maxabs)
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
        
        # Get the position of the axes on the figure
        ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
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
        
        cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                       cbax_width, cbax_height])
        cbar = plt.colorbar(cax=cbaxes)

        cbaxes.tick_params(labelsize=cbar_fs)
        cbar.ax.tick_params(labelsize=cbar_fs)   #font size for the ticks

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
                fontsize=cbar_fs)

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

def plot_quiver(vr, vt, rr, cost, fig=None, ax=None, minmax=None,\
        plotlatlines=False, rvals=[], cbar_fs=10,\
        showplot=False, scale=None, plot_poles=False, nsample_t=20,\
        nsample_r=10, scale_by_mag=True, plotboundary=True):

    ''' Takes (or creates) set of axes with physical aspect ratio 1x2
    and adds a vector (quiver) plot of [vx, vy] in the meridional plane to 
    the axes.'''

    # First things first, make sure Python does not modify any of the 
    # arrays it was passed (shouldn't fucking have to do this)
    vmag = np.sqrt(np.std(vr**2 + vt**2))
    vr = np.copy(vr)/vmag
    vt = np.copy(vt)/vmag
    rr = np.copy(rr)
    cost = np.copy(cost)
    sint = np.sqrt(1. - cost**2.)

    # Derivative grid info
    ri, ro = np.min(rr), np.max(rr)
    nr = len(rr)
    nt = len(cost)

    # Create a default set of figure axes if they weren't already
    # specified by user
    if fig is None or ax is None:
        fig, ax = default_axes_1by2()
        showplot = True # probably in this case the user just
        # ran plot_azav from the command line wanting to view
        # view the plot
   
    # Get the position of the axes on the figure
    pos = ax.get_position().get_points()
    ax_left, ax_bottom = pos[0]
    ax_right, ax_top = pos[1]
    
    # Calculate the grid on which to plot
    rr2d = rr.reshape((1, nr))
    tt = np.arccos(cost)
    cost2d = cost.reshape((nt, 1))
    sint2d = sint.reshape((nt, 1))
    xx = rr2d*sint2d/ro
    zz = rr2d*cost2d/ro

    # Compute the vector field in cylindrical coordinates
    vx = sint2d*vr + cost2d*vt
    vz = cost2d*vr - sint2d*vt

    if not plot_poles:
        it_75_south = np.argmin(np.abs(tt - 11*np.pi/12))
        it_75_north = np.argmin(np.abs(tt - np.pi/12))
        xx = xx[it_75_south:it_75_north]
        zz = zz[it_75_south:it_75_north]
        vx = vx[it_75_south:it_75_north]
        vz = vz[it_75_south:it_75_north]
    if not scale_by_mag:
        v_mag = np.sqrt(vx**2 + vz**2)
        vx /= (20*v_mag)
        vz /= (20*v_mag)

    # Specify linewidths to be used in the meridional plane to plot the 
    # boundaries
    lw = 1.
    contour_lw = 0.2
   
    # Make quiver plot
    plt.sca(ax)
    # Scale down the arrays
    nskip_t, nskip_r = nt//nsample_t, nr//nsample_r
    xx = xx[::nskip_t, ::nskip_r]
    zz = zz[::nskip_t, ::nskip_r]
    vx = vx[::nskip_t, ::nskip_r]
    vz = vz[::nskip_t, ::nskip_r]

    plt.quiver(xx, zz, vx, vz, scale=scale)
   
    # Plot the boundary of the meridional plane
    if plotboundary:
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
    for rval in rvals: 
        plt.sca(ax)
        rval /= ro
        plt.plot(rval*sint, rval*cost, 'k--', linewidth=.7)

    # Set ax ranges to be just outside the boundary lines
    lilbit = 0.01
    ax.set_xlim((-lilbit, 1 + lilbit))
    ax.set_ylim((-1 - lilbit, 1 + lilbit))
    ax.axis('off') 

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
