########################################################################
# This file contains utilities useful for plotting time-yy traces
# of the AZ_Avgs data in Rayleigh  (usually time-latitude or time-radius)
# Author: Loren Matilsky
# Created: 04/06/2020

import numpy as np
import matplotlib as mpl
from matplotlib import ticker
#mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
plt.rcParams['contour.negative_linestyle'] = 'solid'
from common import *
from plotcommon import *

def plot_tl(field, times, yy, fig=None, ax=None, cmap='RdYlBu_r',\
    units='', minmax=None, xminmax=None, posdef=False, logscale=False,\
    symlog=False, plotcontours=False, plotfield=True, nlevs=10,\
    levels=None, plottimes=None, cbar_fs=10, navg=1,\
    lw_scaling=1., showplot=False, plot_cbar=True, linthresh=None,\
    linscale=None, yvals=None, rbcz=None, minmaxrz=None,\
    linthreshrz=None, linscalerz=None, nosci=False, cbar_scaling=1.):

    ''' Takes (or creates) set of axes
    and adds a plot of [field] in time-yy space to the axes,
    with colorbar to the right.'''

    # Work with copy of field (not actual field)
    field = np.copy(field)
    times = np.copy(times)

    # Get dimensions of field
    ntimes = len(times)
    ny = len(yy)

    # First, average field in time
    over2 = navg//2
    field_timeavg = np.zeros((ntimes - navg + 1, ny))

    for i in range(navg):
        field_timeavg += field[i:ntimes - navg + 1 + i]
    field = field_timeavg/navg

    times = times[over2:ntimes - over2]

    if xminmax is None: # By default use all times available
        it1 = 0
        it2 = ntimes - navg
    else:
        it1 = np.argmin(np.abs(times - xminmax[0]))
        it2 = np.argmin(np.abs(times - xminmax[1]))
    times = times[it1:it2+1]
    field = field[it1:it2+1]

    # Make 2D grids from times/yy
    times_2d, yy_2d = np.meshgrid(times, yy, indexing='ij')
    
    if not rbcz is None: # plotting two domains
        irbcz = np.argmin(np.abs(yy - rbcz))
        fieldcz = field[:, :irbcz+1]
        fieldrz = field[:, irbcz+1:]
        times_2d_cz = times_2d[:, :irbcz+1]
        times_2d_rz = times_2d[:, irbcz+1:]
        yy_2d_cz = yy_2d[:, :irbcz+1]
        yy_2d_rz = yy_2d[:, irbcz+1:]

    # If using logscale, you better have positive values!
    if logscale:
        posdef = True
	
    # Get default bounds if not specified
    if rbcz is None:
        if minmax is None:
            trimmed_field = field[:, 1:-1]
            minmax = get_satvals(trimmed_field, posdef=posdef,\
                logscale=logscale, symlog=symlog)
    else:
        if minmax is None:
            trimmed_fieldcz = fieldcz[:, 1:-1]
            minmax = get_satvals(trimmed_fieldcz, posdef=posdef,\
                    logscale=logscale, symlog=symlog)
        if minmaxrz is None:
            # Cut away the data near the domain boundaries
            trimmed_fieldrz = fieldrz[:, 1:-1]
            minmaxrz = get_satvals(trimmed_fieldrz, posdef=posdef,\
                    logscale=logscale, symlog=symlog)

    # Factor out the exponent on the field and put it on the color bar
    # assuming we want scientific notation (turn off by setting nosci=True)
       
    # for the linear-scaled color bars (default and posdef)
    if not (logscale or symlog or nosci) and plotfield:
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

    # Specify linewidths to be used, one for the 
    # boundary (lw) and one for the contours (contour_lw) and 
    # time cuts
    contour_lw = 0.3*lw_scaling

    if plotfield:
        if rbcz is None:
            if logscale:
                log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
                levs = np.logspace(log_min, log_max, 150)
                im = ax.contourf(times_2d, yy_2d, field, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)  
            elif posdef:
                levs = np.linspace(minmax[0], minmax[1], 150)
                im = ax.contourf(times_2d, yy_2d, field, cmap='plasma', levels=levs)
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
                        nlevs_per_interval, endpoint=False)
                levels_mid = np.linspace(-linthresh, linthresh,\
                        nlevs_per_interval, endpoint=False)
                levels_pos = np.logspace(log_thresh, log_max,\
                        nlevs_per_interval)
                levs = np.hstack((levels_neg, levels_mid, levels_pos))
                im = ax.contourf(times_2d, yy_2d, field, cmap='RdYlBu_r',\
                    norm=colors.SymLogNorm(linthresh=linthresh,\
                    linscale=linscale, vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)
            else:
                im = ax.contourf(times_2d, yy_2d, field, cmap='RdYlBu_r',\
                        levels=np.linspace(minmax[0], minmax[1], 150))                
        else:
            if logscale:
                # First plot field in CZ
                log_min, log_max = np.log10(minmax[0]), np.log10(minmax[1])
                levs = np.logspace(log_min, log_max, 150)
                im = ax.contourf(times_2d_cz, yy_2d_cz, fieldcz, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmax[0], vmax=minmax[1]),\
                    levels=levs)  
                # Now plot field in RZ
                log_min, log_max = np.log10(minmaxrz[0]),\
                        np.log10(minmaxrz[1])
                levs = np.logspace(log_min, log_max, 150)
                imrz = ax.contourf(times_2d_rz, yy_2d_rz, fieldrz, cmap='Greys',\
                    norm=colors.LogNorm(vmin=minmaxrz[0],\
                    vmax=minmaxrz[1]), levels=levs)  
            elif posdef:
                # First plot field in CZ
                levs = np.linspace(minmax[0], minmax[1], 150)
                im = ax.contourf(times_2d_cz, yy_2d_cz, fieldcz, cmap='plasma',\
                        levels=levs)
                # Then plot field in RZ
                levs = np.linspace(minmaxrz[0], minmaxrz[1], 150)
                imrz = ax.contourf(times_2d_rz, yy_2d_rz, fieldrz, cmap='plasma',\
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
                im = ax.contourf(times_2d_cz, yy_2d_cz, fieldcz, cmap='RdYlBu_r',\
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

                imrz = ax.contourf(times_2d_rz, yy_2d_rz, fieldrz, cmap='RdYlBu_r',\
                    norm=colors.SymLogNorm(linthresh=linthreshrz,\
                    linscale=linscalerz, vmin=minmaxrz[0],\
                    vmax=minmaxrz[1]), levels=levs)
            else:
                # First plot field in CZ
                im = ax.contourf(times_2d_cz, yy_2d_cz, fieldcz, cmap='RdYlBu_r',\
                        levels=np.linspace(minmax[0], minmax[1], 150)) 
                # Then plot field in RZ
                imrz = ax.contourf(times_2d_rz, yy_2d_rz, fieldrz, cmap='RdYlBu_r',\
                        levels=np.linspace(minmaxrz[0], minmaxrz[1], 150)) 

        if plot_cbar:
            # Get the position of the axes on the figure
            ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
            ax_width = ax_right - ax_left
            ax_height = ax_top - ax_bottom
            fig_width_inches, fig_height_inches = fig.get_size_inches()
            fig_aspect = fig_height_inches/fig_width_inches
          
            # Set the colorbar ax to be just to the right of the subplot,
            # with the same height as the subplot
            cbax_center_y = ax_bottom + ax_height/2.
            cbax_aspect = 20.
            cbax_height = ax_height*cbar_scaling
            cbax_width = cbax_height*fig_aspect/cbax_aspect
            
            cbax_left = ax_right + 1./4./fig_width_inches
            cbax_bottom = ax_bottom
            
            cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                           cbax_width, cbax_height])
            cbar = plt.colorbar(im, cax=cbaxes)
    
            cbaxes.tick_params(labelsize=cbar_fs)
            cbar.ax.tick_params(labelsize=cbar_fs)   #font size for the ticks

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
                ticklabels[0] = sci_format(minmax[0])
                ticklabels[nlog] = sci_format(-linthresh)
                ticklabels[nticks//2] = r'$0$'
                ticklabels[nlog + nlin - 1] = sci_format(linthresh)
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
    
            # Put the units (and possibly the exponent) to left of colorbar
            fig.text(cbax_left - 0.3*cbax_width, cbax_center_y,\
                    cbar_label, ha='right', va='center', rotation=90,\
                    fontsize=cbar_fs)

            if not rbcz is None: # Make a colorbar for the RZ
                cbax_left = ax_right + 1.0/fig_width_inches + cbax_width
                cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                               cbax_width, cbax_height])
                cbar = plt.colorbar(imrz, cax=cbaxes)
        
                cbaxes.tick_params(labelsize=cbar_fs)
                cbar.ax.tick_params(labelsize=cbar_fs)   #font size for the ticks

                if logscale:
                    locator = ticker.LogLocator(subs='all')
                    cbar.set_ticks(locator)
                    cbar_label = units
                elif posdef:
                    cbar_label = (r'$\times10^{%i}\ $' %exp) + units
                    cbar.set_ticks([minmaxrz[0], minmaxrz[1]])
                    cbar.set_ticklabels(['%1.1f' %minmaxrz[0],\
                            '%1.1f' %minmaxrz[1]])
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
                    ticklabels[0] = sci_format(minmaxrz[0])
                    ticklabels[nlog] = sci_format(-linthresh)
                    ticklabels[nticks//2] = r'$0$'
                    ticklabels[nlog + nlin - 1] = sci_format(linthresh)
                    ticklabels[nticks - 1] = sci_format(minmaxrz[1])
                    cbar.set_ticklabels(ticklabels)
                else:
                    cbar_label = (r'$\times10^{%i}\ $' %exp) + units
                    cbar.set_ticks([minmaxrz[0], 0, minmaxrz[1]])
                    cbar.set_ticklabels(['%1.1f' %minmaxrz[0], '0', '%1.1f'\
                            %minmaxrz[1]])
        
                # Put the units (and possibly the exponent) to left of colorbar
                fig.text(cbax_left - 0.3*cbax_width, cbax_center_y,\
                        cbar_label, ha='right', va='center', rotation=90,\
                        fontsize=cbar_fs)

    # Plot contours in time-yy space, if desired
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
            ax.contour(times_2d, yy_2d, field, colors=contour_color,\
                    levels=levels, linewidths=contour_lw)
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
            ax.contour(times_2d_cz, yy_2d_cz, fieldcz, colors=contour_color,\
                    levels=levelscz, linewidths=contour_lw)
            ax.contour(times_2d_rz, yy_2d_rz, fieldrz, colors=contour_color,\
                    levels=levelsrz, linewidths=contour_lw)

    # Plot various times, if desired
    if not plottimes is None:
        for time in plottimes: 
            plt.sca(ax)
            # "dimensional" rval (in dimensions of ro!)
            plt.plot(time + np.zeros(ny), yy, 'k--',\
                    linewidth=0.7*lw_scaling)

    # Plot various yvals, if desired
    # If rbcz has been provided, it should be one of the radii plotted
    if not rbcz is None:
        if yvals is None:
            yvals = []
        yvals.append(rbcz)
    
    if not yvals is None:
        for yval in yvals: 
            plt.sca(ax)
            plt.plot(times, yval + np.zeros_like(times), 'k--',\
                    linewidth=.7*lw_scaling)

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    if showplot:
        plt.show()
    if plotfield:
        return im
