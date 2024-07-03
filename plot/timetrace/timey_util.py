########################################################################
# This file contains utilities useful for plotting time-yy traces
# of the AZ_Avgs data in Rayleigh  (usually time-latitude or time-radius)
# Author: Loren Matilsky
# Created: 04/06/2020

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from common import *
from plotcommon import *

# default timey fig dimensions
timey_fig_dimensions = dict({'sub_width_inches': 7.5, 'sub_height_inches': 2.0, 'margin_top_inches': 1.25, 'sub_margin_left_inches': 5/8, 'sub_margin_right_inches': 7/8})

# plot_timey needs my_contourf/my_pcolormesh args, then some
plot_timey_kwargs_default = dict({'ycut': None, 'yminmax': None, 'ymin': None, 'ymax': None,'minmax2': None, 'xvals': np.array([]), 'yvals': np.array([]), 'navg': None, 'plotboundary': True, 'linestyles1': np.array(['-']), 'linewidths1': np.array([default_lw]), 'linecolors1': np.array(['k']),\
       'linestyles2': np.array(['-']), 'linewidths2': np.array([default_lw]), 'linecolors2': np.array(['k']),\
       'pcolormesh': False})

# need to change a few default my_contourf settings,
# but don't change default kwargs dictionary directly
tmp1 = my_contourf_kwargs_default.copy()
tmp2 = my_pcolormesh_kwargs_default.copy()
tochange = {'plotcontours': False, 'cbar_pos': 'right', 'allticksoff': False}
tmp1.update(tochange)
tmp2.update(tochange)
plot_timey_kwargs_default.update(tmp1)
plot_timey_kwargs_default.update(tmp2)

# plot time "lat or rad"
def plot_timey(field, times, yy, fig, ax, **kwargs):

    find_bad_keys(plot_timey_kwargs_default, kwargs, 'plot_timey')
    kw = update_dict(plot_timey_kwargs_default, kwargs)

    if kw.pcolormesh:
        plotting_func = my_pcolormesh
        kw_plotting_func = update_dict(tmp2, kwargs)
    else:        
        plotting_func = my_contourf
        kw_plotting_func = update_dict(tmp1, kwargs)

    # Work with copy of field (not actual field)
    field_full = np.copy(field)
    times = np.copy(times)

    # Get dimensions of field
    ntimes = len(times)
    ny = len(yy)

    # First, average field in time possibly
    if not kw.navg is None:
        over2 = kw.navg//2
        field_timeavg = np.zeros((ntimes - kw.navg + 1, ny))
        for i in range(kw.navg):
            field_timeavg += field_full[i:ntimes - kw.navg + 1 + i]
        field_full = field_timeavg/kw.navg
        times = times[over2:ntimes - over2]

    # set yminmax if not set by user
    if kw.yminmax is None:
        # set ymin possibly
        if kw.ymin is None:
            kw.ymin = np.min(yy)
        # set ymax possibly
        if kw.ymax is None:
            kw.ymax = np.max(yy)
        kw.yminmax = kw.ymin, kw.ymax

    iymin = np.argmin(np.abs(yy - kw.yminmax[0]))
    iymax = np.argmin(np.abs(yy - kw.yminmax[1]))

    # Make 2D grids from times/yy
    times_2d_full, yy_2d_full = np.meshgrid(times, yy, indexing='ij')

    if kw.ycut is None: # just plotting 1 domain
        field1 = field_full
        times_2d1 = times_2d_full
        yy_2d1 = yy_2d_full
    else:
        iycut = np.argmin(np.abs(yy - kw.ycut))
        field1 = field_full[:, :iycut+1]
        field2 = field_full[:, iycut+1:]
        times_2d1 = times_2d_full[:, :iycut+1]
        times_2d2 = times_2d_full[:, iycut+1:]
        yy_2d1 = yy_2d_full[:, :iycut+1]
        yy_2d2 = yy_2d_full[:, iycut+1:]

    # plot the first field
    if kw.pcolormesh:
        kw_plotting_func.x = times_2d1[:, 0]
        kw_plotting_func.y = yy_2d1[0, :]
        plotting_func(field1, fig, ax, **kw_plotting_func)
    else:
        plotting_func(times_2d1, yy_2d1, field1, fig, ax, **kw_plotting_func)

    if not kw.ycut is None:
        kw_plotting_func.minmax = kw.minmax2
        if kw.posdef: 
            kw_plotting_func.cmap = 'cividis'
        else:
            kw_plotting_func.cmap = 'PuOr_r'    
        kw_plotting_func.cbar_no = 2
        if kw.pcolormesh:
            kw_plotting_func.x = times_2d2[:, 0]
            kw_plotting_func.y = yy_2d2[0, :]
            plotting_func(field2, fig, ax, **kw_plotting_func)
            #ax.set_ylim(np.min(yy), np.max(yy))
        else:
            plotting_func(times_2d2, yy_2d2, field2, fig, ax, **kw_plotting_func)

    # potentially plot coordinate lines
    for ind in [1, 2]:
        if ind == 1:
            vals = make_array(kw.xvals, tolist=True)
            linecolors = kw.linecolors1
            linestyles = kw.linestyles1
            linewidths = kw.linewidths1
        if ind == 2:
            vals = make_array(kw.yvals, tolist=True)
            linecolors = kw.linecolors2
            linestyles = kw.linestyles2
            linewidths = kw.linewidths2

        # make lists from everything that needs to be
        linecolors = make_array(linecolors, tolist=True, length=len(vals))
        linestyles = make_array(linestyles, tolist=True, length=len(vals))
        linewidths = make_array(linewidths, tolist=True, length=len(vals))

        npoints = 100
        xaxis = np.linspace(times[0], times[-1], npoints)
        yaxis = np.linspace(yy[0], yy[-1], npoints)
        for i in range(len(vals)):
            val = vals[i]
            if ind == 1:
                xline, yline = val + np.zeros(npoints), yaxis
            if ind == 2:
                xline, yline = xaxis, val + np.zeros(npoints)
            ax.plot(xline, yline, linewidth=linewidths[i], linestyle=linestyles[i], color=linecolors[i])

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    plt.xticks(fontsize=kw.fontsize)
    plt.yticks(fontsize=kw.fontsize)
