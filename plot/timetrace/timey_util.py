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

# plot_timey needs my_contourf args, then some
plot_timey_kwargs_default = dict({'ycut': None, 'xminmax': None, 'xmin': None, 'xmax': None, 'minmax2': None, 'timevals': np.array([]), 'yvals': np.array([]), 'navg': None, 'plotboundary': True, 'linestyles1': np.array(['-']), 'linewidths1': np.array([default_lw]), 'linecolors1': np.array(['k']),\
       'linestyles2': np.array(['-']), 'linewidths2': np.array([default_lw]), 'linecolors2': np.array(['k'])})

# need to change a few default my_contourf settings
my_contourf_kwargs_default.update(dict({'plotcontours': False, 'cbar_pos': 'right', 'allticksoff': False}))
plot_timey_kwargs_default.update(my_contourf_kwargs_default)

# plot time "lat or rad"
def plot_timey(field, times, yy, fig, ax, **kwargs):
    find_bad_keys(plot_timey_kwargs_default, kwargs, 'plot_timey')
    kw = update_dict(plot_timey_kwargs_default, kwargs)

    kw_my_contourf = update_dict(my_contourf_kwargs_default, kwargs)

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

    # set xminmax if not set by user
    if kw.xminmax is None:
        # set xmin possibly
        if kw.xmin is None:
            kw.xmin = np.min(times)
        # set xmax possibly
        if kw.xmax is None:
            kw.xmax = np.max(times)
        kw.xminmax = kw.xmin, kw.xmax

    ixmin = np.argmin(np.abs(times - kw.xminmax[0]))
    ixmax = np.argmin(np.abs(times - kw.xminmax[1]))

    # Now shorten all the "x" arrays
    times = times[ixmin:ixmax + 1]
    field_full = field_full[ixmin:ixmax + 1]

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
    my_contourf(times_2d1, yy_2d1, field1, fig, ax, **kw_my_contourf)

    if not kw.ycut is None:
        # plot second field first
        kw_my_contourf.minmax = kw.minmax2
        if kw.posdef: 
            kw_my_contourf.cmap = 'cividis'
        else:
            kw_my_contourf.cmap = 'PuOr_r'    
        kw_my_contourf.cbar_no = 2

        # will need to change some kwargs:
        my_contourf(times_2d2, yy_2d2, field2, fig, ax, **kw_my_contourf)

    # potentially plot coordinate lines
    for ind in [1, 2]:
        if ind == 1:
            vals = make_array(kw.timevals, tolist=True)
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
        for i in range(len(vals)):
            val = vals[i]
            if ind == 1:
                xline, yline = val + np.zeros(100), yy
            if ind == 2:
                xline, yline = times, val + np.zeros(100)
            ax.plot(xline, yline, linewidth=linewidths[i], linestyle=linestyles[i], color=linecolors[i])

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
    plt.xticks(fontsize=kw.fontsize)
    plt.yticks(fontsize=kw.fontsize)
