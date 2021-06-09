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

# plot time "lat or rad"
def plot_tlr(field, times, yy, fig, ax, ycut=None, tminmax=None, tmin=None, tmax=None, minmax2=None, timevals=np.array([]), yvals=np.array([]), navg=None, **kwargs_supplied):
    # **kwargs_supplied corresponds to my_contourf
    kwargs_default = {**kwargs_contourf}
    kwargs = dotdict(update_kwargs(kwargs_supplied, kwargs_default))

    # Work with copy of field (not actual field)
    field_full = np.copy(field)
    times = np.copy(times)

    # Get dimensions of field
    ntimes = len(times)
    ny = len(yy)

    # First, average field in time possibly
    if not navg is None:
        over2 = navg//2
        field_timeavg = np.zeros((ntimes - navg + 1, ny))
        for i in range(navg):
            field_timeavg += field_full[i:ntimes - navg + 1 + i]
        field_full = field_timeavg/navg
        times = times[over2:ntimes - over2]

    # set tminmax if not set by user
    if tminmax is None:
        # set xmin possibly
        if tmin is None:
            tmin = np.min(times)
        # set tmax possibly
        if tmax is None:
            tmax = np.max(times)
        tminmax = tmin, tmax

    itmin = np.argmin(np.abs(times - tminmax[0]))
    itmax = np.argmin(np.abs(times - tminmax[1]))

    # Now shorten all the "x" arrays
    times = times[itmin:itmax + 1]
    field = field[itmin:itmax + 1]

    # Make 2D grids from times/yy
    times_2d_full, yy_2d_full = np.meshgrid(times, yy, indexing='ij')

    if ycut is None: # just plotting 1 domain
        kwargs1 = dict(kwargs)
        kwargs1['plotboundary'] = False 
        field1 = field_full
        times_2d1 = times_2d_full
        yy_2d1 = yy_2d_full
    else:
        iycut = np.argmin(np.abs(yy - ycut))
        field1 = field[:, :iycut+1]
        field2 = field[:, iycut+1:]
        times_2d1 = times_2d[:, :iycut+1]
        times_2d2 = times_2d[:, iycut+1:]
        yy_2d1 = yy_2d[:, :iycut+1]
        yy_2d2 = yy_2d[:, iycut+1:]

        # plot second field first
        # will need to change some kwargs:
        kwargs2 = dict(kwargs)
        kwargs2['plotboundary'] = False 
        kwargs2['minmax'] = minmax2
        if kwargs.posdef:
            kwargs2['cmap'] = 'cividis'
        else:
            kwargs2['cmap'] = 'PuOr_r'    
        kwargs2['cbar_no'] = 2
        my_contourf(times_2d2, yy_2d2, field2, fig, ax, **kwargs2)

    # regardless, plot the first field
    my_contourf(times_2d1, yy_2d1, field1, fig, ax, **kwargs1)

    # potentially plot coordinate lines
    my_contourf(times_2d_full, yy_2d_full, field_full, fig, ax, plotfield=False, plotcontours=False, func1=times_2d_full, vals1=timevals, func2=yy_2d_full, vals2=yvals, plotboundary=kwargs['plotboundary'])

    # Get ticks everywhere
    plt.sca(ax)
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')
