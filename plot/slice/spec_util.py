# Author: Loren Matilsky
# Updated: 07/09/2021

from matplotlib import ticker, colors
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from transform_coordinates import mollweide_transform

plot_spec_lm_kwargs_default = dict({'lvals': None, 'mvals': None, 'linewidth': default_lw,\
    'minmax': None, 'lminmax', 'mminmax',\
    # colorbar stuff
    'plotcbar': True, 'cmap': None, 'norm': None, 'linear': False, 'units': '', 'fontsize': default_labelsize})
plot_spec_lm_kwargs_default.update(add_cbar_kwargs_default)

def plot_spec_lm(field, fig, ax, **kwargs):
    kw = update_dict(my_contourf_kwargs_default, kwargs)
    kw_add_cbar = update_dict(add_cbar_kwargs_default, kwargs)
    find_bad_keys(plot_spec_lm_kwargs_default, kwargs, 'plot_spec_lm')

    # make sure Python does not modify any of the arrays it was passed
    field = np.copy(field)

    # full (l, m) grid:
    nell, nm = np.shape(field)
    lvals_all = np.arange(nell)
    mvals_all = np.arange(nm)

    if not kw.lminmax is None:
        il1 = np.argmin(np.abs(lvals_all - kw.lminmax[0]))
        il2 = np.argmin(np.abs(lvals_all - kw.lminmax[1]))
    else:
        il1, il2 = 0, nell - 1

    if not kw.mminmax is None:
        im1 = np.argmin(np.abs(mvals_all - kw.mminmax[0]))
        im2 = np.argmin(np.abs(mvals_all - kw.mminmax[1]))
    else:
        im1, im2 = 0, nm - 1

    # now adjust everything by the (l, m) range we want
    lvals = lvals[il1:il2+1]
    mvals = mvals[im1:im2+1]
    field = field[il1:il2+1, im1:im2+1]

    lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')
    lvals_2d, mvals_2d = xy_grid(lvals_2d, mvals_2d)

    # Get minmax, if not specified
    if kw.minmax is None:
        power_not0 = np.copy(field)
        # power gets wierd (close to 0?) at the two
        # most extreme l-values
        if il1 == 0 and il2 == nell - 1: 
            field_not0 = np.copy(field_loc[1:-1, :])
        if il1 == 0: 
            power_not0 = power_not0[1:, :]
        if il2 == nell - 1: 
            power_not0 = power_not0[:-1, :]
        power_not0 = power_not0[power_not0 != 0.]
        if kw.linear: # NOT the default...
            kw.minmax = kw_add_cbar.minmax =\
                contourf_minmax(power_not0, posdef=True)
        else:
            kw.minmax = kw_add_cbar.minmax =\
                contourf_minmax(power_not0, logscale=True)
   
    if kw.norm is None and not kw.linear: # the default
        kw.norm = colors.LogNorm(vmin=kw.minmax[0], vmax=kw.minmax[1])
    
    if kw.cmap is None:
        kw.cmap = 'jet'
    im = plt.pcolormesh(lvals_2d_new, mvals_2d_new, field, cmap=kw.cmap, norm=kw.norm)  

    # now deal with color bar, if one is desired
    if kw.plotcbar:
        add_cbar(fig, ax, im, **kw_add_cbar)

