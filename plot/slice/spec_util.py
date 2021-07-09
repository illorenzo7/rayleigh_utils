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
    'plotcbar': True, 'cbar_thick': 1/16, 'cbar_aspect': 1/20, 'cbar_prec': 2, 'cbar_no': 1, 'cbar_pos': 'bottom', 'cmap': None, 'norm': None, 'linear': False, 'units': '', 'nosci': False, 'fontsize': default_labelsize})

def plot_spec_lm(field, fig, ax, **kwargs):
    kw = update_dict(my_contourf_kwargs_default, kwargs)
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
    lvals_2d_new, mvals_2d_new = xy_grid(lvals_2d, mvals_2d)

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
            kw.minmax = contourf_minmax(power_not0, posdef=True)
        else:
            kw.minmax = contourf_minmax(power_not0, logscale=True)
   
    if kw.norm is None and not kw.linear: # the default
        kw.norm = colors.LogNorm(vmin=minmax_loc[0], vmax=minmax_loc[1])
    
    if kw.cmap is None:
        kw.cmap = 'jet'
    im = plt.pcolormesh(lvals_2d_new, mvals_2d_new, field, cmap=kw.cmap, norm=kw.norm)  

    # Set up the colorbar "by hand"
    ax_left, ax_right, ax_bottom, ax_top = axis_range(ax)
    ax_width = ax_right - ax_left
    ax_height = ax_top - ax_bottom
    cbax_center_y = ax_bottom + ax_height/2.
    cbax_left = ax_right + 1/8/fig_width_inches
    cbax_width = 1/8/fig_width_inches
    cbax_aspect = 20.
    cbax_height = cbax_width*cbax_aspect/fig_aspect
    cbax_bottom = cbax_center_y - cbax_height/2. 
    
    cbaxes = fig.add_axes([cbax_left, cbax_bottom,\
                   cbax_width, cbax_height])
    cbar = plt.colorbar(im, cax=cbaxes)

    plt.sca(ax)
    # set bounds
#    plt.xlim(0.5, nell - 0.5)
#    plt.ylim(0.5, nm - 0.5)

    # label axes
    plt.xlabel(r'${\rm{spherical\ harmonic\ degree}}\ \ell$', fontsize=fs)
    plt.ylabel(r'${\rm{azimuthal\ order}}\ m$', fontsize=fs)

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    # Get colorbar
#    cbar = plt.colorbar()

    # Make title
    # Compute l_rms and m_rms
    l_rms = np.sum(power_loc*lvals_2d)/np.sum(power_loc)
    m_rms = np.sum(power_loc*mvals_2d)/np.sum(power_loc)

    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    
    
    if rotation:
        time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.1f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'
    else:
        time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
                + time_label + (r'$\ (\Delta t = %.3f\ $'\
                %((t2 - t1)/time_unit)) + time_label + ')'

    # Make title
    title = dirname_stripped +\
        '\n' + r'$\rm{specav\_lm}$' + '     '  + time_string +\
        '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
        '\n' + (r'$\ell_{\rm{rms}} = %.1f$' %l_rms) + '     ' +\
        (r'$m_{\rm{rms}} = %.1f$' %m_rms)
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=fs, **csfont)   

    plt.savefig(plotdir + savename, dpi=300)
    if showplot:
        plt.show()
    plt.close()
