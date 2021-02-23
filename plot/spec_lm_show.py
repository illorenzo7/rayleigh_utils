# Author: Loren Matilsky
# Created: 02/26/2020
# This script shows powerspectrum (amplitude-squared power) vs spherical-
# harmonic modes l,m for all l, m at 
# a given radius, from the Shell_Spectra data in
# the Rayleigh run directory indicated by [dirname]. 
# does not save plot

import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import ticker, colors
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from varprops import texunits, texlabels, var_indices
from common import *
from plotcommon import axis_range, xy_grid
from rayleigh_diagnostics import Shell_Spectra

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Spectra/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
iiter = nfiles - 1 # by default plot the last iteration
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr'
minmax = None
lminmax = None
mminmax = None
fs = 12.
tag = ''

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-var':
        varname = args[i+1]
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-lminmax':
        lminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-mminmax':
        mminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-iter':
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-sec':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='sec')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-day':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='day')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-prot':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='prot')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# File name to read
iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell spectrum
spec = Shell_Spectra(radatadir + fname, '')
rvals = spec.radius/rsun
qv = spec.qv
lut = spec.lut
nr = spec.nr
nell = spec.nell
nm = spec.nm
lvals = np.arange(nell, dtype='float')
mvals = np.arange(nm, dtype='float')

# Get local time (in seconds)
t_loc = spec.time[0]

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(rvals - rval))
rval = rvals[ir] # in any case, this is the actual rvalue we get

if not lminmax is None:
    il1 = np.argmin(np.abs(lvals - lminmax[0]))
    il2 = np.argmin(np.abs(lvals - lminmax[1]))
else:
    il1, il2 = 0, nell - 1

if not mminmax is None:
    im1 = np.argmin(np.abs(mvals - mminmax[0]))
    im2 = np.argmin(np.abs(mvals - mminmax[1]))
else:
    im1, im2 = 0, nm - 1

# This will set the aspect ratio for the plot
nl_used = il2 - il1 + 1
nm_used = im2 - im1 + 1

# Create the plot using subplot axes
if nl_used >= nm_used:
    subplot_width_inches = 5.
    subplot_height_inches = nm_used/nl_used*subplot_width_inches
else:
    subplot_height_inches = 5.
    subplot_width_inches = nl_used/nm_used*subplot_height_inches

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_left_inches = 5./8.
margin_right_inches = 1.
margin_top_inches = 1.
margin_inches = 1./8.

fig_width_inches = subplot_width_inches + margin_left_inches +\
        margin_right_inches
fig_height_inches = margin_bottom_inches + subplot_height_inches +\
    margin_top_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_left = margin_left_inches/fig_width_inches
margin_right = margin_right_inches/fig_width_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig_aspect = fig_height_inches/fig_width_inches

lvals = lvals[il1:il2+1]
mvals = mvals[im1:im2+1]

lvals_2d, mvals_2d = np.meshgrid(lvals, mvals, indexing='ij')
lvals_2d_new, mvals_2d_new = xy_grid(lvals_2d, mvals_2d)

# Get full power
fullpower = np.abs(spec.vals[:, :, :, :, 0])**2


# Get power associated with desired quantity 
# (should really have a contigency
# where the routine exits if the desired quantity isn't present
if not (varname == 'vtot' or varname == 'btot'):
    desired_qv = var_indices[varname]
    iq = np.argmin(np.abs(qv - desired_qv))
    varlabel = texlabels.get(varname, 'qval = ' + varname)
else:
    if varname == 'vtot':
        desired_qv_vals = [1, 2, 3]
        varlabel = r'$|\mathbf{v}|$'
    elif varname == 'btot':
        desired_qv_vals = [801, 802, 803]
        varlabel = r'$|\mathbf{B}|$'
    iq_vals = []
    for desired_qv in desired_qv_vals:
        iq_vals.append(np.argmin(np.abs(qv - desired_qv)))

# Get power in l,m space for desired variable
if varname == 'vtot' or varname == 'btot':
    power = np.abs(fullpower[:, :, :, iq_vals[0]]) +\
            np.abs(fullpower[:, :, :, iq_vals[1]]) +\
            np.abs(fullpower[:, :, :, iq_vals[2]]) 
else:
    power = np.abs(fullpower[:, :, :, iq])

# Now make plots
print('Plotting spec_lm: ' + varname +\
        (', r/rsun = %0.3f (ir = %02i), ' %(rval, ir)) + 'iter ' +\
        fname + ' ...')
power_loc = power[il1:il2+1, im1:im2+1, ir]

# Get minmax, if not specified
if minmax is None:
    if il1 == 0 and il2 == nell - 1: 
        # power gets wierd (close to 0?) at the two
        # most extreme l-values
        power_not0 = np.copy(power_loc[1:-1, :])
    elif il1 == 0: 
        power_not0 = np.copy(power_loc[1:, :])
    elif il2 == nell - 1: 
        power_not0 = np.copy(power_loc[:-1, :])
    else:
        power_not0 = np.copy(power_loc)
    power_not0 = power_not0[power_not0 != 0.]
    minmax_loc = get_satvals(power_not0, logscale=True)
else:
    minmax_loc = minmax

# Make plot
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes([margin_left, margin_bottom, subplot_width,\
        subplot_height])
plt.sca(ax)
im = plt.pcolormesh(lvals_2d_new, mvals_2d_new, power_loc,\
        cmap='jet',\
    norm=colors.LogNorm(vmin=minmax_loc[0], vmax=minmax_loc[1]))  

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
    time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
            ' (1 ' + time_label + (' = %.2f days)'\
            %(time_unit/86400.))
else:
    time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
            ' (1 ' + time_label + (' = %.1f days)'\
            %(time_unit/86400.))

# Make title
title = dirname_stripped +\
    '\n' + r'$\rm{spec\_lm}$' + '     '  + time_string +\
    '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
    '\n' + (r'$\ell_{\rm{rms}} = %.1f$' %l_rms) + '     ' +\
    (r'$m_{\rm{rms}} = %.1f$' %m_rms)
fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
     verticalalignment='bottom', horizontalalignment='center',\
     fontsize=fs, **csfont)   

plt.show()
