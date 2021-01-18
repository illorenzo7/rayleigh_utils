import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import *
from plotcommon import axis_range
from sslice_util import plot_ortho
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
logscale = False
symlog = False
iiter = nfiles - 1 # by default plot the last iteration
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr' # by default plot the radial velocity
clon = 0
clat = 20
bold_patch = None
thickcenter = True

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-ir'):
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif (arg == '-var'):
        varname = args[i+1]
    elif (arg == '-iter'):
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
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-log':
        logscale = True
    elif arg == '-symlog':
        symlog = True
    elif arg == '-patch':
        bold_patch = np.float(args[i+1]), np.float(args[i+2])
    elif arg == '-nothick':
        thickcenter = False

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')
vals = get_sslice(a, varname, dirname=dirname)

# We also need the radius, which we can get from the reference state
eq = get_eq(dirname)
radius = eq.radius

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
field = vals[:, :, ir]
rval = radius[ir] # in any case, this is the actual rvalue we get

# Create the plot using subplot axes
# Offset axes slightly (at the end) to deal with annoying white space cutoff
fig_width_inches = 6.

# General parameters for main axis/color bar
margin_bottom_inches = 1.
margin_top_inches = 1.
margin_inches = 1/8

subplot_width_inches = fig_width_inches - 2*margin_inches
subplot_height_inches = subplot_width_inches
fig_height_inches = margin_bottom_inches + subplot_height_inches +\
    margin_top_inches
fig_aspect = fig_height_inches/fig_width_inches

# "Non-dimensional" figure parameters
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Main axis plot: orthographic projection
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes([margin_x, margin_bottom, subplot_width, subplot_height])

plot_ortho(field, radius, a.costheta, fig=fig, ax=ax, ir=a.inds[ir],\
        minmax=minmax, clon=clon, clat=clat, varname=varname,\
        logscale=logscale, symlog=symlog, bold_patch=bold_patch,\
        thickcenter=thickcenter) 

# Make title
ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
ax_delta_x = ax_xmax - ax_xmin
ax_delta_y = ax_ymax - ax_ymin
ax_center_x = ax_xmin + 0.5*ax_delta_x

varlabel = texlabels.get(varname, 'qval = ' + varname)
title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %(rval/rsun)) +\
            '     ' + ('iter = ' + fname)
fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
     verticalalignment='bottom', horizontalalignment='center',\
     fontsize=14, **csfont)   

fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
    strip_dirname(dirname), ha='center', va='bottom', **csfont,\
    fontsize=14)

plt.show()
