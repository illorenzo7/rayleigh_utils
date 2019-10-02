import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun, reverse_dict,\
        get_satvals, saturate_array
from plotcommon import axis_range, default_axes_1by1
from translate_times import translate_times
from rayleigh_diagnostics import Equatorial_Slices
from varprops import texlabels, var_indices, texunits
from eqslice_util import plot_eqslice

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Equatorial_Slices/'
plotdir = dirname + '/plots/eqslice/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
posdef = False
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
varname = 'all' # by default plot all the fluid variables
plotlonlines = True

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-symlog':
        symlog = True
    elif arg == '-log':
        logscale = True
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
    elif arg == '-nolon':
        plotlonlines = False

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired equatorial slice
print ("Reading " + radatadir + fname + '...')
eq = Equatorial_Slices(radatadir + fname, '')
rr = eq.radius
phi = eq.phi

var_indices_rev = reverse_dict(var_indices)

if varname == 'all': # plot everything
    index_list = eq.qv
else:
    index_list = [var_indices[varname]]

# Figure dimensions
width = 4.
subplot_width_inches = width
subplot_height_inches = width
margin_inches = 1./8.
margin_bottom_inches = 3./4.
margin_top_inches = 3./4.

fig_width_inches = subplot_width_inches + 2.*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

margin_x = margin_inches/fig_width_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches


for var_index in index_list:
    varname = var_indices_rev[var_index]
    
    field = eq.vals[:, :, eq.lut[var_index], 0]
    if var_index in [1, 2, 3]:
        field /= 100.
    elif var_index in [501, 502]: # take out azavg
        field -= np.mean(field, axis=0).reshape((1, eq.nr))
        varname += '_prime'

    texlabel = texlabels[varname]

    # Get default bounds if not specified
    if minmax is None:
        # Cut away the data near the domain boundaries
        minmax = get_satvals(field, posdef=posdef,\
            logscale=logscale, symlog=symlog)

    savename = 'eqslice_' + varname + '_iter' + fname + '.png'
    print('Plotting eqslice: ' + varname + ', iter ' + fname + ' ...')

    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))

    plot_eqslice(field, rr, phi, fig=fig, ax=ax, units=texunits[varname],\
            symlog=symlog, plotlonlines=plotlonlines)

    # Make title
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    

    title = texlabel + '     ' + ('iter = ' + fname)    

    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         va='bottom', ha='center',\
         fontsize=14, **csfont)   
    
    fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
        strip_dirname(dirname), ha='center', va='bottom', **csfont,\
        fontsize=14)
    plt.savefig(plotdir + savename, dpi=300)

    if len(index_list) == 1:
        plt.show()
    plt.close()
