# Created 04/24/2020
# Plots and saves equatorial slices associated for a range of variables at 
# a particular time
# Similar to moll_vars_sample
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
from translate_times import translate_times
from eqslice_util import plot_eqslice
from rayleigh_diagnostics import Equatorial_Slices
from get_eqslice import get_eqslice
from get_parameter import get_parameter
from varprops import texlabels, texunits
from time_scales import compute_Prot, compute_tdt

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Data with Equatorial_Slices
radatadir = dirname + '/Equatorial_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
symlog = False
logscale = False
iiter = nfiles - 1 # by default plot the last iteration
posdef = False
plotcontours = True
plotlonlines = True
varlist = None
nlevs = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-log':
        logscale = True
    elif arg == '-symlog':
        symlog = True
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nolon':
        plotlatlines = True
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
    elif arg == '-vars':
        varlist = args[i+1].split()
    elif arg == '-nlevs':
        nlevs = int(args[i+1])

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Read in desired meridional slice
fname = file_list[iiter]
print ("read Equatorial_Slices/" + fname)
eq = Equatorial_Slices(radatadir + fname, '')

# Get local time (in seconds)
t_loc = eq.time[0]

# Figure dimensions
width = 5.
subplot_width_inches = width
subplot_height_inches = width
margin_inches = 1./8.
margin_bottom_inches = 3./4.
margin_top_inches = 1.

fig_width_inches = subplot_width_inches + 2.*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/eq/vars_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over vars and make plots
if varlist is None:
    varlist = ['vr_prime', 'vt_prime', 'vp_prime', 'omr_prime',\
            'omt_prime', 'omp_prime', 's_prime', 'p_prime']
    magnetism = get_parameter(dirname, 'magnetism')
    if magnetism:
        varlist.extend(['br', 'bt', 'bp'])


for varname in varlist: 
    # Get Tex units and variable label
    units = texunits[varname]
    texlabel = texlabels[varname]

    # Name to save plot
    savename = 'eq_iter' + fname + '_' + varname  + '.png'

    # See if posdef should be true
    if 'sq' in varname:
        posdef = True

    # Get the field
    field = get_eqslice(eq, varname, dirname=dirname)

    # Display at terminal what we are plotting
    print('Plotting eq: ' + varname + ', iter ' + fname)

    # create axes
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
            subplot_height))
    plot_eqslice (field, eq.radius, eq.phi, fig=fig, ax=ax, units=units,\
            minmax=minmax, plotlonlines=plotlonlines, posdef=posdef,\
            symlog=symlog, logscale=logscale, plotcontours=plotcontours,\
            nlevs=nlevs)

    # Make title
    # Get the time label
    if rotation:
        time_string = ('t = %.1f ' %(t_loc/time_unit)) + time_label +\
                '(1 ' + time_label + (' = %.2f days)'\
                %(time_unit/86400.))
    else:
        time_string = ('t = %.3f ' %(t_loc/time_unit)) + time_label +\
                '(1 ' + time_label + (' = %.1f days)'\
                %(time_unit/86400.))

    fsize = 12.
    line_height = 1./4./fig_height_inches
    fig.text(margin_x, 1. - margin_y, dirname_stripped, ha='left',\
            va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1. - margin_y - line_height, 'Equatorial Slice',\
            ha='left', va='top', fontsize=fsize, **csfont)
    fig.text(margin_x, 1. - margin_y - 2*line_height,\
             time_string, ha='left', va='top', fontsize=fsize,\
             **csfont)
    fig.text(margin_x, 1 - margin_y - 3*line_height, texlabel,\
             ha='left', va='top', fontsize=fsize, **csfont)
    plt.savefig(plotdir + savename, dpi=300)
    plt.close() 
