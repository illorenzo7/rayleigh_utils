# Created 04/24/2020
# Plots and saves equatorial slice associated with a given variable
# in a time series
# Similar to moll_times_sample
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
from eqslice_util import plot_eqslice
from rayleigh_diagnostics import Equatorial_Slices
from get_eqslice import get_eqslice
from varprops import texlabels, texunits

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Equatorial_Slices
radatadir = dirname + '/Equatorial_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
symlog = False
logscale = False
varname = 'vr' # by default plot the radial velocity
posdef = False
plotcontours = True
plotlonlines = True
nlevs = None

plotdir = None

args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple

for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-var' or arg == '-qval':
        varname = args[i+1]
    elif arg == '-log':
        logscale = True
    elif arg == '-symlog':
        symlog = True
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-nolon':
        plotlonlines = True
    elif arg == '-nlevs':
        nlevs = int(args[i+1])

# See if posdef should be true
if 'sq' in varname:
    posdef = True

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

# Get tex units and label
units = texunits.get(varname, 'cgs')
texlabel = texlabels.get(varname, 'qval = ' + varname)

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/eq/times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over each desired iteration and make plots
for i in range(index_first, index_last + 1):
    # Read in desired equatorial slice
    fname = file_list[i]
    print ("read Equatorial_Slices/" + fname)
    eq = Equatorial_Slices(radatadir + fname, '')

    # Loop over the slices in the file 
    for j in range(eq.niter):
        # Get field 
        field = get_eqslice(eq, varname, dirname=dirname, j=j)
        # Get local time (in seconds)
        t_loc = eq.time[j]
        iter_loc = eq.iters[j]

        # Make name to save the plot
        savename = 'eq_' + varname + '_iter' + str(iter_loc).zfill(8) +\
                '.png'

        # Display at terminal what we are plotting
        print('Plotting eq: ' + varname + ', iter ' + fname)

        # create axes
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = fig.add_axes((margin_x, margin_bottom, subplot_width,\
                subplot_height))
        plot_eqslice (field, eq.radius, eq.phi, fig=fig, ax=ax,\
                units=units,\
                minmax=minmax, plotlonlines=plotlonlines, posdef=posdef,\
                symlog=symlog, logscale=logscale, plotcontours=plotcontours,\
                nlevs=nlevs)

        # Make title
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
        plt.savefig(plotdir + savename, dpi=300)
        fig.text(margin_x, 1 - margin_y - 3*line_height, texlabel,\
                ha='left', va='top', fontsize=fsize, **csfont)
        plt.close() 
