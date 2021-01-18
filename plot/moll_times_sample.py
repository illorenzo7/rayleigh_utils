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
from sslice_util import plot_moll
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices, GridInfo
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Data with Shell_Slices
radatadir = dirname + '/Shell_Slices/'

# Get list of shell slices at all possible times to plot
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get specific range desired for plotting
args = sys.argv[2:]
nargs = len(args)

the_tuple = get_desired_range(int_file_list, args)
if the_tuple is None: # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  
else:
    index_first, index_last = the_tuple
nfiles = index_last - index_first + 1 # this is the number of files we 
# will plot (modulo the nskip or ntot filters)

# Other defaults
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr' # by default plot the radial velocity
clon = 0.
minmax = None
nskip = 1 # by default don't skip any slices in the range

# Change the defaults using the CLAs
for i in range(nargs):
    arg = args[i]
    if arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nskip':
        nskip = int(args[i+1])
    elif arg == '-ntot':
        ntot = int(args[i+1])
        nskip = nfiles//ntot

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# We need the radius, which we can get from grid_info
gi = GridInfo(dirname + '/grid_info', '')
radius = gi.radius

# Create the plot template
fig_width_inches = 6.

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 5./8.
margin_inches = 1./8.

subplot_width_inches = fig_width_inches - 2*margin_inches
subplot_height_inches = 0.5*subplot_width_inches
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

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/moll/times_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over each desired iteration and make plots
for i in range(index_first, index_last + 1, nskip):
    # Read in desired shell slice
    fname = file_list[i]
    a = Shell_Slices(radatadir + fname, '')

    # Find desired radius (by default ir=0--near outer surface)
    if not rval is None:
        ir = np.argmin(np.abs(a.radius/rsun - rval))
    rval = a.radius[ir]/rsun 
    # in any case, this is the actual rvalue we get

    # Loop over the slices in the file 
    for j in range(a.niter):
        # Get local time (in seconds)
        t_loc = a.time[j]
        iter_loc = a.iters[j]

        # Savename
        savename = 'moll_' + varname + ('_rval%0.3f' %rval) + '_iter' +\
                str(iter_loc).zfill(8) + '.png'
        print('Plotting: ' + savename)
        vals = get_sslice(a, varname, dirname=dirname, j=j)
        field = vals[:, :, ir]
        
        # Make axes and plot the Mollweide projection
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
                subplot_height])
        
        plot_moll(field, a.costheta, fig=fig, ax=ax, clon=clon,\
                varname=varname, minmax=minmax)     

        # Make title
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
        varlabel = texlabels.get(varname, varname)

        title = dirname_stripped +\
            '\n' + r'$\rm{Mollweide}$' + '     '  + time_string +\
            '\n' + varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval)
        fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
             verticalalignment='bottom', horizontalalignment='center',\
             fontsize=10, **csfont)   
        
        plt.savefig(plotdir + savename, dpi=300)
        plt.close()
