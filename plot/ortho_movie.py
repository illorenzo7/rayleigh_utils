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

# Get list of shell slices at all possible times to plot
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get specific range desired for plotting
args = sys.argv[2:]
nargs = len(args)

if '-n' in args or '-range' in args or '-centerrange' in args \
    or '-all' in args or '-iter' in args:        
    index_first, index_last = get_desired_range(int_file_list, args)
else:  # By default plot all available shell_Slices
    index_first, index_last = 0, nfiles - 1  

# Other defaults
minmax = None
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr' # by default plot the radial velocity
clon = 0.
clat = 20.

posdef = False
logscale = False

count = 0

for i in range(nargs):
    arg = args[i]
    if arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-clat':
        clat = float(args[i+1])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-var':
        varname = args[i+1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
        posdef = True
    elif arg == '-posdef':
        posdef = True
    elif arg == '-count':
        count = int(args[i+1])

# The variables with 'sq' in the name are positive-definite
if 'sq' in varname:
    posdef = True

# Create the plot template
fig_width_inches = 4.

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 3./8.
margin_inches = 1./8.

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

# Read in first shell slice for various initial purposes
a0 = Shell_Slices(radatadir + file_list[index_first], '')

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
rval = a0.radius[ir]/rsun # in any case, this is the actual rvalue we get

# Get default minmax value from first shell_slice, if none specified
vals0 = get_sslice(a0, varname, dirname=dirname)
field0 = vals0[:, :, ir]
if minmax is None:
    minmax = get_satvals(field0, posdef, logscale)

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/ortho/movies/' + varname + '/' +\
        ('/rval%0.3f' %rval) + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Calculate the rotation period
om0 = get_parameter(dirname, 'angular_velocity')
prot = 2.*np.pi/om0

# Loop over each desired iteration and make a plot
for iiter in range(index_first, index_last + 1):
    # Read in desired shell slice
    fname = file_list[iiter]
    if iiter > index_first:
        a = Shell_Slices(radatadir + fname, '')
    else:
        a = a0
    time = a.time[0]
    t_prot = time/prot
    
    savename = "img%04i.png" %count
    print('Plotting ortho: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname + (' (count = %04i) ...' %count))
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]
    
    # Make axes and plot the orthographic projection
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])
    
    plot_ortho(field, a.radius, a.costheta, ir=ir, fig=fig, ax=ax,\
            clon=clon, minmax=minmax, varname=varname, posdef=posdef,\
            logscale=logscale, lw_scaling=0.7)     

    # Make title
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    
    
    varlabel = texlabels.get(varname, varname)
    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$t = %.1f\ P_{\rm{rot}}$' %t_prot)    
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=10, **csfont)   
    
    plt.savefig(plotdir + savename, dpi=150)
    plt.close()
    count += 1
