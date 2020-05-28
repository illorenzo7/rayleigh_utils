import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun, get_desired_range,\
        get_satvals
from plotcommon import axis_range
from sslice_util import plot_moll
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices
from varprops import texlabels
from get_parameter import get_parameter

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
minmax_vr = None
minmax_bp = None
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
clon = 0.
skip = 1 # by default plot all files the desired plotting range
    # set skip > 1 to only plot every skip^th file

for i in range(nargs):
    arg = args[i]
    if arg == '-clon':
        clon = float(args[i+1])
    elif arg == '-ir':
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-minmax':
        minmax_vr = float(args[i+1]), float(args[i+2])
        minmax_bp = float(args[i+3]), float(args[i+4])
    elif arg == '-skip':
        skip = int(args[i+1])

# Create the plot template
fig_width_inches = 12.

# General parameters for main axis/color bar
margin_bottom_inches = 1./2.
margin_top_inches = 3./8.
margin_inches = 1./8.

subplot_width_inches = 0.5*(fig_width_inches - 3.*margin_inches)
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

# Read in first shell slice for various initial purposes
a0 = Shell_Slices(radatadir + file_list[index_first], '')

# Find desired radius (by default ir=0--near outer surface)
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
rval = a0.radius[ir]/rsun # in any case, this is the actual rvalue we get

# Get default minmax value from first shell_slice, if none specified
vals_vr0 = get_sslice(a0, 'vr', dirname=dirname)/100.
# convert cm/s --> m/s
vals_bp0 = get_sslice(a0, 'bp', dirname=dirname)
vr0 = vals_vr0[:, :, ir]
bp0 = vals_bp0[:, :, ir]

if minmax_vr is None or minmax_bp is None:
    minmax_vr = get_satvals(vr0, False, False)
    minmax_bp = get_satvals(bp0, False, False)
    # logscale and posef should be "false" for vr and bp

# Make plot (sub-)directory if it doesn't already exist
plotdir = dirname + '/plots/moll/movies/vr_and_bp/' +\
        ('/rval%0.3f' %rval) + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Calculate the rotation period
om0 = get_parameter(dirname, 'angular_velocity')
prot = 2.*np.pi/om0

# Loop over each desired iteration and make a plot
count = 0
for iiter in range(index_first, index_last + 1, skip):
    # Read in desired shell slice
    fname = file_list[iiter]
    if iiter > index_first:
        a = Shell_Slices(radatadir + fname, '')
    else:
        a = a0
    time = a.time[0]
    t_prot = time/prot
    
    savename = "img%04i.png" %count
    print(('Plotting moll: vr_and_bp,  r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname + (' (count = %04i) ...' %count))
    vals_vr = get_sslice(a, 'vr', dirname=dirname)
    vals_bp = get_sslice(a, 'bp', dirname=dirname)
    field_vr = vals_vr[:, :, ir]/100. # change cm/s --> m/s
    field_bp = vals_bp[:, :, ir]
   
    # Make axes and plot the Mollweide projections
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax_left = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])
    ax_right = fig.add_axes([2.*margin_x + subplot_width, margin_bottom,\
            subplot_width, subplot_height])
   
    plot_moll(field_vr, a.costheta, fig=fig, ax=ax_left, clon=clon,\
            minmax=minmax_vr, varname='vr', lw_scaling=0.7)
    plot_moll(field_bp, a.costheta, fig=fig, ax=ax_right, clon=clon,\
            minmax=minmax_bp, varname='bp', lw_scaling=0.7)

    # Make title
    ax_left_xmin, ax_left_xmax, ax_left_ymin, ax_left_ymax =\
            axis_range(ax_left)
    ax_right_xmin, ax_right_xmax, ax_right_ymin, ax_right_ymax =\
            axis_range(ax_right)
    x_between = 0.5*(ax_left_xmax + ax_right_xmin)
    ax_ymax = ax_left_ymax
    ax_delta_y = ax_left_ymax - ax_left_ymin
    
    varlabel_vr = texlabels['vr']
    varlabel_bp = texlabels['bp']

    title = varlabel_vr + ' and ' + varlabel_bp + '     ' + \
            (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + (r'$t = %.1f\ P_{\rm{rot}}$' %t_prot)    
    fig.text(x_between, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=14, **csfont)   
    
    plt.savefig(plotdir + savename, dpi=150)
    plt.close()
    count += 1
