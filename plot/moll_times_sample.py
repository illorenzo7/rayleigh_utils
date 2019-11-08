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
        allthrees_start
from plotcommon import axis_range
from sslice_util import plot_moll
from get_sslice import get_sslice
from rayleigh_diagnostics import Shell_Slices, ReferenceState
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
else:  # By default plot the last 10 Shell_Slices
    index_first, index_last = nfiles - 11, nfiles - 1  

# Other defaults
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
            # to a local radius divided by rsun
varname = 'vr' # by default plot the radial velocity
            
clon = 0.
clat = 20.

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

# We need the radius, which we can get from the reference state
ref = ReferenceState(dirname + '/reference', '')
radius = ref.radius

# Create the plot template
fig_width_inches = 12.

# General parameters for main axis/color bar
margin_bottom_inches = 1.
margin_top_inches = 1.
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
plotdir = dirname + '/plots/moll/times_sample/' + varname + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Loop over each desired iteration and make a plot
for iiter in range(index_first, index_last + 1):
    # Read in desired shell slice
    fname = file_list[iiter]
    a = Shell_Slices(radatadir + fname, '')
    
    # Find desired radius (by default ir=0--near outer surface)
    if not rval is None:
        ir = np.argmin(np.abs(a.radius/rsun - rval))
    rval = a.radius[ir]/rsun # in any case, this is the actual rvalue we get
    
    savename = 'moll_iter' + fname + ('_rval%0.3f_' %(rval/rsun)) +\
            varname  + '.png'
    print('Plotting moll: ' + varname + (', r/rsun = %0.3f (ir = %02i), '\
            %(rval, ir)) + 'iter ' + fname + ' ...')
    vals = get_sslice(a, varname, dirname=dirname)
    field = vals[:, :, ir]
    
    # Make axes and plot the orthographic projection
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width,\
            subplot_height])
    
    plot_moll(field, a.costheta, fig=fig, ax=ax, clon=clon, varname=varname)     

    # Make title
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin
    ax_center_x = ax_xmin + 0.5*ax_delta_x    
    
    varlabel = texlabels[varname]
    prot = 2*np.pi/8.61e-6
    t_loc = a.time[0]/prot - allthrees_start
    title = varlabel + '     ' + (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
            '     ' + ('iter = ' + fname) + '\t' +\
            (r'$t = %.1f\ P_{\rm{rot}}$' %t_loc)
    fig.text(ax_center_x, ax_ymax + 0.02*ax_delta_y, title,\
         verticalalignment='bottom', horizontalalignment='center',\
         fontsize=14, **csfont)   
    
    fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
        strip_dirname(dirname), ha='center', va='bottom', **csfont,\
        fontsize=14)
    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
