import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname
from sslice_util import plot_ortho
from rayleigh_diagnostics import Shell_Slices
from get_parameter import get_parameter

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)

minmax = None
iiter = nfiles - 1 # by default plot the last iteration
idepth = 0 # by default plot just below the surface
clon = 0
clat = 20

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-d'):
        idepth = int(args[i+1])
    elif (arg == '-iter'):
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif (arg == '-clat'):
        clat = float(args[i+1])
    elif (arg == '-clon'):
        clon = float(args[i+1])

iter_val = int_file_list[iiter]
fname = file_list[iiter]

# Read in desired shell slice
a = Shell_Slices(radatadir + fname, '')

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

# Loop over depths and make plots
plotdir = dirname + '/plots/ortho/fluidvars_sample/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    varlist = ['vr', 'vt', 'vp', 's_prime', 'p_prime', 's_prime_sph', 'br', 'bt', 'bp' ]
else:
    varlist = ['vr', 'vt', 'vp', 's_prime', 'p_prime', 's_prime_sph']

count = 0
for varname in varlist: 
    savename = 'ortho_iter' + fname + ('_depth%02i_' %idepth) +\
            str(count).zfill(2) + '-' + varname  + '.png'
    count += 1
    print('Plotting ortho: ' + varname + (', depth %02i, ' %idepth) +\
          'iter ' + fname + ' ...')
    fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
    ax = fig.add_axes([margin_x, margin_bottom, subplot_width, subplot_height])
    
    plot_ortho(fig, ax, a, dirname, varname, idepth=idepth, minmax=minmax,\
                clon=clon, clat=clat) 
    fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
            strip_dirname(dirname), ha='center', va='bottom', **csfont, fontsize=14)
    plt.savefig(plotdir + savename, dpi=300)
    plt.close()
