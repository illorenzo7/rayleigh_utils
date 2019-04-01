# Created: 03/29/2019
# Author: Loren Matilsky
#
# Purpose: create a string Mollweide slices for every shell slice in a
# directory, to blend together in an ffmpeg movie
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
sys.path.append(os.environ['rapp'])
from common import get_file_lists, strip_dirname, rsun
from sslice_util import plot_moll, get_satvals, get_sslice
from rayleigh_diagnostics import Shell_Slices
from get_parameter import get_parameter
from varprops import texlabels

# Get command line arguments
dirname = sys.argv[1]
radatadir = dirname + '/Shell_Slices/'

file_list, int_file_list, nfiles = get_file_lists(radatadir)
Omega0 = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/Omega0

minmax = None
varname = 'vr' # by default plot the radial velocity
ir = 0 # by default plot just below the surface
rval = None # can also find ir by finding the closest point
clon = 0

args = sys.argv[2:]
nargs = len(args)
try:
    index_first, index_last = get_desired_range(int_file_list, args)
except:
    index_first, index_last = 0, nfiles - 1  
    # By default plot all the shell slices
fnames = file_list[index_first:index_last+1] 

for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-var'):
        varname = args[i+1]
    elif (arg == '-ir'):
        ir = int(args[i+1])
    elif arg == '-rval':
        rval = float(args[i+1])
    elif (arg == '-clon'):
        clon = float(args[i+1])

posdef = False
if 'sq' in varname:
    posdef = True

varlabel = texlabels[varname]
count = 0

# Get some useful info from the first sslice
a = Shell_Slices(radatadir + fnames[0], '')
# Figure out ir if user specified rval
if not rval is None:
    ir = np.argmin(np.abs(a.radius/rsun - rval))
# Replace desired rval (or "None") by the actual rval
rval = a.radius[ir]/rsun
# Start time
time0 = a.time[0]

field = get_sslice(a, varname, dirname=dirname)[:, :, ir]
if minmax is None: # Set saturation values by the field from
    # the first slice
    field = get_sslice(a, varname, dirname=dirname)[:, :, ir]
    minmax = get_satvals(field, posdef=posdef)

for fname in fnames:
    # Read in desired shell slice
    a = Shell_Slices(radatadir + fname, '')
    for j in range(a.niter):
        time = a.time[j]
            
        # Create the plot using subplot axes
        # Offset axes slightly (at the end) to deal with annoying white space cutoff
        fig_width_inches = 12.

        # General parameters for main axis/color bar
        margin_bottom_inches = 1.
        margin_top_inches = 1.
        margin_inches = 1/8

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

        # Loop over depths and make plots
        rval = a.radius[ir]/rsun
        plotdir = dirname + '/plots/moll/movies/' + varname + '/' +\
                ('rval%0.3f' %rval) + '/'
        print(plotdir)
#        if not os.path.isdir(plotdir):
#            os.makedirs(plotdir)
            
        savename = 'img' + str(count).zfill(4) + '.png'
        print(savename)
        print('Plotting moll: ' + varname + (', rval = %0.3f, ' %rval) +\
              'iter ' + fname + ' ...')
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = fig.add_axes([margin_x, margin_bottom, subplot_width, subplot_height])
        
        plot_moll(fig, ax, a, dirname, varname, ir=ir, minmax=minmax,\
                    clon=clon, plot_title=False) 
        title = varlabel + '     ' +\
                (r'$r/R_\odot\ =\ %0.3f$' %rval) +\
                '     ' + (r'$t = %02.1f\ P_{\rm{rot}}$' %(time/Prot))
        fig.text(margin_x + 0.5*subplot_width, 1. - 0.5*margin_top,\
                title, ha='center', va='bottom', **csfont, fontsize=14)
#        plt.savefig(plotdir + savename, dpi=300)
#        plt.close()
        plt.show()
        count += 1
