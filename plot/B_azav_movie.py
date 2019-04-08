# Author: Loren Matilsky
# Created: 04/06/2019
# This script plots mean magnetic field components in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. 
# Uses an instantaneous AZ_Avgs file unless told otherwise
# Saves plot in
# [dirname]_B_azav_[iter].npy
# or [dirname]_B_azav_[first iter]_[last iter].npy if a time average was 
# specified

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
from binormalized_cbar import MidpointNormalize
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from azavg_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_file_lists,\
        get_desired_range
from rayleigh_diagnostics import AZ_Avgs
from get_parameter import get_parameter

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

radatadir = dirname + '/AZ_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Get rotation rate
Omega0 = get_parameter(dirname, 'angular_velocity')
Prot = 2*np.pi/Omega0

# Set defaults
save = True
plotcontours = True
my_boundstype = 'manual'
user_specified_minmax = False 
count = 0

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)

# Get desired range of files we wish to plot
if '-range' in args or '-n' in args or '-centerrange' in args or\
        '-all' in args or '-iter' in args:
    index_first, index_last = get_desired_range(int_file_list, args)
else:
    index_first, index_last = 0, nfiles - 1  
    # By default, don't average over any files;
    # just plot the last file

iter1, iter2 = int_file_list[index_first], int_file_list[index_last]

for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        user_specified_minmax = True
        min_br, max_br  = float(args[i+1]), float(args[i+2])
        min_bt, max_bt  = float(args[i+3]), float(args[i+4])
        min_bp, max_bp  = float(args[i+5]), float(args[i+6])
    elif arg == '-nosave':
        save = False
    elif arg == '-start':
        count = int(args[i+1])

# Get grid information from first AZ_Avgs file
az0 = AZ_Avgs(radatadir + file_list[index_first], '')
rr = az0.radius
ri, ro = np.min(rr), np.max(rr)
d = ro - ri
rr_depth = (ro - rr)/d
rr_height = (rr - ri)/d
sint = az0.sintheta
cost = az0.costheta
tt = np.arccos(cost)
tt_lat = (np.pi/2 - tt)*180/np.pi
nr = az0.nr
nt = az0.ntheta

# compute some derivative quantities for the grid
tt_2d, rr_2d = np.meshgrid(tt, rr, indexing='ij')
sint_2d = np.sin(tt_2d); cost_2d = np.cos(tt_2d)
xx = rr_2d*sint_2d
zz = rr_2d*cost_2d

# Get saturation values for B field (if user didn't specify them)
if (not user_specified_minmax):
    br0 = az0.vals[:, :, az0.lut[801], 0]
    bt0 = az0.vals[:, :, az0.lut[802], 0]
    bp0 = az0.vals[:, :, az0.lut[803], 0]
    nstd = 5
    min_br, max_br = -nstd*np.std(br0), nstd*np.std(br0)
    min_bt, max_bt = -nstd*np.std(bt0), nstd*np.std(bt0)
    min_bp, max_bp = -nstd*np.std(bp0), nstd*np.std(bp0)

# Create the save directory if it doesn't already exist
plotdir = dirname + '/plots/Bazav_movie/' + '/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Before plotting, set up figure dimensions
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 3/4 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1 # margin to accommodate just subplot titles
ncol = 3 # put three plots per row
nrow = 1

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*subplot_height_inches + margin_top_inches +\
    (nrow - 1)*margin_subplot_top_inches + margin_inches 
    # Room for titles on each row and a regular margin on the bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

#   Plot desired range, 1 image at a time
print ('Considering AZ_Avgs files %s through %s to plot the B field...'\
       %(file_list[index_first], file_list[index_last]))
for i in range(index_first, index_last + 1):
    print ('Plotting AZ_Avgs/%s ...' %file_list[i])
    if i == index_first:
        az = az0
    else:   
        az = AZ_Avgs(radatadir + file_list[i], '')

    local_ntimes = az.niter
    for j in range(local_ntimes):
        savename = 'img' + str(count).zfill(4) + '.png'
        count += 1
        br = az.vals[:, :, az0.lut[801], j]
        bt = az.vals[:, :, az0.lut[802], j]
        bp = az.vals[:, :, az0.lut[803], j]

        field_components = [br, bt, bp]
        titles = [r'$B_r$', r'$B_\theta$', r'$B_\phi$']
        units = r'$\rm{G}$'
        my_mins = [min_br, min_bt, min_bp]
        my_maxes = [max_br, max_bt, max_bp]

        # Generate a figure of the correct dimensions
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

        for iplot in range(3):
            ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
            ax_bottom = 1 - margin_top - subplot_height - \
                    (iplot//ncol)*(subplot_height + margin_subplot_top)
            ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
            plot_azav (fig, ax, field_components[iplot], rr, cost, sint,\
                   units = units, boundstype = my_boundstype,\
                   caller_minmax = (my_mins[iplot], my_maxes[iplot]),\
                   norm=MidpointNormalize(0), plotcontours=False)
            ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

        # Make the title indicating the simulation time
        fsize = 12
        time = az.time[j]
        title = r'$t = %02.1f\ P_{\rm{rot}}$' %(time/Prot)
        fig.text(margin_x + 0.5*subplot_width + (margin_x + subplot_width),\
                1. - 0.3*margin_top, title, ha='center',\
                va='center', **csfont)
        print("Saving " + plotdir + savename + ' ...')
        plt.savefig(plotdir + savename, dpi=200)
        plt.close()
