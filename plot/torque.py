import numpy as np
from azavg_util import plot_azav
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from binormalized_cbar import MidpointNormalize
from common import get_widest_range_file, get_iters_from_file, strip_dirname

dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

datadir = dirname + '/data/'
plotdir = dirname + '/plots/'

if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Read command-line arguments (CLAs)
my_boundstype = 'manual'
user_specified_minmax = False
showplot = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    if (arg == '-show'):
        showplot = True

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)
   
# Get the torques:
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
iter1, iter2 = get_iters_from_file(AZ_Avgs_file)
print ('Getting torques from ' + datadir + AZ_Avgs_file + ' ...')
vals, qv, counts, iters1, iters2 = np.load(datadir + AZ_Avgs_file)
ind_pp = np.argmin(np.abs(qv - 1801))
ind_mm = np.argmin(np.abs(qv - 1802))
ind_cor = np.argmin(np.abs(qv - 1803))
ind_visc = np.argmin(np.abs(qv - 1804))

torque_rs, torque_mc, torque_visc = -vals[:, :, ind_pp],\
        -vals[:, :, ind_mm] + vals[:, :, ind_cor],\
        vals[:, :, ind_visc]
torque_tot = torque_rs + torque_mc + torque_visc
maxabs = max(np.max(np.abs(torque_rs)), np.max(np.abs(torque_mc)), np.max(np.abs(torque_visc)))

if not user_specified_minmax: 
    my_min, my_max = -maxabs, maxabs

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches (i.e., 8x11.5 paper with 1/2-inch
                # margins)
margin_inches = 0.25 # margin width in inches (for both x and y) and horizontally 
            # in between figures
nplots = 4

subplot_width_inches = (fig_width_inches - 2*margin_inches -\
    (nplots - 1)*margin_inches)/nplots # Make the subplot width so that three subplots
    # fit together side-by-side with margins in between them and at the left
    # and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = subplot_height_inches + 2*margin_inches # Make overall
    # figure height able to accommodate margins on the top and bottom
fig_aspect = fig_height_inches/fig_width_inches

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches

# Location of "main axis" (axis on figure accommodating the margins)
# in figure units
main_axis_left = margin_x
main_axis_bottom = margin_y
main_axis_width = 1 - 2*margin_x
main_axis_height = 1 - 2*margin_y

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

# Set up three sets of axes to hold the subplots
ax1 = fig.add_axes((main_axis_left, main_axis_bottom,\
                    subplot_width, subplot_height))
ax2 = fig.add_axes((main_axis_left + subplot_width + margin_x,\
                    main_axis_bottom, subplot_width, subplot_height))
ax3 = fig.add_axes((main_axis_left + 2*(subplot_width + margin_x),\
                    main_axis_bottom, subplot_width, subplot_height))

ax4 = fig.add_axes((main_axis_left + 3*(subplot_width + margin_x),\
                    main_axis_bottom, subplot_width, subplot_height))


plot_azav (fig, ax1, torque_rs, rr, cost, sint, plotcontours=False, 
           units = r'$\frac{\rm{g}}{\rm{cm}\ \rm{s}^2}$',
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

plot_azav (fig, ax2, torque_mc, rr, cost, sint, plotcontours=False, 
           units = r'$\frac{\rm{g}}{\rm{cm}\ \rm{s}^2}$',
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

plot_azav (fig, ax3, torque_visc, rr, cost, sint, plotcontours=False, 
           units = r'$\frac{\rm{g}}{\rm{cm}\ \rm{s}^2}$',
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

plot_azav (fig, ax4, torque_tot, rr, cost, sint, plotcontours=False, 
           units = r'$\frac{\rm{g}}{\rm{cm}\ \rm{s}^2}$',
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

ax1.set_title(r'$\tau_{\rm{rs}}$', fontsize=16, verticalalignment='top')
ax2.set_title(r'$\tau_{\rm{mc}}$', fontsize=16, verticalalignment='top')
ax3.set_title(r'$\tau_{\rm{v}}$', fontsize=16, verticalalignment='top')
ax4.set_title(r'$\tau_{\rm{tot}}$', fontsize=16, verticalalignment='top')

savefile = plotdir + dirname_stripped + '_torque_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'

print ('Saving torques at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
if (showplot): plt.show()
plt.close()
