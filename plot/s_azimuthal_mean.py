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

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

sph_mean_file = get_widest_range_file(datadir, 's_spherical_mean')
az_mean_file = get_widest_range_file(datadir, 's_azimuthal_mean')
s_sph = np.load(datadir + sph_mean_file)
s_az = np.load(datadir + az_mean_file)

s_fluc = s_az - s_sph.reshape((1,nr))

# Set up the actual figure from scratch
fig_width_inches = 3.5 # TOTAL figure width, in inches (i.e., 8x11.5 paper with 1/2-inch
                # margins)
margin_inches = 0.25 # margin width in inches (for both x and y) and horizontally 
            # in between figures
nplots = 1

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

std_plus = np.sqrt(np.mean(s_fluc[np.where(s_fluc > 0)]**2))
std_minus = np.sqrt(np.mean(s_fluc[np.where(s_fluc < 0)]**2))
my_min, my_max = -3*std_minus, 3*std_plus

plot_azav (fig, ax1, s_fluc, rr, cost, sint, plotcontours=False, 
           units = r'$\frac{\rm{erg}}{\rm{K}\ \rm{g}}$',
           boundstype = 'manual', caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

savename = plotdir + dirname_stripped + '_s_azimuthal_mean.png'

plt.savefig(savename, dpi=300)
plt.close()