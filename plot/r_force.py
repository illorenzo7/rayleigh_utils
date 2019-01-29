# Author: Loren Matilsky
# Created: 01/28/2019
# This script plots the radial forces in the meridional plane (advection, 
# Coriolis, pressure, viscosity, and JxB (if present))
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one assocsiated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_r_force_[first iter]_[last iter].png

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'Times New Roman'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from azavg_util import plot_azav
from common import get_widest_range_file, strip_dirname
from get_parameter import get_parameter
from binormalized_cbar import MidpointNormalize

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
# Read command-line arguments (CLAs)
my_boundstype = 'manual'
user_specified_minmax = False

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

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get AZ_Avgs file
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting r_forces from ' + datadir + AZ_Avgs_file + ' ...')
di = np.load(datadir + AZ_Avgs_file).item()

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
 
ind_adv = lut[1201] # gets minus sign
ind_cor = lut[1219]
ind_prs = lut[1237]
ind_visc = lut[1228]

r_force_adv = -vals[:, :, ind_adv]
r_force_cor = vals[:, :, ind_cor]
r_force_prs = vals[:, :, ind_prs]
r_force_visc = vals[:, :, ind_visc]
r_force_tot = r_force_adv + r_force_cor + r_force_prs +\
    r_force_visc

max_sig = max(np.std(r_force_adv), np.std(r_force_cor),\
              np.std(r_force_prs), np.std(r_force_visc))

if magnetism:
    ind_mag = lut[1248]
    r_force_mag = vals[:, :, ind_mag]       
    max_sig = max(max_sig, np.std(r_force_mag))
    r_force_tot += r_force_mag
    
if not user_specified_minmax: 
    my_min, my_max = -3*max_sig, 3*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 3/8 # wider top margin to accommodate subplot titles
nplots = 5 + magnetism
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = nrow*(subplot_height_inches + margin_top_inches) +\
    margin_inches # Room for titles on each row and a regular margin on the 
                  # bottom
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

r_forces = [r_force_adv, r_force_cor, r_force_prs,\
                r_force_visc, r_force_tot]

titles =\
[r'$(\mathbf{f}_{\rm{adv}})_r$', r'$(\mathbf{f}_{\rm{cor}})_r$',\
 r'$(\mathbf{f}_{\rm{p}})_r$', r'$(\mathbf{f}_{\rm{v}})_r$',\
 r'$(\mathbf{f}_{\rm{tot}})_r$']
units = r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-2}$'

if magnetism:
    r_forces.insert(4, r_force_mag)
    titles.insert(4, r'$(\mathbf{f}_{\rm{mag}})_r$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - ((iplot//ncol) + 1)*(subplot_height + margin_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (fig, ax, r_forces[iplot], rr, cost, sint, plotcontours=False, 
           units = units,
           boundstype = my_boundstype, caller_minmax = (my_min, my_max),\
           norm=MidpointNormalize(0))

    ax.set_title(titles[iplot], verticalalignment='top', **csfont)

savefile = plotdir + dirname_stripped + '_r_forces_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

print ('Saving r_forces at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()