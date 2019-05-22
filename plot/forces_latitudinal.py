# Author: Loren Matilsky
# Created: 01/28/2019
# This script plots the latitudinal forces in the meridional plane (advection, 
# Coriolis, pressure, viscosity, and JxB (if present))
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one assocsiated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_theta_force_[first iter]_[last iter].png

import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
sys.path.append(os.environ['pl'])
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_dict
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
showplot = True
saveplot = True
plotcontours = True
minmax = None
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_boundstype = 'manual'
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif arg == '-noshow':
        showplot = False
    elif arg == '-nosave':
        saveplot = False
    elif arg == '-nocontour':
        plotcontours = False
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# See if magnetism is "on"
try:
    magnetism = get_parameter(dirname, 'magnetism')
except:
    magnetism = False # if magnetism wasn't specified, it must be "off"

# Get AZ_Avgs file
print ('Getting theta_forces from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Get grid info
rr, tt, cost, sint, ri,ro,d = di['rr'], di['tt'], di['cost'],\
        di['sint'], di['ri'], di['ro'], di['d']
nr, nt = len(rr), len(tt)

ind_adv = lut[1202] # gets minus sign
ind_cor = lut[1220]
ind_prs = lut[1238]
ind_visc = lut[1229]

theta_force_adv = -vals[:, :, ind_adv]
theta_force_cor = vals[:, :, ind_cor]
theta_force_prs = vals[:, :, ind_prs]
theta_force_visc = vals[:, :, ind_visc]
theta_force_tot = theta_force_adv + theta_force_cor + theta_force_prs +\
    theta_force_visc

max_sig = max(np.std(theta_force_adv), np.std(theta_force_cor),\
              np.std(theta_force_prs), np.std(theta_force_visc))

if magnetism:
    ind_mag = lut[1249]
    theta_force_mag = vals[:, :, ind_mag]       
    max_sig = max(max_sig, np.std(theta_force_mag))
    theta_force_tot += theta_force_mag
    
if minmax is None:
    minmax = -3.*max_sig, 3.*max_sig

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 2 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1 # margin to accommodate just subplot titles
nplots = 5 + magnetism
ncol = 3 # put three plots per row
nrow = np.int(np.ceil(nplots/3))

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


theta_forces = [theta_force_adv, theta_force_cor, theta_force_prs,\
                theta_force_visc, theta_force_tot]

titles =\
[r'$(\mathbf{f}_{\rm{adv}})_\theta$', r'$(\mathbf{f}_{\rm{cor}})_\theta$',\
 r'$(\mathbf{f}_{\rm{p}})_\theta$', r'$(\mathbf{f}_{\rm{v}})_\theta$',\
 r'$(\mathbf{f}_{\rm{tot}})_\theta$']
units = r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-2}$'

if magnetism:
    theta_forces.insert(4, theta_force_mag)
    titles.insert(4, r'$(\mathbf{f}_{\rm{mag}})_\theta$')

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (theta_forces[iplot], rr, cost, sint, fig=fig, ax=ax,\
        units = units, minmax=minmax, norm=MidpointNormalize(0.),\
        plotcontours=plotcontours)

    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Latitudinal force balance (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_theta_forces_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if saveplot:
    print ('Saving theta_forces at ' + savefile + ' ...')
    plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
plt.close()
