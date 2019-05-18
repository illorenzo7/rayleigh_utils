# Author: Loren Matilsky
# Created: 01/29/2019
# This script plots the kinetic energy in the meridional plane (KE of
# diff. rot., KE of mer. circ., convective KE of r, theta, phi -- 
# 5 components)
# ...for the Rayleigh run directory indicated by [dirname]. 
# To use an AZ_Avgs file
# different than the one assocsiated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_r_force_[first iter]_[last iter].png

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
from azavg_util import plot_azav
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

# Get AZ_Avgs file
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
print ('Getting ke terms from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt = di['tt']
tt_lat = di['tt_lat']

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
 
ind_rke_tot = lut[402] 
ind_tke_tot = lut[403] 
ind_pke_tot = lut[404] 

ind_rke_fluc = lut[410] 
ind_tke_fluc = lut[411] 
ind_pke_fluc = lut[412] 

rke_tot = vals[:, :, ind_rke_tot]
tke_tot = vals[:, :, ind_tke_tot]
pke_tot = vals[:, :, ind_pke_tot]

rke_fluc = vals[:, :, ind_rke_fluc]
tke_fluc = vals[:, :, ind_tke_fluc]
pke_fluc = vals[:, :, ind_pke_fluc]

rke_mean = rke_tot - rke_fluc
tke_mean = tke_tot - tke_fluc
pke_mean = pke_tot - pke_fluc

mer_ke = rke_mean + tke_mean
ke_mean = mer_ke + pke_mean 
ke_fluc = rke_fluc + tke_fluc + pke_fluc
ke = ke_mean + ke_fluc

#max_sig = max(np.std(r_force_adv), np.std(r_force_cor),\
#              np.std(r_force_prs), np.std(r_force_buoy),\
#              np.std(r_force_visc))

if not user_specified_minmax: 
#    my_min, my_max = -3*max_sig, 3*max_sig
    my_min, my_max = 0., np.max(ke)

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_top_inches = 2 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1 # margin to accommodate just subplot titles
nplots = 8
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

ke_terms = [pke_mean, mer_ke, rke_fluc, tke_fluc, pke_fluc, ke_mean,\
        ke_fluc, ke]

titles =\
[r'$\overline{\rm{KE}}_{\rm{DR}}$', r'$\overline{\rm{KE}}_{\rm{MC}}$',\
    r'$\rm{KE}^\prime_r$', r'$\rm{KE}^\prime_\theta$',\
    r'$\rm{KE}^\prime_\phi$', r'$\overline{\rm{KE}}_{\rm{tot}}$',\
    r'$\rm{KE}^\prime_{\rm{tot}}$', r'$\rm{KE}_{\rm{tot}}$']

units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-2}$'

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(nplots):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - \
            (iplot//ncol)*(subplot_height + margin_subplot_top)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (fig, ax, ke_terms[iplot], rr, cost, sint, 
           units = units,
           caller_minmax = (my_min, my_max), mycmap='Greens')

    ax.set_title(titles[iplot], va='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Kinetic energies (zonally averaged)',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_ke_merplane_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

print ('Saving KE terms at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
