# Author: Loren Matilsky
# Created: 05/14/2018
# This script plots the meridional circulation cells the meridional plane 
# for the Rayleigh run directory indicated by [dirname], using the AZ_Avgs
# data. To use an AZ_Avgs file different than the one associated with the 
# longest averaging range, run with option
# -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'Times New Roman'}
from binormalized_cbar import MidpointNormalize
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from azavg_util import plot_azav, streamfunction
from common import get_widest_range_file, strip_dirname
from rayleigh_diagnostics import ReferenceState

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Get density
ref = ReferenceState(dirname + '/reference', '')
rho = ref.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults
my_boundstype = 'manual'
user_specified_minmax = False 
my_nlevs = 20
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif (arg == '-nlevs'):
        my_nlevs = int(args[i+1])
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = np.load(datadir + AZ_Avgs_file).item()
#vals, qv, counts, iters1, iters2 = np.load(datadir + AZ_Avgs_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost = di['cost']
sint = di['sint']

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Compute the mass flux
rhovm = rho*np.sqrt(vr_av**2 + vt_av**2)

# Compute the streamfunction
psi = streamfunction(vr_av, vt_av, rr, cost)

# Make CCW negative and CW positive
rhovm *= np.sign(psi)

# Maximum/minimum mass flux over the meridional plane (excluding polar regions)
it15, it75 = np.argmin(np.abs(tt - 11*np.pi/12)),\
    np.argmin(np.abs(tt - np.pi/12)) # ignore problematic poles 
global_min, global_max = np.min(rhovm[it15:it75, :]), np.max(rhovm[it15:it75, :])
maxabs = np.max((np.abs(global_min), np.abs(global_max)))

if (not user_specified_minmax):
    my_min, my_max = -maxabs, maxabs

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 5/8 # larger top margin to make room for titles

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches + margin_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_y, subplot_width, subplot_height))

# Plot mass flux
plot_azav (fig, ax, rhovm, rr, cost, sint,\
    units = r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-1}$',
    plotcontours=False, norm=MidpointNormalize(0),\
    boundstype=my_boundstype, caller_minmax = (my_min, my_max))

# Plot streamfunction contours
lilbit = 0.01
maxabs = np.max(np.abs(psi))
plot_azav (fig, ax, psi, rr, cost, sint, plotfield=False,\
    norm=MidpointNormalize(0),\
    levels=(-maxabs/2, -maxabs/4,  -lilbit*maxabs, 0,\
            lilbit*maxabs, maxabs/4, maxabs/2),\
    boundstype=my_boundstype, caller_minmax = (my_min, my_max))

# Make title
fsize = 10
fig.text(margin_x, 1 - 1/8*margin_top,\
         r'$|\langle\overline{\rho}\mathbf{v}_m\rangle|$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 3/8*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 5/8*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)

savefile = plotdir + dirname_stripped + '_mercirc_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()