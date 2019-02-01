# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
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
from common import get_widest_range_file, strip_dirname

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

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
di = np.load(datadir + AZ_Avgs_file, encoding='latin1').item()

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Get necessary grid info
rr = di['rr']
cost = di['cost']
sint = di['sint']
tt_lat = di['tt_lat']
xx = di['xx']

# Get differential rotation in the rotating frame. 
Om = vp_av/xx
diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

# Maximum differential rotation over whole meridional plane
it15, it75 = np.argmin(np.abs(tt_lat - 15)), np.argmin(np.abs(tt_lat - 75))
     # ignore problematic poles 
global_min, global_max = np.min(diffrot[it15:it75, :]),\
     np.max(diffrot[it15:it75, :])
Delta_Om = global_max - global_min
maxabs = np.max((np.abs(global_min), np.abs(global_max)))

if (not user_specified_minmax):
    my_min, my_max = -maxabs, maxabs

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 3/2 # larger top margin to make room for titles

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

plot_azav (fig, ax, diffrot, rr, cost, sint, units = 'nHz', nlevs=my_nlevs,
        norm=MidpointNormalize(0),\
        boundstype=my_boundstype, caller_minmax = (my_min, my_max))

# Make title + label diff. rot. contrast and no. contours
fsize = 12
fig.text(margin_x, 1 - 0.05*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.2*margin_top, r'$\Omega - \Omega_0$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.35*margin_top,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.5*margin_top,\
         r'$\Delta\Omega_{\rm{tot}} = %.1f\ nHz$' %Delta_Om,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.65*margin_top,\
         'nlevs = %i' %my_nlevs,
         ha='left', va='top', fontsize=fsize, **csfont)
savefile = plotdir + dirname_stripped + '_diffrot_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
