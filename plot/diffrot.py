# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. To use a vavg file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from binormalized_cbar import MidpointNormalize
import sys
import os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from azavg_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_iters_from_file

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
showplot = False
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
    elif (arg == '-show'):
        showplot = True
    elif (arg == '-nlevs'):
        my_nlevs = int(args[i+1])
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
mydict = (np.load(datadir + AZ_Avgs_file)).item()

#iter1, iter2 = get_iters_from_file(AZ_Avgs_file)
iter1, iter2 = mydict['iter1'], mydict['iter2']
vals = mydict['vals']
lut = mydict['lut']
#ind_vr, ind_vt, ind_vp = np.argmin(np.abs(qv - 1)), np.argmin(np.abs(qv - 2)),\
#    np.argmin(np.abs(qv - 3))

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Get grid info
rr = mydict['rr']
nr = mydict['nr']
tt = mydict['tt']
nt = mydict['nt']
cost = mydict['cost']
sint = mydict['sint']
ri = mydict['ri']
ro = mydict['ro']
d = mydict['d']

#rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
#nr, nt = len(rr), len(tt)
rr_2d = rr.reshape((1,nr))
sint_2d = sint.reshape((nt, 1))
rsint = rr_2d*sint_2d


# Get differential rotation in the rotating frame. 
Om = vp_av/rsint
diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

# Maximum differential rotation over whole meridional plane
it15, it75 = np.argmin(np.abs(tt - 11*np.pi/12)),\
    np.argmin(np.abs(tt - np.pi/12)) # ignore problematic poles 
global_min, global_max = np.min(diffrot[it15:it75, :]), np.max(diffrot[it15:it75, :])
Delta_Om = global_max - global_min
maxabs = np.max((np.abs(global_min), np.abs(global_max)))

if (not user_specified_minmax):
    my_min, my_max = -maxabs, maxabs

# Create plot
fig_width, fig_height = 4, 7
margin_inches = .5
aspect = fig_height/fig_width
margin_x = margin_inches/fig_width
margin_y = margin_inches/fig_height

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
plt.subplots_adjust(left=margin_x, right=1-margin_x, bottom=margin_y, top=1-margin_y)
plot_width, plot_height = 1 - 2*margin_x, 1 - 2*margin_y

plot_azav (fig, ax, diffrot, rr, cost, sint, units = 'nHz', nlevs=my_nlevs,
        norm=MidpointNormalize(0),\
        boundstype=my_boundstype, caller_minmax = (my_min, my_max))

# Make title
plt.title(dirname_stripped + '\n ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8))
plt.text(1 - 2*margin_x, 1 - 2*margin_y,  r'$\Delta\Omega_{\rm{tot}} = %.1f$' %Delta_Om)
plt.text(1 - 2*margin_x, 1 - 3*margin_y, 'nlevs = %i' %my_nlevs)
savefile = plotdir + dirname_stripped + '_diffrot_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
else:
    plt.close()
