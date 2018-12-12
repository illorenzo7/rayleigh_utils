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
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from azavg_util import plot_azav, streamfunction
from common import get_widest_range_file, strip_dirname, get_iters_from_file
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

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)
rr_2d = rr.reshape((1,nr))
sint_2d = sint.reshape((nt, 1))
rsint = rr_2d*sint_2d

# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
vals, qv, counts, iters1, iters2 = np.load(datadir + AZ_Avgs_file)
iter1, iter2 = get_iters_from_file(AZ_Avgs_file)
ind_vr, ind_vt, ind_vp = np.argmin(np.abs(qv - 1)), np.argmin(np.abs(qv - 2)),\
    np.argmin(np.abs(qv - 3))

vr_av, vt_av, vp_av = vals[:, :, ind_vr], vals[:, :, ind_vt],\
        vals[:, :, ind_vp]

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
fig_width, fig_height = 4, 7
margin_inches = .5
aspect = fig_height/fig_width
margin_x = margin_inches/fig_width
margin_y = margin_inches/fig_height

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
plt.subplots_adjust(left=margin_x, right=1-margin_x, bottom=margin_y, top=1-margin_y)
plot_width, plot_height = 1 - 2*margin_x, 1 - 2*margin_y

plot_azav (fig, ax, rhovm, rr, cost, sint, units = r'$\rm{g}\ \rm{cm}^{-1}\ \rm{s}^{-1}$',
        plotcontours=False, norm=MidpointNormalize(0),\
        boundstype=my_boundstype, caller_minmax = (my_min, my_max))

lilbit = 0.02
maxabs = np.max(np.abs(psi))
plot_azav (fig, ax, psi, rr, cost, sint, units = r'$g\ cm^{-2}\ s^{-1}$', plotfield=False,
        norm=MidpointNormalize(0), levels=(-maxabs/2, -maxabs/4,  -lilbit*maxabs, 0,\
                lilbit*maxabs, maxabs/4, maxabs/2),\
        boundstype=my_boundstype, caller_minmax = (my_min, my_max))

# Make title
plt.title(dirname_stripped + '\n ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8))
#plt.text(1 - 2*margin_x, 1 - 2*margin_y,  r'$\Delta\Omega_{\rm{tot}} = %.1f$' %Delta_Om)
#plt.text(1 - 2*margin_x, 1 - 3*margin_y, 'nlevs = %i' %my_nlevs)
savefile = plotdir + dirname_stripped + '_mercirc_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
if showplot:
    plt.show()
else:
    plt.close()
