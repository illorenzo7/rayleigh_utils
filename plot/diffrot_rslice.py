# Author: Loren Matilsky
# Created: 05/14/2018
# Last Modified: 12/10/2018
# This script generates differential rotation plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_rslice_[first iter]_[last iter].npy

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import sys, os
from get_parameter import get_parameter
from common import strip_dirname, get_widest_range_file, get_iters_from_file

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
showplot = False
lats = [0, 15, 30, 45, 60, 75]
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-show'):
        showplot = True
    elif (arg == '-lats'):
        lats_str = args[i+1].split()
        lats = []
        for j in range(len(lats_str)):
            lats.append(float(lats_str[j]))
    elif (arg == '-dimple'):
        dimple = True
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90 - lats
theta_vals = colats*np.pi/180

# Get grid info
rr, tt, cost, sint, rr_depth, ri, ro, d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)
rr_2d = rr.reshape((1, nr))
sint_2d = sint.reshape((nt, 1))
rsint = rr_2d*sint_2d

# Read in vavg data
print ('Reading AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
vals, qv, counts, iters1, iters2 = np.load(datadir + AZ_Avgs_file)
ind_vr, ind_vt, ind_vp = np.argmin(np.abs(qv - 1)), np.argmin(np.abs(qv - 2)),\
    np.argmin(np.abs(qv - 3))

vr_av, vt_av, vp_av = vals[:, :, ind_vr], vals[:, :, ind_vt],\
        vals[:, :, ind_vp]

iter1, iter2 = get_iters_from_file(AZ_Avgs_file)

# Get frame rate rotation and compute differential rotation in the 
# lab frame. 
Om0 = get_parameter(dirname, 'angular_velocity')
Om = vp_av/rsint + Om0
Om *= 1.0e9/2/np.pi # convert from rad/s --> nHz

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value

# Plot rotation vs radius at the desired latitudes
for theta_val in theta_vals:
    diffs = np.abs(tt - theta_val)
    index = np.argmin(diffs)
    latitude = 90 - theta_val*180/np.pi
    ax.plot(rr/ro, Om[index,:], linewidth=2., 
            label = r'$\rm{%2.1f}$' %latitude + r'$^\circ$')
    maxes.append(np.max(Om[index,:]))
    mins.append(np.min(Om[index,:]))

# Global extrema
mmax = np.max(maxes)
mmin = np.min(mins)
difference = mmax - mmin
buffer = difference*0.2 # "Guard" the yrange of the plot with whitespace

# Label the axes
plt.xlabel(r'$r/r_o$',fontsize=12)
plt.ylabel(r'$\Omega/2\pi \ \rm{(nHz)}$',fontsize=12)

# Set the axis limits
xmin, xmax = ri/ro, 1
ymin, ymax = mmin - buffer, mmax + buffer
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
delta_x = xmax - xmin
delta_y = ymax - ymin
xvals = np.linspace(xmin, xmax, 100)
yvals = np.linspace(ymin, ymax, 100)

# Create a title    
plt.title(dirname_stripped + ', Diff. Rot., radial slice, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8))
plt.legend(title='latitude')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='on', right='on', direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_diffrot_rslice_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.pdf'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile)
if (showplot): plt.show()
plt.close()
