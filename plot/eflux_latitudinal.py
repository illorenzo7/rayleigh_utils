###############################################
# Author: Loren Matilsky
# Date created: 11/18
# Last modified: 11/18/2018
#
# This script computes the "cone"-averages latitudinal energy flux in the 
# spherical domain as a function of radius. 
# Plots various average energy fluxes, integrated over cones of opening angle theta
# Since Rayleigh has no "cone"-averages, we assume symmetry in phi and use the AZ_Avgs

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from common import get_widest_range_file, strip_dirname, get_iters_from_file

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Get command-line arguments to adjust the interval of averaging files
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Make the plot name, labelling the first/last iterations we average over
iter1, iter2 = get_iters_from_file(AZ_Avgs_file)
savename = dirname_stripped + '_eflux_latitudinal_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
nr, nt = len(rr), len(tt)

#Create the plot
lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting AZ_Avgs data from %s ...' %AZ_Avgs_file)
vals, lut = np.load(datadir + AZ_Avgs_file)
qindex_eflux = lut[1456]
qindex_cflux = lut[1471]
qindex_kflux = lut[1924]
qindex_vflux = lut[1936]

# Get the fluxes in the whole meridional plane
eflux = vals[:, :, qindex_eflux]
cflux = vals[:, :, qindex_cflux]
kflux = vals[:, :, qindex_kflux]
vflux = vals[:, :, qindex_vflux]
tflux = eflux + cflux + kflux + vflux # compute the total flux

# Compute the integrated fluxes
# At each point in the meridional plane we associate a "ring" of width dr and circumference 2 pi r sin(theta)
dr = np.zeros_like(rr)
dr[1:nr-1] = 0.5*(rr[:nr-2] - rr[2:])
dr[0] = dr[1]; dr[-1] = dr[-2]
areas = 2*np.pi*sint.reshape((nt, 1))*rr.reshape((1, nr))*dr.reshape((1, nr))

eflux_int = np.sum(eflux*areas, axis=1)
cflux_int = np.sum(cflux*areas, axis=1)
kflux_int = np.sum(kflux*areas, axis=1)
vflux_int = np.sum(vflux*areas, axis=1)
tflux_int = np.sum(tflux*areas, axis=1)

# Create the plot; start with plotting all the energy fluxes
solar_lum = 3.846e33 # Normalize by the solar luminosity
lats = 180*(np.pi/2 - tt)/np.pi
plt.plot(lats, eflux_int/solar_lum, 'm', label = r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(lats, cflux_int/solar_lum, label = r'$\rm{F}_{cond}$', linewidth=lw)
plt.plot(lats, kflux_int/solar_lum, label = r'$\rm{F}_{KE}$', linewidth=lw)
plt.plot(lats, vflux_int/solar_lum, label = r'$\rm{F}_{visc}$', linewidth=lw)
plt.plot(lats, tflux_int/solar_lum, label=r'$\rm{F}_{total}$',\
        linewidth=lw, color='black')

# Get the y-axis in scientific notation
#plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = -90, 90
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits 
maxabs = max((np.max(np.abs(eflux_int)), np.max(np.abs(cflux_int)),\
        np.max(np.abs(kflux_int)), np.max(np.abs(vflux_int)), np.max(np.abs(tflux_int))))

ymin, ymax = -1.2*maxabs/solar_lum, 1.2*maxabs/solar_lum
delta_y = ymax - ymin
plt.ylim(ymin, ymax)

# Label the axes
plt.xlabel(r'$\rm{Latitude} \ (^\circ)$', fontsize=12)
plt.ylabel('(Integrated Energy Flux)' + r'$/L_\odot$',\
        fontsize=12)

# Create a see-through legend
leg=plt.legend(loc='lower left',shadow=True, ncol=3,fontsize=10)
leg.get_frame().set_alpha(0.3)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
