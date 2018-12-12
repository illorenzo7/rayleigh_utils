###############################################
# Author: Loren Matilsky
# Date created: 02/14/2018
# Last modified: 11/24/2018
#
#  This script computes the hell-averages radial energy flux in the 
#  spherical domain as a function of radius. 
#  Plots various spherically intebrated energy fluxes 

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
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
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Get command-line arguments to adjust the interval of averaging files
user_specified_minmax = False
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-usefile'):
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = eflux_file.split('/')[-1]
    elif (arg == '-minmax'):
        user_specified_minmax = True
        my_min, my_max = float(args[i+1]), float(args[i+2])
        
# Make the plot name, labelling the first/last iterations we average over
iter1, iter2 = get_iters_from_file(Shell_Avgs_file)
savename = dirname_stripped + '_eflux_radial_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Get grid info
rr,tt,cost,sint,rr_depth,ri,ro,d = np.load(datadir + 'grid_info.npy')
rr_n = rr/1e8 # rr "normalized" by a 1 Mm

#Create the plot
lw = 1.5 # Bit thicker lines

# Read in the flux data
vals, lut = np.load(datadir + Shell_Avgs_file)
qindex_hflux = lut[1433]
qindex_eflux = lut[1455]
qindex_cflux = lut[1470]
qindex_kflux = lut[1923]
qindex_vflux = lut[1935]

hflux = vals[:, qindex_hflux]
eflux = vals[:, qindex_eflux]
cflux = vals[:, qindex_cflux]
kflux = vals[:, qindex_kflux]
vflux = -vals[:, qindex_vflux]
tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux


# Compute the integrated fluxes
fpr = 4*np.pi*rr**2
hflux_int = hflux*fpr
eflux_int = eflux*fpr
cflux_int = cflux*fpr
kflux_int = kflux*fpr
vflux_int = vflux*fpr
tflux_int = tflux*fpr

# Create the plot; start with plotting all the energy fluxes
solar_lum = 3.846e33 # Normalize by the solar luminosity
plt.plot(rr_n, hflux_int/solar_lum, label = r'$\rm{F}_{heat}$', linewidth=lw)
plt.plot(rr_n, eflux_int/solar_lum, 'm', label = r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(rr_n, cflux_int/solar_lum, label = r'$\rm{F}_{cond}$', linewidth=lw)
plt.plot(rr_n, kflux_int/solar_lum, label = r'$\rm{F}_{KE}$', linewidth=lw)
plt.plot(rr_n, vflux_int/solar_lum, label = r'$\rm{F}_{visc}$', linewidth=lw)
plt.plot(rr_n, tflux_int/solar_lum, label=r'$\rm{F}_{total}$',\
        linewidth=lw, color='black')

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top='True', right='True', direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 
ymin, ymax = -0.7, 1.3
if user_specified_minmax:
    ymin, ymax = my_min, my_max
delta_y = ymax - ymin
plt.ylim(ymin, ymax)

# Label the axes
plt.xlabel(r'$\rm{Radius} \ (Mm)$', fontsize=12)
plt.ylabel(r'$4\pi r^2\ \rm{\times \ (Energy \ Flux)}\ /\ L_\odot$',\
        fontsize=12)

# Create a see-through legend
plt.legend(loc='lower left',shadow=True, ncol=3,fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
