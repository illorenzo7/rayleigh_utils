###############################################
# Author: Loren Matilsky
# Date created: 02/14/2018
#
# This script plots the radial energy fluxes as functions of
# radius using from the Shell_Avgs data

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
user_specified_rnorm = False
magnetism = False
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
    elif (arg == '-mag'):
        magnetism = True
    elif (arg == '-rnorm'):
        user_specified_rnorm = True
        user_supplied_rnorm = float(args[i+1])

#Create the plot
lw = 1.5 # Bit thicker lines

# Read in the flux data
di = np.load(datadir + Shell_Avgs_file).item()
vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_eflux_radial_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.pdf'

qindex_hflux = lut[1433]
qindex_eflux = lut[1455]
qindex_cflux = lut[1470]
qindex_kflux = lut[1923]
qindex_vflux = lut[1935]
if magnetism:
    qindex_mflux = lut[2001]
    mflux = vals[:, qindex_mflux]

hflux = vals[:, qindex_hflux]
eflux = vals[:, qindex_eflux]
cflux = vals[:, qindex_cflux]
kflux = vals[:, qindex_kflux]
vflux = -vals[:, qindex_vflux]

tflux = hflux + eflux + cflux + kflux + vflux # compute the total flux
if magnetism:
    tflux += mflux
    
# Compute the integrated fluxes
fpr = 4*np.pi*rr**2
hflux_int = hflux*fpr
eflux_int = eflux*fpr
cflux_int = cflux*fpr
kflux_int = kflux*fpr
vflux_int = vflux*fpr
tflux_int = tflux*fpr
if magnetism:
    mflux_int = mflux*fpr

# Create the plot; start with plotting all the energy fluxes
Lsun = 3.846e33 # Normalize the energy flux by the solar luminosity
Rsun = 6.955e10 # Normalize the radius by the solar radius

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if not user_specified_rnorm:
    rr_n = rr/Rsun
else:
    rr_n = rr/user_supplied_rnorm                                           

plt.plot(rr_n, hflux_int/Lsun, label=r'$\rm{F}_{heat}$', linewidth=lw)
plt.plot(rr_n, eflux_int/Lsun, 'm', label = r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(rr_n, cflux_int/Lsun, label = r'$\rm{F}_{cond}$', linewidth=lw)
plt.plot(rr_n, kflux_int/Lsun, label = r'$\rm{F}_{KE}$', linewidth=lw)
plt.plot(rr_n, vflux_int/Lsun, label = r'$\rm{F}_{visc}$', linewidth=lw)
plt.plot(rr_n, tflux_int/Lsun, label=r'$\rm{F}_{total}$',\
        linewidth=lw, color='black')
if magnetism:
    plt.plot(rr_n, mflux_int/Lsun, label=r'$\rm{F}_{Poynting}$',\
        linewidth=lw)

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
if not user_specified_rnorm:
    plt.xlabel(r'$r/R_\odot$',fontsize=12)
else:
    plt.xlabel(r'r/(%.1e cm)' %user_supplied_rnorm, fontsize=12)
plt.ylabel(r'$4\pi r^2\ \rm{\times \ (energy \ flux)}\ /\ L_\odot$',\
        fontsize=12)

# Make title
plt.title(dirname_stripped + '\n' + 'energy flux, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8))

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the eflux plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename)

# Show the plot
plt.show()
