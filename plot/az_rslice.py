# Author: Loren Matilsky
# Created: 06/18/2019
# This script generate generic quantities (speicified at the command line
# by providing a list of quantity codes) plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# For viewing purposes only (does not save a plot)

# Import relevant modules
import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
from get_parameter import get_parameter
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun

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
rnorm = None
lats = [0., 15., 30., 45., 60., 75.]
qvals = [1, 2, 3]
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
logscale = False

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-lats':
        lats_str = args[i+1].split()
        lats = []
        for j in range(len(lats_str)):
            lats.append(float(lats_str[j]))
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-qvals':
        qvals = [] 
        qvals_str = args[i+1].split()
        for j in range(len(qvals_str)):
            qvals.append(int(qvals_str[j]))
    elif arg == '-log':
        logscale = True

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90. - lats
theta_vals = colats*np.pi/180.

# Read in AZ_Avg data
print ('Reading AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost, sint = di['cost'], di['sint']
xx = di['xx']
ri = di['ri']

nq = len(qvals)
ncol = 3
nrow = int(np.ceil(nq/ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*3., nrow*3.), sharex=True)
if nrow == 1:
    axs = axs.reshape((1, ncol))

iplot = 0
for qval in qvals:
    irow = iplot//ncol
    icol = iplot - ncol*irow
    ax = axs[irow, icol]
    qty = vals[:, :, lut[qval]]
                           
    # User can specify what to normalize the radius by
    # By default, normalize by the solar radius
    if rnorm is None:
        rr_n = rr/rsun
    else:
        rr_n = rr/rnorm

    # Plot rotation vs radius at the desired latitudes
    for theta_val in theta_vals:
        diffs = np.abs(tt - theta_val)
        index = np.argmin(diffs)
        latitude = 90. - theta_val*180./np.pi
        ax.plot(rr_n, qty[index,:],\
                label = r'$\rm{%2.1f}$' %latitude + r'$^\circ$')

    # Label the axes
    if rnorm is None:
        plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
    else:
        plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12,\
                **csfont)
    ax.set_title('qval = %i' %qval)

    if logscale:
        ax.set_yscale('symlog')

    # Set the axis limits
    xmin, xmax = np.min(rr_n), np.max(rr_n)

    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    iplot += 1

plt.xlim((xmin, xmax))
plt.tight_layout()
plt.legend(title='latitude')
plt.show()
