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
sys.path.append(os.environ['rapp'])
from get_parameter import get_parameter
from common import get_file_lists, rsun
from rayleigh_diagnostics import Meridional_Slices

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]

radatadir = dirname + '/Meridional_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
rnorm = None
lats = [0., 15., 30., 45., 60., 75.]
qvals = [1, 2, 3]
logscale = False
iiter = nfiles - 1 # by default, plot the last data file in the list
iphi = 0 # first longitude in slice array, by default
rvals_to_plot = None # user can specify r-values to mark by vertical lines

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
    elif (arg == '-iter'):
        desired_iter = int(args[i+1])
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-sec':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='sec')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-day':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='day')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-prot':
        time = float(args[i+1])
        di_trans = translate_times(time, dirname, translate_from='prot')
        desired_iter = di_trans['val_iter']
        iiter = np.argmin(np.abs(int_file_list - desired_iter))
    elif arg == '-iphi':
        iphi = int(args[i+1])
    elif arg == '-rvals':
        rvals_to_plot = []
        rvals_str = args[i+1].split()
        for j in range(len(rvals_str)):
            rvals_to_plot.append(float(rvals_str[j]))

iter_val = int_file_list[iiter]
fname = file_list[iiter]


# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90. - lats
theta_vals = colats*np.pi/180.

# Read in Meridional_Slices data
print ('Reading ' + radatadir + fname + ' ...')
mer = Meridional_Slices(radatadir + fname, '')

phival = mer.phi[iphi]*180./np.pi
vals = mer.vals[iphi, :, :, :, 0] 
lut = mer.lut
# only consider the first item of the record
# (good idea to keep nrec = 1 for Meridional_Slices in general)

rr = mer.radius
tt = np.arccos(mer.costheta)

nq = len(qvals)
ncol = 3
nrow = int(np.ceil(nq/ncol))

subplot_side = 4. # 4 inch subplot
fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*subplot_side,\
        nrow*subplot_side), sharex=True)
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
        latitude = 90. - tt[index]*180./np.pi
        ax.plot(rr_n, qty[index,:],\
                label = r'$\rm{%2.2f}$' %latitude + r'$^\circ$',\
                linewidth=0.5)
    if not rvals_to_plot is None:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax)
        yvals = np.linspace(ymin, ymax, 100)
        for rval in rvals_to_plot:
            ax.plot(np.zeros(100) + rval/rsun, yvals, 'k--', \
                    linewidth=0.5)

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
