# Author: Loren Matilsky
# Created: 06/18/2019
# This script plots generic quantities (speicified at the command line
# by providing a list of quantity codes) with respect to latitude for
# the Rayleigh run directory indicated by [dirname]. 
# For viewing purposes only (does not save a plot)
# Plots at various r-values; change the default with "-rvals"

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
sys.path.append(os.environ['raco'])
from get_parameter import get_parameter
from common import get_file_lists, rsun
from rayleigh_diagnostics import Meridional_Slices
from get_eq import get_eq

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
# Get radial geometry from ReferenceState
eq = get_eq(dirname)
ri, ro = np.min(eq.radius), np.max(eq.radius)

radatadir = dirname + '/Meridional_Slices/'
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
lats = [0., 15., 30., 45., 60., 75.]
rvals = np.linspace(ri/rsun, ro/rsun, 7)
qvals = [1, 2, 3]
logscale = False
iiter = nfiles - 1 # by default, plot the last data file in the list
iphi = 0 # first longitude in slice array, by default

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for j in range(len(rvals_str)):
            rvals.append(float(rvals_str[j]))
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
lats = 90. - tt*180./np.pi

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
                           
    # Plot rotation vs radius at the desired latitudes
    for rval in rvals:
        diffs = np.abs(rr/rsun - rval)
        ir = np.argmin(diffs)
        rval = rr[ir]/rsun
        ax.plot(lats, qty[:, ir],\
                label = r'$\rm{%0.3f}$' %rval, linewidth=0.5)

    # Label the axes
    plt.xlabel('latitude ' + r'$(^\circ)$', fontsize=12, **csfont)
    ax.set_title('qval = %i' %qval)

    if logscale:
        ax.set_yscale('symlog')


    # Get ticks everywhere
    plt.minorticks_on()
    plt.tick_params(top=True, right=True, direction='in', which='both')

    iplot += 1

# Set the axis limits
xmin, xmax = -90., 90.
plt.xlim((xmin, xmax))
plt.tight_layout()
plt.legend(title='r/rsun')
plt.show()
