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
from get_parameter import get_parameter
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun, get_file_lists
from rayleigh_diagnostics import ReferenceState

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
# Get radial geometry from ReferenceState
ref = ReferenceState(dirname + '/reference')
ri, ro = np.min(ref.radius), np.max(ref.radius)

# Set defaults
rvals = np.linspace(ri/rsun, ro/rsun, 7)
qvals = [1, 2, 3]
datadir = dirname + '/data/'
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
logscale = False

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
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-qvals':
        qvals = [] 
        qvals_str = args[i+1].split()
        for j in range(len(qvals_str)):
            qvals.append(int(qvals_str[j]))
    elif arg == '-log':
        logscale = True

# Read in AZ_Avg data
print ('Reading AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost, sint = di['cost'], di['sint']
lats = di['tt_lat']
xx = di['xx']
ri = di['ri']

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
