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
sys.path.append(os.environ['raco'])
from common import *
        get_iters_from_file, get_dict, rsun
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

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
    elif arg == '-rvals':
        rvals_to_plot = []
        rvals_str = args[i+1].split()
        for j in range(len(rvals_str)):
            rvals_to_plot.append(float(rvals_str[j]))

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

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

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

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Put some metadata in upper left
fig_width_inches, fig_height_inches = fig.get_size_inches()
margin_x = 1./8./fig_width_inches
margin_y = 1./8./fig_height_inches
fsize = 12
line_height = 1./4./fig_height_inches
fig.text(margin_x, 1. - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - margin_y - line_height,\
        time_string, ha='left', va='top', fontsize=fsize, **csfont)

plt.xlim((xmin, xmax))
plt.tight_layout()
plt.subplots_adjust(top=0.7)
plt.legend(title='latitude')
plt.show()
