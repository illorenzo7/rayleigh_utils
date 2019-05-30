# Author: Loren Matilsky
# Date created: 06/08/2017
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys, os
from subprocess import call
from common import get_file_lists, get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]

# Data and plot directories
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)
dirname_stripped = strip_dirname(dirname)

# Find the etrace file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the trace.
trace_G_Avgs_file = get_widest_range_file(datadir, 'trace_G_Avgs')

# Set defaults
xiter = False
notfrom0 = False
magnetism = False
minmax = None
xminmax = None
plotall = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-xiter': # plot w.r.t. iterations
        xiter = True
    elif arg == '-usefile':
        trace_G_Avgs_file = args[i+1]
        trace_G_Avgs_file = trace_G_Avgs_file.split('/')[-1]
    elif arg == '-notfrom0':
        notfrom0 = True
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-plotall':
        plotall = True


# Tag the plot by whether or not the x axis is in "time" or "iteration"
if (xiter):
    tag = '_xiter'
else:
    tag = '_xtime'

# Read in the KE data (dictionary form)
di = get_dict(datadir + trace_G_Avgs_file)

vals = di['vals']
lut = di['lut']
times = di['times']
iters = di['iters']
iter1 = di['iter1']
iter2 = di['iter2']

# Get global rotation rate, if present
rotation = get_parameter(dirname, 'rotation')
if rotation:
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity
    tnorm = Prot # normalize the time axis by rotation period if applicable
else:
    ktop = get_parameter(dirname, 'kappa_top')
    try:
        rmin = get_parameter(dirname, 'rmin')
        rmax = get_parameter(dirname, 'rmax')
    except: # two domains stitched together
        domain_bounds = get_parameter(dirname, 'domain_bounds')
        rmin = np.min(domain_bounds)
        rmax = np.max(domain_bounds)
    depth = rmax - rmin
    tdt = d**2/ktop
    tnorm = tdt

# Make appropriate file name to save
savename = dirname_stripped + '_Ltrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

Lz = vals[lut[1819]]
if plotall:
    Lx = vals[lut[1820]]
    Ly = vals[lut[1821]]

Lpz = vals[lut[1822]]
if plotall:
    Lpx = vals[lut[1823]]
    Lpy = vals[lut[1824]]

Lmz = vals[lut[1825]]
if plotall:
    Lmx = vals[lut[1826]]
    Lmy = vals[lut[1827]]

# Get global min/max vals
varlist = [Lz, Lpz, Lmz]
if plotall:
    varlist.extend([Lx, Ly, Lpx, Lpy, Lmx, Lmy])
mmax = -np.inf
mmin = np.inf
for var in varlist:
    mmax = max(mmax, np.max(var))
    mmin = min(mmin, np.min(var))

if (not xiter):
    xaxis = times/tnorm
else:
    xaxis = iters

if (notfrom0):
    x_min = np.min(xaxis)
else:
    x_min = 0

# create figure with  3 panels in a row (total, mean and fluctuating amom)
fig, axs = plt.subplots(3, 1, figsize=(5, 10), sharex=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# Make thin lines to see structure of variation
lw = 0.5

# first plot: average angular momentum density trace      
ax1.set_title(dirname_stripped + '\n ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) +\
          '\n\nangular momentum')
ax1.plot(xaxis, Lz, 'k', linewidth=lw, label=r'$\mathcal{L}_z$')
if plotall:
    ax1.plot(xaxis, Lx, 'r', linewidth=lw, label=r'$\mathcal{L}_x$')
    ax1.plot(xaxis, Ly, 'g', linewidth=lw, label=r'$\mathcal{L}_y$')

# If we're looking for machine-precision variations, just determine
# min and max ranges for y-axes automatically in Python...
if minmax is None:
    pass
#    diff = mmax - mmin
#    buff = 0.05*diff
#    ax1.set_ylim(mmin - buff, mmax + buff)
else:
    ax1.set_ylim(minmax)
    
# Set x limits  
if xminmax is None:
    ax1.set_xlim((np.min(xaxis), np.max(xaxis)))
else:
    ax1.set_xlim(xminmax)

# legend
ax1.legend(ncol=2, fontsize=8)

# Make axes use scientific notation
ax1.ticklabel_format(scilimits = (-3,4), useMathText=True)
ax2.ticklabel_format(scilimits = (-3,4), useMathText=True)
ax3.ticklabel_format(scilimits = (-3,4), useMathText=True)

# Make the second plot (angular momentum of fluctuating motions)
ax2.plot(xaxis, Lpz, 'k', linewidth=lw, label=r'$\mathcal{L}_z^\prime$')
if plotall:
    ax2.plot(xaxis, Lpx, 'r', linewidth=lw, label=r'$\mathcal{L}_x^\prime$')
    ax2.plot(xaxis, Lpy, 'g', linewidth=lw, label=r'$\mathcal{L}_y^\prime$')

# Title and axis label
ax2.set_title('convection amom')
# Put the y-label on the middle plot
ax2.set_ylabel(r'$\rm{angular\ momentum\ density\ (g\ cm^{-1}\ s^{-1})}$')

# Third plot: angular momentum of mean energies
ax3.plot(xaxis, Lz, 'k', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_z$')
if plotall:
    ax3.plot(xaxis, Lx, 'r', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_x$')
    ax3.plot(xaxis, Ly, 'g', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_y$')

# title and x-axis label
ax3.set_title('mean-motion amom')

# Put the x-axis label on the bottom
if (xiter):
    ax3.set_xlabel('iteration #')
else:
    if rotation:
        ax3.set_xlabel(r'$\rm{t\ (P_{rot})}$')
    else:
        ax3.set_xlabel(r'$\rm{t\ (T_{diff})}$')

# Get ticks everywhere
plt.sca(ax1)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(ax2)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

plt.sca(ax3)
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Space the subplots to make them look pretty
plt.tight_layout
plt.subplots_adjust(left=0.15, bottom=0.08, top=0.85, wspace=0.4)

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
