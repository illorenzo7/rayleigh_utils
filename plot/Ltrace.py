# Author: Loren Matilsky
# Date created: 06/08/2017
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pickle
import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapp'])
from subprocess import call
from common import get_file_lists, get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from time_scales import compute_Prot, compute_tdt

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
from0 = False
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
    elif arg == '-from0':
        from0 = True
    elif arg == '-minmax':
        try: # first see if user wants to set all three ranges set
            # If "plotall", then use the first range for the total amom,
            # second for the convective, third for the mean
            minmax = float(args[i+1]), float(args[i+2]),\
                    float(args[i+3]), float(args[i+4]),\
                    float(args[i+5]), float(args[i+6])
        except: # if not, they want the same pair used for each subplot
            minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-plotall':
        plotall = True


# Tag the plot by whether or not the x axis is in "time" or "iteration"
if xiter:
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

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

# Get the first and last indices in the time direction, depending on 
# xminmax
ntimes = len(times)
if xminmax is None:
    it1, it2 = 0, ntimes - 1 # Plot over whole range by default
else:
    if xiter:
        it1 = np.argmin(np.abs(iters - xminmax[0]))
        it2 = np.argmin(np.abs(iters - xminmax[1]))
    else:
        it1 = np.argmin(np.abs(times/time_unit - xminmax[0]))
        it2 = np.argmin(np.abs(times/time_unit - xminmax[1]))

# Make appropriate file name to save
savename = dirname_stripped + '_Ltrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'

# Get data in apropriate time range
times = times[it1:it2+1]
t1 = np.min(times)
t2 = np.max(times)
iters = iters[it1:it2+1]
Lz = vals[lut[1819], it1:it2+1]
if plotall:
    Lx = vals[lut[1820], it1:it2+1]
    Ly = vals[lut[1821], it1:it2+1]

Lpz = vals[lut[1822], it1:it2+1]
if plotall:
    Lpx = vals[lut[1823], it1:it2+1]
    Lpy = vals[lut[1824], it1:it2+1]

Lmz = vals[lut[1825], it1:it2+1]
if plotall:
    Lmx = vals[lut[1826], it1:it2+1]
    Lmy = vals[lut[1827], it1:it2+1]

# Get global min/max vals
varlist = [Lz, Lpz, Lmz]
if plotall:
    varlist.extend([Lx, Ly, Lpx, Lpy, Lmx, Lmy])
mmax = -np.inf
mmin = np.inf
for var in varlist:
    mmax = max(mmax, np.max(var))
    mmin = min(mmin, np.min(var))

if not xiter:
    xaxis = times/time_unit
else:
    xaxis = iters

if from0:
    xmin = 0.
else:
    xmin = np.min(xaxis)

if xminmax is None:
    xminmax = xmin, np.max(xaxis)

# create figure with  3 panels in a row (total, mean and fluctuating amom)
fig, axs = plt.subplots(3, 1, figsize=(5, 10), sharex=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# Make thin lines to see structure of variation
lw = 0.5

# first plot: average angular momentum density trace      
# Include title

# Time string showing trace interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
ax1.set_title(dirname_stripped + '\n ' + time_string +\
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
    if len(minmax) == 2:
        ax1.set_ylim(minmax)
        ax2.set_ylim(minmax)
        ax3.set_ylim(minmax)
    elif len(minmax) == 6:
        ax1.set_ylim(minmax[0], minmax[1])
        ax1.set_ylim(minmax[2], minmax[3])
        ax1.set_ylim(minmax[4], minmax[5])
 
# Set x limits  
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
ax3.plot(xaxis, Lmz, 'k', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_z$')
if plotall:
    ax3.plot(xaxis, Lmx, 'r', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_x$')
    ax3.plot(xaxis, Lmy, 'g', linewidth=lw, label=r'$\langle\mathcal{L}\rangle_y$')

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
