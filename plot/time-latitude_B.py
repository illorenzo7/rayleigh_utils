# Author: Loren Matilsky
# Date created: 06/08/2017
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.environ['co'])
from common import get_file_lists, get_widest_range_file, strip_dirname
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
time_latitude_file = get_widest_range_file(datadir, 'time-latitude')

user_specified_minmax = False
desired_depth = 0.5 # by default, plot time-latitude diagram for fields 
    # at mid-CZ

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        minmax = args[i + 1].split()
        my_min, my_max = float(minmax[0]), float(minmax[1])
    elif (arg == '-usefile'):
        trace_G_Avgs_file = args[i+1]
        trace_G_Avgs_file = trace_G_Avgs_file.split('/')[-1]
    elif (arg == '-depth'):
        desired_depth = float(args[i+1])

# Read in the time-latitude data (dictionary form)
print ('Reading in time-latitude trace from ' + datadir +\
       time_latitude_file + ' ...')
di = np.load(datadir + time_latitude_file, encoding='latin1').item()
vals = di['vals']

times = di['times']
iters = di['iters']
tt = di['tt']
tt_lat = di['tt_lat']
depths = np.array(di['depths'])
rinds = np.array(di['rinds'])
qvals = np.array(di['qvals'])

niter = di['niter']
ntheta = di['ntheta']
ndepths = di['ndepths']
nq = di['nq']

iter1 = di['iter1']
iter2 = di['iter2']
rr_depth = di['rr_depth']

# Get global rotation rate, if present
rotation = get_parameter(dirname, 'rotation')
if rotation:
    angular_velocity = get_parameter(dirname, 'angular_velocity')
    Prot = 2*np.pi/angular_velocity
    tnorm = Prot # normalize the time axis by rotation period if applicable
else:
    tnorm = 86400. # normalize by "days"

br_index = np.argmin(np.abs(qvals - 801))
bt_index = np.argmin(np.abs(qvals - 802))
bp_index = np.argmin(np.abs(qvals - 803))

idepth = np.argmin(np.abs(depths - desired_depth))
actual_depth = rr_depth[rinds[idepth]]
# Get the magnetic energies if they are available

br = vals[:, :, idepth, br_index]
bt = vals[:, :, idepth, bt_index]
bp = vals[:, :, idepth, bp_index]

# Make meshgrid of time/latitude
times2, tt_lat2 = np.meshgrid(times/tnorm, tt_lat, indexing='ij')

# Make appropriate file name to save
savename = dirname_stripped + '_time-latitude_B_' +\
    ('depth%0.2f_' %actual_depth) + str(iter1).zfill(8) + '_' +\
    str(iter2).zfill(8) + '.png'
    
# Create figure with  3 panels in a row (time-latitude plots of br, btheta,\
    # and bphi)
fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# first plot: total kinetic energy trace      
ax1.set_title(dirname_stripped + '\n ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) +\
          '\n\ntotal energy')
ax1.plot(xaxis, ke, 'k', linewidth=lw, label=r'$\rm{KE_{tot}}$')
ax1.plot(xaxis, rke, 'r', linewidth=lw, label=r'$\rm{KE}_r$')
ax1.plot(xaxis, tke, 'g', linewidth=lw, label=r'$\rm{KE}_\theta$')
ax1.plot(xaxis, pke, 'b', linewidth=lw, label=r'$\rm{KE}_\phi$')

# If magnetic, plot magnetic energies!
if magnetism:
    ax1.plot(xaxis, me, 'k--', linewidth=lw, label=r'$\rm{ME_{tot}}$')
    ax1.plot(xaxis, rme, 'r--', linewidth=lw, label=r'$\rm{ME}_r$')
    ax1.plot(xaxis, tme, 'g--', linewidth=lw, label=r'$\rm{ME}_\theta$')
    ax1.plot(xaxis, pme, 'b--', linewidth=lw, label=r'$\rm{ME}_\phi$')

if not ylog:
    ax1.set_ylim(0, mmax*1.05)
else:
    ax1.set_yscale('log')
    ax1.set_ylim((mmin/3.0, mmax*3.0))
    
# Set x limits  
ax1.set_xlim((x_min, np.max(xaxis)))

# legend
ax1.legend(ncol=2, fontsize=8)

# Make axes use scientific notation
# Only x-axis if on log scale
if ylog:
    ax1.ticklabel_format(axis='x', scilimits = (-3,4), useMathText=True)
else:
    ax1.ticklabel_format(scilimits = (-3,4), useMathText=True)

# Make the second plot (kinetic energy of the mean motions)
ax2.plot(xaxis, mke, 'k', linewidth=lw)
ax2.plot(xaxis, mrke, 'r', linewidth=lw)
ax2.plot(xaxis, mtke, 'g', linewidth=lw)
ax2.plot(xaxis, mpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    ax2.plot(xaxis, mme, 'k--', linewidth=lw)
    ax2.plot(xaxis, mrme, 'r--', linewidth=lw)
    ax2.plot(xaxis, mtme, 'g--', linewidth=lw)
    ax2.plot(xaxis, mpme, 'b--', linewidth=lw)

# Title and axis label
ax2.set_title('mean energy')
# Put the y-label on the middle plot
ax2.set_ylabel(r'$\rm{energy\ density\ (erg}\ cm^{-3})$')

# Third plot: fluctuating energy
ax3.plot(xaxis, fke, 'k', linewidth=lw)
ax3.plot(xaxis, frke, 'r', linewidth=lw)
ax3.plot(xaxis, ftke, 'g', linewidth=lw)
ax3.plot(xaxis, fpke, 'b', linewidth=lw)

# If magnetic, plot magnetic energies!
if magnetism:
    ax3.plot(xaxis, fme, 'k--', linewidth=lw)
    ax3.plot(xaxis, frme, 'r--', linewidth=lw)
    ax3.plot(xaxis, ftme, 'g--', linewidth=lw)
    ax3.plot(xaxis, fpme, 'b--', linewidth=lw)

# title and x-axis label
ax3.set_title('fluctuating energy')

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
