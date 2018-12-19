# Author: Loren Matilsky
# Date created: 06/08/2017
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from subprocess import call
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
trace_G_Avgs_file = get_widest_range_file(datadir, 'trace_G_Avgs')

# Set defaults
xiter = False
notfrom0 = False
magnetism = False
ylog = False

# Get command-line arguments
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-xiter'): # plot w.r.t. iterations
        xiter = True
    elif (arg == '-usefile'):
        trace_G_Avgs_file = args[i+1]
        trace_G_Avgs_file = trace_G_Avgs_file.split('/')[-1]
    elif (arg == '-notfrom0'):
        notfrom0 = True
    elif (arg == '-mag'):
        magnetism = True
    elif (arg == '-log'):
        ylog = True

# Tag the plot by whether or not the x axis is in "time" or "iteration"
if (xiter):
    tag = '_xiter'
else:
    tag = '_xtime'

# Read in the KE data (dictionary form)
print ('Reading in KE components data from ' + datadir + trace_G_Avgs_file + ' ...')
di = np.load(datadir + trace_G_Avgs_file).item()
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
savename = dirname_stripped + '_etrace_' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.pdf'

ke = vals[lut[401]]
rke = vals[lut[402]]
tke = vals[lut[403]]
pke = vals[lut[404]]

mke = vals[lut[405]]
mrke = vals[lut[406]]
mtke = vals[lut[407]]
mpke = vals[lut[408]]

fke = vals[lut[409]]
frke = vals[lut[410]]
ftke = vals[lut[411]]
fpke = vals[lut[412]]

# Get the magnetic energies if they are available
if (magnetism):
    me = vals[lut[1101]]
    rme = vals[lut[1102]]
    tme = vals[lut[1103]]
    pme = vals[lut[1104]]

    mme = vals[lut[1105]]
    mrme = vals[lut[1106]]
    mtme = vals[lut[1107]]
    mpme = vals[lut[1108]]

    fme = vals[lut[1109]]
    frme = vals[lut[1110]]
    ftme = vals[lut[1111]]
    fpme = vals[lut[1112]]

# Get global min/max vals
if magnetism:
    mmax = np.max((np.max(ke), np.max(me)))
    mmin = np.min((np.min(mrke), np.min(mtke), np.min(mpke), np.min(frke), np.min(ftke), np.min(fpke),\
     np.min(mrme), np.min(mtme), np.min(mpme), np.min(frke), np.min(ftme), np.min(fpme)))
else:
    mmax = np.max(ke)
    mmin = np.min((np.min(mrke), np.min(mtke), np.min(mpke), np.min(frke), np.min(ftke), np.min(fpke)))   

if (not xiter):
    xaxis = times/tnorm
else:
    xaxis = iters

if (notfrom0):
    x_min = np.min(xaxis)
else:
    x_min = 0

# create figure with  3 panels in a row (total, mean and fluctuating energy)
fig, axs = plt.subplots(3, 1, figsize=(5, 10), sharex=True, sharey=True)
ax1 = axs[0]; ax2 = axs[1]; ax3 = axs[2]

# first plot: total kinetic energy trace      
ax1.set_title('total energy')
ax1.plot(xaxis, ke, 'k', label=r'$\rm{KE_{tot}}$')
ax1.plot(xaxis, rke, 'r', label=r'$\rm{KE}_r$')
ax1.plot(xaxis, tke, 'g', label=r'$\rm{KE}_\theta$')
ax1.plot(xaxis, pke, 'b', label=r'$\rm{KE}_\phi$')

# If magnetic, plot magnetic energies!
if magnetism:
    ax1.plot(xaxis, me, 'k--', label=r'$\rm{ME_{tot}}$')
    ax1.plot(xaxis, rme, 'r--', label=r'$\rm{ME}_r$')
    ax1.plot(xaxis, tme, 'g--', label=r'$\rm{ME}_\theta$')
    ax1.plot(xaxis, pme, 'b--', label=r'$\rm{ME}_\phi$')

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
ax2.plot(xaxis, mke, 'k')
ax2.plot(xaxis, mrke, 'r')
ax2.plot(xaxis, mtke, 'g')
ax2.plot(xaxis, mpke, 'b')

# If magnetic, plot magnetic energies!
if magnetism:
    ax2.plot(xaxis, mme, 'k--')
    ax2.plot(xaxis, mrme, 'r--')
    ax2.plot(xaxis, mtme, 'g--')
    ax2.plot(xaxis, mpme, 'b--')

# Title and axis label
ax2.set_title('mean energy')
# Put the y-label on the middle plot
ax2.set_ylabel(r'$\rm{energy\ density\ (erg}\ cm^{-3})$')

# Third plot: fluctuating energy
ax3.plot(xaxis, fke, 'k')
ax3.plot(xaxis, frke, 'r')
ax3.plot(xaxis, ftke, 'g')
ax3.plot(xaxis, fpke, 'b')

# If magnetic, plot magnetic energies!
if magnetism:
    ax3.plot(xaxis, fme, 'k--')
    ax3.plot(xaxis, frme, 'r--')
    ax3.plot(xaxis, ftme, 'g--')
    ax3.plot(xaxis, fpme, 'b--')

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
plt.subplots_adjust(left=0.15, bottom=0.08, top=0.95, wspace=0.4)

# Save the plot
print ('Saving the etrace plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename)

# Show the plot
plt.show()