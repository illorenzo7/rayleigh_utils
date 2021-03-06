###############################################
# Author: Loren Matilsky
# Date created: 06/26/2019
#
# This script plots various radial profiles associated with a CZ-RZ
# system with tachocline boundary layer.
# Uses data from the Shell_Avgs

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = [] # user can specify radii to mark by vertical lines

plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if (arg == '-usefile'):
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-rnorm'):
        rnorm = float(args[i+1])
    elif arg == '-depths':
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
            rvals.append(rval)
    elif arg == '-depthscz':
        rm = domain_bounds[1]
        dcz = ro - rm
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*dcz
            rvals.append(rval)
    elif arg == '-depthsrz':
        rm = domain_bounds[1]
        drz = rm - ri
        strings = args[i+1].split()
        for st in strings:
            rval = rm - float(st)*drz
            rvals.append(rval)
    elif arg == '-rvals':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)*rsun
            rvals.append(rval)
    elif arg == '-rvalscm':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = float(st)
            rvals.append(rval)

#Create the plot
lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in Shell_Avgs data
print ('Getting data from ' + datadir + Shell_Avgs_file + ' ...')
di_sh = get_dict(datadir + Shell_Avgs_file)
vals_sh = di_sh['vals']
lut_sh = di_sh['lut']
rr = di_sh['rr']
eflux = vals_sh[:, lut_sh[1455]]
iter1, iter2 = get_iters_from_file(Shell_Avgs_file)

# get convective velocities
vr = np.sqrt(vals_sh[:, lut_sh[422]])
vh = np.sqrt(vals_sh[:, lut_sh[423]] + vals_sh[:, lut_sh[424]])

print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di_az = get_dict(datadir + AZ_Avgs_file)
vals_az = di_az['vals']
lut_az = di_az['lut']
xx = di_az['xx'] # moment arm r*sintheta
nr, nt = di_az['nr'], di_az['nt']

# Compute rotation rate
om = vals_az[:, :, lut_az[3]]/xx
om_eq = om[nt//2, :]

# Get dsdr
eq = get_eq(dirname)
dsdr = eq.dsdr

# Get nu profile
nu = eq.nu

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_tach_rprofiles_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Create the plot; start with plotting all the energy fluxes

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
# Set the x limits
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

plt.plot(rr_n, dsdr/np.max(dsdr), label=r'$ds/dr$',\
        linewidth=lw)
plt.plot(rr_n, nu/np.max(nu), label=r'$\nu(r)$',\
        linewidth=lw)
plt.plot(rr_n, eflux/np.max(eflux), label=r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(rr_n, om_eq/np.max(om_eq), label=r'$\Omega_{\rm{eq}}$',\
        linewidth=lw)
plt.plot(rr_n, vr/np.max(vr), label=r'$v_r^\prime(r)$',\
        linewidth=lw)
plt.plot(rr_n, vh/np.max(vh), label=r'$v_h^\prime(r)$',\
        linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 
if minmax is None:
    minmax = plt.gca().get_ylim()
plt.ylim(minmax)

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

# Mark radii if desired
if not rvals is None:
    yvals = np.linspace(minmax[0], minmax[1], 100)
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
#        plt.ylim(ymin, ymax)
        plt.plot(rval_n + np.zeros(100), yvals, 'k--')

# Make title
plt.title(dirname_stripped + '\n' + 'tach. radial profiles ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)

# Create a see-through legend
plt.legend(shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the tach_rprofiles plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
