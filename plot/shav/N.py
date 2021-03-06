# Author: Loren Matilsky
# Created: 04/08/2020
# This script generates the Brunt-Vaisala frequency N
# for the Rayleigh run directory indicated by [dirname]. 
# Saves plot in
# [dirname]_Brunt.png

# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with plots, make the plotting directory if it doesn't
# already exist    

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = [] # user can specify radii to mark by vertical lines
tag = ''

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
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

# Get reference state info
eq = get_eq(dirname)
rr = eq.radius
g = eq.gravity
dsdr = eq.dsdr
N2 = g/c_P*dsdr

# Derivative grid info
nr = len(rr)
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value
                                               
# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plot Re vs radius
ax.plot(rr_n, N2)

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12)

plt.ylabel(r'${\rm{N^2}} = g\ (d\overline{S}/dr)/c_p\ [s^{-2}]$',\
        fontsize=12)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute max/min of N 
val_min = np.min(N2)
val_max = np.max(N2)

if minmax is None:
    if logscale:
        ratio = val_max/val_min
        ybuffer = 0.2*ratio
        ymin = val_min/ybuffer
        ymax = val_max*ybuffer
    else:
        difference = val_max - val_min
        ybuffer = 0.2*difference
        ymin, ymax = val_min - ybuffer, val_max + ybuffer
else:
    ymin, ymax = minmax

plt.ylim((ymin, ymax))
if logscale:
    plt.yscale('log')

xvals = np.linspace(xmin, xmax, 100)
yvals = np.linspace(ymin, ymax, 100)

# Mark radii if desired
if not rvals is None:
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
        plt.plot(rval_n + np.zeros(100), yvals, 'k--')

# Create a title    
# Put max frequency (and time-step limit) in title
N_max = np.sqrt(val_max)
dtlimit = 1./N_max
plt.title(dirname_stripped + '\n' +'Brunt-Vaisala Frequency' +\
        '\n' + r'$N_{\rm{max}} = %1.1e\ s^{-1},\ (\Delta t)_N = %1.1e\ s$'\
        %(N_max, dtlimit))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()
plotdir = my_mkdir(clas0['plotdir'])

savefile = plotdir + dirname_stripped + '_Brunt' + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
