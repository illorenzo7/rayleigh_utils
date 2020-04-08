# Author: Loren Matilsky
# Created: 04/08/2020
# This script generates the differential-rotation Reynolds number (Re)
# vs radius, using the velocity and length scale associated with the 
# average (horizontal) shear between the equator and 60 deg
# Analyzes the Rayleigh run directory indicated by [dirname].
# To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_Re_DR_[first iter]_[last iter].png

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from get_eq import get_eq
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = None # user can specify radii to mark by vertical lines
rvals = None
tag = ''
latrange = 0., 60. # By default compute average shear between 
    # equator and 60 deg
the_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-latrange':
        latrange = float(args[i+1]), float(args[i+2])
    elif arg == '-usefile':
        the_file = args[i+1].split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for j in range(len(rvals_str)):
            rvals.append(float(rvals_str[j]))

# Read in vavg data
print ('Reading AZ_Avgs data from ' + datadir + the_file + ' ...')
di = get_dict(datadir + the_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
tt_lat = di['tt_lat']
rr_depth = di['rr_depth']

# Mean rotation velocity amplitudes
mean_vp = vals[:, :, lut[3]]

# Get indices associated with latrange
it1 = np.argmin(np.abs(tt_lat - latrange[0]))
it2 = np.argmin(np.abs(tt_lat - latrange[1]))

# Compute differential velocity over latrange
Dvp = np.abs(mean_vp[it2] - mean_vp[it1])

# ...and the lengthscale
Dtheta = tt[it1] - tt[it2]
length_scale = Dtheta*rr

# Get molecular diffusivity 
eq = get_eq(dirname)
nu = eq.nu

# Compute differential-rotation Reynolds number
Re = Dvp*length_scale/nu

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plot Re vs radius
ax.plot(rr_n, Re)

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel(r'${\rm{Re_{DR}}} = \Delta\langle v_\phi\rangle r\Delta\theta \nu^{-1}$',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
# of the shell to avoid extreme values associated with boundary conditions
ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
        np.argmin(np.abs(rr_depth - 0.95))
Re_min = np.min(Re[ir1:ir2])
Re_max = np.max(Re[ir1:ir2])

if minmax is None:
    if logscale:
        ratio = Re_max/Re_min
        ybuffer = 0.2*ratio
        ymin = Re_min/ybuffer
        ymax = Re_max*ybuffer
    else:
        difference = Re_max - Re_min
        ybuffer = 0.2*difference
        ymin, ymax = Re_min - ybuffer, Re_max + ybuffer
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
plt.title(dirname_stripped + '\n' +'DR Reynolds number, ' +\
        ('%.1f deg to %.1f deg' %(tt_lat[it1], tt_lat[it2])) +\
          '\n' + str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_Re_DR_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
