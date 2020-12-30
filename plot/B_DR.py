# Author: Loren Matilsky
# Created: 04/08/2020
# This script generates the buoyancy parameter, using large scale shear
# (between equator and 60 deg by default) for the velocity and length scale
# Analayzes Rayleigh run directory indicated by [dirname]. 
# Saves plot in
# [dirname]_Brunt.png

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
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict, rsun, c_P
from get_eq import get_eq

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
rvals = [] # user can specify radii to mark by vertical lines
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
        for rval_str in rvals_str:
            rvals.append(float(rval_str))

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

# Get reference state info
eq = get_eq(dirname)
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

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

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

# Compute Buoyancy parameter
B = N2*length_scale**2/Dvp**2

# Plot Re vs radius
ax.plot(rr_n, B)

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel(r'${\rm{B}} = N^2\ (r\Delta\theta)^2 (\Delta\langle v_\phi\rangle)^{-2}$',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute max/min of N 
val_min = np.min(B)
val_max = np.max(B)

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
plt.title(dirname_stripped + '\n' +'Buoyancy Parameter for DR', **csfont)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_buoyancy_parameter' + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
