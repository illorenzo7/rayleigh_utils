# plot the Alvfen speed
# Author: Loren Matilsky
# Created: 10/31/2019
# This script generates the spherically averaged convective Reynolds 
# number (Re), plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. 
# Computes contributions from v_r, v_theta, and v_phi separately 
# To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_Re_rslice_[first iter]_[last iter].png

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
from common import *

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# get rho
eq = get_eq(dirname)
rho = eq.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Set defaults
rnorm = None
minmax = None
logscale = False
rvals = [] # user can specify radii to mark by vertical lines
tag = ''
use_hrho = False
the_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-hrho':
        use_hrho = True
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

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Read in vavg data
print ('Reading Shell_Avgs data from ' +  the_file)
di = get_dict( the_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = get_iters_from_file(the_file)

di_grid = get_grid_info(dirname)
rr = di_grid['rr']

# Derivative grid info
nr = len(rr)
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri

# Convective velocity amplitudes, get these from KE
# use convective velocities to ignore whopping DR
frke = vals[:, 0, lut[410]]
ftke = vals[:, 0, lut[411]]
fpke = vals[:, 0, lut[412]]

vsq_r = 2.*frke/rho
vsq_t = 2.*ftke/rho
vsq_p = 2.*fpke/rho
vsq = vsq_r + vsq_t + vsq_p
vamp = np.sqrt(vsq)

# B amplitudes, get these from ME
# use full amplitude, since both lon. mean and fluc 
# are likely relevant (i.e., m = 0, 1, 2, could all play a role...)
rme = vals[:, 0, lut[1102]]
tme = vals[:, 0, lut[1103]]

eightpi = 8.*np.pi

bsq_r = rme*eightpi
bsq_t = tme*eightpi

bsq_pol = bsq_r + bsq_t
v_A = np.sqrt(bsq_pol/4./np.pi/rho)

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
ax.plot(rr_n, v_A)

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel('Alfven speed (cgs)',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
# of the shell to avoid extreme values associated with boundary conditions
rr_depth = (ro - rr)/shell_depth
ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
        np.argmin(np.abs(rr_depth - 0.95))
amp_min = np.min((v_A)[ir1:ir2])
amp_max = np.max((v_A)[ir1:ir2])
if minmax is None:
    fact = 0.2
    if logscale:
        ratio = amp_max/amp_min
        ybuffer = ratio**fact
        ymin = amp_min/ybuffer
        ymax = amp_max*ybuffer
        print (ymin, ymax)
    else:
        difference = amp_max - amp_min
        ybuffer = fact*difference
        ymin, ymax = amp_min - ybuffer, amp_max + ybuffer
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
plt.title(dirname_stripped + '\n' +'v_Alfven (pol.) / v_convection ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_v_alfven_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
