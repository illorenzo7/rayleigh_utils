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
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-usefile':
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
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

# Read in vavg data
print ('Reading Shell_Avgs data from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']

# Derivative grid info
nr = len(rr)
ri, ro = np.min(rr), np.max(rr)
shell_depth = ro - ri

# Convective velocity amplitudes, get these from KE
frke = vals[:, lut[410]]
ftke = vals[:, lut[411]]
fpke = vals[:, lut[412]]

vsq_r = frke/rho
vsq_t = ftke/rho
vsq_p = fpke/rho
vsq = vsq_r + vsq_t + vsq_p

amp_v = np.sqrt(vsq)/100.
amp_vr = np.sqrt(vsq_r)/100.
amp_vt = np.sqrt(vsq_t)/100.
amp_vp = np.sqrt(vsq_p)/100.

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
ax.plot(rr_n, amp_v, label= r'$(v^\prime)_{\rm{rms}}$')
ax.plot(rr_n, amp_vr, label= r'$(v^\prime_r)_{\rm{rms}}$')
ax.plot(rr_n, amp_vt, label = r'$(v^\prime_\theta)_{\rm{rms}}$')
ax.plot(rr_n, amp_vp, label = r'$(v^\prime_\phi)_{\rm{rms}}$')

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel('velocity (m/s)',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
# of the shell to avoid extreme values associated with boundary conditions
rr_depth = (ro - rr)/shell_depth
ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
        np.argmin(np.abs(rr_depth - 0.95))
amp_min = min(np.min(amp_vr[ir1:ir2]), np.min(amp_vt[ir1:ir2]),\
        np.min(amp_vp[ir1:ir2]))
amp_max = np.max(amp_v[ir1:ir2])
print (amp_min, amp_max)
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
plt.title(dirname_stripped + '\n' +'Velocity Amplitudes, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_vamp_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
