# Author: Loren Matilsky
# Created: 08/14/2019
# This script generates the spherically averaged magnetic 
# Reynolds number <Re_m> for the mean flows.
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
from rayleigh_diagnostics import TransportCoeffs, ReferenceState
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
tag = ''
use_hrho = False
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
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
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))

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

# Mean velocity amplitudes...
vsq_r = vals[:, lut[414]] - vals[:, lut[422]]
vsq_t = vals[:, lut[415]] - vals[:, lut[423]]
vsq_p = vals[:, lut[416]] - vals[:, lut[424]]

vsq = vsq_r + vsq_t + vsq_p

# Get molecular diffusivity from 'transport' file
t = TransportCoeffs(dirname + '/transport')
eta = t.eta

# Compute convective Reynolds number
if use_hrho:
    ref = ReferenceState(dirname + '/reference')
    hrho = -1./ref.dlnrho
    Rm = np.sqrt(vsq)*hrho/nu
    Rm_r = np.sqrt(vsq_r)*hrho/eta
    Rm_t = np.sqrt(vsq_t)*hrho/eta
    Rm_p = np.sqrt(vsq_p)*hrho/eta
else:
    Rm = np.sqrt(vsq)*shell_depth/eta
    Rm_r = np.sqrt(vsq_r)*shell_depth/eta
    Rm_t = np.sqrt(vsq_t)*shell_depth/eta
    Rm_p = np.sqrt(vsq_p)*shell_depth/eta

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
ax.plot(rr_n, Rm, label=r'$\rm{Rm}$')
ax.plot(rr_n, Rm_r, label=r'${\rm{Rm}}_r$')
ax.plot(rr_n, Rm_t, label = r'${\rm{Rm}}_\theta$')
ax.plot(rr_n, Rm_p, label = r'${\rm{Rm}}_\phi$')

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

if use_hrho:
    plt.ylabel(r'${\rm{Rm}} = H_\rho \langle v\rangle_{\rm{rms}} \eta^{-1}$',fontsize=12,\
        **csfont)
else:
    plt.ylabel(r'${\rm{Rm}} = (r_o-r_i)\langle v\rangle_{\rm{rms}}\eta^{-1}$',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

# Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
# of the shell to avoid extreme values associated with boundary conditions
rr_depth = (ro - rr)/shell_depth
ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
        np.argmin(np.abs(rr_depth - 0.95))
Rm_min = min(np.min(Rm_r[ir1:ir2]), np.min(Rm_t[ir1:ir2]),\
        np.min(Rm_p[ir1:ir2]))
Rm_max = np.max(Rm[ir1:ir2])

if minmax is None:
    if logscale:
        ratio = Rm_max/Rm_min
        ybuffer = 0.2*ratio
        ymin = Rm_min/ybuffer
        ymax = Rm_max*ybuffer
    else:
        difference = Rm_max - Rm_min
        ybuffer = 0.2*difference
        ymin, ymax = Rm_min - ybuffer, Rm_max + ybuffer
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
plt.title(dirname_stripped + '\n' +'Mean-Flow Reynolds number, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_Rm_mean_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
