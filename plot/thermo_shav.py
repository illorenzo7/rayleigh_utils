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

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Set defaults
rnorm = None
minmax = None
subtop = False  # if True, subtract the top value of the entropy before\
        # computing the fluctuations
rvals = [] # user can specify radii to mark by vertical lines
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
    elif arg == '-hrho':
        use_hrho = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-sub':
        subtop = True
        print ("subtop = True")
        print ("subtracting top value of entropy from radial profile")

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


# Read reference state
eq = get_eq(dirname)
prs_spec_heat = 3.5e8
ref_rho = eq.density
ref_prs = eq.pressure
ref_temp = eq.temperature

# Read in entropy and pressure, nond
entropy = vals[:, lut[501]]/prs_spec_heat
if subtop:
    entropy -= entropy[0]

prs = vals[:, lut[502]]/ref_prs

# Calculate temp. from EOS
poly_n = 1.5
temp = prs/(poly_n + 1.) + entropy

# Calculate  density from Ideal Gas Law
rho = prs - temp

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

# Plot thermo vars vs radius
ax.plot(rr_n, entropy, label=r'$\langle S\rangle_{\rm{sph}}/c_{\rm{p}}$')
ax.plot(rr_n, prs, label=r'$\langle P\rangle_{\rm{sph}}/\overline{P}(r)$')
ax.plot(rr_n, temp, label=r'$\langle T\rangle_{\rm{sph}}/\overline{T}(r)$')
ax.plot(rr_n, rho, label=r'$\langle \rho\rangle_{\rm{sph}}/\overline{\rho}(r)$')

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

plt.ylabel('thermo. pert.',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))

if minmax is None:
    # Compute maximum/minimum thermo. pert. (ignore the upper/lower 5%
    # of the shell to avoid extreme values 
    rr_depth = (ro - rr)/shell_depth
    #mindepth, maxdepth = 0.05, 0.95
    mindepth, maxdepth = 0.0, 1.0
    ir1, ir2 = np.argmin(np.abs(rr_depth - mindepth)),\
            np.argmin(np.abs(rr_depth - maxdepth))

    mmin = min(np.min(entropy[ir1:ir2+1]), np.min(prs[ir1:ir2+1]),\
            np.min(temp[ir1:ir2+1]), np.min(rho[ir1:ir2+1]))
    mmax = max(np.max(entropy[ir1:ir2+1]), np.max(prs[ir1:ir2+1]),\
            np.max(temp[ir1:ir2+1]), np.max(rho[ir1:ir2+1]))
    difference = mmax - mmin
    ybuffer = 0.2*difference
    ymin, ymax = mmin - ybuffer, mmax + ybuffer
else:
    ymin, ymax = minmax

plt.ylim((ymin, ymax))

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
plt.title(dirname_stripped + '\n' +'Thermodynamic Perturbations ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_thermo_shellav_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile)
plt.savefig(savefile, dpi=300)
plt.show()
