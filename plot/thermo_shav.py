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
from rayleigh_diagnostics import TransportCoeffs, ReferenceState
from reference_tools import equation_coefficients
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


# Read reference state
prs_spec_heat = 3.5e8
try:
    ref = ReferenceState(dirname + '/reference', '')
    ref_rho = ref.density
    ref_prs = ref.pressure
    ref_temp = ref.temperature
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    ref_rho = eq.functions[0]
    ref_temp = eq.functions[3]
    gam = 5.0/3.0
    gas_R = (gam - 1.0)/gam*prs_spec_heat
    ref_prs = gas_R*ref_rho*ref_temp

# Read in entropy and pressure, nond
entropy = vals[:, lut[501]]/prs_spec_heat
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

# Compute maximum/minimum Reynolds numbers (ignore the upper/lower 5%
# of the shell to avoid extreme values associated with boundary conditions
rr_depth = (ro - rr)/shell_depth
ir1, ir2 = np.argmin(np.abs(rr_depth - 0.05)),\
        np.argmin(np.abs(rr_depth - 0.95))
my_min = min(np.min(entropy[ir1:ir2]), np.min(prs[ir1:ir2]),\
        np.min(temp[ir1:ir2]))
my_max = np.max(rho[ir1:ir2])

if minmax is None:
    if logscale:
        ratio = my_max/my_min
        ybuffer = 0.2*ratio
        ymin = my_min/ybuffer
        ymax = my_max*ybuffer
    else:
        difference = my_max - my_min
        ybuffer = 0.2*difference
        ymin, ymax = my_min - ybuffer, my_max + ybuffer
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
plt.title(dirname_stripped + '\n' +'Thermodynamic Perturbations ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend()

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_thermo_shellav_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()