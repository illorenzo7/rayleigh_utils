# Plot location of Chebyshev collocation points along with various radial
# locations (see how many grid points are inside different zones)
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from compute_grid_info import compute_grid_info
from common import *
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

rnorm = None
rval = None
ncheby = None
domain_bounds = None
nt = None

fname = 'grid_info'
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        fname = args[i+1]
    elif arg == '-rval':
        rval = float(args[i+1])
    elif arg == '-nr':
        ncheby = (int(args[i+1]),)
    elif arg == '-nt':
        nt = int(args[i+1])
    elif arg == '-ncheby':
        ncheby_str = args[i+1].split()
        ncheby = []
        for val_str in ncheby_str:
            ncheby.append(int(val_str))
        ncheby = tuple(ncheby)
    elif arg == '-rminmax':
        domain_bounds = (float(args[i+1]), float(args[i+2]))
    elif arg == '-dombounds':
        dombounds_str = args[i+1].split()
        domain_bounds = []
        for val_str in dombounds_str:
            domain_bounds.append(float(val_str))
        domain_bounds = tuple(domain_bounds)
    elif arg == '-rnorm':
        rnorm = float(args[i+1])

# Get relevant info from main_input file
if nt is None:
    print("Getting nt from main_input")
    nt = get_parameter(dirname, 'n_theta')
else:
    print("Using user-specified nt")
print ("nt = %i" %nt)

#use_extrema = get_parameter(dirname, 'use_extrema')
use_extrema = False # also, dummy; may change this if I ever start using
    # the "correct" Chebyshev weights like Connor

if ncheby is None:
    print("Getting ncheby from main_input")
    ncheby, dummy = get_domain_bounds(dirname)
else:
    print("Using user-specified ncheby")
print_tuple(ncheby, "%i", prepend="ncheby = ")

if domain_bounds is None:
    print("Getting domain_bounds from main_input")
    dummy, domain_bounds = get_domain_bounds(dirname)
else:
    print("Using user-specified domain_bounds")
print_tuple(domain_bounds, "%1.3e", prepend="domain_bounds = ")

nr, nt, nphi, rr, rw, tt, cost, sint, tw, phi, dphi =\
        compute_grid_info(domain_bounds, ncheby, nt,\
        use_extrema=use_extrema)

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rnorm = rsun
if rval is None:
    rval = np.max(rr)
ir0 = np.argmin(np.abs(rr - rval))
rr /= rnorm                                           
r0 = rr[ir0]
rval /= rnorm
dr = np.zeros(nr)
dr[:-1] = rr[:-1] - rr[1:]
dr[-1] = dr[-2]

# Generate the plot (1 row of 3 subplots)
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(12,4))

# Plot delta_r
size = 0.5
axs[0].scatter(rr, dr, color='k', s=size)
# Make the "rval" location red and big
axs[0].scatter(r0, dr[ir0], color='r', s=3*size)

# Plot r0 * delta_theta
tt_lat = 180./np.pi*(np.pi/2. - tt)
dt = np.zeros(nt)
dt[:-1] = tt[:-1] - tt[1:]
dt[-1] = dt[-2]
axs[1].scatter(tt_lat, r0*dt, color='k', s=size)

# Plot r0 * sin(theta) * dphi
dphi = 2*np.pi/(2*nt)
axs[2].scatter(tt_lat, r0*sint*dphi, color='k', s=size)

# Set x-limits for axs[0]
plt.sca(axs[0])
ri, ro = np.min(rr), np.max(rr)
Dr = ro - ri
plt.xlim(ri - 0.01*Dr, ro + 0.01*Dr)
plt.xlabel(r'$r/({\rm{rnorm}} = %.1e\ \rm{cm})$' %rnorm, fontsize=12, **csfont)
plt.ylabel(r'$\delta r/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# ... for axs[1]
plt.sca(axs[1])
plt.xlim(-90., 90.)
plt.xlabel('latitude', fontsize=12, **csfont)
plt.ylabel(r'$r_0\delta \theta/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# ... for axs[2]
plt.sca(axs[2])
plt.xlim(-90., 90.)
plt.xlabel('latitude', fontsize=12, **csfont)
plt.ylabel(r'$r_0\sin\theta\delta \phi/\rm{rnorm}$')
# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Make title(s)
axs[0].set_title(dirname_stripped + '\n' + 'Grid Resolution ', **csfont)
axs[1].set_title((r'$r_0/{\rm{rnorm}}=%.3f$' %r0) + '\n' +\
        (r'$\delta r_0/{\rm{rnorm}}=%1.3e$' %dr[ir0]))
axs[2].set_title((r'$r_0{\rm{max}}(\delta\theta)/{\rm{rnorm}}=%1.3e$' %(r0*np.max(dt))) + '\n' +\
        (r'$r_0{\rm{max}}(\sin\theta\delta\phi)/{\rm{rnorm}}=%1.3e$' %(r0*np.max(sint)*dphi)))

# Display the plot
plt.show()
