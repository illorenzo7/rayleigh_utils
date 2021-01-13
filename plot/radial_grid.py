# Plot location of Chebyshev collocation points along with various radial
# locations (see how many grid points are inside different zones)
# revised 11/30/2020
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
from get_parameter import get_parameter
from get_domain_bounds import get_domain_bounds
from common import *
dirname = sys.argv[1]
rnorm = None
rvals = []
ncheby = None
domain_bounds = None
xminmax = None

fname = 'grid_info'
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        fname = args[i+1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-nr':
        ncheby = (int(args[i+1]),)
    elif arg == '-ncheby':
        ncheby_str = args[i+1].split()
        ncheby = []
        for val_str in ncheby_str:
            ncheby.append(int(val_str))
        ncheby = tuple(ncheby)
    elif arg == '-xminmax':
        xminmax = np.float(args[i+1]), np.float(args[i+2])
    elif arg == '-dombounds':
        dombounds_str = args[i+1].split()
        domain_bounds = []
        for val_str in dombounds_str:
            domain_bounds.append(float(val_str))
        domain_bounds = tuple(domain_bounds)
    elif arg == '-rnorm':
        rnorm = float(args[i+1])

# Get relevant info from main_input file
#nt = get_parameter(dirname, 'n_theta')
nt = 96 # dummy variable
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
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plot the Chebyshev collocation points
zero = np.zeros(nr)
plt.scatter(rr_n, zero, color='k', s=.3)


# Set the y-limits (to arbitrary (-1,1))
ymin, ymax = -1, 1
plt.ylim(ymin, ymax)
yvals = np.linspace(ymin, ymax, 100)

# Plot desired radial locations, if any
# Also set x-limits
if not rvals is None:
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
        plt.plot(rval_n + np.zeros(100), yvals, 'k')
    min_rval = np.min(rvals)
    max_rval = np.max(rvals)
    rbuffer = max_rval - min_rval
    frac = 0.2
    xmin, xmax = min_rval - frac*rbuffer, max_rval + frac*rbuffer
else:
    xmin, xmax = np.min(rr), np.max(rr)

# Set x-limits
if rnorm is None:
    xmin /= rsun
    xmax /= rsun
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    xmin /= rnorm
    xmax /= rnorm
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

# Set the x limits
if xminmax is None:
    plt.xlim(xmin, xmax)
else:
    plt.xlim(xminmax[0], xminmax[1])

# Display the plot
plt.show()
