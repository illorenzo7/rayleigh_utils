###############################################
# Author: Loren Matilsky
# Date created: 11/05/2019
#
# This script plots the volume heating terms as functions of
# radius using from the Shell_Avgs data

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import get_widest_range_file, strip_dirname,\
        get_iters_from_file, get_dict, rsun, c_P, rms
from get_parameter import get_parameter
from rayleigh_diagnostics import GridInfo
from get_eq import get_eq

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
the_file = get_widest_range_file(datadir, 'Shell_Avgs')

# Get command-line arguments to adjust the interval of averaging files
xminmax = None
minmax = None
rnorm = None
rvals = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))

lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting l=0 forcing terms from ' + datadir + the_file)
di = get_dict(datadir + the_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
nr = di['nr']

# Get the rho*T
eq = get_eq(dirname)
rhot = eq.density*eq.temperature

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

savename = dirname_stripped + '_hydro_balance_' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

eq = get_eq(dirname)
dlnrho = eq.dlnrho
rho = eq.density
rhog = rho*eq.gravity
rr = eq.radius

#f_prs_deriv = -vals[:, lut[508]]
#f_prs = vals[:, lut[502]]*dlnrho
f_prs = -vals[:, lut[508]] + vals[:, lut[502]]*dlnrho
f_buoy = rhog*vals[:,lut[501]]/c_P
f_uphi = rho*vals[:,lut[404]]/(rho/2.)/rr
f_utheta = rho*vals[:,lut[403]]/(rho/2.)/rr
#f_tot = f_prs_deriv + f_prs + f_buoy + f_uphi + f_utheta
f_tot = f_prs + f_buoy

#print ("rms(f_prs_deriv) = %8.2e" %rms(f_prs_deriv))
print ("rms(f_prs) = %8.2e" %rms(f_prs))
print ("rms(f_buoy) = %8.2e" %rms(f_buoy))
#print ("rms(f_uphi) = %8.2e" %rms(f_uphi))
#print ("rms(f_utheta) = %8.2e" %rms(f_utheta))
print ("rms(f_tot) = %8.2e" %rms(f_tot))

# Create the plot

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

#plt.plot(rr_n, f_prs_deriv, 'b', label='prs deriv', linewidth=lw)
plt.plot(rr_n, f_prs, 'c', label='prs', linewidth=lw)
plt.plot(rr_n, f_buoy, 'r', label='buoy', linewidth=lw)
#plt.plot(rr_n, f_uphi, 'g', label='rho u_phi^2/r', linewidth=lw)
#plt.plot(rr_n, f_utheta, 'm', label='rho u_theta^2/r', linewidth=lw)
plt.plot(rr_n, f_tot, 'k', label='tot', linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
if xminmax is None:
    xminmax = np.min(rr_n), np.max(rr_n)
delta_x = xminmax[1] - xminmax[0]
plt.xlim(xminmax[0], xminmax[1])

# Set the y-limits

if not minmax is None:
    plt.ylim(minmax[0], minmax[1])

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

# Mark radii if desired
if not rvals is None:
    yvals = np.linspace(minmax[0], minmax[1], 100)
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
#        plt.ylim(ymin, ymax)
        plt.plot(rval_n + np.zeros(100), yvals, 'k--')

plt.title("l=0 force balance", **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=2, fontsize=8)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the energy eqn. plot at ' + plotdir + savename)
#plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
