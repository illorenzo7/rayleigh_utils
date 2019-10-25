###############################################
# Author: Loren Matilsky
# Date created: 06/26/2019
#
# This script plots various radial profiles associated with a CZ-RZ
# system with tachocline boundary layer.
# Uses data from the Shell_Avgs, reference state, and 
# transport coefficients

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
        get_iters_from_file, get_dict, rsun, lsun
from get_parameter import get_parameter
from rayleigh_diagnostics import ReferenceState, TransportCoeffs
from reference_tools import equation_coefficients
# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Find the Shell_Avgs file(s) in the data directory. If there are multiple, by
# default choose the one with widest range in the average
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = None # user can specify radii to mark by vertical lines

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-usefile'):
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif (arg == '-minmax'):
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-rnorm'):
        rnorm = float(args[i+1])
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))

#Create the plot
lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in Shell_Avgs data
print ('Getting data from ' + datadir + Shell_Avgs_file + ' ...')
di_sh = get_dict(datadir + Shell_Avgs_file)
vals_sh = di_sh['vals']
lut_sh = di_sh['lut']
rr = di_sh['rr']
eflux = vals_sh[:, lut_sh[1455]]
iter1, iter2 = get_iters_from_file(Shell_Avgs_file)

# get convective velocities
vr = np.sqrt(vals_sh[:, lut_sh[422]])
vh = np.sqrt(vals_sh[:, lut_sh[423]] + vals_sh[:, lut_sh[424]])

print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di_az = get_dict(datadir + AZ_Avgs_file)
vals_az = di_az['vals']
lut_az = di_az['lut']
xx = di_az['xx'] # moment arm r*sintheta
nr, nt = di_az['nr'], di_az['nt']

# Compute rotation rate
om = vals_az[:, :, lut_az[3]]/xx
om_eq = om[nt//2, :]

# Get dsdr
try:
    ref = ReferenceState(dirname + '/reference')
    dsdr = ref.dsdr
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    dsdr = eq.functions[13]

# Get nu profile
try:
    t = TransportCoeffs(dirname + '/transport')
    nu = t.nu
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    # nu(r) = c_5 * f_3
    nu = eq.constants[4]*eq.functions[2]

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_tach_rprofiles_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Create the plot; start with plotting all the energy fluxes

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
# Set the x limits
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

plt.plot(rr_n, dsdr/np.max(dsdr), label=r'$ds/dr$',\
        linewidth=lw)
plt.plot(rr_n, nu/np.max(nu), label=r'$\nu(r)$',\
        linewidth=lw)
plt.plot(rr_n, eflux/np.max(eflux), label=r'$\rm{F}_{enth}$',\
        linewidth=lw)
plt.plot(rr_n, om_eq/np.max(om_eq), label=r'$\Omega_{\rm{eq}}$',\
        linewidth=lw)
plt.plot(rr_n, vr/np.max(vr), label=r'$v_r^\prime(r)$',\
        linewidth=lw)
plt.plot(rr_n, vh/np.max(vh), label=r'$v_h^\prime(r)$',\
        linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 
if minmax is None:
    minmax = plt.gca().get_ylim()
plt.ylim(minmax)

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

# Make title
plt.title(dirname_stripped + '\n' + 'tach. radial profiles ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)

# Create a see-through legend
plt.legend(shadow=True, ncol=3, fontsize=10)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the tach_rprofiles plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
