###############################################
# Author: Loren Matilsky
# Date created: 12/19/2019
#
# This script plots the volume heating term I think might be "missing"
# from our flux equation

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import numpy as np
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from common import *
from rayleigh_diagnostics import GridInfo, Meridional_Slices

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

radatadir = dirname + '/Meridional_Slices/'
# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)


print ("Reading grid_info")
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights
tw = gi.tweights
rr = gi.radius
shell_volume = 4.0*np.pi/3.0*(np.max(rr)**3.0 - np.min(rr)**3.0)

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = []

# Get rho*T ds dr
eq = get_eq(dirname)
rhotdsdr = eq.density*eq.temperature*eq.dsdr

# Defaults and CLAs
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-depths':
        rvals = []
        strings = args[i+1].split()
        for st in strings:
            rval = ro - float(st)*d
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

if nargs == 0:
    index = -1 # By default plot the last Mer slice
else:
    index, dummy = get_desired_range(int_file_list, args)

print ("Reading Meridional_Slices/" + file_list[index])
mer = Meridional_Slices(radatadir + file_list[index], '')

lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Get the v_r*S product, averaging over the Meridional Slice longitudes
vrs = np.mean(mer.vals[:, :, :, mer.lut[1], 0]*\
        mer.vals[:, :, :, mer.lut[501], 0], axis=0)

# Average over the sphere
nt, nr = np.shape(vrs)
vrs = np.sum(vrs*tw.reshape((nt, 1)), axis=0)

cp = 3.5e8
missing_term = rhotdsdr/cp*vrs

tot_missing = shell_volume*np.sum(missing_term*rw)

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

plt.plot(rr_n, missing_term, 'm', label='missing heat', linewidth=lw)

# Get the y-axis in scientific notation
plt.ticklabel_format(useMathText=True, axis='y', scilimits=(0,0))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits

if minmax is None:
    ymin = np.min(missing_term)    
    ymax = np.max(missing_term)
    delta_y = ymax - ymin
    ybuffer = 0.1*delta_y
    minmax = ymin - 3*ybuffer, ymax + ybuffer
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

plt.ylabel('erg/s/cm^3', fontsize=12, **csfont)

# Make title
basetitle = r'$\overline{\rho}\overline{T}(S/c_p)v_r(d\overline{S}/dr)$' 

lum = 3.846e33

title = dirname_stripped + '\n' + file_list[index] +\
    '\nintegrated: %1.3e erg/cm^3/s\n = %1.3e lsun'\
          %(tot_missing/shell_volume, tot_missing/lum)
plt.title(title, **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=2, fontsize=8)

# Last command
plt.tight_layout()

# Save the plot
savename = dirname_stripped + '_missing_heat_term_vs_r_' +\
        file_list[index] + '.png'
print ('Saving the missing heating plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
