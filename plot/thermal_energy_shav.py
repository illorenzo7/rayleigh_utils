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
        get_iters_from_file, get_dict, rsun
from get_parameter import get_parameter
from rayleigh_diagnostics import ReferenceState, GridInfo
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

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = None
entropy_equation = False

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-usefile':
        Shell_Avgs_file = args[i+1]
        Shell_Avgs_file = Shell_Avgs_file.split('/')[-1]
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for rval_str in rvals_str:
            rvals.append(float(rval_str))
    elif arg == '-s':
        entropy_equation = True

lw = 1. # regular lines
#lw = 1.5 # Bit thicker lines

# Read in the flux data
print ('Getting heating terms from ' + datadir + Shell_Avgs_file + ' ...')
di = get_dict(datadir + Shell_Avgs_file)
vals = di['vals']
lut = di['lut']
nq = di['nq']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
nr = di['nr']

# Get the rho*T
try:
    ref = ReferenceState(dirname + '/reference')
    rhot = ref.density*ref.temperature
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    rhot = eq.functions[0]*eq.functions[3]

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

# Make the plot name, labelling the first/last iterations we average over
if entropy_equation:
    basename = '_entropy_equation_shav_'
else:
    basename = '_thermal_energy_shav_'

savename = dirname_stripped + basename +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

advec_tot = -vals[:, lut[1401]]
advec_fluc = -vals[:, lut[1402]]
advec_mean = advec_tot - advec_fluc
advec_vr = -vals[:, lut[1406]]
cond_heating = vals[:, lut[1421]]
int_heating = vals[:, lut[1434]]
visc_heating = vals[:, lut[1435]]
tot_heating = advec_tot + cond_heating + int_heating + visc_heating
if magnetism:
    joule_heating_tot = vals[:, lut[1436]]
    joule_heating_fluc = vals[:, lut[1437]]
    joule_heating_mean = joule_heating_tot - joule_heating_fluc
    tot_heating += joule_heating_tot

# Compute the INTEGRATED total heating
gi = GridInfo(dirname + '/grid_info')
rw = gi.rweights
shell_volume = 4.0*np.pi/3.0*(np.max(rr)**3.0 - np.min(rr)**3.0)
tot_heating_integrated = shell_volume*np.sum(tot_heating*rw)

if entropy_equation:
    advec_tot /= rhot
    advec_fluc /= rhot
    advec_mean /= rhot
    advec_vr /= rhot
    cond_heating /= rhot
    int_heating /= rhot
    visc_heating /= rhot
    tot_heating /= rhot
    if magnetism:
        joule_heating_tot /= rhot
        joule_heating_fluc /= rhot
        joule_heating_mean /= rhot

# Create the plot

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

plt.plot(rr_n, advec_tot, 'm', label='advection tot', linewidth=lw)
plt.plot(rr_n, advec_fluc, 'm--', label='advection fluc', linewidth=lw)
plt.plot(rr_n, advec_mean, 'm:', label='advection mean', linewidth=lw)
plt.plot(rr_n, cond_heating, 'r', label='conductive heating', linewidth=lw)
plt.plot(rr_n, int_heating, 'g', label='internal heating', linewidth=lw)
plt.plot(rr_n, visc_heating, 'c', label='viscous heating', linewidth=lw)
if magnetism:
    plt.plot(rr_n, joule_heating_tot, 'b', label='Joule heating tot',\
            linewidth=lw)
    plt.plot(rr_n, joule_heating_fluc, 'b--', label='Joule heating fluc',\
            linewidth=lw)
    plt.plot(rr_n, joule_heating_mean, 'b:', label='Joule heating mean',\
            linewidth=lw)
plt.plot(rr_n, tot_heating, 'k', label='total heating')

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
    ymin = min(np.min(advec_tot), np.min(advec_fluc),\
            np.min(advec_mean), np.min(cond_heating), np.min(int_heating),\
            np.min(visc_heating), np.min(tot_heating)),
    ymax = max(np.max(advec_tot), np.max(advec_fluc),\
            np.max(advec_mean), np.max(cond_heating), np.max(int_heating),\
            np.max(visc_heating), np.max(tot_heating))
    if magnetism:
        ymin = min(ymin, np.min(joule_heating_tot),\
                np.min(joule_heating_fluc), np.min(joule_heating_mean))
        ymax = max(ymax, np.max(joule_heating_tot),\
                np.max(joule_heating_fluc), np.max(joule_heating_mean))

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

if entropy_equation:
    plt.ylabel(r'$\partial S/\partial t\ \rm{(erg\ g^{-1}\ K^{-1}\ s^{-1})}$',\
            fontsize=12, **csfont)
else:
    plt.ylabel('heating (' + r'$\rm{erg\ cm^{-3}\ s^{-1}}$' + ')',\
            fontsize=12, **csfont)

# Make title
if entropy_equation:
    basetitle = 'entropy eqn., '
else:
    basetitle = 'thermal energy eqn., ' 

lum = 3.846e33

title = dirname_stripped + '\n' + basetitle +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) +\
          '\nintegrated total heating: %1.3e erg/cm^3/s\n = %1.3e lsun'\
          %(tot_heating_integrated/shell_volume, tot_heating_integrated/lum)
plt.title(title, **csfont)

# Create a see-through legend
plt.legend(loc='lower left', shadow=True, ncol=2, fontsize=8)

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the energy eqn. plot at ' + plotdir + savename + ' ...')
#plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
