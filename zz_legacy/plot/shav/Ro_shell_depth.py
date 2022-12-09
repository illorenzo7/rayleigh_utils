###############################################
# Author: Loren Matilsky
# Date created: 04/22/2020
#
# This script plots Rossby numbers using the vorticity-based length_scale
# i.e. Ro = sqrt(vort)/(2 Omega_0),
# where vort takes on different components of the vorticity
# functions of radius, using data from the Shell_Avgs/Shell_Spectra

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

# Get the run directory on which to perform the analysis
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Get command-line arguments to adjust the interval of averaging files
minmax = None
rnorm = None
rvals = []
log = False

plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
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
    elif arg == '-log':
        log = True

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'

# Get the rotation rate
Om0 = 2.*np.pi/compute_Prot(dirname)

# Get the lengthscales
di = get_length_scales(dirname)
rr = di['rr']
nr = di['nr']
iter1, iter2 = di['iter1'], di['iter2']
shell_depth = di['shell_depth']

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
rotation = get_parameter(dirname, 'rotation')
if rotation:
    time_unit = compute_Prot(dirname)
    time_label = r'$\rm{P_{rot}}$'
else:
    time_unit = compute_tdt(dirname)
    time_label = r'$\rm{TDT}$'

if plotdir is None:
    plotdir = dirname + '/plots/'
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

# Get the velocity amplitude
the_file = get_widest_range_file(datadir, 'Shell_Avgs')
print ("Ro_vort: Getting convective velocity amplitudes from ", the_file)
di_sh = get_dict(datadir + the_file)
vals = di_sh['vals']
lut = di_sh['lut']
vramp = np.sqrt(vals[:, lut[422]])
vhamp = np.sqrt(vals[:, lut[423]] + vals[:, lut[424]])
vamp = np.sqrt(vals[:, lut[422]] + vals[:, lut[423]] + vals[:, lut[424]])

# Compute the shell-depth Rossby number
Ror_sd = vramp/(2.*Om0*shell_depth)
Roh_sd = vhamp/(2.*Om0*shell_depth)
Ro_sd = vamp/(2.*Om0*shell_depth)

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

vars_to_plot = [Ror_sd, Roh_sd, Ro_sd]
tex_names = [r'$v_r^\prime/2\Omega_0(r_o - r_i)$',\
        r'$v_h^\prime/2\Omega_0(r_o - r_i)$',\
        r'$v^\prime/2\Omega_0(r_o - r_i)$']

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_Ro_shell_depth_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Loop through and plot length scales
count = 0
for var in vars_to_plot:
    plt.plot(rr_n, var, label=tex_names[count])
    count += 1

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits if desired
if not minmax is None:
    plt.ylim(minmax[0], minmax[1])

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)
plt.ylabel('Shell-depth Rossby number', fontsize=12, **csfont)

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

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + ' ' + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
plt.title(dirname_stripped + '\n' + 'Shell-depth Rossby number\n' +\
          time_string, **csfont)

# Create a see-through legend
plt.legend(shadow=True, fontsize=14, framealpha=0.5)

if log:
    plt.yscale('log')

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the Rossby number plot at ' + plotdir + savename)
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()
