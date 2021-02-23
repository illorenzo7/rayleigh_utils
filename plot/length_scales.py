###############################################
# Author: Loren Matilsky
# Date created: 04/14/2020
#
# This script plots varius length_scales as functions of
# radius using from the Shell_Avgs/Shell_Spectra data

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
ynorm = None
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
    elif arg == '-ynorm':
        ynorm = float(args[i+1])
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

# Directory with plots, make the plotting directory if it doesn't
# already exist    

di = get_length_scales(dirname)
rr = di['rr']
rr_spec = di['rr_spec']
nr = di['nr']
nr_spec = di['nr_spec']
iter1, iter2 = di['iter1'], di['iter2']

# Get grid spacing
dr = np.zeros(nr)
dr[:-1] = rr[:-1] - rr[1:]
dr[-1] = dr[-2]

# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
    rr_spec_n = rr_spec/rsun
else:
    rr_n = rr/rnorm                                           
    rr_spec_n = rr_spec/rnorm

if ynorm is None:
    ynorm = rsun

scale_names = ['L_rho', 'L_omr', 'L_omh', 'L_om', 'L_vr', 'L_vh', 'L_v']
tex_names = [r'$L_\rho$', r'$L_{\omega_r^\prime}$', r'$L_{\omega_h^\prime}$', r'$L_{\omega^\prime}$',  r'$L_{v_r^\prime}$', r'$L_{v_h^\prime}$',  r'$L_{v^\prime}$']

magnetism = get_parameter(dirname, 'magnetism')
if magnetism:
    scale_names += ['L_J', 'L_Bp', 'L_Bm', 'L_B']
    tex_names += [r'$L_{J^\prime}$', r'$L_{B_\phi^\prime}$', r'$L_{B_m^\prime}$',  r'$L_{B^\prime}$']

# Make the plot name, labelling the first/last iterations we average over
savename = dirname_stripped + '_length_scales_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

# Loop through and plot length scales
count = 0
for scale_name in scale_names:
    length_scale = di[scale_name]
    if len(length_scale) == nr: # plot a line
        plt.plot(rr_n, length_scale/rsun, label=tex_names[count])
    elif len(length_scale) == nr_spec:
        plt.scatter(rr_spec_n, length_scale/rsun,\
                label=tex_names[count], s=10.)
    count += 1

# Plot the grid resolution (will hopefully appear underneath everything
# else!
plt.plot(rr_n, dr/rsun, 'k', label=r'$\delta r$')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set the y-limits (the following values seem to "work well" for my models
# so far...perhaps adjust this in the future. 

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

plt.ylabel('length/ynorm (default rsun)',\
        fontsize=12, **csfont)

# Make title
plt.title(dirname_stripped + '\n' + 'length scales, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)

# Create a see-through legend
plt.legend(shadow=True, ncol=3, fontsize=14, framealpha=0.5)

if log:
    plt.yscale('log')

# Last command
plt.tight_layout()

# Save the plot
print ('Saving the length_scales plot at ' + plotdir + savename + ' ...')
plt.savefig(plotdir + savename, dpi=300)

# Show the plot
plt.show()

plt.show()
