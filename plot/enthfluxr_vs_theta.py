# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. To use an AZ_Avgs file
# different than the one associated with the longest averaging range, use
# -usefile [complete name of desired AZ_Avgs file]
# Saves plot in
# [dirname]_diffrot_[first iter]_[last iter].npy

import numpy as np
import pickle
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

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

# Set defaults
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
minmax = None

# Read in CLAs (if any) to change default variable ranges and other options
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)
iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
eflux_azav = vals[:, :, lut[1455]]

# Get necessary grid info
rr = di['rr']
ri = di['ri']
ro = di['ro']
rr_depth = di['rr_depth']

tt_lat = di['tt_lat']

nr, nt = di['nr'], di['nt']

# Plot radial conductive flux
# 0, 5, 10, 15, and 25 per cent from bottom
rvals_desired = np.linspace(0, 1, 10)
nrvals = len(rvals_desired)
for i in range(nrvals):
    rval_desired = rvals_desired[i]
    ir_to_plot = np.argmin(np.abs(rr_depth - rval_desired))
    plt.plot(tt_lat, eflux_azav[:, ir_to_plot],\
            label=r'$r/r_o=%0.3f$' %(rr[ir_to_plot]/ro))

# Plot the luminosity that must be driven through (heat) as a flux:
lum = get_parameter(dirname, 'luminosity')
Flux_out = lum/4/np.pi/ro**2
plt.plot(tt_lat, Flux_out*np.ones(nt), 'k--', alpha=0.3, label=r'$L_*/4\pi r_o^2$')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = -90, 90
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set y limits if user wanted you to
if not minmax is None:
    plt.ylim(minmax[0], minmax[1])

# Create a see-through legend
leg=plt.legend(shadow=True,fontsize=8, loc=2)
leg.get_frame().set_alpha(0.3)

# Label the axes
plt.xlabel(r'$\rm{Latitude} \ (^\circ)$', fontsize=12)
plt.ylabel(r'$\mathcal{F}_{{\rm{enth}},r} = -\overline{\rho}c_p\langle Tv_r\rangle_\phi\ [\rm{erg}\ \rm{cm}^{-2}\ \rm{s}^{-1}]$',\
        fontsize=12)

# Make title
plt.title(dirname_stripped + '     ' + str(iter1).zfill(8) + ' to ' +\
        str(iter2).zfill(8) + '\n' +\
        'radial enthalpy flux, different radii', **csfont)

# Last command
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_enthfluxr_vs_theta_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
