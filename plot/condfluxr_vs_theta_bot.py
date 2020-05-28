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
from common import get_widest_range_file, strip_dirname, get_dict
from get_parameter import get_parameter
from get_eq import get_eq
from time_scales import compute_Prot, compute_tdt
from translate_times import translate_times

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

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
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)
iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
cflux_azav = vals[:, :, lut[1470]]

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


# Get necessary grid info
rr = di['rr']
ri = di['ri']
ro = di['ro']
rr_depth = di['rr_depth']

tt_lat = di['tt_lat']

nr, nt = di['nr'], di['nt']

# Plot radial conductive flux
# 0, 5, 10, 15, and 25 per cent from bottom
rvals_desired = np.array([0.75, 0.85, 0.9, 0.95, 1.])
nrvals = len(rvals_desired)
for i in range(nrvals):
    rval_desired = rvals_desired[i]
    ir_to_plot = np.argmin(np.abs(rr_depth - rval_desired))
    plt.plot(tt_lat, cflux_azav[:, ir_to_plot],\
            label=r'$r/r_o=%0.3f$' %(rr[ir_to_plot]/ro))

# Plot the luminosity that must be driven through (heat) as a flux:
eq = get_eq(dirname)
lum = eq.lum 
Flux_in = lum/4/np.pi/ri**2
plt.plot(tt_lat, Flux_in*np.ones(nt), 'k--', alpha=0.3, label=r'$L_*/4\pi r_i^2$')

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
leg=plt.legend(shadow=True,fontsize=10)
leg.get_frame().set_alpha(0.3)

# Label the axes
plt.xlabel(r'$\rm{Latitude} \ (^\circ)$', fontsize=12)
plt.ylabel(r'$\mathcal{F}_{{\rm{cond}},r} = -\kappa \overline{\rho}\overline{T}\partial\langle S\rangle_\phi/\partial r\ [\rm{erg}\ \rm{cm}^{-2}\ \rm{s}^{-1}]$',\
        fontsize=12)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label  + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
plt.title(dirname_stripped + '\nradial conductive flux, bottom\n' +\
        time_string, **csfont)

# Last command
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_condfluxr_vs_theta_bot_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
