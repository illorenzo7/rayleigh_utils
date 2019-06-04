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
from binormalized_cbar import MidpointNormalize
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['co'])
from common import get_widest_range_file, strip_dirname, get_dict

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
my_boundstype = 'manual'
user_specified_minmax = False
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-minmax'):
        my_min, my_max = float(args[i+1]), float(args[i+2])
        user_specified_minmax = True
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)
iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']
S_azav = vals[:, :, lut[501]]

# Read in Shell_Avgs data
print ('Getting Shell_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
Shell_Avgs_file = get_widest_range_file(datadir, 'Shell_Avgs')
di_sh = get_dict(datadir + Shell_Avgs_file)
vals_sh = di_sh['vals']
lut_sh = di_sh['lut']
S_shav = vals_sh[:, lut_sh[501]]

# Get necessary grid info
rr = di['rr']
ro = di['ro']
rr_depth = di['rr_depth']

tt_lat = di['tt_lat']

nr, nt = di['nr'], di['nt']

# Get deviation of azimuthal average of entropy from the spherical average
S_prime = S_azav - S_shav.reshape((1, nr))

# Plot the entropy deviation vs. latitude, at several depths: 
# 0, 5, 10, 15, and 25 per cent down
rvals_desired = np.array([1, 0.95, 0.9, 0.85, 0.75])
nrvals = len(rvals_desired)
ir_to_plot = np.zeros(nrvals, dtype=int)
for i in range(nrvals):
    rval_desired = rvals_desired[i]
    ir_to_plot[i] = np.argmin(np.abs(rr_depth - rval_desired))
    plt.plot(tt_lat, S_prime[:, ir_to_plot[i]],\
            label=r'$r/r_o=%0.3f$' %(rr[ir_to_plot[i]]/ro))

# Get the zero line
plt.plot(tt_lat, np.zeros(nt), 'k', linewidth=0.7)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = -90, 90
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set y limits if user wanted you to
if user_specified_minmax:
    plt.ylim(my_min, my_max)

# Create a see-through legend
leg=plt.legend(shadow=True,fontsize=10)
leg.get_frame().set_alpha(0.3)

# Label the axes
plt.xlabel(r'$\rm{Latitude} \ (^\circ)$', fontsize=12)
plt.ylabel(r'$\langle S\rangle_\phi - \langle S\rangle_{\theta,\phi}\ [\rm{erg}\ \rm{g}^{-1}\ \rm{K}^{-1}]$',\
        fontsize=12)

# Make title
plt.title(dirname_stripped + '     ' + str(iter1).zfill(8) + ' to ' +\
        str(iter2).zfill(8) + '\n' +\
        'entropy (deviation from l=0), bottom', **csfont)

# Last command
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_entropy_vs_theta_bot_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
