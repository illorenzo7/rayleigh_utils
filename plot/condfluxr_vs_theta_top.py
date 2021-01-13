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
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)

# Set defaults
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read in CLAs (if any) to change default variable ranges and other options
minmax = None
xminmax = None
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-xminmax':
        xminmax = float(args[i+1]), float(args[i+2])
    elif arg == '-usefile':
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
ro = di['ro']
rr_depth = di['rr_depth']
tt_lat = di['tt_lat']
nr, nt = di['nr'], di['nt']

# compute the x limits
if xminmax is None:
    xmin, xmax = -75., 75.
else:
    xmin, xmax = xminmax
 
# Get minmax lat info
ilat1 = np.argmin(np.abs(tt_lat - xmin))
ilat2 = np.argmin(np.abs(tt_lat - xmax))

# Compute luminosity that must be carried out as a flux
eq = get_eq(dirname)
lum = eq.lum 
eq_flux = lum/4/np.pi/ro**2

# Plot the entropy deviation vs. latitude, at several depths: 
# 0, 5, 10, 15, and 25 per cent down
depths_desired = np.array([0., 0.05, 0.1, 0.15, 0.25])
ndepths = len(depths_desired)
vmins = []
vmaxes = []
for i in range(ndepths):
    depth_desired = depths_desired[i]
    ir_to_plot = np.argmin(np.abs(rr_depth - depth_desired))
    plt.plot(tt_lat[ilat1:ilat2+1], cflux_azav[ilat1:ilat2+1,\
            ir_to_plot]/eq_flux,\
            label=r'$r/r_o=%0.3f$' %(rr[ir_to_plot]/ro))
    vmins.append(np.min(cflux_azav[ilat1:ilat2+1, ir_to_plot]))
    vmaxes.append(np.max(cflux_azav[ilat1:ilat2+1, ir_to_plot]))
vmins = np.array(vmins)/eq_flux
vmaxes = np.array(vmaxes)/eq_flux
vmin, vmax = np.min(vmins), np.max(vmaxes)

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
plt.xlim(xmin, xmax)

# Set y limits
if minmax is None:
    ybuffer = 0.2*(vmax - vmin)
    ymin, ymax = vmin - ybuffer, vmax + ybuffer
else:
    ymin, ymax = minmax
plt.ylim(ymin, ymax)

# Plot reference line at y = 1
plt.plot(tt_lat[ilat1:ilat2+1], np.ones(ilat2 - ilat1 + 1), 'k--')

# Plot lat = 0 (equator)
yvals = np.linspace(ymin, ymax, 100)
plt.plot(np.zeros(100), yvals, 'k--')

# Create a see-through legend
leg=plt.legend(shadow=True,fontsize=10)
leg.get_frame().set_alpha(0.3)

# Label the axes
plt.xlabel(r'$\rm{Latitude} \ (deg)$', fontsize=12)
plt.ylabel(r'$\mathcal{F}_{{\rm{cond}},\ r}/(L_*/4\pi r_{\rm{o}}^2)$',
        fontsize=12)

# Label averaging interval
if rotation:
    time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.1f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'
else:
    time_string = ('t = %.3f to %.3f ' %(t1/time_unit, t2/time_unit))\
            + time_label + (r'$\ (\Delta t = %.3f\ $'\
            %((t2 - t1)/time_unit)) + time_label + ')'

# Make title
plt.title(dirname_stripped + '\nradial conductive flux, top\n' +\
        time_string, **csfont)

# Last command
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_condfluxr_vs_theta_top_' + str(iter1).zfill(8) +\
    '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
