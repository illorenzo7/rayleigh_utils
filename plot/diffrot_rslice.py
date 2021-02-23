# Author: Loren Matilsky
# Created: 05/14/2018
# This script generates differential rotation plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_rslice_[first iter]_[last iter].npy

# Import relevant modules
import numpy as np
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
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

# Set defaults
rnorm = None
minmax = None
lats = [0., 15., 30., 45., 60., 75.]
the_file = get_widest_range_file(datadir, 'AZ_Avgs')
rvals = []

# Read command-line arguments (CLAs)
plotdir = None

args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-plotdir':
        plotdir = args[i+1]
    if arg == '-lats':
        lats_str = args[i+1].split()
        lats = []
        for j in range(len(lats_str)):
            lats.append(float(lats_str[j]))
    elif arg == '-usefile':
        the_file = args[i+1]
        the_file = the_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
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

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90. - lats
theta_vals = colats*np.pi/180.

# Read in vavg data
print ('Reading AZ_Avgs data from ' + datadir + the_file)
di = get_dict(datadir + the_file)

vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost, sint = di['cost'], di['sint']
xx = di['xx']
ri = di['ri']

vr_av, vt_av, vp_av = vals[:, :, lut[1]], vals[:, :, lut[2]],\
        vals[:, :, lut[3]]

# Get the time range in sec
t1 = translate_times(iter1, dirname, translate_from='iter')['val_sec']
t2 = translate_times(iter2, dirname, translate_from='iter')['val_sec']

# Get the baseline time unit
time_unit = compute_Prot(dirname)
time_label = r'$\rm{P_{rot}}$'

# Get frame rate rotation and compute differential rotation in the 
# lab frame. 
Om0 = 2.*np.pi/time_unit
Om = vp_av/xx + Om0
Om *= 1.0e9/2/np.pi # convert from rad/s --> nHz
Om0 *= 1.0e9/2/np.pi # convert from rad/s --> nHz

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value
                                               
# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rnorm = rsun
    xlabel = r'$r/R_\odot$'
else:
    xlabel = r'r/(%.1e cm)' %rnorm
rr_n = rr/rnorm

# Plot rotation vs radius at the desired latitudes
for theta_val in theta_vals:
    diffs = np.abs(tt - theta_val)
    index = np.argmin(diffs)
    latitude = 90 - theta_val*180/np.pi
    ax.plot(rr_n, Om[index,:],\
            label=r'$\rm{%2.1f}$' %latitude + r'$^\circ$')
    maxes.append(np.max(Om[index,:]))
    mins.append(np.min(Om[index,:]))

# show the frame rotation rate
plt.plot(rr_n, Om0 + np.zeros_like(rr_n), 'k--', label=r'$\Omega_0=%.1f\ \rm{nHz}$' %Om0)

# Global extrema
mmax = np.max(maxes)
mmin = np.min(mins)

# Label the axes
plt.xlabel(xlabel, fontsize=12, **csfont)
plt.ylabel(r'$\Omega/2\pi \ \rm{(nHz)}$',fontsize=12, **csfont)

# Set the axis limits
# x axis
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))
# y axis
if minmax is None:
    ydiff = mmax - mmin
    ybuffer = ydiff*0.05 # "Guard" the yrange of the plot with whitespace
    minmax = mmin - ybuffer, mmax + ybuffer
plt.ylim((minmax[0], minmax[1]))

# Mark radii if desired
if not rvals is None:
    yvals = np.linspace(minmax[0], minmax[1], 100)
    for rval in rvals:
        rval_n = rval/rnorm
        plt.plot(rval_n + np.zeros(100), yvals, 'k--')

# Create a title    
# Label averaging interval
time_string = ('t = %.1f to %.1f ' %(t1/time_unit, t2/time_unit))\
        + time_label + '\n' + (r'$\ (\Delta t = %.1f\ $'\
        %((t2 - t1)/time_unit)) + time_label + ')'

plt.title(dirname_stripped + '\n' +\
        r'$\Omega(r,\theta),\ \rm{rslice},\ $' +\
          time_string, **csfont)
plt.legend(title='latitude')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_diffrot_rslice_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
