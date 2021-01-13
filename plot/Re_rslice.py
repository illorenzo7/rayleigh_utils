# Author: Loren Matilsky
# Created: 04/01/2018
# This script generates the convective Reynolds number (Re), plotted along 
# radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_Re_rslice_[first iter]_[last iter].png

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from rayleigh_diagnostics import TransportCoeffs
from reference_tools import equation_coefficients
from common import *
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
rnorm = None
minmax = None
logscale = False
rvals = [] # user can specify radii to mark by vertical lines
just_vr = False
rvals = []
tag = ''
lats = [0., 15., 30., 45., 60., 75.]
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if arg == '-lats':
        lats_str = args[i+1].split()
        lats = []
        for j in range(len(lats_str)):
            lats.append(float(lats_str[j]))
    elif arg == '-usefile':
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif arg == '-rnorm':
        rnorm = float(args[i+1])
    elif arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-log':
        logscale = True
    elif arg == '-vr':
        just_vr = True
    elif arg == '-tag':
        tag = '_' + args[i+1]
    elif arg == '-rvals':
        rvals_str = args[i+1].split()
        rvals = []
        for j in range(len(rvals_str)):
            rvals.append(float(rvals_str[j]))

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90. - lats
theta_vals = colats*np.pi/180.

# Read in vavg data
print ('Reading AZ_Avgs data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)
vals = di['vals']
lut = di['lut']
iter1, iter2 = di['iter1'], di['iter2']
rr = di['rr']
tt = di['tt']
cost, sint = di['cost'], di['sint']
xx = di['xx']
ri = di['ri']

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

# Convective velocity amplitudes...
vsq_r, vsq_t, vsq_p = vals[:, :, lut[422]], vals[:, :, lut[423]],\
    vals[:, :, lut[424]], 
vsq = vsq_r
if not just_vr:
    print ("adding all velocity components together (not just vr)")
    vsq += vsq_t + vsq_p
d = di['d']

# Get molecular diffusivity from 'transport' file or equation_coefficients
try:
    t = TransportCoeffs(dirname + '/transport')
    nu = t.nu
except:
    eq = equation_coefficients()
    eq.read(dirname + '/equation_coefficients')
    # nu(r) = c_5 * f_3
    nu = eq.constants[4]*eq.functions[2]

# Compute magnetic Reynolds number
Re = np.sqrt(vsq)*d/nu

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value
                                               
# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if rnorm is None:
    rr_n = rr/rsun
else:
    rr_n = rr/rnorm                                           

# Plot rotation vs radius at the desired latitudes
for theta_val in theta_vals:
    diffs = np.abs(tt - theta_val)
    index = np.argmin(diffs)
    latitude = 90. - theta_val*180./np.pi
    ax.plot(rr_n, Re[index,:],\
            label = r'$\rm{%2.1f}$' %latitude + r'$^\circ$')

    maxes.append(np.max(Re[index,:]))
    mins.append(np.min(Re[index,:]))

# Global extrema
mmax = np.max(maxes)
mmin = np.min(mins)
difference = mmax - mmin
buffer = difference*0.2 # "Guard" the yrange of the plot with whitespace

# Label the axes
if rnorm is None:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %rnorm, fontsize=12, **csfont)

if just_vr:
    plt.ylabel(r'${\rm{Re}} = (r_o-r_i)v_r^\prime \nu^{-1}$',fontsize=12,\
        **csfont)
else:
    plt.ylabel(r'${\rm{Re}} = (r_o-r_i)v^\prime \nu^{-1}$',fontsize=12,\
        **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))
if not logscale or not minmax is None:
    if minmax is None:
        ymin, ymax = mmin - buffer, mmax + buffer
    else:
        ymin, ymax = minmax
    plt.ylim((ymin, ymax))
if logscale:
    plt.yscale('log')

ymin, ymax = ax.get_ylim()
xvals = np.linspace(xmin, xmax, 100)
yvals = np.linspace(ymin, ymax, 100)

# Mark radii if desired
if not rvals is None:
    for rval in rvals:
        if rnorm is None:
            rval_n = rval/rsun
        else:
            rval_n = rval/rnorm
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

# Create a title    
plt.title(dirname_stripped + '\n' +'Convective Reynolds number\n' +\
          time_string, **csfont)
plt.legend(title='latitude')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_Re_rslice_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + tag + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
