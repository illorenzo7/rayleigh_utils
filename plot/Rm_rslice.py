# Author: Loren Matilsky
# Created: 04/01/2018
# This script generates magnetic Reynolds # plotted along radial lines for
# the Rayleigh run directory indicated by [dirname]. To use  time-averaged 
# AZ_Avgs file different than the one associated with the longest averaging 
# range, use -usefile [complete name of desired vavg file]
# Saves plot in
# [dirname]_diffrot_rslice_[first iter]_[last iter].npy

# Import relevant modules
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
from rayleigh_diagnostics import TransportCoeffs
from common import strip_dirname, get_widest_range_file,\
        get_iters_from_file, get_dict

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
user_specified_rnorm = False
lats = [0, 15, 30, 45, 60, 75]
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read command-line arguments (CLAs)
args = sys.argv[2:]
nargs = len(args)
for i in range(nargs):
    arg = args[i]
    if (arg == '-lats'):
        lats_str = args[i+1].split()
        lats = []
        for j in range(len(lats_str)):
            lats.append(float(lats_str[j]))
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
    elif (arg == '-rnorm'):
        user_specified_rnorm = True
        user_supplied_rnorm = float(args[i+1])

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90 - lats
theta_vals = colats*np.pi/180

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

vsq_r, vsq_t, vsq_p = vals[:, :, lut[422]], vals[:, :, lut[423]],\
    vals[:, :, lut[424]], 
vsq = vsq_r + vsq_t + vsq_p
d = di['d']

# Get magnetic diffusivity from 'transport' file
t = TransportCoeffs(dirname + '/transport')
eta = t.eta

# Compute magnetic Reynolds number
Rm = np.sqrt(vsq)*d/eta

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value
                                               
# User can specify what to normalize the radius by
# By default, normalize by the solar radius
Rsun = 6.955e10
if not user_specified_rnorm:
    rr_n = rr/Rsun
else:
    rr_n = rr/user_supplied_rnorm                                           

# Plot rotation vs radius at the desired latitudes
for theta_val in theta_vals:
    diffs = np.abs(tt - theta_val)
    index = np.argmin(diffs)
    latitude = 90 - theta_val*180/np.pi
    ax.plot(rr_n, Rm[index,:],\
            label = r'$\rm{%2.1f}$' %latitude + r'$^\circ$')
    maxes.append(np.max(Rm[index,:]))
    mins.append(np.min(Rm[index,:]))

# Global extrema
mmax = np.max(maxes)
mmin = np.min(mins)
difference = mmax - mmin
buffer = difference*0.2 # "Guard" the yrange of the plot with whitespace

# Label the axes
if not user_specified_rnorm:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %user_supplied_rnorm, fontsize=12, **csfont)
    
plt.ylabel(r'${\rm{Rm}} = (r_o-r_i)v^\prime \eta^{-1}$',fontsize=12, **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
ymin, ymax = mmin - buffer, mmax + buffer
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
delta_x = xmax - xmin
delta_y = ymax - ymin
xvals = np.linspace(xmin, xmax, 100)
yvals = np.linspace(ymin, ymax, 100)

# Create a title    
plt.title(dirname_stripped + '\n' +'Magnetic Reynolds number, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend(title='latitude')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_Rm_rslice_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
