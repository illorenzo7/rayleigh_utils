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
from common import strip_dirname, get_widest_range_file,\
    get_iters_from_file, rsun, get_dict
from get_parameter import get_parameter

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
user_specified_minmax = False
lats = [0, 15, 30, 45, 60, 75]
enstrophy_file = get_widest_range_file(datadir, 'enstrophy_from_mer')
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
        enstrophy_file = args[i+1]
        enstrophy_file = enstrophy_file.split('/')[-1]
    elif arg == '-rnorm':
        user_specified_rnorm = True
        user_supplied_rnorm = float(args[i+1])
    elif arg == '-minmax':
        user_specified_minmax = True
        my_min, my_max = float(args[i+1]), float(args[i+2])

# Get angular velocity
Om0 = get_parameter(dirname, 'angular_velocity')

# Get the spherical theta values associated with [lats]       
lats = np.array(lats)
colats = 90 - lats
theta_vals = colats*np.pi/180

# Read in enstrophy data
# First, try to get all the data from the AZ_Avgs file

# Get enstrophy from the AVERAGE vorticities (I almost always output these)
# and basic grid info
print("Getting enstrophy from the averaged vorticities and grid info from " +\
        AZ_Avgs_file + ' ...')
di_az = get_dict(datadir + AZ_Avgs_file)
vals_az = di_az['vals']
lut_az = di_az['lut']
nq_az = di_az['nq']

rr = di_az['rr']
tt = di_az['tt']
cost, sint = di_az['cost'], di_az['sint']
xx = di_az['xx']

ens_r_mean = vals_az[:, :, lut_az[301]]**2
ens_theta_mean = vals_az[:, :, lut_az[302]]**2
ens_phi_mean = vals_az[:, :, lut_az[303]]**2
ens_mean = ens_r_mean + ens_theta_mean + ens_phi_mean

# Try now to also get total enstrophy from AZ_Avgs data:
print("Trying to also get total enstrophy from " +\
        AZ_Avgs_file + ' ...')
if lut_az[314] < nq_az and lut_az[315] < nq_az and lut_az[316] < nq_az:
    ens_r = vals_az[:, :, lut_az[314]]
    ens_theta = vals_az[:, :, lut_az[315]]
    ens_phi = vals_az[:, :, lut_az[316]]

else:
    print("Could not get total enstrophy from AZ_Avgs file. Trying to get it from 'enstrophy_from_mer' data ...")
    if enstrophy_file == '':
        print("No 'enstrophy_from_mer' data found. Please run\n     python $co/timeavg/enstrophy_from_mer.py [OPTIONS]'\nand try plotting again.")
        sys.exit()
    else:
        di = get_dict(datadir + enstrophy_file)
        vals = di['vals']
        lut = di['lut']
        iter1, iter2 = di['iter1'], di['iter2']
        ens_r = np.mean(vals[:, :, :, 0], axis=0)
        ens_theta = np.mean(vals[:, :, :, 1], axis=0)
        ens_phi = np.mean(vals[:, :, :, 2], axis=0)

# Compute total enstrophy from component enstrophies
ens = ens_r + ens_theta + ens_phi

# Compute vorticity-based Rossby number
Ro = np.sqrt(ens - ens_mean)/Om0

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Get extrema values for diff. rot.
maxes = [] # Get the max-value of Omega for plotting purposes
mins = []  # ditto for the min-value
                                               
# User can specify what to normalize the radius by
# By default, normalize by the solar radius
if not user_specified_rnorm:
    rr_n = rr/rsun
else:
    rr_n = rr/user_supplied_rnorm                                           

# Plot rotation vs radius at the desired latitudes
for theta_val in theta_vals:
    diffs = np.abs(tt - theta_val)
    index = np.argmin(diffs)
    latitude = 90 - theta_val*180/np.pi
    ax.plot(rr_n, Ro[index,:],\
            label = r'$\rm{%2.1f}$' %latitude + r'$^\circ$')
    maxes.append(np.max(Ro[index,:]))
    mins.append(np.min(Ro[index,:]))

if not user_specified_minmax:
    # Global extrema
    mmax = np.max(maxes)
    mmin = np.min(mins)
    difference = mmax - mmin
    ybuffer = difference*0.2 
    # "Guard" the yrange of the plot with whitespace
    ymin, ymax = mmin - ybuffer, mmax + ybuffer
else:
    ymin, ymax = my_min, my_max

# Label the axes
if not user_specified_rnorm:
    plt.xlabel(r'$r/R_\odot$',fontsize=12, **csfont)
else:
    plt.xlabel(r'r/(%.1e cm)' %user_supplied_rnorm, fontsize=12, **csfont)
    
plt.ylabel(r'${\rm{Ro}}_{\rm{vort}} \equiv \omega^\prime \Omega_0^{-1}$',fontsize=12, **csfont)

# Set the axis limits
xmin, xmax = np.min(rr_n), np.max(rr_n)
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
delta_x = xmax - xmin
delta_y = ymax - ymin
xvals = np.linspace(xmin, xmax, 100)
yvals = np.linspace(ymin, ymax, 100)

# Create a title    
plt.title(dirname_stripped + '\n' +'"Vorticity" Rossby number, ' +\
          str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8), **csfont)
plt.legend(title='latitude')

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_Ro_vort_rslice_' +\
    str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print('Saving plot at ' + savefile + ' ...')
plt.savefig(savefile, dpi=300)
plt.show()
