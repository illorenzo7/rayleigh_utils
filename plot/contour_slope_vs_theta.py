# Author: Loren Matilsky
# Created: 01/15/2020
# This script plots the normalized slope angle (tilt) of angular velocity
# contours in the meridional plane:
# T = (2/pi)*arctan(abs((dOmega/dz)/(dOmega/dlambda)))
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
csfont = {'fontname':'DejaVu Serif'}
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_file_lists,\
        get_desired_range, get_dict
from get_parameter import get_parameter
from derivs import drad, dth

# Get directory name and stripped_dirname for plotting purposes
dirname = sys.argv[1]
dirname_stripped = strip_dirname(dirname)

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    
datadir = dirname + '/data/'
plotdir = dirname + '/plots/'
if (not os.path.isdir(plotdir)):
    os.makedirs(plotdir)

radatadir = dirname + '/AZ_Avgs/'

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
save = True
plotcontours = True
my_nlevs = 20
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')

# Read in CLAs (if any) to change default variable ranges and other options
minmax = None

args = sys.argv[2:]
nargs = len(args)

# Change other defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-nosave':
        save = False
    elif arg == '-nlevs':
        my_nlevs = int(args[i+1])
    elif arg == '-nocontour':
        plotcontours = False
    elif (arg == '-usefile'):
        AZ_Avgs_file = args[i+1]
        AZ_Avgs_file = AZ_Avgs_file.split('/')[-1]
        
# Read in AZ_Avgs data
print ('Getting data from ' + datadir + AZ_Avgs_file + ' ...')
di = get_dict(datadir + AZ_Avgs_file)

iter1, iter2 = di['iter1'], di['iter2']
vals = di['vals']
lut = di['lut']

# Grid info
rr = di['rr']
rr_depth = di['rr_depth']
ri, ro = di['ri'], di['ro']
tt = di['tt']
tt_lat = di['tt_lat']
cost = di['cost']
sint = di['sint']
cost_2d = di['cost_2d']
sint_2d = di['sint_2d']
rr_2d = di['rr_2d']
nr, nt = di['nr'], di['nt']
xx, zz = di['xx'], di['zz']

# Coriolis term:
Om0 = get_parameter(dirname, 'angular_velocity')
vp = vals[:, :, lut[3]]
Om = vp/xx # rotation rate

# Compute the finite difference z/lambda derivatives of the rotation rate
dOmdr = drad(Om, rr)
dOmdt = dth(Om, tt)
dOmdz = cost_2d*dOmdr - sint_2d*dOmdt/rr_2d
dOmdl = sint_2d*dOmdr + cost_2d*dOmdt/rr_2d
# Get the normalized tilt angle in degrees
T = -180/np.pi*np.arctan(dOmdz/dOmdl)

# Plot the tilt angle
# 0, 5, 10, 15, and 25 per cent from bottom
rvals_desired = np.linspace(0, 1, 10)
nrvals = len(rvals_desired)
for i in range(nrvals):
    rval_desired = rvals_desired[i]
    ir_to_plot = np.argmin(np.abs(rr_depth - rval_desired))
    plt.plot(tt_lat, T[:, ir_to_plot],\
            label=r'$r/r_o=%0.3f$' %(rr[ir_to_plot]/ro))

# Get ticks everywhere
plt.minorticks_on()
plt.tick_params(top=True, right=True, direction='in', which='both')

# Set the x limits
xmin, xmax = -90, 90
delta_x = xmax - xmin
plt.xlim(xmin, xmax)

# Set y limits if user wanted you to
if not minmax is None:
    plt.ylim(ymin[0], ymin[1])

# Create a see-through legend
leg=plt.legend(shadow=True,fontsize=8, loc=2)
leg.get_frame().set_alpha(0.3)

# Label the axes
plt.xlabel(r'$\rm{Latitude\ (deg)}$', fontsize=12)
plt.ylabel('DR tilt angle (deg)',\
        fontsize=12)

# Make title
plt.title(dirname_stripped + '     ' + str(iter1).zfill(8) + ' to ' +\
        str(iter2).zfill(8) + '\n' +\
        'Normalized D.R. tilt angle, different radii', **csfont)

# Last command
plt.tight_layout()

savefile = plotdir + dirname_stripped + '_contour_slope_vs_theta_' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
