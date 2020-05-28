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
rbcz = None

# Get all the file names in datadir and their integer counterparts
file_list, int_file_list, nfiles = get_file_lists(radatadir)

# Set defaults
save = True
plotcontours = True
my_nlevs = 20
AZ_Avgs_file = get_widest_range_file(datadir, 'AZ_Avgs')
minmax = None

# Read in CLAs (if any) to change default variable ranges and other options
args = sys.argv[2:]
nargs = len(args)

# Change other defaults
for i in range(nargs):
    arg = args[i]
    if arg == '-minmax':
        minmax = float(args[i+1]), float(args[i+2])
    elif arg == '-rbcz':
        rbcz = float(args[i+1])
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
tt = di['tt']
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
# Get the normalized tilt angle
T = -2/np.pi*np.arctan(dOmdz/dOmdl)

# Create plot
subplot_width_inches = 2.5
subplot_height_inches = 5.
margin_inches = 1/8
margin_top_inches = 3/4 # larger top margin to make room for titles
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)

fig_width_inches = subplot_width_inches + 2*margin_inches
fig_height_inches = subplot_height_inches + margin_top_inches +\
        margin_bottom_inches

fig_aspect = fig_height_inches/fig_width_inches
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
ax = fig.add_axes((margin_x, margin_bottom, subplot_width, subplot_height))

plot_azav (T, rr, cost, fig=fig, ax=ax, posdef=False, minmax=minmax)

# Make title + label diff. rot. contrast and no. contours
fsize = 12
fig.text(margin_x, 1 - 1/16/fig_height_inches, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 5/16/fig_height_inches, r'$-(2/\pi){\rm{tan}}^{-1}(\partial\Omega/\partial z)/(\partial\Omega/\partial\lambda))$',\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 9/16/fig_height_inches,\
         str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8),\
         ha='left', va='top', fontsize=fsize, **csfont)
#fig.text(margin_x, 1 - 0.5*margin_top,\
#         r'$\Delta\Omega_{\rm{tot}} = %.1f\ nHz$' %Delta_Om,\
#         ha='left', va='top', fontsize=fsize, **csfont)
savefile = plotdir + dirname_stripped + '_contour_slope_' +\
        str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'
print ('Saving plot at %s ...' %savefile)
plt.savefig(savefile, dpi=300)
plt.show()
