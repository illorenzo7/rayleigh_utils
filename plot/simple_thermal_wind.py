# Author: Loren Matilsky
# Created: 04/06/2019
# This script plots mean velocity components in the meridional plane 
# for the Rayleigh run directory indicated by [dirname]. 
# Uses an instantaneous AZ_Avgs file unless told otherwise
# Saves plot in
# [dirname]_v_azav_[iter].npy
# or [dirname]_v_azav_[first iter]_[last iter].npy if a time average was 
# specified

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
from azav_util import plot_azav
from common import get_widest_range_file, strip_dirname, get_file_lists,\
        get_desired_range, get_dict
from get_eq import get_eq
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
rbcz = None

# Read in CLAs (if any) to change default variable ranges and other options
minmax = None

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

# Coriolis term:
Om0 = get_parameter(dirname, 'angular_velocity')
vp = vals[:, :, lut[3]]
# Compute the finite difference axial derivative of vp:
dvpdr = drad(vp, rr)
dvpdt = dth(vp, tt)
dvpdz = cost_2d*dvpdr - sint_2d*dvpdt/rr_2d
T1 = 2*Om0*dvpdz

# Baroclinic term:
eq = get_eq(dirname)

g = eq.gravity
cp = get_parameter(dirname, 'pressure_specific_heat')

rho_bar = eq.density
t_bar = eq.temperature
kappa = eq.kappa
cond_flux_theta = vals[:, :, lut[1471]]
dsdrt = -cond_flux_theta/(rho_bar*t_bar*kappa).reshape((1, nr))
T2 = -g.reshape((1, nr))*dsdrt/cp

# Set up the actual figure from scratch
fig_width_inches = 7 # TOTAL figure width, in inches
    # (i.e., 8x11.5 paper with 1/2-inch margins)
margin_inches = 1/8 # margin width in inches (for both x and y) and 
    # horizontally in between figures
margin_bottom_inches = 0.75*(2 - (rbcz is None)) 
    # larger bottom margin to make room for colorbar(s)
margin_top_inches = 1 # wider top margin to accommodate subplot titles AND metadata
margin_subplot_top_inches = 1/4 # margin to accommodate just subplot titles
ncol = 3 # put three plots per row: T1, T2, and T1 + T2
nrow = 1

subplot_width_inches = (fig_width_inches - (ncol + 1)*margin_inches)/ncol
    # Make the subplot width so that ncol subplots fit together side-by-side
    # with margins in between them and at the left and right.
subplot_height_inches = 2*subplot_width_inches # Each subplot should have an
    # aspect ratio of y/x = 2/1 to accommodate meridional planes. 
fig_height_inches = margin_top_inches + nrow*(subplot_height_inches +\
        margin_subplot_top_inches + margin_bottom_inches)
fig_aspect = fig_height_inches/fig_width_inches

# "Margin" in "figure units"; figure units extend from 0 to 1 in BOTH 
# directions, so unitless dimensions of margin will be different in x and y
# to force an equal physical margin
margin_x = margin_inches/fig_width_inches
margin_y = margin_inches/fig_height_inches
margin_bottom = margin_bottom_inches/fig_height_inches
margin_top = margin_top_inches/fig_height_inches
margin_subplot_top = margin_subplot_top_inches/fig_height_inches

# Subplot dimensions in figure units
subplot_width = subplot_width_inches/fig_width_inches
subplot_height = subplot_height_inches/fig_height_inches

field_components = [T1, T2, T1+T2]
titles = [r'$T_1\equiv2\Omega_0\partial\langle v_\phi\rangle/\partial z$',\
        r'$T_2\equiv -(g/rc_p)\partial\langle S\rangle/\partial \theta$',\
        r'$T_1+T_2$']
units = r'$\rm{s}^{-2}$'

# Generate the actual figure of the correct dimensions
fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

for iplot in range(3):
    ax_left = margin_x + (iplot%ncol)*(subplot_width + margin_x)
    ax_bottom = 1 - margin_top - subplot_height - margin_subplot_top -\
            (iplot//ncol)*(subplot_height + margin_subplot_top +\
            margin_bottom)
    ax = fig.add_axes((ax_left, ax_bottom, subplot_width, subplot_height))
    plot_azav (field_components[iplot], rr, cost, fig=fig, ax=ax,\
           units=units, nlevs=my_nlevs, minmax=minmax,\
           plotcontours=plotcontours)
    ax.set_title(titles[iplot], verticalalignment='bottom', **csfont)

# Put some metadata in upper left
fsize = 12
fig.text(margin_x, 1 - 0.1*margin_top, dirname_stripped,\
         ha='left', va='top', fontsize=fsize, **csfont)
fig.text(margin_x, 1 - 0.3*margin_top, 'Simple thermal wind balance',\
         ha='left', va='top', fontsize=fsize, **csfont)
iter1, iter2 = di['iter1'], di['iter2']
iter_string = str(iter1).zfill(8) + ' to ' + str(iter2).zfill(8) 

fig.text(margin_x, 1 - 0.5*margin_top, iter_string,\
         ha='left', va='top', fontsize=fsize, **csfont)

savename = dirname_stripped + '_simple_tw_' + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + '.png'
if save:
    plt.savefig(plotdir + savename, dpi=100)
    print ("Saving first two terms of TW balance in " + plotdir + savename + ' ...')
plt.show()
