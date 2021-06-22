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
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import plot_azav_grid
from common import *
from plotcommon import *
from cla_util import *

# Read command-line arguments (CLAs)
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], dataname)

print ('Getting quantities from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
if dataname == 'AZ_Avgs':
    lut = di['lut']

# see if the user wants a separate plot of lat. averaged quantities
if 'shav' in clas:
    shav = True
else:
    shav = False

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

# save the figure
plotdir = make_plotdir(dirname, clas['plotdir'], '/plots/azav/')
savefile = plotdir + clas['routinename'] + clas['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas['showplot']:
    plt.show()
plt.close()

savename = dirname_stripped + '_simple_tw_' + str(iter1).zfill(8) + '_' +\
        str(iter2).zfill(8) + '.png'
if save:
    plt.savefig(plotdir + savename, dpi=100)
    print ("Saving first two terms of TW balance in " + plotdir + savename + ' ...')
plt.show()
