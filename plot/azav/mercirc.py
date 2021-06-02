# Author: Loren Matilsky
# Created: 05/14/2018
# This script plots the meridional circulation cells the meridional plane 
# for the Rayleigh run directory indicated by [dirname], using the AZ_Avgs
# data. To use an AZ_Avgs file different than the one associated with the 
# longest averaging range, run with option
# -usefile [complete name of desired vavg file]
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
from azav_util import plot_azav, streamfunction
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname, wrap=True)

# domain bounds
ncheby, domain_bounds = get_domain_bounds(dirname)
ri = np.min(domain_bounds)
ro = np.max(domain_bounds)
d = ro - ri

# Get density
eq = get_eq(dirname)
rho = eq.density

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
lut = di['lut']

# Get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
cost = di_grid['cost']

# get the meridional circulation flows
vr_av, vt_av = vals[:, :, lut[1]], vals[:, :, lut[2]]

# Compute the mass flux
rhovm = rho*np.sqrt(vr_av**2 + vt_av**2)

# Compute the streamfunction
psi = streamfunction(rho*vr_av, rho*vt_av, rr, cost)

# Make CCW negative and CW positive
rhovm *= np.sign(psi)

# Create plot
width_inches = 3.25
sub_aspect = 2
margin_top_inches = 1.25 # larger top margin to make room for titles
margin_bottom_inches = 0.7
# larger bottom margin to make room for colorbar(s)
if 'rbcz' in clas:
    margin_bottom_inches *= 2

fig, axs, fpar = make_figure(nplots=1, sub_aspect=sub_aspect, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches, width_inches=width_inches)

# Plot mass flux
plot_azav (rhovm, rr, cost, fig=fig, ax=axs[0,0], cbar_prec=1, nosci=True,\
    units=r'$\rm{g}\ \rm{cm}^{-2}\ \rm{s}^{-1}$', plotcontours=False,\
    **clas)

# Plot streamfunction contours
lilbit = 0.01
maxabs = np.max(np.abs(psi))
contourlevels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0.,\
        lilbit*maxabs, maxabs/4., maxabs/2.)
plot_azav (psi, rr, cost, fig=fig, ax=axs[0,0], plotfield=False,\
    contourlevels=contourlevels, **clas)

# make title 
iter1, iter2 = get_iters_from_file(the_file) 
time_string = get_time_info(dirname, iter1, iter2)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
line_height = 1/4/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, dirname_stripped,\
         ha='left', va='top', fontsize=default_titlesize, **csfont)
fig.text(margin_x, 1 - margin_y - 2*line_height, 'Mass flux (circulation)',\
         ha='left', va='top', fontsize=default_titlesize, **csfont)
fig.text(margin_x, 1 - margin_y - 3*line_height, time_string,\
         ha='left', va='top', fontsize=default_titlesize, **csfont)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
