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
from azav_util import plot_azav
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

# get data
if 'the_file' in clas: 
    the_file = clas['the_file']
else:
    the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + the_file)
di = get_dict(the_file)
vals = di['vals']
lut = di['lut']
vp_av = vals[:, :, lut[3]]

# Get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
cost = di_grid['cost']
tt_lat = di_grid['tt_lat']
xx = di_grid['xx']

# Get differential rotation in the rotating frame. 
Om = vp_av/xx
diffrot = Om*1.0e9/2/np.pi # rad/s --> nHz

# DR contrast between 0 and 60 degrees
it0, it60 = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60))
Delta_Om = diffrot[it0, 0] - diffrot[it60, 0]

# create plot
nplots = 1
sub_width_inches = 2.
sub_aspect = 2
margin_top_inches = 1 # larger top margin to make room for titles
margin_bottom_inches = 1/2
# larger bottom margin to make room for colorbar(s)
if 'rbcz' in clas:
    margin_bottom_inches *= 2

# make plot
fig, axs, fpar = make_figure(nplots=nplots, sub_width_inches=sub_width_inches, sub_aspect=sub_aspect, margin_top_inches=margin_top_inches, margin_bottom_inches=margin_bottom_inches)
ax = axs[0, 0]

plot_azav (diffrot, rr, cost, fig, axs[0, 0], units='nHz', plotlatlines=False, nosci=True, cbar_prec=1, **clas)
        
# make title 
iter1, iter2 = get_iters_from_file(the_file)
time_string = get_time_info(dirname, iter1, iter2) 
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
the_title = dirname_stripped + '\n' +  r'$\Omega - \Omega_0$' + '\n' + time_string + '\n' + (r'$\Delta\Omega_{\rm{60}} = %.1f\ \rm{nHz}$' %Delta_Om)
fig.text(margin_x, 1 - margin_y, the_title,\
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
