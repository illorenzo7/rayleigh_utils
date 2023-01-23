# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot rotation-rate contours in meridional plane

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
from azav_util import *
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname, wrap=True)

# allowed args + defaults
# key unique to this script
kwargs_default = dict({'the_file': None})

# also need make figure kwargs
azav_fig_dimensions['margin_top_inches'] += 0.25
make_figure_kwargs_default.update(azav_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# and of course need plot_azav kwargs
plot_azav_kwargs_default['plotlatlines'] = False
kwargs_default.update(plot_azav_kwargs_default)

# overwrite defaults, first main kwargs
kw = update_dict(kwargs_default, clas)
kw_plot_azav = update_dict(plot_azav_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)
if not kw.rcut is None:  
    # need room for two colorbars and line up top stating rcut 
    kw_make_figure.margin_top_inches += 1/4
    kw_make_figure.sub_margin_bottom_inches *= 2

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
vp_av = vals[:, :, lut[3]]

# get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
cost = di_grid['cost']
tt_lat = di_grid['tt_lat']
xx = di_grid['xx']

# frame rate
eq = get_eq(dirname)
Om0 = 2*np.pi/eq.prot

# differential rotation in the rotating frame. 
Om = vp_av/xx

# DR contrast between 0 and 60 degrees
it0, it60_N, it60_S = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat + 60))
Delta_Om = Om[it0, 0] - (Om[it60_N, 0] + Om[it60_S, 0])/2

# make plot
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]

plot_azav (Om/Om0, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2, threelines=True) 
maintitle = dirname_stripped + '\n' +  r'$\Omega/\Omega_0 - 1$' + '\n' + time_string + '\n' + r'$\Delta\Omega_{\rm{60}}/\Omega_0$' + (' = %1.2e' %(Delta_Om/Om0))
if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut
    
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, maintitle,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
