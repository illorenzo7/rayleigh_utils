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
make_figure_kwargs_default.update(azav_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# of course, plot_azav kwargs, but need to change a few
plot_azav_kwargs_default['plotlatlines'] = False
plot_azav_kwargs_default['units'] = 'nHz'
plot_azav_kwargs_default['nosci'] = True
plot_azav_kwargs_default['cbar_prec'] = 1
kwargs_default.update(plot_azav_kwargs_default)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_azav = update_dict(plot_azav_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)
if not kw.rbcz is None:  # need room for two colorbars
    kw_make_figure.margin_bottom_inches *= 2

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
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
Om *= 1.0e9/2/np.pi # rad/s --> nHz

# DR contrast between 0 and 60 degrees
it0, it60_N, it60_S = np.argmin(np.abs(tt_lat)), np.argmin(np.abs(tt_lat - 60)), np.argmin(np.abs(tt_lat + 60))
Delta_Om = Om[it0, 0] - (Om[it60_N, 0] + Om[it60_S, 0])/2

# make plot
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]

plot_azav (Om, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
the_title = dirname_stripped + '\n' +  r'$\Omega - \Omega_0$' + '\n' + time_string + '\n' + r'$\Delta\Omega_{\rm{60}}$' + (' = %.1f nHz' %Delta_Om)
fig.text(margin_x, 1 - margin_y, the_title,\
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
