# Author: Loren Matilsky
# Created: 05/14/2018
# Updated: 2022
# Plots mass flux (stream lines and magnitude)
# savename: mercirc-[first iter]_[last iter].png
 
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
    kw_make_figure.sub_margin_bottom_inches *= 2
    print(kw_make_figure.sub_margin_bottom_inches)

# Get density
eq = get_eq(dirname)
rho = eq.rho

# Directory with data and plots, make the plotting directory if it doesn't
# already exist    

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
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

# make plot

fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0, 0]

# Plot mass flux
kw_plot_azav.plotcontours = False
plot_azav (rhovm, rr, cost, fig, ax, **kw_plot_azav)

# Plot streamfunction contours
lilbit = 0.01
maxabs = np.max(np.abs(psi))
contourlevels = (-maxabs/2., -maxabs/4., -lilbit*maxabs, 0.,\
        lilbit*maxabs, maxabs/4., maxabs/2.)
kw_plot_azav.plotcontours = kw.plotcontours
kw_plot_azav.plotfield = False
kw_plot_azav.contourlevels = contourlevels
plot_azav (psi, rr, cost, fig, ax, **kw_plot_azav)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file) 
time_string = get_time_string(dirname, iter1, iter2, threelines=True)
margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
the_title = dirname_stripped + '\n' + 'Mass flux (circulation)\n' + time_string
ax.set_title(the_title, fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
