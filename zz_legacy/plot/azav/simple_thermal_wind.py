# Author: Loren Matilsky
# Created: 05/14/2018
# Updated: 2021
# prototype for AZ_Avgs scripts. Savename:
# diffrot-[first iter]_[last iter].npy

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
rr_2d = di_grid['rr_2d']
nr = di_grid['nr']
cost = di_grid['cost']
cost_2d = di_grid['cost_2d']
sint_2d = di_grid['sint_2d']
tt = di_grid['tt']

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
ref_rho = eq.density
ref_temp = eq.temperature
kappa = eq.kappa
cond_flux_theta = vals[:, :, lut[1471]]
dsdrt = -cond_flux_theta/(ref_rho*ref_temp*kappa).reshape((1, nr))
T2 = -g.reshape((1, nr))*dsdrt/c_P

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Basic thermal wind balance' + '\n' +\
        time_string

# terms to plot and sub-titles
terms = [T1, T2, T1 + T2]
titles = [r'$T_1\equiv2\Omega_0\partial\langle v_\phi\rangle/\partial z$',\
        r'$T_2\equiv -(g/rc_p)\partial\langle S\rangle/\partial \theta$',\
        r'$T_1+T_2$']

# make figure using usual routine
fig = plot_azav_grid (terms, rr, cost, maintitle=maintitle, titles=titles, **kw_plot_azav)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + 'azav_simple_tw' + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
