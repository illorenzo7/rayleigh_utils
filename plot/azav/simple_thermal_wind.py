# Author: Loren Matilsky
# Created: 12/19/2022
#
# Description: Script to plot the "simple" thermal-wind equation in 
# the meridional plane:
# d(Om^2)/dz = (g/r^2 sin(th))*(d/dth)(S'/c_p)
# without LBR this should also include 
# (1/r^2 sin(th)) * (d(S_bar)/dr)/c_p * (d/dth)(P'/rho_bar)

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
#azav_fig_dimensions['margin_top_inches'] += 1.5 
# make room for subplot labels
#plot_azav_grid_kwargs_default.update(azav_fig_dimensions)
kwargs_default.update(plot_azav_grid_kwargs_default)

# overwrite defaults, first main kwargs
kw = update_dict(kwargs_default, clas)
kw_plot_azav_grid = update_dict(plot_azav_grid_kwargs_default, clas)

# need a bit extra room for subplot labels
kw_plot_azav_grid.sub_margin_top_inches += 1/4

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

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
rr_2d = di_grid['rr_2d']
nr = di_grid['nr']
cost = di_grid['cost']
cost_2d = di_grid['cost_2d']
sint_2d = di_grid['sint_2d']
tt = di_grid['tt']
xx = di_grid['xx']

# Coriolis term:
Om0 = get_parameter(dirname, 'angular_velocity')
vp = vals[:, :, lut[3]]
Om = Om0 + vp/xx

# Compute the finite difference axial derivative of Om^2
dom2dr = drad(Om**2, rr)
dom2dt = dth(Om**2, tt)
dom2dz = cost_2d*dom2dr - sint_2d*dom2dt/rr_2d
T1 = xx*dom2dz

# Baroclinic term:
eq = get_eq(dirname)
cond_flux_theta = vals[:, :, lut[1471]]
dsdrt = -cond_flux_theta/(eq.rho*eq.tmp*eq.kappa).reshape((1, nr))
T2 = -eq.grav.reshape((1, nr))*dsdrt/eq.c_p

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
kw_plot_azav_grid.maintitle = dirname_stripped + '\n' +\
        'basic thermal wind balance' + '\n' +\
        time_string

# terms to plot and sub-titles
terms = [T1, T2, T1 + T2]
kw_plot_azav_grid.titles = [r'$T_1\equiv2\Omega_0\partial\langle v_\phi\rangle/\partial z$',\
        r'$T_2\equiv -(g/rc_p)\partial\langle S\rangle/\partial \theta$',\
        r'$T_1+T_2$']

# make figure using usual routine
fig = plot_azav_grid (terms, rr, cost, **kw_plot_azav_grid)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + 'azav_simple_tw' + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
