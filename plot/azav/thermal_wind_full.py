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

# Compute the finite difference curls of the r/theta forces
force_r_adv = -vals[:, :, lut[1201]]
force_r_adv_rs = -vals[:, :, lut[1210]]
force_r_cor = -vals[:, :, lut[1219]]
force_r_adv_mm = force_r_adv - force_r_adv_rs + force_r_cor
force_r_buoy = vals[:, :, lut[1216]]
force_r_visc = vals[:, :, lut[1228]]
if clas0['magnetism']:
    force_r_mag_mm = vals[:, :, lut[1248]]
    force_r_mag_ms = vals[:, :, lut[1260]]

force_t_adv = -vals[:, :, lut[1202]]
force_t_adv_rs = -vals[:, :, lut[1211]]
force_t_cor = -vals[:, :, lut[1220]]
force_t_adv_mm = force_t_adv - force_t_adv_rs + force_t_cor
force_t_visc = vals[:, :, lut[1229]]
if clas0['magnetism']:
    force_r_mag_mm = vals[:, :, lut[1249]]
    force_r_mag_ms = vals[:, :, lut[1261]]

# take the curl (divide by rho first)
eq = get_eq(dirname)
rho_2d = eq.rho.reshape((1, nr))

svort_adv_rs = curlphi(force_r_adv_rs/rho_2d, force_t_adv_rs/rho_2d, rr, tt)
svort_adv_mm = curlphi(force_r_adv_mm/rho_2d, force_t_adv_mm/rho_2d, rr, tt)
svort_visc = curlphi(force_r_visc/rho_2d, force_t_visc/rho_2d, rr, tt)
svort_buoy = curlphi(force_r_visc/rho_2d, force_t_visc/rho_2d, rr, tt)
if clas0['magnetism']:
    svort_mag_mm = curlphi(force_r_mag_mm/rho_2d, force_t_mag_mm/rho_2d, rr, tt)
    svort_mag_ms = curlphi(force_r_mag_ms/rho_2d, force_t_mag_ms/rho_2d, rr, tt)

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
kw_plot_azav_grid.maintitle = dirname_stripped + '\n' +\
        'full thermal wind balance' + '\n' +\
        time_string

# terms to plot and sub-titles
terms = [svort_buoy, svort_adv_mm, svort_adv_rs, svort_visc]
kw_plot_azav_grid.titles = ['svort_buoy', 'svort_adv_mm', 'svort_adv_rs', 'svort_visc']
if clas0['magnetism']:
    terms.append(svort_mag_mm)
    terms.append(svort_mag_ms)
    kw_plot_azav_grid.titles.append('svort_mag_mm')
    kw_plot_azav_grid.titles.append('svort_mag_ms')

kw_plot_azav_grid.totsig = np.ones(len(terms))

# make figure using usual routine
fig = plot_azav_grid (terms, rr, cost, **kw_plot_azav_grid)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
