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
kw_default = dict({'the_file': None, 'simple': False})

# also need make figure kw
#azav_fig_dimensions['margin_top_inches'] += 1.5 
# make room for subplot labels
#kw_plot_azav_grid_default.update(azav_fig_dimensions)
kw_default.update(kw_plot_azav_grid_default)

# overwrite defaults, first main kw
kw = update_dict(kw_default, clas)
kw_plot_azav_grid = update_dict(kw_plot_azav_grid_default, clas)

# need a bit extra room for subplot labels
kw_plot_azav_grid.sub_margin_top_inches += 1/4

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')

print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# Get necessary grid info
ntheta = np.shape(vals)[0]
di_grid = get_grid_info(dirname, ntheta=ntheta)
rr = di_grid['rr']
rr_2d = di_grid['rr_2d']
nr = di_grid['nr']
cost = di_grid['cost']
cost_2d = di_grid['cost_2d']
sint_2d = di_grid['sint_2d']
tt = di_grid['tt']
xx = di_grid['xx']

# Coriolis term:
eq = get_eq(dirname)
vp = vals[:, :, lut[3]]

# get the zonal vorticity (ish, don't forget about the rho)
vort_phi = vals[:, :, lut[303]]
u_theta = vals[:, :, lut[2]]
lhs = eq.rho*(eq.dlnrho*u_theta + vort_phi)

# Compute the finite difference curls of the r/theta forces
force_r_adv = -vals[:, :, lut[1201]]
force_r_adv_rs = -vals[:, :, lut[1210]]
force_r_cor = vals[:, :, lut[1219]]
force_r_adv_mm = force_r_adv - force_r_adv_rs + force_r_cor
force_r_buoy = vals[:, :, lut[1216]]
force_r_visc = vals[:, :, lut[1228]]
if clas0['magnetism']:
    force_r_mag = vals[:, :, lut[1248]]
    force_r_mag_ms = vals[:, :, lut[1260]]
    force_r_mag_mm = force_r_mag - force_r_mag_ms

force_t_adv = -vals[:, :, lut[1202]]
force_t_adv_rs = -vals[:, :, lut[1211]]
force_t_cor = vals[:, :, lut[1220]]
force_t_adv_mm = force_t_adv - force_t_adv_rs + force_t_cor
force_t_visc = vals[:, :, lut[1229]]
if clas0['magnetism']:
    force_t_mag = vals[:, :, lut[1249]]
    force_t_mag_ms = vals[:, :, lut[1261]]
    force_t_mag_mm = force_t_mag - force_t_mag_ms

# take the curl (divide by rho first)
rho_2d = eq.rho.reshape((1, nr))

svort_adv_rs = curlphi(force_r_adv_rs/rho_2d, force_t_adv_rs/rho_2d, rr, tt)
svort_adv_mm = curlphi(force_r_adv_mm/rho_2d, force_t_adv_mm/rho_2d, rr, tt)
svort_visc = curlphi(force_r_visc/rho_2d, force_t_visc/rho_2d, rr, tt)
svort_buoy = curlphi(force_r_buoy/rho_2d, np.zeros_like(force_r_buoy), rr, tt)
if clas0['magnetism']:
    svort_mag_mm = curlphi(force_r_mag_mm/rho_2d, force_t_mag_mm/rho_2d, rr, tt)
    svort_mag_ms = curlphi(force_r_mag_ms/rho_2d, force_t_mag_ms/rho_2d, rr, tt)

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
if kw.simple:
    simple_or_full = 'simple'
else:
    simple_or_full = 'full'

kw_plot_azav_grid.maintitle = dirname_stripped + ('\n%s thermal wind balance\n' %simple_or_full) + time_string

# terms to plot and sub-titles
terms = [lhs, svort_buoy, svort_adv_mm, svort_buoy + svort_adv_mm, svort_adv_rs, svort_visc]
titles = ['curl(rho u)_phi', 'dS/dtheta', 'dOm/dz', '(b) + (c)', 'RS term', 'visc term']
if kw.simple:
    terms = terms[:4]
    titles = titles[:4]
if clas0['magnetism'] and not kw.simple:
    terms.append(svort_mag_mm)
    terms.append(svort_mag_ms)
    titles.append('mean mag.')
    titles.append('MS')
kw_plot_azav_grid.titles = np.array(titles)
totsig = np.ones(len(terms))
totsig[0] = totsig[3] = 0
kw_plot_azav_grid.totsig = totsig

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
