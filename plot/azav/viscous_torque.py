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
ntheta = np.shape(vals)[0]
gi = get_grid_info(dirname, ntheta=ntheta)

# Coriolis term:
eq = get_eq(dirname)
vp = vals[:, :, lut[3]]
Om = vp/gi.xx

# get angular momentum viscous fluxes, and thus derivatives of <v_phi>
amom_visc_r = vals[:, :, lut[1813]]
amom_visc_t = vals[:, :, lut[1814]]

eq = get_eq(dirname)
nu = eq.nu.reshape((1, gi.nr))
rho = eq.rho.reshape((1, gi.nr))
dlnu = eq.dlnu.reshape((1, gi.nr))
dlnrho = eq.dlnrho.reshape((1, gi.nr))
prefactor = -1./(rho*nu*gi.xx**2.)

# get diffrot and its derivs
dOmdr = prefactor*amom_visc_r
dOmdt = prefactor*amom_visc_t
dOmdl = gi.cost_2d*dOmdt + gi.sint_2d*dOmdr

# get the full viscous torque
taunu = vals[:, :, lut[1804]]

# now get the two first-derivative terms
drhonu = nu*rho*(dlnrho + dlnu)
term1 = drhonu*gi.xx**2*dOmdr
term2 = 2*gi.xx*rho*nu*dOmdl

# try to get Laplacian via numerical differentiation
# 
#Om_lap = drad(dOmdr, gi.rr) + 2./gi.rr*dOmdr +\
#        1./gi.rr*(dth(dOmdt, gi.tt) + gi.cott_2d*dOmdt)
#term3 = rho*nu*gi.xx**2.*Om_lap

# get the main second derivative term we think might be playing a role
Om_dr2 = drad(dOmdr, gi.rr)
Om_dt2 = 1./gi.rr*dth(dOmdt, gi.tt)
term3 = rho*nu*gi.xx**2*Om_dr2
term4 = rho*nu*gi.xx**2*Om_dt2

# terms to plot and sub-titles
terms = [term1, term2, taunu - term1 - term2, term3, term4, taunu]
#terms = [term1, term2, term3 -term2, taunu - term1 - term2, taunu]
kw_plot_azav_grid.nrow = 1
kw_plot_azav_grid.titles =\
        [r'$[(\partial/\partial r)(\tilde{\rho}\tilde{\nu})]\lambda^2\partial\Omega/\partial r$',\
        r'$2\lambda\tilde{\rho}\tilde{\nu}\partial\Omega/\partial\lambda$',\
        r'$\tilde{\rho}\tilde{\nu}\lambda^2\nabla^2\Omega$',\
        r'$\tilde{\rho}\tilde{\nu}\lambda^2\partial^2\Omega/\partial r^2$',\
        r'$\tilde{\rho}\tilde{\nu}\sin^2\theta\partial^2\Omega/\partial \theta^2$',\
        r'$\tau_v=\nabla\cdot(\tilde{\rho}\tilde{\nu}\lambda^2\nabla\Omega)$']

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
kw_plot_azav_grid.maintitle = dirname_stripped + ('\n%viscous torque decomposition\n') + time_string

# make figure using usual routine
fig = plot_azav_grid (terms, gi.rr, gi.cost, **kw_plot_azav_grid)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
