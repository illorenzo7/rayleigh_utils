# Author: Loren Matilsky
# Created: 12/19/2022
# Updated: 02/10/2025
#
# Description: Script to plot thermal <P>, <S>, and <T>
#
# to normalize thermal variables by the background state, run with
#      --rel
# to normalize using non-default dissipation number, run with:
#      --nrho
#      --beta
#      --gamma
# To subtract the spherical mean from the thermal variables, run with
#      --sub

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapp'])
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['raco'] + '/reference_state')
from arbitrary_atmosphere import compute_Di_v
from azav_util import *
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
rotation = clas0['rotation']
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
# key unique to this script
kw_default = dict({'the_file': None, 'rel': False, 'sub': False, 'nrho': 3, 'beta': 0.759, 'gamma': gamma_ideal})

# also need make figure kw

# of course, plot_azav_grid kw
kw_plot_azav_grid_default['margin_top_inches'] = 1
kw_default.update(kw_plot_azav_grid_default)

# overwrite defaults
kw = update_dict(kw_default, clas)
kw_plot_azav_grid = update_dict(kw_plot_azav_grid_default, clas)

# check for bad keys
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)
if not kw.rcut is None:  
    # need room for two colorbars and line up top stating rcut 
    kw_plot_azav_grid.margin_top_inches += 1/4
    kw_plot_azav_grid.sub_margin_bottom_inches *= 2

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'AZ_Avgs')
print ('reading ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

# Get necessary grid info
gi = get_grid_info(dirname)

# Compute the zonally averaged thermo. vars
entr = vals[:, :, lut[501]]
prs = vals[:, :, lut[502]]

# reference state variables
eq = get_eq(dirname)
reftype = eq.reference_type
rho_2d = (eq.rho).reshape((1, gi.nr))
prs_2d = (eq.prs).reshape((1, gi.nr))
tmp_2d = (eq.tmp).reshape((1, gi.nr))

# to get relative thermal perturbations (and thus temp.)
# need to do some special stuff here
# get the conversion constants, k_s and k_p
if reftype == 2:
    c_p = get_parameter(dirname, 'pressure_specific_heat')
    k_s = 1./c_p
    k_p = k_T = 1.
elif reftype in [4, 5]: 
    c2 = eq.constants[1] # Ra/Pr or Ro_c^2
    diss = compute_Di_v(kw.gamma, kw.beta, kw.nrho)
    if rotation and reftype == 4: 
        # assume my custom states are nondimensionalized by rotational time
        k = 1.19e-5
    else: # assume nondimensionalized by viscous time
        k = 6.38e-12
    k_s = k_T = c2*k
    k_p = kw.gamma/(kw.gamma-1.)*diss*k

# calculate relative thermal perturbations
entr_rel = k_s*entr
prs_rel = k_p*prs/prs_2d

# Calculate temperature from EOS
tmp_rel = (kw.gamma - 1)/kw.gamma * prs_rel + entr_rel
tmp = tmp_rel/k_T # this is kind of a cluge and unphysical,
# but clearly k_S is the logical conversion factor for T as well

if kw.sub:
    # compute the spherically averaged thermo. vars
    entr_sph = np.sum(entr*gi.tw_2d, axis=0)
    prs_sph = np.sum(prs*gi.tw_2d, axis=0)
    tmp_sph = np.sum(tmp*gi.tw_2d, axis=0)
    # subtract the spherical mean from the zonal mean
    entr -= entr_sph.reshape((1, gi.nr))
    prs -= prs_sph.reshape((1, gi.nr))
    tmp -= tmp_sph.reshape((1, gi.nr))

# set the plot name (base of it) here
basename = 'thermo'
if kw.rel:
    terms = [entr_rel, prs_rel, tmp_rel]
    titles = [r'$\hat{s}/c_p$', r'$\hat{p}/\tilde{p}$', r'$\hat{T}/\tilde{T}$']
    basename += '_rel'
    if kw.sub:
        titletag = '(relative, sub. sph. mean)'
    else:
        titletag = '(relative, full field)'
else:
    terms = [entr, prs, tmp]
    titles = [r'$\hat{s}$', r'$\hat{p}$', r'$\hat{T}$']
    basename += '_dim'
    if kw.sub:
        titletag = '(dimensional, sub. sph. mean)'
    else:
        titletag = '(dimensional, full field)'

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'thermal variables' + '\n' + titletag +\
        '\n' + time_string
if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut
kw_plot_azav_grid.maintitle = maintitle
kw_plot_azav_grid.titles = titles

# make figure using usual routine
fig = plot_azav_grid (terms, gi.rr, gi.cost, **kw_plot_azav_grid)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + basename + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
