# Author: Loren Matilsky
# Created: 02/10/2025
#
# Description: Script to plot diff rot, mc, and <S> simultaneously, 
# like Nick's old routine

# to print more info about the differential rotation, run with
#      --diffrot
# to normalize thermal variables by the background state, run with
#      --nond
# To subtract the spherical mean from the thermal variables, run with
#      --sub

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
dirname_stripped = strip_dirname(dirname)

# get desired shells to average over for DR numbers
rvals = clas.rvals
if rvals is None:
    rvals = interpret_rvals(dirname, ['rmin', 'rmax'])
nshells = len(rvals) - 1

# allowed args + defaults
# key unique to this script
kwargs_default = dict({'the_file': None, 'the_file2': None, 'nond': False, 'sub': False, 'verbose': False})

# also need make figure kwargs
plot_azav_grid_kwargs_default['margin_top_inches'] = 1
azav_fig_dimensions['margin_top_inches'] += 0.5*nshells
nlines = get_num_lines(clas0.dirname_label)
azav_fig_dimensions['margin_top_inches'] += (nlines-1)*default_line_height
make_figure_kwargs_default.update(azav_fig_dimensions)


# of course, plot_azav_grid kwargs
kwargs_default.update(plot_azav_grid_kwargs_default)

# overwrite defaults
kw = update_dict(kwargs_default, clas)
kw_plot_azav_grid = update_dict(plot_azav_grid_kwargs_default, clas)

# check for bad keys
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)
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
rr = gi.rr
nr = gi.nr
cost = gi.cost

# reference state variables
eq = get_eq(dirname)
rho_2d = (eq.rho).reshape((1, nr))
prs_2d = (eq.prs).reshape((1, nr))
tmp_2d = (eq.tmp).reshape((1, nr))

# Compute the NOND zonally averaged thermo. vars
ent_az = vals[:, :, lut[501]]

# need some more info, if it's nonD

# Calculate temp from EOS
eq.gamma=5/3
tmp_az = (eq.gamma - 1)/eq.gamma*prs_az + ent_az 

# Compute the spherically averaged thermo. vars
ent_sph = (vals_sph[:, 0, lut_sph[501]]).reshape((1, nr))
prs_sph = (vals_sph[:, 0, lut_sph[502]]/eq.prs).reshape((1, nr))
tmp_sph = (eq.gamma - 1)/eq.gamma*prs_sph + ent_sph

if kw.sub:
    # Now subtract the spherical mean from the zonal mean
    ent_az -= ent_sph
    prs_az -= prs_sph
    tmp_az -= tmp_sph

# set the plot name (base of it) here
basename = 'thermo'
if kw.nond:
    terms = [ent_az, prs_az, tmp_az]
    titles = [r'$S/c_P$', r'$P/\overline{P}$', r'$T/\overline{T}$']
    basename += '_nond'
    if kw.sub:
        titletag = '(nondimensional, sub. sph.)'
    else:
        titletag = '(nondimensional, full field)'
else:
    terms = [ent_az, prs_az*prs_2d, tmp_az*tmp_2d]
    titles = ['S', 'P', 'T']
    basename += '_dim'
    if kw.sub:
        titletag = '(dimensional, sub. sph.)'
    else:
        titletag = '(dimensional, full field)'

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'thermal variables: AZ_Avgs - Shell_Avgs' + '\n' + titletag +\
        '\n' + time_string
if not kw.rcut is None:
    maintitle += '\nrcut = %1.3e' %kw.rcut
kw_plot_azav_grid.maintitle = maintitle
kw_plot_azav_grid.titles = titles

# make figure using usual routine

# for diff rot.
#plot_azav_kwargs_default['plotlatlines'] = False
fig = plot_azav_grid (terms, rr, cost, **kw_plot_azav_grid)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + basename + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
