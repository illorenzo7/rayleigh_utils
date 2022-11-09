# Author: Loren Matilsky
# Created: 01/29/2019
# plots <P>, <S> (spherical average subtracted) 
# and derives <rho>, <T>
# to normalize by the background reference state, select
# --nond

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

# allowed args + defaults
# key unique to this script
kwargs_default = dict({'the_file': None, 'the_file2': None, 'nond': False})

# also need make figure kwargs
make_figure_kwargs_default.update(azav_fig_dimensions)
kwargs_default.update(make_figure_kwargs_default)

# of course, plot_azav kwargs
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
print ('reading ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']

if kw.the_file2 is None:
    kw.the_file2 = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
    
print ('reading ' + kw.the_file2)
di_sph = get_dict(kw.the_file2)
vals_sph = di_sph['vals']
lut_sph = di_sph['lut']

# Get necessary grid info
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
nr = di_grid['nr']
cost = di_grid['cost']

# reference state variables
eq = get_eq(dirname)
rho_2d = (eq.rho).reshape((1, nr))
prs_2d = (eq.prs).reshape((1, nr))
tmp_2d = (eq.tmp).reshape((1, nr))

# Compute the NOND zonally averaged thermo. vars
ent_az = vals[:, :, lut[501]]/eq.c_p
prs_az = vals[:, :, lut[502]]/prs_2d

# Calculate temp from EOS
tmp_az = (eq.gamma - 1)/eq.gamma*prs_az + ent_az 

# Compute the spherically averaged thermo. vars
ent_sph = (vals_sph[:, 0, lut_sph[501]]/eq.c_p).reshape((1, nr))
prs_sph = (vals_sph[:, 0, lut_sph[502]]/eq.prs).reshape((1, nr))
tmp_sph = (eq.gamma - 1)/eq.gamma*prs_sph + ent_sph

# Now subtract the spherical mean from the zonal mean
ent = ent_az - ent_sph
prs = prs_az - prs_sph
tmp = tmp_az - tmp_sph

# set the plot name (base of it) here
basename = 'thermo'
if kw.nond:
    terms = [ent, prs, tmp]
    titles = [r'$S/c_P$', r'$P/\overline{P}$', r'$T/\overline{T}$']
    basename += '_nond'
    titletag = '(nondimensional)'
else:
    terms = [ent*eq.c_p, prs*prs_2d, tmp*tmp_2d]
    titles = ['S', 'P', 'T']
    basename += '_dim'
    titletag = '(dimensional)'

# make the main title
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2)
maintitle = dirname_stripped + '\n' +\
        'Thermal variables: Az. Avg. - Sph. Avg.' + '\n' + titletag +\
        '\n' + time_string

# make figure using usual routine
fig = plot_azav_grid (terms, rr, cost, maintitle=maintitle, titles=titles, **kw_plot_azav)

# save the figure
plotdir = my_mkdir(clas0['plotdir'] + 'azav/')
savefile = plotdir + basename + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
plt.close()
