# Author: Loren Matilsky
# Created: 03/09/2023
# This script plots the spherically aveaged terms in the thermal energy 
# (or entropy) equation as functions of radius using the Shell_Avgs data

import numpy as np
import sys, os
sys.path.append(os.environ['raco'])
sys.path.append(os.environ['rapl'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# equation coefficients
eq = get_eq(dirname)

# allowed args + defaults
lineplot_kwargs_default['legfrac'] = 0.3
lineplot_kwargs_default['plotleg'] = True
make_figure_kwargs_default.update(lineplot_fig_dimensions)

kwargs_default = dict({'the_file': None,  'entropy': False})
kwargs_default.update(make_figure_kwargs_default)

kw = update_dict(kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)

if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
if kw.mark_bcz: # make room for the bcz label
    kw_make_figure.margin_top_inches += default_line_height
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get data
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
print (buff_line)
print ('running ' + clas0['routinename'])
print (buff_line)
print ('reading ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
qv = di['qv']
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
nr = di_grid['nr']

# get rho*T
rhot = eq.rho*eq.tmp

# Determine the simulation is magnetic
magnetism = get_parameter(dirname, 'magnetism')

# Make the plot name, labelling the first/last iterations we average over
if kw.entropy:
    basename = 'entropy_equation_shav'
    baselabel = 'entropy equation'
else:
    basename = 'thermal_energy_shav'
    baselabel = 'thermal energy equation'

prs_work = -vals[:, 0, lut[1901]]
buoy_work = vals[:, 0, lut[1904]]
visc_work = vals[:, 0, lut[1907]]
advec_work = vals[:, 0, lut[1910]]

profiles = [prs_work, buoy_work, visc_work, advec_work]
kw_lineplot.labels = ['pressure work', 'buoy. work', 'visc. work', 'advec. work']
tot_heating = advec_tot + cond_heating + int_heating + visc_heating

if magnetism:
    joule_heating = vals[:, 0, lut[1436]]
    profiles.append(joule_heating)
    tot_heating += joule_heating
    kw_lineplot.labels.append('Joule')

profiles.append(tot_heating)
kw_lineplot.labels.append('total')

# create the plot; start with plotting all the heating terms
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = 'heating per vol.'

# change some things if plotting entropy
if kw.entropy:
    for i in range(len(profiles)):
        profiles[i] /= rhot
        if kw_lineplot.labels[i] == 'Q(r)':
            kw_lineplot.labels[i] == 'Q(r)/(rho*T)'
    kw_lineplot.ylabel = 'dS/dt'

lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  baselabel + '\n' + time_string
if kw.mark_bcz:
    the_title += ('\n' + r'$r_{BCZ} = %1.3e$' %rbcz_est)
    the_title += ('\n' + r'$r_{os} = %1.3e$' %rov_est)

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + basename + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print (buff_line)
    print ('saving figure at:')
    print(savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
print (buff_line)
