# Author: Loren Matilsky
# Created: 01/27/2023
# This script plots the radial energy fluxes as functions of radius
# using the Shell_Avgs data

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
#eq = get_eq(dirname)
#c4 = eq.constants[3]
#c8 = eq.constants[7]
#c10 = eq.constants[9]

# allowed args + defaults
kw_lineplot_default['legfrac'] = 0.3
kw_lineplot_default['plotleg'] = True
kw_make_figure_default.update(lineplot_fig_dimensions)

kw_default = dict({'the_file': None})
kw_default.update(kw_make_figure_default)

kw = update_dict(kw_default, clas)
kw_lineplot = update_dict(kw_lineplot_default, clas)
kw_make_figure = update_dict(kw_kw_make_figure_default, clas)

if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
if kw.mark_bcz: # make room for the bcz label
    kw_make_figure.margin_top_inches += default_line_height
find_bad_keys(kw_default, clas, clas0['routinename'], justwarn=True)

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

#  get all the fluxes
F_RS = vals[:, 0, lut[1807]]
F_MC = vals[:, 0, lut[1809]] + vals[:, 0, lut[1811]]
F_V = vals[:, 0, lut[1813]] # remember the negative HAS BEEN included here

F_tot = F_RS + F_MC + F_V # compute the total flux

profiles = [F_RS, F_MC, F_V]
kw_lineplot.labels = ['Reyn. Stress', 'Mer. Circ.', 'Visc.']

if clas0['magnetism']:
    # A Space Oddysey is actually (-c_4) TIMES the correct Poynting flux
    F_MS = vals[:, 0, lut[1815]]
    F_MM = vals[:, 0, lut[1817]]
    profiles += [F_MS, F_MM]
    kw_lineplot.labels += ['Max. Str.', 'Mean Mag.']
    F_tot += (F_MS + F_MM)

profiles.append(F_tot)
kw_lineplot.labels.append('total')

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel  = 'ang. mom. flux'
#kw_lineplot.ylabel = r'$4\pi r^2$' + '(flux)/L'

lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'radial L flux' + '\n' + time_string

margin_x = fpar['margin_left'] + fpar['sub_margin_left']
margin_y = default_margin/fpar['height_inches']
fig.text(margin_x, 1 - margin_y, the_title,\
         ha='left', va='top', fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print (buff_line)
    print ('saving figure at:')
    print(savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
print (buff_line)
