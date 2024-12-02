# Author: Loren Matilsky
# Created: 10/31/2019

# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.environ['rapl'])
sys.path.append(os.environ['raco'])
from common import *
from plotcommon import *
from cla_util import *

# Get CLAs
args = sys.argv
clas0, clas = read_clas(args)
dirname = clas0['dirname']
dirname_stripped = strip_dirname(dirname)

# allowed args + defaults
kwargs_default = dict({'the_file': None})
kwargs_default.update(make_figure_kwargs_default)
lineplot_kwargs_default['legfrac'] = 0.25
lineplot_kwargs_default['buff_ignore'] = buff_frac # ignore nastiness 
# at endpoints
kwargs_default.update(lineplot_kwargs_default)
kw = update_dict(kwargs_default, clas)
kw_make_figure = update_dict(make_figure_kwargs_default, clas)
kw_lineplot = update_dict(lineplot_kwargs_default, clas)
if not kw.xcut is None: # make room for label on right
    kw_make_figure.sub_margin_right_inches = default_margin_xlabel
    kw_make_figure.margin_top_inches += default_line_height
find_bad_keys(kwargs_default, clas, clas0['routinename'], justwarn=True)

# get data and grid
if kw.the_file is None:
    kw.the_file = get_widest_range_file(clas0['datadir'], 'Shell_Avgs')
print ('Getting data from ' + kw.the_file)
di = get_dict(kw.the_file)
vals = di['vals']
lut = di['lut']
di_grid = get_grid_info(dirname)
rr = di_grid['rr']
nr = di_grid['nr']

# Create the plot; start with plotting all the energy fluxes
kw_make_figure.ncol = 3

# row of three figures
fig, axs, fpar = make_figure(**kw_make_figure)

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = '4 pi r^2 rho <L>_sph'
kw_lineplot.plotleg = True

# plot the integrated angular momentum (in rotating frame)
ax = axs[0,0]
L_shell = 4*np.pi*rr**2 * vals[:, 0, lut[1819]]
profiles = [L_shell]
kw_lineplot.labels = ['horiz. integrated L (rot. frame)']
lineplot(rr, profiles, ax, **kw_lineplot)

# plot the solid-body rate of the frame
# since everything is nondimensional, set \Omega_0 = 1/2

# get rho
eq = get_eq(dirname)
rho = eq.rho

# compute solid body angular momentum
L_solid = 4*np.pi/3*rho*rr**4
ax = axs[0,1]
profiles = [L_solid]
kw_lineplot.labels = ['solid-body angular momentum']
lineplot(rr, profiles, ax, **kw_lineplot)

# finally, do the relative variation 
ax = axs[0,2]
profiles = [L_shell/L_solid]
kw_lineplot.labels = ['L_shell/L_solid']
lineplot(rr, profiles, ax, **kw_lineplot)


#
# set up profiles and labels
# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'horizontally integrated angular momentum' +\
        '\n' + time_string
axs[0,0].set_title(the_title, fontsize=default_titlesize)

# mark zero lines
for ax in axs.flatten():
    mark_axis_vals(ax, 'x')
    mark_axis_vals(ax, 'y')

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
