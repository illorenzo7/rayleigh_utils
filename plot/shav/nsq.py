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
kwargs_default = dict({'the_file': None, 'mark_bcz': False})
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
if kw.mark_bcz: # make room for the bcz label
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

# get rho
eq = get_eq(dirname)
reftype = eq.reference_type

# mean entropy gradient
dsdr = vals[:, 0, lut[507]]

# compute prefactor, based on reference type
if reftype == 5:
    ek = 1./eq.constants[0]
    ra_over_pr = eq.constants[1]
    reduced_ra = ra_over_pr*ek**2
    nsq = reduced_ra*eq.functions[1]*dsdr

profiles = [nsq]
#kw_lineplot.labels = [r'$u_r0$', r'$u_\theta$', r'$u_\phi$', r'$\mathbf{u}$']

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = r'$N^2/4\Omega_0^2$'
kw_lineplot.plotleg = True
lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  ' buoyancy frequency' + '\n' + time_string
ax.set_title(the_title, fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
