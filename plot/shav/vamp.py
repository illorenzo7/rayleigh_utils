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
rho = eq.density

# Convective velocity amplitudes, get these from KE
frke = vals[:, 0, lut[410]]
ftke = vals[:, 0, lut[411]]
fpke = vals[:, 0, lut[412]]

vsq_r = frke/rho
vsq_t = ftke/rho
vsq_p = fpke/rho
vsq = vsq_r + vsq_t + vsq_p

amp_v = np.sqrt(vsq)/100.
amp_vr = np.sqrt(vsq_r)/100.
amp_vt = np.sqrt(vsq_t)/100.
amp_vp = np.sqrt(vsq_p)/100.

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111)

profiles = [amp_vr, amp_vt, amp_vp, amp_v]
kw_lineplot.labels = ['r', r'$\theta$', r'$\phi$', 'tot']

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = r'$r/R_\odot$'
kw_lineplot.ylabel = 'velocity (m/s)'
lineplot(rr/rsun, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  'radial energy flux' + '\n' + time_string
ax.set_title(the_title, fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
