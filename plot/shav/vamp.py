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
kwargs_default = dict({'the_file': None, 'type': 'fluc'})
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

# get rho
eq = get_eq(dirname)

# Convective velocity amplitudes, get these from KE
rke_fluc = vals[:, 0, lut[410]]
tke_fluc = vals[:, 0, lut[411]]
pke_fluc = vals[:, 0, lut[412]]

# Total velocity amplitudes, get these from KE
rke = vals[:, 0, lut[402]]
tke = vals[:, 0, lut[403]]
pke = vals[:, 0, lut[404]]

# mean velocity amplitudes from tot - fluc
rke_mean = rke - rke_fluc
tke_mean = tke - tke_fluc
pke_mean = pke - pke_fluc

if kw.type == 'fluc':
    rke_use = rke_fluc
    tke_use = tke_fluc
    pke_use = pke_fluc
elif kw.type == 'tot':
    rke_use = rke
    tke_use = tke
    pke_use = pke
elif kw.type == 'mean':
    rke_use = rke_mean
    tke_use = tke_mean
    pke_use = pke_mean
elif kw.type == 'rat': # compute a bunch of ratios
    vsq = (rke + tke + pke)/eq.rho # total
    vsq_fluc = (rke_fluc + tke_fluc + pke_fluc)/eq.rho # convect.
    vsq_fluc_r = rke_fluc/eq.rho # vertical convect.
    vsq_fluc_h = (tke_fluc + pke_fluc)/eq.rho # horizontal convect.
    vsq_dr = (pke_mean)/eq.rho # diff. rot.
    vsq_mc = (rke_mean + tke_mean)/eq.rho # mer. circ.

    amp_v = np.sqrt(vsq)
    amp_fluc = np.sqrt(vsq_fluc)
    amp_fluc_r = np.sqrt(vsq_fluc_r)
    amp_fluc_h = np.sqrt(vsq_fluc_h)
    amp_dr = np.sqrt(vsq_dr)
    amp_mc = np.sqrt(vsq_mc)

else:
    print("type must be fluc, mean, tot, or rat")
    print("exiting now")
    sys.exit()

# compute desired profiles
if kw.type == 'rat':
    profiles = [amp_fluc/amp_v, amp_dr/amp_v, amp_mc/amp_v, amp_fluc_r/amp_fluc_h]
    kw_lineplot.labels = [r'$|\mathbf{u}^\prime|/|\mathbf{u}|$', r'$|\mathbf{u}_{DR}|/|\mathbf{u}|$', r'$|\mathbf{u}_{MC}|/|\mathbf{u}|$', r'$|u_r^\prime|/|\mathbf{u}_h^\prime|$']
else:
    vsq_r = rke_use/eq.rho
    vsq_t = tke_use/eq.rho
    vsq_p = pke_use/eq.rho
    vsq = vsq_r + vsq_t + vsq_p

    amp_v = np.sqrt(vsq)
    amp_vr = np.sqrt(vsq_r)
    amp_vt = np.sqrt(vsq_t)
    amp_vp = np.sqrt(vsq_p)

    profiles = [amp_vr, amp_vt, amp_vp, amp_v]
    kw_lineplot.labels = [r'$u_r$', r'$u_\theta$', r'$u_\phi$', r'$\mathbf{u}$']

# Create the plot; start with plotting all the energy fluxes
fig, axs, fpar = make_figure(**kw_make_figure)
ax = axs[0,0]

# x and y labels
kw_lineplot.xlabel = 'radius'
kw_lineplot.ylabel = 'rms velocity'
kw_lineplot.plotleg = True
lineplot(rr, profiles, ax, **kw_lineplot)

# make title 
iter1, iter2 = get_iters_from_file(kw.the_file)
time_string = get_time_string(dirname, iter1, iter2) 
the_title = dirname_stripped + '\n' +  ' velocity amplitudes (' +\
        kw.type + ')' + '\n' + time_string
ax.set_title(the_title, fontsize=default_titlesize)

# save the figure
plotdir = my_mkdir(clas0['plotdir'])
savefile = plotdir + clas0['routinename'] + '_' + kw.type + clas0['tag'] + '-' + str(iter1).zfill(8) + '_' + str(iter2).zfill(8) + '.png'

if clas0['saveplot']:
    print ('saving figure at ' + savefile)
    plt.savefig(savefile, dpi=300)
if clas0['showplot']:
    plt.show()
